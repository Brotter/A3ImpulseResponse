"""
I split up Ben's signal chain stuff and antenna stuff for usability and modded both for a4
"""
import os as os
import numpy as np
import pylab as lab

from scipy import signal
from scipy.interpolate import interp1d, Akima1DInterpolator

from glob import glob

import tfUtils as tfu
import transferFunctionFinal as tff
import transferFunctionA4add as a4

import cmath #for complex sqrt

import time #for drawing with a pause

try:
    import seaborn as sns #just a pretty plotting package, it isn't important
except:
    print("You don't have seaborn but thats okay!")


#==============================================================================
# Directories which are referred to in this script.
#==============================================================================

# Directory where data saved in GitHub reposity for A3ImpulseResponse.
A3Dir = "./"

# Directories with antenna pulse data.
localDir = "../calibrationLinks/"
#localDir = "/Volumes/ANITA3Data/"

# Where the cable info is stored.
cablesBaseDir = localDir + "antarctica14/S21ExternalRFChainForSingleChannelCallibration/"
#cablesBaseDir = "/Volumes/ANITA3Data/antarctica14/S21ExternalRFChainForSingleChannelCallibration/"
palCablesBaseDir = os.environ['PALESTINE_DATA']
#palCablesBaseDir = localDir + "palestine14/SeaveyAntennas/S21s/"
palAntS11sDir = localDir + "palestine14/SeaveyAntennas/S11s/"


# Signal chain data (57dB seems to be a happy middle power).
waveformDir = localDir + "waveforms/"
#waveformDir = "/Users/brotter/benCode/impulseResponse/integratedTF/waveforms_57dB/"

roofDir = "Rooftop_Seevey_Antenna/Antenna_Impulse_Testing/"
currLocalRoofDir = localDir + roofDir + "Attempt_2_02July2012/Hpol/"
remoteRoofDir = "/storageb/ANITA-3/" + roofDir
currRemoteRoofDir = remoteRoofDir + "Attempt_2_02July2012/Hpol/"

chamberDir = "Anechoic_Chamber_New_Horn_Impulse/ANITA2Antennas/HPol_Pulse"
remoteChamberDir = "/storage/Krauss_Backup/ANITA/" + chamberDir

chamberIdentDir = "Anechoic_Chamber_New_Horn_Impulse/IdenticalAntennas/Ant2Blue_to_Ant0/HPol_Pulse/HPol_Pulse_New_Cable/"
chamberRefPulse = localDir + chamberIdentDir + "Impulse_NewRedCable_23May2012.dat"

chamberBoresightPulse = localDir + chamberIdentDir + "Impulse_p180Deg_t090Deg_22May2012_H_Two.dat"




#A global constant that I don't want to add in a bunch of places because I should use classes
suffix="_SimpVoltCalib"


#==============================================================================
# Importing things
#==============================================================================


def importSurf(chan):
    fileName = "averageSignalChainGraphs/" + chan[2:] + str(int(chan[:2])) + "wf.txt"
    graphT,graphV = np.loadtxt(fileName, delimiter = " ").T
    
#    dataY /= 1000 #mv->V  #not mV!  It is in ADC counts so I shouldn't touch this...

    #make the time range start at zero
    graphT -= graphT[0]

    dataF,dataFFT = tfu.genFFT(graphT,graphV)
        
    return graphT, graphV, dataF, dataFFT
    

def importScope():
    #Calibration Scope
    fileName = "scope_data/fastImpulseAvg017.csv"
    dataX, dataY = np.loadtxt(fileName, delimiter=",", usecols=(3,4)).T

    dataX *= 1e9 #s ->ns

    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tfu.genFFT(dataX,dataY)
        
    return dataX, dataY, dataF, dataFFT


def getRegeneratedCables(fileName,dF=5./512.,length=513):
    """
    New getRegeneratedCables():

    imports the data (correctly with the s2pParser fix), then resamples the gain and phase
    to meet the specified dF and length (which should stay constant because the time domain I choose defines
    them)
    phase -> just a constant group delay (and I don't care about absolute phase)
    gain -> interpolated inside the range you can, zero outside? (sure why not)

    some of this was moved to getRegeneratedCables

    Lets skip causality for now (that seems to be a bigger problem at the end anyway)

    """

    cableFreq,cableFFT = getCables(fileName)

    cableNewFreq,cableNewFFT = tfu.regenerateCable(cableFreq,cableFFT,dF=dF,length=length)
    

    print("getRegeneratedCables(): dFin=:",cableFreq[1]-cableFreq[0])
    print("getRegeneratedCables(): dFout=:",cableNewFreq[1]-cableNewFreq[0])
    print("getRegeneratedCables(): len(f)in=:",len(cableFreq))
    print("getRegeneratedCables(): len(f)out:",len(cableNewFreq))



    return cableNewFreq, cableNewFFT


def getCables(fileName):
    """
    basically just import the cables.

    *network analyzer cable measurements (not all are)*
    
    also transforms it from measured log magnitude to linear magnitude

    """
    cableFreq,cableGainLog,cablePhase = s2pParser(cablesBaseDir + fileName)
    cableGainLin = 10.**(cableGainLog/20.)

    cableFFT = tfu.gainAndPhaseToComplex(cableGainLin,cablePhase)

    return cableFreq,cableFFT

def getCablesForChannel(chan):
    cableStart = 0
    if chan[3] == "H":
        cableStart = 8
    cableIndex = int(chan[:2])%8 + cableStart
    if cableIndex < 10:
        fname = "16way_data/16WAY_CH0{0}.csv".format(cableIndex)
    else:
        fname = "16way_data/16WAY_CH{0}.csv".format(cableIndex)

    cableFreq,cableGainLog,cableGD =  s2pParser(fname)
    cableGainLin = 10.**(cableGainLog/20.)
    return cableFreq,cableGainLin,cableGD

def getExtraCable(refOrTest):
    gdname = "ImpulseTestingUpload/IMPULSE{0}MEASGD.csv".format(refOrTest)
    magname = "ImpulseTestingUpload/IMPULSE{0}MEASGD.csv".format(refOrTest)
    if(refOrTest == "ATTEN"):
        gdname = "ImpulseTestingUpload/AttenSetup-GD.csv"
        magname = "ImpulseTestingUpload/AttenSetup-S.csv"
    if(refOrTest == "EXTRAATTEN"):
        gdname = "ImpulseTestingUpload/AttenSetupExtradB-GD.csv"
        magname = "ImpulseTestingUpload/AttenSetupExtradB-S.csv"
    if(refOrTest == "Z"):
        gdname = "ImpulseTestingUpload/ZGD.csv"
        magname = "ImpulseTestingUpload/Z.csv"

    freq = []
    mag = []
    gd = []
    record = 0
    with open(magname) as inFile:
        for line in inFile:
            split = line.strip().split(',')
            if len(split) == 2:
                if split[1] == "S21":
                    record = 1
            if(line[0] == 'E' and record == 1):
                break
            if(line[0] != '!' and line[0] != 'B' and record == 1):
                freq.append(float(split[0])/1e9)
                mag.append(float(split[1].strip()))
    record = 0
    with open(gdname) as inFile:
        for line in inFile:
            split = line.strip().split(',')
            if len(split) == 2:
                if split[1] == "S21":
                    record = 1
            if(line[0] == 'E' and record == 1):
                break
            if(line[0] != '!' and line[0] != 'B' and record == 1):
                gd.append(float(split[1].strip()))
    npmag = np.array(mag)
    cableGainLin = 10.**(npmag/20.)

    return np.array(freq),cableGainLin,np.array(gd)

def regenerateCablesFromGroupDelay(cableF, cableG, cableGD):
    cableP = phaseFromGroupDelay(cableF, cableGD)
    cableFFT = tfu.gainAndPhaseToComplex(cableG,cableP)
    newCableF,newCableFFT = tfu.regenerateCable(cableF,cableFFT,dF=5./512.,length=513)
    return newCableF, newCableFFT

def s2pParser(fileName):
    """
      reads in s21 files
    """

    freq = []
    mag = []
    gd = []
    record = 0
    with open(fileName) as inFile:
        for line in inFile:
            split = line.strip().split(',')
            if len(split) == 2:
                if split[1] == "S21":
                    record += 1
            if(line[0] == 'E' and record == 1):
                record = 2 
            if(record == 0):
                continue
            if(line[0] != '!' and line[0] != 'B' and record == 1):
                freq.append(float(split[0])/1e9)
                mag.append(float(split[1].strip()))
            if(line[0] != '!' and line[0] != 'B' and line[0] != 'E' and record == 3):
                gd.append(float(split[1].strip())*1e9)

    return np.array(freq), np.array(mag), np.array(gd)




def s11Parser(fileName):

    data = np.genfromtxt(fileName,delimiter=",",comments="!",skip_footer=1,skip_header=18).T
    
    f = data[0]
    fft = data[1] + data[2]*1j

    return f,fft

#==============================================================================
# Processing
#==============================================================================

P2SF = False
P2SFFT = False
P2AF = False
P2AFFT = False

def doSigChainWithCables(chan,showPlots=False,savePlots=False,writeFiles=False):
#    if savePlots or showPlots:
#        lab.close("all")

    #phase shifts are likely pointless!
    #they were: scopeFFT=1.26
    #           phaseShift=1.698

    #get scope data from waveforms (extracted from raw data from another script)
    #NOTE: removed the kSlice stuff for now
    scopeRawX,scopeRawY,scopeRawF,scopeRawFFT = importScope()
    scopeX,scopeY = tff.processWaveform(scopeRawX,scopeRawY,"calScope")
    scopeF,scopeFFT = tfu.genFFT(scopeX,scopeY)
            
    #1.26 is the max coherance phase?
#    scopeFFT = tfu.fftPhaseShift(scopeFFT,1.26)

    #Get the cable's H(f) transfer function for PULSER TO SCOPE
#    global cableScopeF #once upon a time regenerating the cables took awhile
#    global cableScopeFFT
#    if type(cableScopeF) != np.ndarray:
    print( "Getting Pulser to Scope Cables...")
    ###### NOTE: not sure these are the correct cables ######
#        cableScopeF,cableScopeFFT= getCables\\\\\\\("A-B_PULSER-SCOPE.s2p",tZero=17,resample="interp")
    cableScopeRawF, cableScopeRawG, cableScopeRawGD = getExtraCable("Z")
    cableScopeF, cableScopeFFT = regenerateCablesFromGroupDelay(cableScopeRawF, cableScopeRawG, cableScopeRawGD)
    cableScopeX,cableScopeY = tfu.genTimeSeries(cableScopeF,cableScopeFFT)

    #Get the cable's H(f) transfer function for PULSER TO AMPA
#    global P2AF
#    global P2AFFT
#    if type(P2AF) != np.ndarray:
    print( "Getting Pulser to Ampa Cables...")
#        P2AF,P2AFFT= getCables("A-C_PULSER-TEST_66DB.s2p",tZero=True,hanning=True,resample="fftFit")
    cableAmpaRawF,cableAmpaRawG,cableAmpaRawGD = getCablesForChannel(chan)
    cableAmpaF,cableAmpaFFT = regenerateCablesFromGroupDelay(cableAmpaRawF, cableAmpaRawG, cableAmpaRawGD)
    cableAmpaX,cableAmpaY = tfu.genTimeSeries(cableAmpaF,cableAmpaFFT)

    #Get the other cables (ref and test)
    if(chan == "01BV"):
        extraCableF,extraCableG,extraCableGD = getExtraCable("REF")
    else:
        extraCableF,extraCableG,extraCableGD = getExtraCable("TEST")

    extraCableF,extraCableFFT = regenerateCablesFromGroupDelay(extraCableF, extraCableG, extraCableGD)
    extraCableX,extraCableY = tfu.genTimeSeries(extraCableF,extraCableFFT)

    if(chan == "01BV"):
        attenCableF, attenCableG, attenCableGD = getExtraCable("EXTRAATTEN")
    else:
        attenCableF, attenCableG, attenCableGD = getExtraCable("ATTEN")
    attenCableF, attenCableFFT = regenerateCablesFromGroupDelay(attenCableF, attenCableG, attenCableGD)
    #PULSER TO SCOPE
    #deconvolve cable pulse * cable = scope -> pulse = scope/cable
    #finds just the pulser impulse
    pulseDeconvFFT = scopeFFT/cableScopeFFT
    pulseDeconvFFT[np.isposinf(np.real(pulseDeconvFFT))] = 0
    pulseDeconvFFT[np.isneginf(np.real(pulseDeconvFFT))] = 0
    pulseDeconvFFT = np.nan_to_num(pulseDeconvFFT)
    pulseDeconvFFT[pulseDeconvFFT > 1e10] = 0
    pulseDeconvFFT[pulseDeconvFFT < -1e10] = 0
    pulseDeconvX,pulseDeconvY = tfu.genTimeSeries(scopeF,pulseDeconvFFT)
    pulseDeconvY -= np.mean(pulseDeconvY)

    #PULSER TO AMPA
    #convolve pulse with cables transfer function to get pulse at AMPA
    if(chan == "01BV"):
        ampaInputFFT = extraCableFFT*attenCableFFT*pulseDeconvFFT
    else:
        ampaInputFFT = extraCableFFT*cableAmpaFFT*attenCableFFT*pulseDeconvFFT

    ampaInputFFT = np.nan_to_num(ampaInputFFT)
    ampaInputF = scopeF
    #do a time domain thing
    ampaInputX,ampaInputY = tfu.genTimeSeries(ampaInputF,ampaInputFFT)
    ampaInputY = np.roll(ampaInputY,100-np.argmax(ampaInputY))
#    ampaInputY = tfu.hanningTail(ampaInputY,np.argmax(ampaInputY), 20)
#    ampaInputY = tfu.hanningTailNeg(ampaInputY,np.argmax(ampaInputY),20)
    ampaInputY -= np.mean(ampaInputY)

#    return ampaInputY,ampaInputY2

    #then back to fft again
    ampaInputF,ampaInputFFT2 = tfu.genFFT(ampaInputX,ampaInputY)
    #this is a stupid thing with complex division in python and inf/nan values
    ampaInputFFT2[np.absolute(ampaInputFFT2) < 1e-18] = 0
#    ampaInputFFT2[0] = 0


    ampaInputFFT = ampaInputFFT2

#    return ampaInputFFT,ampaInputFFT2
    

    if (showPlots or savePlots) and (chan == "01BH"):
        
        fig2,ax2 = lab.subplots(3,2,figsize=(11,8.5))
        ax2[0][0].set_title("Cables")


        #phase
        ax2[0][0].set_ylabel("phase (radians)")
        
        ax2[0][0].plot(cableScopeRawF,cableScopeRawGD,label="Pulser To Scope Raw",color="green")
        ax2[0][0].plot(cableScopeRawF,cableScopeRawGD,label="Pulser To Scope Processed",color="blue")
        ax2[0][0].legend()

        ax2[0][1].set_ylabel("Group Delay (ns)")
        ax2[0][1].plot(cableAmpaRawF,cableAmpaRawGD,label="Pulser To Ampa Raw",color="green")
        ax2[0][1].plot(cableAmpaRawF,cableAmpaRawGD,label="Pulser To Ampa Processed",color="red")
        ax2[0][1].legend()
        
        #gain
        ax2[1][0].set_ylabel("gain (dB)")
        ax2[1][0].set_xlabel("freq (GHZ)")
        ax2[1][0].plot(cableScopeRawF,cableScopeRawG,label="Pulser To Scope Raw",color="green")
        ax2[1][0].plot(cableScopeF,tfu.calcLogMag(cableScopeF,cableScopeFFT),label="Pulser To Scope Processed",color="blue")
        ax2[1][0].set_ylim([-100,0])
        ax2[1][0].legend()

        ax2[1][1].set_ylabel("gain (dB)")
        ax2[1][1].set_xlabel("freq (GHZ)")
        ax2[1][1].plot(cableAmpaRawF,cableAmpaRawG,label="Pulser To Ampa Raw",color="green")
        ax2[1][1].plot(cableAmpaF,tfu.calcLogMag(cableAmpaF,cableAmpaFFT),label="Pulser To Ampa",color="red")
        ax2[1][1].set_ylim([-100,0])
        ax2[1][1].legend()

        #impulse

        ax2[2][0].set_xlabel("Time (ns)")
        ax2[2][0].set_ylabel("Voltage (V)")
        ax2[2][0].plot(cableScopeX,cableScopeY,label="Cable: Pulser to Scope",color="blue")
        ax2[2][0].plot(cableAmpaX,cableAmpaY,label="Cable: Pulser to AMPA",color="red")
        ax2[2][0].legend()

        ax2[2][1].set_xlabel("Time (ns)")
        ax2[2][1].plot(ampaInputX,ampaInputY,label="Corrected Pulse into AMPA",color="purple")
#        ax2[0][2].set_xlim([0,30])
        ax2[2][1].legend()

        if showPlots:
            fig2.show()
        if savePlots:
            fig2.savefig("plots/doSigChainWithCables_Cables.png")
            



    if savePlots or showPlots:
        fig,ax = lab.subplots(2,2,figsize=(11,8.5))
        fig.set_tight_layout(True)
        ax[0][0].set_title("Raw Signal Chain Calibration Pulser")
        ax[0][0].plot(scopeRawX,scopeRawY,label="raw scope pulse",color="red")
        ax[0][0].plot(scopeX,scopeY,label="processed scope pulse",color="blue")
#        ax[0][0].plot(pulseDeconvX,pulseDeconvY,label="processed scope pulse",color="blue")
        ax[0][0].set_ylabel("Voltage(V)")
        ax[0][0].set_xlabel("Time (ns)")
        ax[0][0].set_xlim([0,200])
        ax[0][0].legend()

        ax[0][1].plot(ampaInputX,ampaInputY,label="cable corrected pulse into ampa",color="purple")
        ax[0][1].set_ylabel("Voltage (mV)")
        ax[0][1].set_xlabel("Time (ns)")
        ax[0][1].legend()

        ax[1][0].plot(scopeRawF,tfu.calcSpecMag(scopeRawF,scopeRawFFT),label="raw scope pulse",color="red")
        ax[1][0].plot(scopeF,tfu.calcSpecMag(scopeF,scopeFFT),label="processed scope pulse",color="blue")
        ax[1][0].set_xlabel("Frequency (GHz)")
        ax[1][0].set_ylabel("Spectral Power (dBm/Hz)")
        ax[1][0].set_ylim([-100,0])
        ax[1][0].legend()

        ax[1][1].plot(ampaInputF,tfu.calcSpecMag(ampaInputF,ampaInputFFT),label="cable corrected pulse into ampa",color="purple")
        ax[1][1].set_xlabel("Frequency (GHz)")
        ax[1][1].set_ylim([-200,0])
#        ax[1].set_xlim([0,4])
        ax[1][1].legend()

        if savePlots:
            fig.savefig("plots/doSigChainWithCables_Input"+chan+".png")
        if showPlots:
            fig.show()




    #get the surf (extracted from ROOT data from another script)
    surfRawX,surfRawY,surfRawF,surfRawFFT = importSurf(chan)
    print( " /|\ /|\ /|\ /|\ ",len(surfRawX),np.diff(surfRawX)[0])
    surfX,surfY = tff.processWaveform(surfRawX,surfRawY,"surf")
    print( " /|\ /|\ /|\ /|\ ",len(surfX),np.diff(surfX)[0])
#    surfY = tfu.lowPass(surfX, surfY, 1.3)
#    surfY = tfu.hanningTail(surfY,np.argmax(surfY)+300,100)
#    surfY = tfu.nyquistLimit(surfY, 3)
    surfY -= np.mean(surfY)

    surfY *= (1./1000) #lets go to V for a sec #TODO: this is now ADC counts/1000.
    surfF,surfFFT = tfu.genFFT(surfX,surfY)
#    surfFFT *= 10**(3.53/20.) #I think there's a missing splitter



    if savePlots or showPlots:
        figS,axS = lab.subplots(2,figsize=(11,8.5))
        axS[0].set_title("Upsampled, correlated, and averaged SURF output ("+chan+")")
        max = np.argmax(surfRawY)
        max = 100
        axS[0].plot(surfRawX,np.roll(surfRawY,100-max),label="Raw Surf Waveform",color="red")
        axS0 = axS[0].twinx()
        axS0.plot(surfX,np.roll(surfY*1000.,100-max),'+',label="Processed Surf Waveform",color="blue")
        axS[0].set_xlabel("Time (ns)")
        axS[0].set_ylabel("ADC counts")
        axS0.set_ylabel("Voltage (mV)")
#        axS[0].set_ylim([-0.2,0.2])
        axS[0].legend(loc="upper left")
        axS0.legend()
        axS[1].plot(surfRawF,tfu.calcSpecMag(surfRawF,surfRawFFT),label="Raw Surf Waveform",color="red")
        axS[1].plot(surfF,tfu.calcSpecMag(surfF,surfFFT),'+',label="Processed Surf Waveform",color="blue")
        axS[1].set_xlabel("Frequency (GHz)")
        axS[1].set_ylabel("\"Spectral Power\" (~dBm/Hz)")
        axS[1].set_ylim([-110,50])
        axS[1].set_xlim([0,3])
        axS[1].legend()
        if savePlots:
            figS.savefig("plots/doSigChainWithCables_SURF"+chan+".png")
        if showPlots:
            figS.show()


    print( "*************",len(surfF),":",np.diff(surfF)[0],len(ampaInputF),":",np.diff(ampaInputF)[0])
    print( "*************",len(surfFFT),len(ampaInputFFT))

    #deconvolve signal chain transfer function out of that!
    tfFFT = surfFFT/ampaInputFFT

    
    #fix the infinities and nans
    tfFFT[np.isposinf(np.real(tfFFT))] = 0
    tfFFT[np.isneginf(np.real(tfFFT))] = 0
    tfFFT = np.nan_to_num(tfFFT)
    tfFFT[tfFFT > 1e10] = 0
    tfFFT[tfFFT < -1e10] = 0

#    return surfFFT,ampaInputFFT,tfFFT

    
    #change it back to time domain
#    for i in range(0,len(surfF)):
#        if surfF[i] > 1.25:
#            tfFFT[i] /= 1e9
    tfX,tfY = tfu.genTimeSeries(surfF,tfFFT)

    #clean up the tail and make it start at the beginning
    tfY = np.roll(tfY,30-np.argmax(tfY))
    tfY = tfu.hanningTail(tfY,320,100)
    tfY = tfu.nyquistLimit(tfY,5.)
    tfY -= np.mean(tfY)



    #vvvvvvvvvvvvv   IMPORTANT
    #CLUDGE for normalization.  The RMS and y-factor analyses both suggest that this process is a factor of two low
#    tfY *= 2
    #^^^^^^^^^^^ 
    
    tfF,tfFFT = tfu.genFFT(tfX,tfY)
    surfY -= np.mean(surfY)
    
    ampaInputY *= 1000
    ampaInputY -= np.mean(ampaInputY)

    if showPlots or savePlots:
        fig3,ax3 = lab.subplots(2,2,figsize=(11,8.5))
        fig3.set_tight_layout(True)

#        fig3.suptitle("Input, Output, and Transfer Function Signals")
        ax30 = ax3[0][0].twinx()
        ax30.plot(ampaInputX,ampaInputY,label="Input Pulse, Cables Convolved",color="red")
        ax30.set_ylabel("Voltage (mV)")
        ax30.set_ylim([-0.2,1])
        ax30.legend()
#        surfYMax = np.argmax(surfY)
#        surfYMax = 200
        ax3[0][0].plot(surfX,surfY,label="SURF waveform",color="blue")
        ax3[0][0].set_ylabel("voltage (mV)")
        ax3[0][0].set_ylim([-0.2,0.2])
        ax3[0][0].set_xlabel("time (ns)")
        ax3[0][0].legend(loc="upper left")
#        tfYMax = np.argmax(tfY)
        tfYMax = 100
        ax3[0][1].plot(tfX,np.roll(tfY,100-tfYMax),label="Signal Chain Transfer Function",color="black")
        ax3[0][1].set_ylim([-80,100])
#        ax3[0].set_xlim([20,80])
        ax3[0][1].set_ylabel("gain normalized voltage")
        ax3[0][1].set_xlabel("time (ns)")
        
        ax3[0][1].legend(loc="upper left")

        ax3[1][0].set_ylabel("Spectral Power Magnitude (dBm/Hz)")
        ax3[1][0].set_xlabel("Frequency (GHz)")
        ax3[1][0].set_xlim([0,3])
        ax3[1][0].plot(ampaInputF,tfu.calcSpecMag(ampaInputF,ampaInputFFT),label="Input Pulse, Cables Removed",color="red")
        ax3[1][0].plot(surfRawF,tfu.calcSpecMag(surfF,surfFFT),label="SURF waveform",color="blue")
        ax3[1][0].legend()
        ax3[1][0].set_ylim([-150,20])

        ax3[1][1].set_ylabel("Gain (dB)")
        ax3[1][1].set_xlabel("Frequency (GHz)")
        ax3[1][1].plot(surfF,tfu.calcLogMag(surfF,tfFFT),label="Signal Chain Transfer Function",color="black")
        ax3[1][1].set_xlim([0,3])
        ax3[1][1].legend()
        ax3[1][1].set_ylim([0,100])




        if showPlots:
            fig3.show()
        if savePlots:
            fig3.savefig("plots/doSigChainWithCables_TF"+chan+".png")

#    for i in range(0,len(surfF)):
#        if surfF[i] > 1.3:
#            tfFFT[i] /= 1e9


    if writeFiles:
        tff.writeOne([tfX,tfY],chan,"sigChain")

    return tfX,tfY,surfF,tfFFT 



def doSigChainAndAntenna(chan,showPlots=False,savePlots=False,writeFiles=False):
    #get antenna
#    antX,antY,antF,antFFT = doRoofAntWithCables()
    channelsWithPhaseInfo = ["02TH", "07TH", "12MH", "13MH", "02TV", "07TV", "12MV", "13MV"]
    if chan in channelsWithPhaseInfo:
        antX,antY,antF,antFFT = a4.doPalAnt(chan,showPlots=showPlots,savePlots=savePlots,writeFiles=writeFiles)
    else:
        antX,antY,antF,antFFT = a4.doPalAntWithAveragePhase(chan,showPlots=showPlots,savePlots=savePlots,writeFiles=writeFiles)
    
    #get sig chain
    print("doSigChainAndAntenna:getSigChain")
    sigChainX,sigChainY,sigChainF,sigChainFFT = doSigChainWithCables(chan,showPlots=showPlots,savePlots=savePlots,writeFiles=writeFiles)

    #convolve the two (full ANITA3 transfer function!)
    a4F = sigChainF
    a4FFT = sigChainFFT * antFFT

    #also there is a scaling from free space (377ohms) to system impedance (50ohms) [ sqrt(50/377) ]
    a4FFT *= np.sqrt(50./377)


    a4X,a4Y = tfu.genTimeSeries(a4F,a4FFT)
    
    if (showPlots or savePlots):
#        lab.close("all")
        fig,ax = lab.subplots(2,4,figsize=(11,8.5))
        ax[0][0].set_title(chan)
        ax[0][0].plot(antX,np.roll(antY,100),label="Antenna")
        ax[0][0].legend()
        ax[0][0].set_ylabel("Antenna Height (V/m)")

        ax[1][0].set_title(chan)
        ax[1][0].plot(antF,tfu.calcAntGain(antF,antFFT))
        ax[1][0].set_ylabel("Antenna Gain (dBi)")
        ax[1][0].set_xlim([0,3])

        ax[0][1].plot(sigChainX,np.roll(sigChainY,100),label="sigChain",color="red")
        ax[0][1].set_ylabel("Gain (unitless)")

        ax[1][1].plot(sigChainF,tfu.calcLogMag(sigChainF,sigChainFFT),color="red")
        ax[1][1].set_ylim([30,80])
        ax[1][1].set_ylabel("Gain (dB)")
        ax[1][1].set_xlim([0,3])

        a4YPeak = np.argmax(a4Y)
        ax[0][2].plot(a4X,np.roll(a4Y,100-a4YPeak),label="Convolution",color="black")
        ax[0][2].legend()
        ax[0][2].set_xlabel("Time (ns)")
        ax[0][2].set_ylabel("Instrument Response \n (V/m)")
        
        ax[1][2].plot(a4F,tfu.calcLogMag(a4F,a4FFT),color="black")
        ax[1][2].set_ylim([0,70])
        ax[1][2].set_xlim([0.1,1.5])
        ax[1][2].set_xlabel("Frequency (GHz)")
        ax[1][2].set_ylabel("Instrument Height (V/m)")
        ax[1][2].set_xlim([0,3])
        
        brX,brY,brF,brFFT = a4.importAverageBRotter() 
        
        ax[0][3].plot(brX,np.roll(brY,100-np.argmax(brY)),label="BRotter average",color="black")
        ax[0][3].legend()
        ax[0][3].set_xlabel("Time (ns)")
        ax[0][3].set_ylabel("Instrument Response \n (V/m)")
        
        ax[1][3].plot(brF,tfu.calcLogMag(brF,brFFT),color="black")
        ax[1][3].set_ylim([0,70])
        ax[1][3].set_xlim([0.1,1.5])
        ax[1][3].set_xlabel("Frequency (GHz)")
        ax[1][3].set_ylabel("Instrument Height (V/m)")
        ax[1][3].set_xlim([0,3])

#        ax[4].plot(a4F,tfu.calcPhase(tfu.minimizeGroupDelayFromFFT(a4F,a4FFT,highLim=256)),color="black")
#        ax[4].set_xlim([0,2])
#        ax[4].set_xlabel("Frequency (GHz)")
#        ax[4].set_ylabel("Group Delay (ns)")
#        ax[4].set_ylim([-20,20])

        if showPlots:
            fig.show()
        if savePlots:
            fig.savefig("plots/doSigChainAndAntenna_TF"+chan+".png")

    #clean it up a bit
    #something was adding a ton of 1.25GHz noise
#    for i in range(0,len(a4FFT)):
#        if a4F[i] > 1.1:
#            a4FFT[i] /= 1e4

#    a4FFT = np.concatenate((a4FFT[:171],np.zeros(342))) #this isn't causal...

    #make it look nice and cut off the junk at the end
    a4Y = np.roll(a4Y,40-np.argmax(a4Y))
    a4Y = tfu.hanningTail(a4Y,300,200)
#    a4X,a4Y = tfu.zeroPadEqual(a4X,a4Y,1024)

    if writeFiles:
        tff.writeOne([a4X,a4Y],chan,"fullTF")
        
    
    return a4X,a4Y



def alignWaveforms(allChans,showPlots=False):
    
    #They need to all be aligned and have the same group delay!
    #I try two ways:
    #1) use shiftToAlign, which does the weird fourier linear regression that messes with the phase
    #2) (THE DEFAULT) fit a gaussian to the peak of the correlation plot and use that as the fractional offset, then resample


    outChans = {}
    desiredDelay = 0
    max = np.argmax(allChans["01BH"][1])
    outChans["01BH"] = allChans["01BH"][0],np.roll(allChans["01BH"][1],max-desiredDelay)
    for chan in allChans:
        if chan == "01BH":
            continue
        outChans[chan] = allChans[chan][0],tfu.phaseFitCorrelate(outChans["01BH"],allChans[chan])[0];


    outChans2 = {}
    for chan in allChans:
        print(chan)
        outChans2[chan] = tfu.shiftToAlign(allChans["01BH"],allChans[chan])


    correlations = {}
    for chan in allChans:
        correlations[chan] = tfu.correlation(outChans["01BH"][1],allChans[chan][1])



    if showPlots:
        highWeird = []
        lowWeird = []

        fig,ax = lab.subplots(4)
        for chan in allChans:
            ax[0].text(0.65, 0.9, "Unaligned", transform=ax[0].transAxes, fontsize=14)
            ax[0].plot(allChans[chan][0],allChans[chan][1])
            

            ax[1].text(0.65, 0.9, "linear regression phase subtraction", transform=ax[1].transAxes, fontsize=14)
            ax[1].plot(outChans[chan][0],outChans[chan][1])
                
            ax[2].text(0.65, 0.9, "fit correlation peak fractional offset", transform=ax[2].transAxes, fontsize=14)
            ax[2].plot(outChans2[chan][0],outChans2[chan][1])
            if outChans2[chan][1][552] > 0.03:
                highWeird.append(chan)
            else:
                lowWeird.append(chan)

            ax[3].text(0.65, 0.9,"Correlation with 01BH" , transform=ax[3].transAxes, fontsize=14)
            ax[3].plot(correlations[chan])


    
        highWeird.sort()
        lowWeird.sort()
        print(len(highWeird),highWeird)
        print(len(lowWeird),lowWeird)

        fig.show()

    return outChans2


    


def doAllInstrumentHeight(showPlots=False,savePlots=False,writeFiles=False):
    #renamed from doTheWholeShebang because that isn't very descriptive
    # Generate the transfer function for all the channels!
    chans = np.loadtxt("chanList.txt",dtype=str)
    allChans = {}
#    lab.close("all")
    for chan in chans:
        try:
            allChans[chan] = doSigChainAndAntenna(chan,savePlots=savePlots,showPlots=showPlots,writeFiles=writeFiles)
        except:
            print(chan+" FAILED")


    allChansAligned = alignWaveforms(allChans)
    rollAllDict(allChansAligned,-500)

    return allChansAligned


def doAllSigChains(showPlots=False,savePlots=False):
    # just does all the signal chains without the antennas
    chans = np.loadtxt("chanList.txt",dtype=str)
    allChans = {}
    for chan in chans:
        try:
            allChans[chan] = doSigChainWithCables(chan,savePlots=savePlots,showPlots=showPlots)[:2]
        except:
            print(chan+" FAILED************************************************")

    return allChans



def rollAllDict(allChans,roll):

    """
    
    If you want to roll all the channels in an "allChans" type object

    """

    for chan in allChans.keys():
        allChans[chan][1] = np.roll(allChans[chan][1],roll)
                

def importAttenInfo(showPlots=False):
    
    #Chan | iRFCM attenuator | Split->Short attenuator | phi sector

    attenInfo = np.loadtxt("attenuatorInfo/Map_sorted_by_Antenna-Table 1.csv",dtype=str,delimiter=",",usecols=(0,2,3,4))[:-2].T

    outInfo = {}
    for chan in attenInfo.T:
        outInfo[chan[0]] = chan[1:]

    print(outInfo)

    if showPlots:
        fig,ax = lab.subplots()
        ax.hist(np.array(attenInfo[2],dtype=float),np.arange(0,7)-0.5,color="blue",edgecolor="black")
        fig.show()


    return outInfo

def phaseFromGroupDelay(f, gd):
    phase = []
    phase.append(gd[0]);
    dF = f[1] - f[0]
    for i in range(len(gd)-1):
        phase.append((gd[i+1] + phase[i]))
    npphase = np.array(phase)
    npphase *= -2*np.pi*dF
    return npphase





