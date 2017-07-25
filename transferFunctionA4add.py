"""
Script written by John Russell which includes functions specific to ANITA-4 over
ANITA-3. Takes some functions from Ben Rotter's "transferFunctionFinal.py" and
modifies them for ANTIA-4 usage.
"""


import numpy as np
import tfUtils as tf
import transferFunctionFinal as tff  #  This is to reference what Ben has already done.
from scipy.interpolate import Akima1DInterpolator
import pylab as lab


#==============================================================================
# Directories which are potentially referred to in this script.
#==============================================================================


# Directory where data saved in GitHub reposity for A3ImpulseResponse.
A4Dir = "./"

# Directories with antenna pulse data.
localDir = "../calibrationLinks/"

# Where the cable info is stored.
cablesBaseDir = localDir + "antarctica14/S21ExternalRFChainForSingleChannelCallibration/"

# Signal chain data (57dB seems to be a happy middle power).
waveformDir = localDir + "waveforms/"

roofDir = "Rooftop_Seevey_Antenna/Antenna_Impulse_Testing/"
currLocalRoofDir = localDir + roofDir + "Attempt_2_02July2012/Hpol/"
remoteRoofDir = "/storageb/ANITA-3/" + roofDir
currRemoteRoofDir = remoteRoofDir + "Attempt_2_02July2012/Hpol/"

chamberDir = "Anechoic_Chamber_New_Horn_Impulse/ANITA2Antennas/HPol_Pulse"
remoteChamberDir = "/storage/Krauss_Backup/ANITA/" + chamberDir

chamberIdentDir = "Anechoic_Chamber_New_Horn_Impulse/IdenticalAntennas/Ant2Blue_to_Ant0/HPol_Pulse/HPol_Pulse_New_Cable/"
chamberRefPulse = localDir + chamberIdentDir + "Impulse_NewRedCable_23May2012.dat"

chamberBoresightPulse = localDir + chamberIdentDir + "Impulse_p180Deg_t090Deg_22May2012_H_Two.dat"


#==============================================================================
# Importing things
#==============================================================================


def importANITA1():
    a1 = np.loadtxt(A4Dir + "ANITA1ImpulseResponse.dat")
#    a1 = np.loadtxt("ANITA1ImpulseResponse.dat")
    return a1.T[0]*1e9, a1.T[1]


def findPalestineAntennaFile(chan, inOrOut):
    antennaInfo = np.loadtxt(A4Dir + "ANITA4-Antenna-assignments.csv",
                             dtype = str, delimiter = ",")

    antennaSN = str(antennaInfo.T[3][np.argwhere(antennaInfo.T[2]==chan[:3])[0][0]])

    dir = localDir + "CSBF16/OFF-BS-PLSR/"
    if inOrOut == "in":
        fileName = dir + "CALIBRATION/06_22-INPUT_PULSE_20dB_DOWN-waveform.csv"
    if inOrOut == "out":
        fileName = dir + "SN" + antennaSN + "/V-TRANS/waveform/06_23-el_0-az_0-V-C-waveform.csv"
 
    return fileName


def importPalAntIn(chan):
    #Palestine Antennas
        # there is only one waveform :( 
        # so I can't do any averaging and I just import the raw data I guess
        # Also this means I probably need to link the whole dir (since they are far away)

    fileName = findPalestineAntennaFile(chan,"in")

    dataX,dataY = np.loadtxt(fileName, delimiter = ",", skiprows = 6, usecols=(3,4)).T
    dataX *= 1e9
    #make the time range start at zero
    dataX -= dataX[0]

    dataY = np.roll(dataY,30-np.argmax(dataY))
    #also I always want to have the fourier spaces of these things
    dataF, dataFFT = tf.genFFT(dataX,dataY)

    return dataX, dataY, dataF, dataFFT


def importPalAntOut(chan):

    fileName = findPalestineAntennaFile(chan, "out")

    dataX,dataY = np.loadtxt(fileName, delimiter = ",", skiprows = 6, usecols=(3,4)).T

    dataX *= 1e9 #s ->ns
    dataY *= -1 #also the polarity is flipped!

    #make the time range start at zero
    dataX -= dataX[0]

    dataY = np.roll(dataY,30-np.argmax(dataY))

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)

    return dataX, dataY, dataF, dataFFT


def importChamberAnts(fileName = localDir + roofDir + "HPulse_10degCCW_148p_090t_10dBatt.dat",numEvents=-1,channel=0):
    """
    I always need to write more parsers, always be parsing (ABP)
    """
    print "Importing locally -> " + fileName

    events = []
    event = ""
    file = open(fileName)
    for line in file:
        if line[0] == '#':
            if event != "":
                events.append(np.array(event))
            event = [[],[]]
            if ( numEvents != -1 and len(events) > numEvents):
                return events
        else:
            event[0].append(float(line.split()[0+channel*2]))
            event[1].append(float(line.split()[1+channel*2]))
    
    events.append(event)

    return events


def importRoofAntIn():
    #Hawaii Rooftop Antenna Input Pulse (FID)
    fileName = A4Dir + "avgFid.txt"
#    fileName = "avgFid.txt"
    dataX,dataY = np.loadtxt(fileName).T

    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)
        
    return dataX, dataY, dataF, dataFFT
                

def importRoofAntOut():
    #Hawaii Rooftop Antenna Input Pulse (Ant)
    fileName = A4Dir + "avgAnt.txt"
#    fileName = "avgAnt.txt"
    dataX,dataY = np.loadtxt(fileName).T

    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX, dataY)
        
    return dataX, dataY, dataF, dataFFT


def importSurf(chan):
    #LAB chip from the SURF
    fileName = waveformDir + str(chan) + "_avgSurfWaveform.txt"
    dataX,dataY = np.loadtxt(fileName).T
    
#    dataY /= 1000 #mv->V  #not mV!  It is in ADC counts so I shouldn't touch this...

    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)
        
    return dataX, dataY, dataF, dataFFT
    

def importScope(chan):
    #Calibration Scope
    fileName = waveformDir + str(chan) + "_avgScopeWaveform.txt"
    dataX,dataY = np.loadtxt(fileName).T

    dataX *= 1e9 #s ->ns

    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)
        
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

    cableNewFreq,cableNewFFT = tf.regenerateCable(cableFreq,cableFFT,dF=dF,length=length)
    
    print "getRegeneratedCables(): dFin=:",cableFreq[1]-cableFreq[0]
    print "getRegeneratedCables(): dFout=:",cableNewFreq[1]-cableNewFreq[0]
    print "getRegeneratedCables(): len(f)in=:",len(cableFreq)
    print "getRegeneratedCables(): len(f)out:",len(cableNewFreq)

    return cableNewFreq, cableNewFFT


def getCables(fileName):
    """
    basically just import the cables.
    
    also transforms it from measured log magnitude to linear magnitude

    """
    cableFreq,cableGainLog,cablePhase = s2pParser(cablesBaseDir + fileName)
    cableGainLin = 10.**(cableGainLog/20.)

    cableFFT = tf.gainAndPhaseToComplex(cableGainLin,cablePhase)

    return cableFreq,cableFFT


def getCablesOLD(fileName,tZero=False,hanning=False,resample=False):
    """
    I'm going to depreciate this because it doesn't really work: use getCables() instead

    A wrapper for getting the network analyzer data for the cables

    if tZero is specified, it will try to phase shift it until it has peak causality (most power after tZero)

    if hanning is specified, it will do a hanning window around the peak point
    
    Three different ways to resample it to the correct window and time binning:
    False = nothing
    if interp is specified, it will do a time series akima spline interpolation to get the correct
       sampling rate and window length for the rest of the waveforms (on by default)
    if fftZeroPad it will do an akima spline in the frequency domain and zero pad everything outside the range
    if fftFit it will fit a poly4 to the magnitude and attempt to regenerate the signal like that
          it also zero pads outside the range
    """

    cableFreq,cableGainLin,cablePhase = s2pParser(cablesBaseDir + fileName)

    #cables from network analyzer were stored in gain and phase...
    cableFFT = tf.gainAndPhaseToComplex(cableGainLin,cablePhase)
    
    fOut = cableFreq
    fftOut = cableFFT

    #For the long cables, there are two time domain pulses for some stupid reason
    #Maybe it wraps around?  Definitely an artifact at least... Windowing it out seems fair
    if hanning==True:
        print "Doing hanning window"
        hanningWidth = 100
        x,y = tf.genTimeSeries(fOut,fftOut)
        startLength = len(x)
        x,y = tf.hanningWindow(x,y,np.argmax(y),totalWidth=hanningWidth,slope=50)
        x,y = tf.zeroPadEqual(x,y,startLength)
        if tZero != False:
            tZero = len(x)/2
        fOut,fftOut = tf.genFFT(x,y)

    #make it causal so it is all pretty :) 
    if tZero != False:
        print "Making Causal with tZero=",tZero
        fftOut  = tf.makeCausalFFT(fftOut,tZero)

    """========================================================================
      Resampling
    ===========================================================================
    
      this might do bad things...
      the goal of it was to get the time binning and window length the same
      for an irfft, but this just worstens the fractional phase offset error at first...
      -> Yeah even for a "phase aligned" waveform it messes it up
    """
    
    if resample=="fftZeroPad":
        fOut,fftOut = tf.complexZeroPadAndResample(fOut,fftOut)

    #I just want 0.1ns sampling and 1024 samples in the time domain so I'll just do that
    #This ALSO introduces a bunch of weird frequency response
    elif resample=="interp":
        print "Interpolating..."
        x,y = tf.genTimeSeries(fOut,fftOut)
        peak = np.argmax(y)
        yInterp = Akima1DInterpolator(x,y) #has a weird frequency dependence
        newX = (np.arange(-512,512)*0.1)+x[peak]
        if newX[0]<0:
            newX -= newX[0]
        newY = np.nan_to_num(yInterp(newX))
        fOut,fftOut = tf.genFFT(newX,newY[::-1])

    #Last thing to try is fitting and reinterpolating it
    elif resample=="fftFit":
        fOut,fftOut = tf.fitAndRegenFFT(fOut,fftOut)

    return fOut,fftOut


def s2pParser(fileName):

    freq = []
    gainLog = []
    phase = []
    with open(fileName) as inFile:
        for line in inFile:
            if (line[0] !='!' and line[0] != '#'):
                split = line.split()
                freq.append(float(split[0])/1e9)
                phase.append(float(split[4])*(np.pi/180.))
                gainLog.append(split[3])
                               
    return np.array(freq),np.array(gainLog,dtype=float),np.array(phase)


#==============================================================================
# Processing.
#==============================================================================


def processWaveform(graphT,graphV,source):
    """
    Single functions to do lots of things!  Input what the "source" is and it
    will do the type of processing that you want :)

    This is also a good way to make a sort of "table" to deterimine
    what exactly is being done for each different type of source waveform
    """

    #I just care what type of thing this is, and source can have a filename too
    if (type(source) == list):
        source = source[0]

    # Notes on what type of "sources" I am looking for and 

    #Raw = what format and dimensions the data was originally taken in
    #Now = after being read and averaged and modified by the c++ code
    #
    #End = what everything should be after running through processing
    #      (10GS/s  1024 samples for all)
    
    #Processing alghorithms things that you can do:
    kPreZeroMean   = False
    kBandPassFilter= False
    kDownsample    = False
    kZeroPad       = False
    kHanningWindow = False
    kPostZeroMean  = False
    kSlice         = False
    kZeroStart     = True #I always want to make them start at zero :|
    
    #          SURF measured by LAB
    #Raw: record length of 259 sampling at 2.5GS/s  
    #Now: recoLen=1024, sampling 10GS/s
    if source.lower() == "surf":
        kPreZeroMean   = True
        kZeroPad       = True
#        kHanningWindow = True
        kPostZeroMean  = True
   
    #          SURF input pulse measured by scope
    #Raw: record length of 5000 | sampling rate of 10GS/s
    #Now: same
    if source.lower() == "calscope":
        kPreZeroMean   = True
#        kHanningWindow = True  #this cuts off a reflection pulse that is only like 35ns late
        kSlice         = True
        kPostZeroMean  = True
        kZeroPad       = True
   
    #          Rooftop antenna testing input (FID)
    #Raw: record length of 1000 | sampling rate of 40GS/s
    #Now: same
    elif source.lower() == "roofantin" or source.lower() == "fid":
        kPreZeroMean   = True
        kBandPassFilter= False #True adds power out of band in transfer function
        kDownsample    = True
        kZeroPad       = True
        kHanningWindow = True
        kPostZeroMean  = True

    #          Rooftop antenna testing output
    #Raw: record length of 5000 | sampling rate of 20GS/s
    #Now: same
    elif source.lower() == "roofantout" or source.lower() == "ant":
        kPreZeroMean   = True
        kBandPassFilter= True
        kDownsample    = True
        kZeroPad       = True
        kHanningWindow = True
        kPostZeroMean  = True

    #          Palestine antenna testing input
    #Raw: record length of 8000 (7994?) | sampling rate of 10GS/s
    #Now: same (only one waveform, so read right into python)
    elif source.lower() == "palin":
        kPreZeroMean   = True
        kZeroPad       = True
        kHanningWindow = False
        kPostZeroMean  = False #created problems with zero division at f=0

    #          Palestine antenna testing output
    #Raw: record length of 8000 (7994?) | samping rate of 10GS/s
    #Now: same (only one waveform, so read right into python)
    elif source.lower() == "palout":
        kPreZeroMean   = True
        kBandPassFilter= False
        kZeroPad       = True
        kHanningWindow = False
        kPostZeroMean  = False

    print "processWaveform(): Starting on " + source

    #print out what you got in initially
    print "Initial: length="+str(len(graphT))+" dT="+str(graphT[1]-graphT[0])+" T="+str((graphT[1]-graphT[0])*len(graphT))

    #pre-processing zero mean it
    if kPreZeroMean:
        initialMean = np.mean(graphV)
        graphV -= initialMean
        print "PreZeroMean: "+str(initialMean)+" -> "+str(np.mean(graphV))

    #find frequency binning and do butterworth IIR filter
    if kBandPassFilter:
        tTotal = (graphT[1]-graphT[0])*len(graphT)
        dF = 1/tTotal
        fMax = (len(graphT)/2)*dF
        graphV = tf.nyquistLimit(graphV,fMax)
        print "Filter: length="+str(len(graphT))+" dT="+str(graphT[1]-graphT[0])+" T="+str((graphT[1]-graphT[0])*len(graphT))
        
    #downsample to 10GS/s [start:stop:step]
    #average together avgN points to get one point
    if kDownsample:
        dT = graphT[1]-graphT[0]
        avgN = np.int(np.round(0.1/dT)) #0.1ns is goal (needs to be int)
        avgY = np.zeros(len(graphV)/avgN)
        for i in range(0,avgN):
            avgY += graphV[i::avgN]/avgN
        avgX = graphT[::avgN]+((avgN-1)*dT/2.)
        graphV = avgY
        graphT = avgX
        print "Downsample: length="+str(len(graphT))+" dT="+str(graphT[1]-graphT[0])+" T="+str((graphT[1]-graphT[0])*len(graphT))

    #find max and do hanning window (cuts waveform to desired length too)
    if kHanningWindow:
        #I have to zero pad out so I can select the center of the window
        graphT,graphV = tf.zeroPadEqual(graphT,graphV,1024*3)
        graphMaxV = np.argmax(graphV)
        graphT,graphV = tf.hanningWindow(graphT,graphV,graphMaxV,512,slope=100)
        print "Hanning: length="+str(len(graphT))+" dT="+str(graphT[1]-graphT[0])+" T="+str((graphT[1]-graphT[0])*len(graphT))

    #find max and do hanning window with SUPER SHARP EDGES(cuts waveform to desired length too)
    if kSlice:
        #I have to zero pad out so I can select the center of the window
        graphT,graphV = tf.zeroPadEqual(graphT, graphV, 1024 * 3)
        graphMaxV = np.argmax(graphV)
        graphT,graphV = tf.hanningWindow(graphT, graphV, graphMaxV + 240, 512, slope = 1)
        print "Slice: length="+str(len(graphT))+" dT="+str(graphT[1]-graphT[0])+" T="+str((graphT[1]-graphT[0])*len(graphT))

    #zero pad back to 1024 because I'm making my hanning window smaller (~50ns wide)
    if kZeroPad:
        graphT,graphV = tf.zeroPadEqual(graphT,graphV,1024)
    
    #zero mean it again at the end
    if kPostZeroMean:
        initialMean = np.mean(graphV)
        graphV -= initialMean
        print "PostZeroMean="+str(initialMean)+" -> "+str(np.mean(graphV))

    #make it start at zero (since all the other crap probably made it move)
    if kZeroStart:
        graphT -= graphT[0]

    return graphT,graphV


def doSigChainWithCables(chan,savePlots=False,showPlots=False,writeFiles=False):
    if savePlots or showPlots:
        lab.close("all")

    #get scope data from waveforms (extracted from raw data from another script)
    if chan=="15BH":
        chan = "15BV"
    scopeRawX,scopeRawY,scopeRawF,scopeRawFFT = importScope(chan)
    scopeX,scopeY = processWaveform(scopeRawX,scopeRawY,"calScope")
    scopeF,scopeFFT = tf.genFFT(scopeX,scopeY)
            
    #1.26 is the max coherance phase?
#    scopeFFT = tf.fftPhaseShift(scopeFFT,1.26)

    #Get the cable's H(f) transfer function for PULSER TO SCOPE
#    global cableScopeF #once upon a time regenerating the cables took awhile
#    global cableScopeFFT
#    if type(cableScopeF) != np.ndarray:
    print "Getting Pulser to Scope Cables..."
#        cableScopeF,cableScopeFFT= getCables\\\\\\\("A-B_PULSER-SCOPE.s2p",tZero=17,resample="interp")
    cableScopeRawF,cableScopeRawG,cableScopeRawP = s2pParser(cablesBaseDir+"A-B_PULSER-SCOPE.s2p")
    cableScopeRawGD = tf.calcGroupDelayFromPhase(np.unwrap(cableScopeRawP),dF=np.diff(cableScopeRawF)[0])
    cableScopeF,cableScopeFFT= getRegeneratedCables("A-B_PULSER-SCOPE.s2p")
    cableScopeGD = tf.calcGroupDelay(cableScopeFFT,inputF=cableScopeF)
    cableScopeX,cableScopeY = tf.genTimeSeries(cableScopeF,cableScopeFFT)

    #Get the cable's H(f) transfer function for PULSER TO AMPA
#    global P2AF
#    global P2AFFT
#    if type(P2AF) != np.ndarray:
    print "Getting Pulser to Ampa Cables..."
#        P2AF,P2AFFT= getCables("A-C_PULSER-TEST_66DB.s2p",tZero=True,hanning=True,resample="fftFit")
    cableAmpaFile = "A-C_PULSER-TEST_56DB.s2p"
    cableAmpaRawF,cableAmpaRawG,cableAmpaRawP = s2pParser(cablesBaseDir+cableAmpaFile)
    cableAmpaRawGD = tf.calcGroupDelayFromPhase(np.unwrap(cableAmpaRawP),dF=np.diff(cableAmpaRawF)[0])
    cableAmpaF,cableAmpaFFT= getRegeneratedCables(cableAmpaFile)
    cableAmpaGD = tf.calcGroupDelay(cableAmpaFFT,inputF=cableAmpaF)
    cableAmpaX,cableAmpaY = tf.genTimeSeries(cableAmpaF,cableAmpaFFT)

    #PULSER TO SCOPE
    #deconvolve cable pulse * cable = scope -> pulse = scope/cable
    #finds just the pulser impulse
    pulseDeconvFFT = scopeFFT/cableScopeFFT

    #PULSER TO AMPA
    #convolve it with that transfer function to get pulse at AMPA
    ampaInputFFT = cableAmpaFFT*pulseDeconvFFT
    ampaInputFFT = np.nan_to_num(ampaInputFFT)
    ampaInputF = scopeF
    #do a time domain thing
    ampaInputX,ampaInputY = tf.genTimeSeries(ampaInputF,ampaInputFFT)
    ampaInputY = np.roll(ampaInputY,100-np.argmax(ampaInputY))

    #then back to fft again
    ampaInputF,ampaInputFFT2 = tf.genFFT(ampaInputX,ampaInputY)
    #this is a stupid thing with complex division in python and inf/nan values
    ampaInputFFT2[np.absolute(ampaInputFFT2) < 1e-18] = 0

    ampaInputFFT = ampaInputFFT2    

    if (showPlots or savePlots) and (chan == "01BH"):
        
        fig2,ax2 = lab.subplots(3,2,figsize=(11,8.5))
        ax2[0][0].set_title("Cables")


        #phase
        ax2[0][0].set_ylabel("phase (radians)")
        
        ax2[0][0].plot(cableScopeRawF,cableScopeRawGD,label="Pulser To Scope Raw",color="green")
        ax2[0][0].plot(cableScopeF,cableScopeGD,label="Pulser To Scope Processed",color="blue")
        ax2[0][0].legend()

        ax2[0][1].set_ylabel("Group Delay (ns)")
        ax2[0][1].plot(cableAmpaRawF,cableAmpaRawGD,label="Pulser To Ampa Raw",color="green")
        ax2[0][1].plot(cableAmpaF,cableAmpaGD,label="Pulser To Ampa Processed",color="red")
        ax2[0][1].legend()
        
        #gain
        ax2[1][0].set_ylabel("gain (dB)")
        ax2[1][0].set_xlabel("freq (GHZ)")
        ax2[1][0].plot(cableScopeRawF,cableScopeRawG,label="Pulser To Scope Raw",color="green")
        ax2[1][0].plot(cableScopeF,tf.calcLogMag(cableScopeF,cableScopeFFT),label="Pulser To Scope Processed",color="blue")
        ax2[1][0].set_ylim([-100,0])
        ax2[1][0].legend()

        ax2[1][1].set_ylabel("gain (dB)")
        ax2[1][1].set_xlabel("freq (GHZ)")
        ax2[1][1].plot(cableAmpaRawF,cableAmpaRawG,label="Pulser To Ampa Raw",color="green")
        ax2[1][1].plot(cableAmpaF,tf.calcLogMag(cableAmpaF,cableAmpaFFT),label="Pulser To Ampa",color="red")
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
        ax[0][0].set_ylabel("Voltage(V)")
        ax[0][0].set_xlabel("Time (ns)")
        ax[0][0].set_xlim([0,200])
        ax[0][0].legend()

        ax[0][1].plot(ampaInputX,ampaInputY,label="cable corrected pulse into ampa",color="purple")
        ax[0][1].set_ylabel("Voltage (mV)")
        ax[0][1].set_xlabel("Time (ns)")
        ax[0][1].legend()

        ax[1][0].plot(scopeRawF,tf.calcSpecMag(scopeRawF,scopeRawFFT),label="raw scope pulse",color="red")
        ax[1][0].plot(scopeF,tf.calcSpecMag(scopeF,scopeFFT),label="processed scope pulse",color="blue")
        ax[1][0].set_xlabel("Frequency (GHz)")
        ax[1][0].set_ylabel("Spectral Power (dBm/Hz)")
        ax[1][0].set_ylim([-100,0])
        ax[1][0].legend()

        ax[1][1].plot(ampaInputF,tf.calcSpecMag(ampaInputF,ampaInputFFT),label="cable corrected pulse into ampa",color="purple")
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
    print " /|\ /|\ /|\ /|\ ",len(surfRawX),np.diff(surfRawX)[0]
    surfX,surfY = processWaveform(surfRawX,surfRawY,"surf")
    print " /|\ /|\ /|\ /|\ ",len(surfX),np.diff(surfX)[0]

    surfY *= (1./1000) #lets go to V for a sec
    surfF,surfFFT = tf.genFFT(surfX,surfY)

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
        axS[1].plot(surfRawF,tf.calcSpecMag(surfRawF,surfRawFFT),label="Raw Surf Waveform",color="red")
        axS[1].plot(surfF,tf.calcSpecMag(surfF,surfFFT),'+',label="Processed Surf Waveform",color="blue")
        axS[1].set_xlabel("Frequency (GHz)")
        axS[1].set_ylabel("\"Spectral Power\" (~dBm/Hz)")
        axS[1].set_ylim([-110,50])
        axS[1].set_xlim([0,3])
        axS[1].legend()
        if savePlots:
            figS.savefig("plots/doSigChainWithCables_SURF"+chan+".png")
        if showPlots:
            figS.show()

    print "*************",len(surfF),":",np.diff(surfF)[0],len(ampaInputF),":",np.diff(ampaInputF)[0]
    print "*************",len(surfFFT),len(ampaInputFFT)

    #deconvolve signal chain transfer function out of that!
    tfFFT = surfFFT/ampaInputFFT
    
    #fix the infinities and nans
    tfFFT[np.isposinf(np.real(tfFFT))] = 0
    tfFFT[np.isneginf(np.real(tfFFT))] = 0
    tfFFT = np.nan_to_num(tfFFT)
    tfFFT[tfFFT > 1e10] = 0
    tfFFT[tfFFT < -1e10] = 0
    
    #change it back to time domain
    tfX,tfY = tf.genTimeSeries(surfF,tfFFT)

    #clean up the tail and make it start at the beginning
    tfY = np.roll(tfY,30-np.argmax(tfY))
    tfY = tf.hanningTail(tfY,320,100)
    tfY = tf.nyquistLimit(tfY,5.)

    tfF,tfFFT = tf.genFFT(tfX,tfY)

    if showPlots or savePlots:
        fig3,ax3 = lab.subplots(2,2,figsize=(11,8.5))
        fig3.set_tight_layout(True)

#        fig3.suptitle("Input, Output, and Transfer Function Signals")
        ax30 = ax3[0][0].twinx()
        ax30.plot(ampaInputX,1000*ampaInputY,label="Input Pulse, Cables Convolved",color="red")
        ax30.set_ylabel("Voltage (mV)")
        ax30.set_ylim([-0.2,1])
        ax30.legend()
        surfYMax = np.argmax(surfY)
        surfYMax = 200
        ax3[0][0].plot(surfX,np.roll(surfY,200-surfYMax),label="SURF waveform",color="blue")
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
        ax3[1][0].plot(ampaInputF,tf.calcSpecMag(ampaInputF,ampaInputFFT),label="Input Pulse, Cables Removed",color="red")
        ax3[1][0].plot(surfRawF,tf.calcSpecMag(surfRawF,surfRawFFT),label="raw SURF waveform",color="blue")
        ax3[1][0].legend()
        ax3[1][0].set_ylim([-150,20])

        ax3[1][1].set_ylabel("Gain (dB)")
        ax3[1][1].set_xlabel("Frequency (GHz)")
        ax3[1][1].plot(surfF,tf.calcLogMag(surfF,tfFFT),label="Signal Chain Transfer Function",color="black")
        ax3[1][1].set_xlim([0,3])
        ax3[1][1].legend()
        ax3[1][1].set_ylim([0,100])

        if showPlots:
            fig3.show()
        if savePlots:
            fig3.savefig("plots/doSigChainWithCables_TF"+chan+".png")

    if writeFiles:
        writeOne([tfX,tfY],chan,"sigChain")

    return tfX,tfY,surfF,tfFFT 


def doPalAnt(chan,savePlots=False,showPlots=False,writeFiles=False):
    """
    
    Develop a transfer function for an antenna measured in palestine 2014, mapped to chan
    
    Also known as the Antenna Height.  See elog 15 (wow its been awhile)

    Normalized to V/m (saved in V - ns fomat)

    """

    if showPlots or savePlots:
        lab.close("all")

    #input pulse
    inRawX,inRawY,inRawF,inRawFFT = importPalAntIn(chan)

    inX,inY = processWaveform(inRawX,inRawY,"palin")
    inY = np.roll(inY,30-np.argmax(inY))

    inF,inFFT = tf.genFFT(inX,inY)
    # Increase by 20dB (from coupler)
    inFFT *= 10
#    #increase by 10dB (from coupler)
##    inFFT *= 10**(10./20)
    inX,inY = tf.genTimeSeries(inF,inFFT)

#    inFFT = tf.minimizeGroupDelayFromFFT(inF,inFFT)


    #output pulse
    outRawX,outRawY,outRawF,outRawFFT = importPalAntOut(chan)
    outX,outY = processWaveform(outRawX,outRawY,"palout")
    outY = np.roll(outY,50-np.argmax(outY))
#    outY = tf.hanningTail(outY,130,30)
##    outX = outX[:1024]
##    outY = outY[:1024]
#    outY = tf.highPass(outX,outY)
    outF,outFFT = tf.genFFT(outX,outY)

#    outFFT = tf.minimizeGroupDelayFromFFT(outF,outFFT)

    if savePlots or showPlots:
        fig0,ax0 = lab.subplots(3,2)

        ax0[0][0].plot(inRawX,inRawY,label="Raw Input Pulse",color="red")
        ax0[0][0].plot(inX,inY,label="Processed Pulse",color="blue")
        ax0[0][0].set_xlim([0,20])
        ax0[0][0].set_xlabel("Time (ns)")
        ax0[0][0].set_ylabel("Voltage (V)")
        ax0[0][0].legend(loc="lower right")

        ax0[0][1].plot(outRawX,outRawY,label="Raw Output Pulse",color="red")
        ax0[0][1].plot(outX,outY,label="Processed Pulse",color="blue")
        ax0[0][1].set_xlim([0,40])
        ax0[0][1].set_xlabel("Time (ns)")
        ax0[0][1].set_ylabel("Voltage (V)")
        ax0[0][1].legend()

        ax0[1][0].plot(inRawF,tf.calcSpecMag(inRawF,inRawFFT),label="Raw Input Pulse",color="red")
        ax0[1][0].plot(inF,tf.calcSpecMag(inF,inFFT),label="Processed Pulse",color="blue")
##        ax0[1][0].set_ylim([-70,-20])
        ax0[1][0].set_ylabel("Spectral Magnitude \n [dBm/Hz]")
        ax0[1][0].legend()

        ax0[2][0].plot(inRawF,tf.calcGroupDelay(inRawFFT,inputF=inRawF),label="Raw Input Pulse",color="red")
        ax0[2][0].plot(inF,tf.calcGroupDelay(inFFT,inputF=inF),label="Processed Pulse",color="blue")
##        ax0[2][0].set_ylim([-25,25])
        ax0[2][0].set_ylabel("Group Delay (ns)")
        ax0[2][0].legend()

        ax0[1][1].plot(outRawF,tf.calcSpecMag(outRawF,outRawFFT),label="Raw Output Pulse",color="red")
        ax0[1][1].plot(outF,tf.calcSpecMag(outF,outFFT),label="Processed Pulse",color="blue")
##        ax0[1][1].set_ylim([-100,-40])
        ax0[1][1].set_ylabel("Gain (dB)")
        ax0[1][1].set_xlabel("Frequency (GHz)")
        ax0[1][1].legend()

        ax0[2][1].plot(outRawF,tf.calcGroupDelay(outRawFFT,inputF=outRawF),label="Raw Output Pulse",color="red")
        ax0[2][1].plot(outF,tf.calcGroupDelay(outFFT,inputF=outF),label="Processed Pulse",color="blue")
##        ax0[2][1].set_ylim([-30,30])
        ax0[2][1].set_ylabel("Group Delay (ns)")
        ax0[2][1].set_xlabel("Frequency (GHz)")
        ax0[2][1].legend()

        if showPlots:
            fig0.show()
            fig0.tight_layout()
            fig0.canvas.draw()
        if savePlots:
            fig0.savefig("plots/doPalAnt_inAndOut"+chan+".png")

    #computeTF?  hmm maybe just do a normal one
    """
    This does the (F_rec/F_src)*(c/if) calculation
    ends up being in units of "meters"
    """

    antTFFFT = np.divide(outFFT,inFFT)*(0.3/(1j*inF))
    antTFFFT[0] = 0

    antTFF = inF
#    antTFX,antTFY,antTFF,antTFFFT = computeTF(inF,inFFT,outF,outFFT)

    """
      one antenna should be the square root of the response
      cmath can't operate on numpy matricies so I have to iterate by hand
      this also doesn't work... maybe I can just do it in group delay and
      magnitude?
      basically since square roots are poorly defined in complex space, you
      NEED to do it with the sines and cosines, which I added to tfUtils!
      
      Make causal?  no probably not, that seems like a "last thing you do" sort
      of technique
    """
    
    #need to take out 1/R flight distance for absolute gain
    # According to notes, 8.89m face to face flight distance
    # Antennas have a depth of 22" = 0.5588m, so r(f) should have two points, (0.180,8.89) and (1.2,10)
    # makes it in units of "meters squared"
    distM = (0.5588*2)/(1.2-0.18)
    distYint = 8.89 - distM*0.18
    dist = antTFF*distM + distYint
    for i in range(0,len(dist)):
        if dist[i] > 10.5:
            dist[i] = 10.5
    
    antTFFFT *= dist

    print "length antTFFFT:",len(antTFFFT)

    if showPlots and 1:
        figDist,axDist = lab.subplots()
        axDist.plot(antTFF,dist)
        axDist.set_xlabel("Frequency (GHz)")
        axDist.set_ylabel("Distance (m)")
        figDist.show()
    
    #take the square root to get the normalized complex height
    # units of "meters again!"
    antTFFFT = tf.sqrtOfFFT2(antTFF,antTFFFT)

    #causality conditions aren't good for measured signals (artifically distort phase)!
#    antTFFFT = tf.makeCausalFFT(antTFFFT,np.argmax(tf.fftw.irfft(antTFFFT)))

    #I have to get the time domain again after that
    # units of "meters / nanosecond"
    antTFX,antTFY = tf.genTimeSeries(antTFF,antTFFFT)

    #remove all the stupid power above 1.3GHz since ANITA can't measure it (probably should be lower)
    #butterworth filter @ 1.3
    antTFY = tf.highPass(antTFX,antTFY)
    antTFY = tf.nyquistLimit(antTFY,5)

    #10BH is upside down?
    #A bunch are actually!  Nothing here constrains absolute polarity
    if np.max(antTFY) < -np.min(antTFY):
        print "Flipped!"
        antTFY *= -1
    
    antTFF,antTFFFT = tf.genFFT(antTFX,antTFY)

    #Produce a plot for the gain normalized to dBi (hence the 4pi)
    if showPlots or savePlots:
        gainTF = ((4*np.pi*antTFF**2)/(0.3**2))*np.abs(antTFFFT)**2
        figGain,axGain = lab.subplots()
        figGain.suptitle(chan)
        axGain.plot(antTFF,10*np.log10(gainTF)) #10 instead of 20 because... its a power gain already I htink

        axGain.set_ylabel("Antenna Gain (dBi)")
        axGain.set_xlabel("Frequency (GHz)")
##        axGain.set_ylim([-20,20])
        axGain.set_xlim([0,3])

        if showPlots:
            figGain.show()
        if savePlots:
            figGain.savefig("plots/doPalAnt_AntGain"+chan+".png")
        
    #clean up the tail and make it start at the beginning
#    antTFY = np.roll(antTFY,20-np.argmax(antTFY))
#    antTFY = tf.hanningTail(antTFY,20,20)
#    antTFFFT = tf.fftw.rfft(antTFY)

    if savePlots or showPlots:
        fig,ax = lab.subplots(3,3)
        
        ax[0][0].plot(inX,inY,label="input Pulse",color="red")
        ax[0][1].plot(outX,outY,label="output Pulse",color="blue")
        ax[0][2].plot(antTFX,antTFY,label="transfer Function",color="black")
        ax[0][0].set_xlim([0,25])
        ax[0][1].set_xlim([0,25])
        ax[0][2].set_xlim([0,25])
        ax[0][0].legend()
        ax[0][1].legend()
        ax[0][2].legend()

        ax[1][0].plot(inF,tf.calcSpecMag(inF,inFFT),label="input Pulse",color="red")
        ax[1][1].plot(outF,tf.calcSpecMag(outF,outFFT),label="output Pulse",color="blue")
        ax[1][2].plot(antTFF,tf.calcLogMag(antTFF,antTFFFT),label="transfer Function",color="black")
        ax[1][0].legend()
        ax[1][1].legend()
        ax[1][2].legend()
       
        ax[2][0].plot(inF,tf.calcGroupDelay(inFFT,inputF=inF),label="input Pulse",color="red")
        ax[2][1].plot(outF,tf.calcGroupDelay(outFFT,inputF=inF),label="output Pulse",color="blue")
        ax[2][2].plot(antTFF,tf.calcGroupDelay(antTFFFT,inputF=antTFF),label="transfer Function",color="black")
        ax[2][0].legend()
        ax[2][1].legend()
        ax[2][2].legend()

        ax[0][0].set_xlabel("Time (ns)")
        ax[0][1].set_xlabel("Time (ns)")
        ax[0][2].set_xlabel("Time (ns)")

        ax[0][0].set_ylabel("Voltage (~V)")
        ax[1][0].set_ylabel("Gain (dB)")
        ax[2][0].set_ylabel("Group Delay (ns)")
        
        ax[2][0].set_xlabel("Frequency (GHz)")
        ax[2][1].set_xlabel("Frequency (GHz)")
        ax[2][2].set_xlabel("Frequency (GHz)")
        
        ax[0][0].set_title("input Pulse")
        ax[0][1].set_title("output Pulse")
        ax[0][2].set_title("transfer Function")
        
        if showPlots:
            fig.show()
            fig.tight_layout()
            fig.canvas.draw()
        if savePlots:
            fig.savefig("plots/doPalAnt_TF"+chan+".png")
        
        if writeFiles:
            writeOne([antTFX,antTFY],chan,"palAnt")
        
        return antTFX,antTFY,antTFF,antTFFFT


def doSigChainAndAntenna(chan,showPlots=False,savePlots=False,writeFiles=False):
    #get antenna
#    antX,antY,antF,antFFT = doRoofAntWithCables()
    antX,antY,antF,antFFT = doPalAnt(chan,showPlots=showPlots,savePlots=savePlots,writeFiles=writeFiles)
    
    #get sig chain
    print "doSigChainAndAntenna:getSigChain"
    sigChainX,sigChainY,sigChainF,sigChainFFT = doSigChainWithCables(chan,showPlots=showPlots,savePlots=savePlots,writeFiles=writeFiles)

    #convolve the two (full ANITA3 transfer function!)
    a3F = sigChainF
    a3FFT = sigChainFFT * antFFT

    #also there is a scaling from free space (377ohms) to system impedance (50ohms) [ sqrt(50/377) ]
    a3FFT *= np.sqrt(50./377)


    a3X,a3Y = tf.genTimeSeries(a3F,a3FFT)
    
    if (showPlots or savePlots):
#        lab.close("all")
        fig,ax = lab.subplots(5,figsize=(11,8.5))
        ax[0].set_title(chan)
        ax[0].plot(antX,np.roll(antY,100),label="Antenna")
        ax[0].legend()
        ax[0].set_ylabel("Antenna Height (V/m)")
        ax[1].plot(sigChainX,np.roll(sigChainY,100),label="sigChain",color="red")
        ax[1].legend()
#        ax[1].set_ylim([-0.2,0.2])
        ax[1].set_ylabel("Gain (unitless)")
        a3YPeak = np.argmax(a3Y)
        ax[2].plot(a3X,np.roll(a3Y,100-a3YPeak),label="Convolution",color="black")
        ax[2].legend()
#        ax[2].set_ylim([-0.6,0.6])
        ax[2].set_xlabel("Time (ns)")
        ax[2].set_ylabel("Instrument Response \n (V/m)")
        
        ax[3].plot(a3F,tf.calcLogMag(a3F,a3FFT),color="black")
#        ax[3].set_ylim([-100,-30])
        ax[3].set_xlim([0.1,1.5])
        ax[3].set_xlabel("Frequency (GHz)")
        ax[3].set_ylabel("Gain (dBm)")

        ax[4].plot(a3F,tf.calcPhase(tf.minimizeGroupDelayFromFFT(a3F,a3FFT,highLim=256)),color="black")
        ax[4].set_xlim([0,2])
#        ax[4].set_ylim([-20,20])

        if showPlots:
            fig.show()
        if savePlots:
            fig.savefig("plots/doSigChainAndAntenna_TF"+chan+".png")

    #make it look nice and cut off the junk at the end
    a3Y = np.roll(a3Y,40-np.argmax(a3Y))
    a3Y = tf.hanningTail(a3Y,300,200)
#    a3X,a3Y = tf.zeroPadEqual(a3X,a3Y,1024)

    if writeFiles:
        writeOne([a3X,a3Y],chan,"fullTF")
    
    return a3X,a3Y

"""
  Evaluating the transfer function chain for a TUFF. Assumes input frequency
  array (f) is in Hertz. By setting each notch's variable capacitors
  (varCap1, varCap2, varCap3) to being negative, the default behavior is
  for the notches to be off.
"""
def doTUFF(f, varCap1 = -1, varCap2 = -1, varCap3 = -1):
    R = 6  #  Parasitic notch resistance (Ohms).
    C = 6e-13  #  Parasitic notch capacitance (Farads).
    L = 56e-9  #  Notch inductance (Henries).
    
    # Calculate total impedance from parallel components.
    Yparallel = 1 / 50  #  Y = 1 / Z called admittance. Start of total notch admittance.
    def Ynotch(capVal):  #  Input admittance function for notch filters.
        Znotch = R + 1j * 2 * np.pi * f * L + (1j * 2 * np.pi * f * (C + capVal))**-1
        return Znotch**-1
    """
      Starting here, we apply additional impedances when notches switched on.
      Chose negative input values to correspond to notches being switched off.
    """
    if (varCap1 >= 0): Yparallel += Ynotch(1.8e-12 + varCap1)
    if (varCap2 >= 0): Yparallel += Ynotch(12e-12 * varCap2 / (12e-12 + varCap2))
    if (varCap3 >= 0): Yparallel += Ynotch(1.5e-12 * varCap3 / (1.5e-12 + varCap3))
    
    #  Calculate complex TUFF gain.
    GTUFF = (1 + 50 * Yparallel)**-1
    
    #  Find transfer function in dB and its phase.
    GdBTUFF = 40 + 20 * np.log10(np.abs(GTUFF))  #  Accounting for 40 dB amplification.
    phaseTUFF = np.angle(GTUFF, deg = True)  #  Phase in degrees.
    return GdBTUFF, phaseTUFF

def doTheWholeShebang(savePlots=False,showPlots=False,writeFiles=False):
    # Generate the transfer function for all the channels!
    chans = np.loadtxt("chanList.txt",dtype=str)
    allChans = {}
#    lab.close("all")
    for chan in chans:
        try:
            allChans[chan] = doSigChainAndAntenna(chan,savePlots=savePlots,showPlots=showPlots,writeFiles=writeFiles)
        except:
            print chan+" FAILED"

    return allChans


#==============================================================================
# Plotting and writing txt files
#==============================================================================


def writeOne(wave,chan,prefix):
    """
      Writes a single file to disk in a space seperated variable format so that
      I can post it
    """
    file = open("transferFunctions/"+prefix+"_"+chan+".txt","w")
    for i in range(0,len(wave[0])):
        file.write(str(wave[0][i])+" "+str(wave[1][i])+"\n")
    file.close()
    print "wrote "+chan+" to file"

    return
