"""
Ben Rotter
University of Hawaii at Manoa
ANITA3
February 2016 - through September 2016

So this is a ton of my code all glommed together into one giant library.  This
should make it easy to generate the full transfer function from all the
different waveforms, and for all the channels easily

It uses tfUtils (which I might turn into a auto-linked brotterlib at some
point) which has a lot of FFT analysis and fun tools to mess with

It uses data from avgSigChainWaveform and avgAntWaveform, which imports,
correlates and averages the waveforms using eventReaderRoot and libFftwWrapper,
things that I don't want to recreate.

This is a selection of a much larger code base that does all the work, which is
available upon request if you want to see a ton of debugging and testing
functions, in addition to some classes that make plotting easier.

Edited by John Russell.
"""

import numpy as np
import pylab as lab

from scipy import signal
from scipy.interpolate import interp1d, Akima1DInterpolator

from glob import glob

import tfUtils as tf

import cmath #for complex sqrt

import time #for drawing with a pause

try:
    import seaborn as sns #just a pretty plotting package, it isn't important
except:
    print "You don't have seaborn but thats okay!"


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


#==============================================================================
# Importing things
#==============================================================================

def importANITA1():
    a1 = np.loadtxt(A3Dir + "ANITA1ImpulseResponse.dat")
#    a1 = np.loadtxt("ANITA1ImpulseResponse.dat")
    return a1.T[0]*1e9, a1.T[1]


def findPalestineAntennaFile(chan, inOrOut):
    antennaInfo = np.loadtxt(A3Dir + "antennaInformation_cut.csv",
                             dtype = str, delimiter = ",")    
#    antennaInfo = np.loadtxt("antennaInformation_cut.csv",dtype=str,delimiter=",")

    antennaNumber = int(antennaInfo.T[1][np.argwhere(antennaInfo.T[2]==chan[:3])[0][0]][1:])

    print chan + " maps to antenna number " + str(antennaNumber)
        
    if antennaNumber == 52:
        print "Antenna 52 doesn't have any calibration associated with it :(, using 51 instead"
        antennaNumber = 51

    dir = localDir + "palestine14/SeaveyAntennas/S21s/"
#    dir = "/Volumes/ANITA3Data/palestine14/SeaveyAntennas/S21s/"
        
    fileName = dir + chan[-1].lower() + "pol_ezLinks/rxp" + str(antennaNumber).zfill(2)
    if inOrOut == "in":
        fileName += "_inputPulse.csv"
    if inOrOut == "out":
        fileName += "_coPol.csv"
 
    return fileName


def importPalAntIn(chan):
    #Palestine Antennas
        # there is only one waveform :( 
        # so I can't do any averaging and I just import the raw data I guess
        # Also this means I probably need to link the whole dir (since they are far away)

    fileName = findPalestineAntennaFile(chan,"in")

    dataX,dataY = np.loadtxt(fileName,comments="\"",delimiter=",",usecols=(3,4)).T

    dataX *= 1e9
    #make the time range start at zero
    dataX -= dataX[0]

    dataY = np.roll(dataY,30-np.argmax(dataY))

    #also I always want to have the fourier spaces of these things
    dataF, dataFFT = tf.genFFT(dataX,dataY)

    return dataX, dataY, dataF, dataFFT


def importPalAntOut(chan):

    fileName = findPalestineAntennaFile(chan, "out")

    dataX,dataY = np.loadtxt(fileName,comments="\"",delimiter=",",usecols=(3,4)).T

    dataX *= 1e9 #s ->ns
#    dataY *= 100 #there is something super wrong with the scaling, 20mV Vpp, 15uV resolution?  Weird...
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
    fileName = A3Dir + "avgFid.txt"
#    fileName = "avgFid.txt"
    dataX,dataY = np.loadtxt(fileName).T

    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)
        
    return dataX, dataY, dataF, dataFFT
                

def importRoofAntOut():
    #Hawaii Rooftop Antenna Input Pulse (Ant)
    fileName = A3Dir + "avgAnt.txt"
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
#        yInterp = interp1d(x,y,kind="quadratic") #is weird in a different way
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
#                gainLin.append(10.**(float(split[3])/10.))
                gainLog.append(split[3])
               
                
    return np.array(freq),np.array(gainLog,dtype=float),np.array(phase)


#==============================================================================
# Processing
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


"""============================================================================
  Generating impulse responses

  getting the cables is really resource intensive, so I'd like to do it once
  and then be done with it
  so I guess I can save it within the module?  Initialize them as False and
  then fill them if they aren't an numpy ndarray (which a boolean isn't)
============================================================================"""

P2SF = False
P2SFFT = False
P2AF = False
P2AFFT = False

def doSigChainWithCables(chan,savePlots=False,showPlots=False,writeFiles=False):
    if savePlots or showPlots:
        lab.close("all")

    #phase shifts are likely pointless!
    #they were: scopeFFT=1.26
    #           phaseShift=1.698

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
#    ampaInputY = tf.hanningTail(ampaInputY,625,100)
#    ampaInputY = tf.hanningTailNeg(ampaInputY,200,100)
#    ampaInputY -= np.mean(ampaInputY)

#    return ampaInputY,ampaInputY2

    #then back to fft again
    ampaInputF,ampaInputFFT2 = tf.genFFT(ampaInputX,ampaInputY)
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

#    return surfFFT,ampaInputFFT,tfFFT

    
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

    #zero out everything above 1.3GHz because that's our nyquist limit
#    for i in range(0,len(surfF)):
#        if surfF[i] > 1.3:
#            tfFFT[i] /= 1e9


    if writeFiles:
        writeOne([tfX,tfY],chan,"sigChain")


    return tfX,tfY,surfF,tfFFT 


def doRoofAntWithCables():

    #get waveforms and process them
    antInX,antInY,antInF,antInFFT = importRoofAntIn()
    antInX,antInY = processWaveform(antInX,antInY,"roofantin")
    antInF,antInFFT = tf.genFFT(antInX,antInY)

    antOutX,antOutY,antOutF,antOutFFT = importRoofAntOut()
    antOutX,antOutY = processWaveform(antOutX,antOutY,"roofantout")
    antOutF,antOutFFT = tf.genFFT(antOutX,antOutY)

    #compute the transfer function (for chan pair system!)
    antTFX,antTFY,antTFF,antTFFFT = computeTF(antInF,antInFFT,antOutF,antOutFFT)

    #factor out 47m flight free space path loss (inverse square law)
#    FSPL = ((4.*np.pi)/47.)**2
#    print "FSPL="+str(FSPL)
#    antTFFFT /= FSPL

    #one antenna should be the square root of the response
    #cmath can't operate on numpy matricies so I have to iterate by hand
#    for i in range(0,len(antTFFFT)):
#        antTFFFT[i] = cmath.sqrt(antTFFFT[i])
    
#    fig0,ax = lab.subplots(3)
#
#    ax[0].plot(antInX,antInY)
#    ax[0].plot(antOutX,antOutY)
#    ax[1].plot(antTFX,antTFY)
#    ax[2].plot(antTFF,antTFFFT)
#    fig0.show()

    return antTFX,antTFY,antTFF,antTFFFT


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
#    inY = tf.hanningTail(inY,130,20)
    inX = inX[:1024]
    inY = inY[:1024]
    inF,inFFT = tf.genFFT(inX,inY)
    #increase by 10dB (from coupler)
    inFFT *= 10**(10./20)
    inX,inY = tf.genTimeSeries(inF,inFFT)

#    inFFT = tf.minimizeGroupDelayFromFFT(inF,inFFT)


    #output pulse
    outRawX,outRawY,outRawF,outRawFFT = importPalAntOut(chan)
    outX,outY = processWaveform(outRawX,outRawY,"palout")
    outY = np.roll(outY,50-np.argmax(outY))
#    outY = tf.hanningTail(outY,130,30)
    outX = outX[:1024]
    outY = outY[:1024]
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
        ax0[1][0].set_ylim([-70,-20])
        ax0[1][0].set_ylabel("Spectral Magnitude \n [dBm/Hz]")
        ax0[1][0].legend()

        ax0[2][0].plot(inRawF,tf.calcGroupDelay(inRawFFT,inputF=inRawF),label="Raw Input Pulse",color="red")
        ax0[2][0].plot(inF,tf.calcGroupDelay(inFFT,inputF=inF),label="Processed Pulse",color="blue")
        ax0[2][0].set_ylim([-25,25])
        ax0[2][0].set_ylabel("Group Delay (ns)")
        ax0[2][0].legend()

        ax0[1][1].plot(outRawF,tf.calcSpecMag(outRawF,outRawFFT),label="Raw Output Pulse",color="red")
        ax0[1][1].plot(outF,tf.calcSpecMag(outF,outFFT),label="Processed Pulse",color="blue")
        ax0[1][1].set_ylim([-100,-40])
        ax0[1][1].set_ylabel("Gain (dB)")
        ax0[1][1].set_xlabel("Frequency (GHz)")
        ax0[1][1].legend()

        ax0[2][1].plot(outRawF,tf.calcGroupDelay(outRawFFT,inputF=outRawF),label="Raw Output Pulse",color="red")
        ax0[2][1].plot(outF,tf.calcGroupDelay(outFFT,inputF=outF),label="Processed Pulse",color="blue")
        ax0[2][1].set_ylim([-30,30])
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
        axGain.set_ylim([-20,20])
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

        ax[1][0].set_ylim([-60,-35])
        ax[1][1].set_ylim([-100,-60])
        ax[1][2].set_ylim([-60,0])

        ax[2][0].set_ylim([-10,10])
        ax[2][1].set_ylim([-20,20])
        ax[2][2].set_ylim([-20,20])
        
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


def palAntCables():
    """
    

    """




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

    #clean it up a bit
    #something was adding a ton of 1.25GHz noise
#    for i in range(0,len(a3FFT)):
#        if a3F[i] > 1.1:
#            a3FFT[i] /= 1e4

#    a3FFT = np.concatenate((a3FFT[:171],np.zeros(342))) #this isn't causal...

    #make it look nice and cut off the junk at the end
    a3Y = np.roll(a3Y,40-np.argmax(a3Y))
    a3Y = tf.hanningTail(a3Y,300,200)
#    a3X,a3Y = tf.zeroPadEqual(a3X,a3Y,1024)

    if writeFiles:
        writeOne([a3X,a3Y],chan,"fullTF")
        
    
    return a3X,a3Y


def computeTF(inF,inFFT,outF,outFFT):
    """
      Generating the transfer function is always the same, so it gets its own class
      Comparable to "deconvwnr()" in matlab (nearly identical)
    
      It has a bunch of other things you can uncomment to use, but most end up being dumb
    """
        #then generate the transfer function!
        #In(f)*H(f) = Out(f) -> H(f) = Out(f)/In(f)
    tfFFT = np.divide(outFFT,inFFT)*(0.3/(1j*inF))

    tfFFT[0] = 0
#    print tfFFT
#    tfFFT = np.divide(outFFT,inFFT)
    tfF = inF
    
    tfY = tf.fftw.irfft(tfFFT)        

#    print tfY

    tfX = np.arange(0,len(tfY))*(1./(2.*tfF[-1]))

    """
      We can also try "cleaning it up" if we want, these lettered things
      basically do that.
      
      A) the group delay offset of this fft is going to be weird, it is the
      addition of the phase offset for both of them (I think?)
      so I have to roll it back a bit or it will wrap around the window,
      probably by the mean of the group delay!
    """

#    tfGrpDlyMean = np.mean(tf.calcGroupDelay(tfFFT,inputF=tfF))
#    dT = tfX[1]-tfX[0]
#    print "tfGrpDlyMeanPt=" + str(tfGrpDlyMean)

        #roll it backwards, 1/8th off window (the mean is at the middle I think)
#    try:
#        tfY = np.roll(tfY,-(int(tfGrpDlyMean/dT)-(len(tfY)*3)/8)) 
#    except:
#        print "Warning in computeTF: coulnd't roll because the group delay is totally borked"

        #B) maybe high pass it (if there is a lot of power there)
    #    tfY = tf.highPass(tfX,tfY)    

        #Bb) also low pass it, what the hell
    #    tfY = tf.lowPass(tfX,tfY)

        #C) or we can take a hanning window
    #    center = np.argmax(tfY)
    #    print "tf center: "+str(center)
    #    tfX,tfY = tf.hanningWindow(tfX,tfY,center,1024,slope=10)

        #X) Thats it, so convert it back and make sure the length is correct
    tfF,tfFFT = tf.genFFT(tfX,tfY)
    
    return tfX,tfY,tfF,tfFFT


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


def alignWaveforms(allChans,showPlots=False):
    
    #They need to all be aligned and have the same group delay!
    #I try two ways:
    #1) use shiftToAlign, which does the weird fourier linear regression that messes with the phase
    #2) (THE DEFAULT) fit a gaussian to the peak of the correlation plot and use that as the fractional offset, then resample


    outChans = {}
    desiredDelay = 40
    max = np.argmax(allChans["01BH"][1])
    outChans["01BH"] = allChans["01BH"][0],np.roll(allChans["01BH"][1],max-desiredDelay)
    for chan in allChans:
        if chan == "01BH":
            continue
        outChans[chan] = allChans[chan][0],tf.phaseFitCorrelate(outChans["01BH"],allChans[chan])[0];


    outChans2 = {}
    for chan in allChans:
        print chan
        outChans2[chan] = tf.shiftToAlign(allChans["01BH"],allChans[chan])


    correlations = {}
    for chan in allChans:
        correlations[chan] = tf.correlation(outChans["01BH"][1],allChans[chan][1])



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
        print len(highWeird),highWeird
        print len(lowWeird),lowWeird

        fig.show()

    return outChans2




def doAllSigChains(savePlots=False,showPlots=False):
    # just does all the signal chains without the antennas
    chans = np.loadtxt("chanList.txt",dtype=str)
    allChans = {}
    for chan in chans:
        try:
            allChans[chan] = doSigChainWithCables(chan,savePlots=savePlots,showPlots=showPlots)[:2]
        except:
            print chan+" FAILED************************************************"

    return allChans


def doAllAntennas(savePlots=False,showPlots=False):
    # just does all the palestine antennas
    chans = np.loadtxt("chanList.txt",dtype=str)
    allChans = {}
    for chan in chans:
        try:
            allChans[chan] = doPalAnt(chan,showPlots=showPlots,savePlots=savePlots)[:2]
        except:
            print chan+" FAILED************************************************"

    return allChans


def plotAllAntennas(allAnts):
    fig,ax = lab.subplots(2)
    for chan in allAnts:
        if chan[-1] == "H":
            ax[0].plot(allAnts[chan][0],allAnts[chan][1])
        else:
            ax[1].plot(allAnts[chan][0],allAnts[chan][1])

    fig.show()


def findSignalToNoise():
    """
      I need to find the signal to noise ratio of the signal to do a Weiner
      Deconvolution
      Basically, out of band our deconvolutions introduce noise
      (small/small = big possibly)
      so we need to take the SNR into account.
      
      The signal in our transfer function is actually the signal in the input
      pulse! It is super well averaged (so there is basically no noise)
      
      The noise in our transfer function then will be whatever the noise
      spectrum of the terminated SURF input is?  Basically doing an average of
      all the FFTs of a bunch of waveforms with no signal...
    """
    return


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

def writeAll(allChans,prefix):
    """
      Writes all the files to disk in a space seperated variable format so that
      I can post them
    """
    for chan in allChans.keys():
        try:
            writeOne(allChans[chan],chan,prefix)
        except:
            print chan+" FAILED"


def saveAllNicePlots(allChans):
    """
      Saves a bunch of plots that I made so I can find channels that didn't
      work correctly
    """

    lab.rcParams['figure.figsize'] = [16.0,11.0]
    fig,ax = lab.subplots(3)

    #get the ANITA1 transfer function
    a1X,a1Y = importANITA1()
    a1Y *= 1000 #V to mV (and a nudge upwards for some reason) (also is it phase flipped(or maybe I am)?)
    
    a1F,a1FFT = tf.genFFT(a1X,a1Y)

    for chan in allChans:
        print chan

        a3X = allChans[chan][0]
        a3Y = allChans[chan][1]

        a3F,a3FFT = tf.genFFT(a3X,a3Y)

        a1YInterp = tf.interp(a3X,a1X,a1Y)
    
        #scale them
        a1MaxVal = np.max(a1Y)
        a3MaxVal = np.max(a3Y)
        scale = a1MaxVal/a3MaxVal
        
        #make them start at the beginning
        argMax = np.argmax(a3Y)
        a3Y = np.roll(a3Y,-argMax+100)

        #line them up
        a1Min = a1X[np.argmin(a1Y)]
        a3Min = a3X[np.argmin(a3Y)]
        maxDiff = a1Min-a3Min
        maxDiff += 0.3 #just a litte nudge

        #plot it
        ax[0].cla()
        ax[0].set_title(chan,fontsize = 30)
        ax[0].plot(a3X,a3Y*scale,'.-',lw=3,label="ANITA3",color="red")
        ax[0].plot(a1X-maxDiff,a1Y,'.',label="ANITA1",color="green")
        ax[0].plot(a3X-maxDiff,a1YInterp,'-',color="green")
        ax[0].legend()
        ax[0].set_xlabel("time (ns)")
        ax[0].set_ylabel("voltage (arbitrary)")
#        ax[0].set_xlim([-20,60])
        
        ax[1].cla()
        ax[1].plot(a3F,tf.calcLogMag(a3F,a3FFT),label="ANITA3",color="red")
        ax[1].plot(a1F,tf.calcLogMag(a1F,a1FFT)-50,label="ANITA1",color="green") #line em up with nudge
        ax[1].legend()
        ax[1].set_xlabel("frequency (GHz)")
        ax[1].set_ylabel("gain (dB)")
        ax[1].set_xlim([0,1.5])
        ax[1].set_ylim([-100,0])
#        ax[1].set_autoscale_on(False)
    
        ax[2].cla()
        ax[2].plot(a3F,tf.calcGroupDelay(a3FFT),label="ANITA3",color="red")
        ax[2].plot(a1F,tf.calcGroupDelay(a1FFT),label="ANITA1",color="green")
#        ax[2].plot(a3F[1:],tf.calcGroupDelay(a3FFT),label="ANITA3",color="red")
#        ax[2].plot(a1F[1:],tf.calcGroupDelay(a1FFT),label="ANITA1",color="green")
        ax[2].set_xlabel("Frequency (GHz)")
        ax[2].set_ylabel("Group Delay (ns)")
        ax[2].set_xlim([0,1.5])

        fig.savefig("autoPlots/"+chan+".png")


def plotCompare(allChans,savePlots=False,ant=False,sig=False):
    """
      Plot a comparison of all the channels in a way that is sort of easy to see
      Also has a good way to exclude channels that you think are "bad" (so you
      can fix them later of course)
    """

    try:
        sns.set_palette(sns.color_palette("husl", n_colors=16)[::-1])
    except:
        print "couldn't set palette"
        pass

    lab.close("all")

    figV,axV = lab.subplots(3,sharex=True,figsize=[11,8.5])
    figH,axH = lab.subplots(3,sharex=True,figsize=[11,8.5])

    figFFTH,axFFTH = lab.subplots(3,sharex=True,figsize=[11,8.5])
    figFFTV,axFFTV = lab.subplots(3,sharex=True,figsize=[11,8.5])

    axFFTH[2].set_xlabel("Frequency (GHz)")
    axFFTH[1].set_ylabel("Gain (dB)")

    axFFTV[2].set_xlabel("Frequency (GHz)")
    axFFTV[1].set_ylabel("Gain (dB)")

    chans = allChans.keys()
    chans.sort()

    """
      why are these bad?????
      calibration runs don't seem to have pulses in them!
      for 67dB runs
      badChans = ["01MH","02TH","03MH","04MV","04TV",
                  "05MH","05TH","05TV","08MH","11TV",
                  "12BV","12TV","13MH","15TH","16TH"]
      
      for 47dB runs (which have a tail that runs off the end of the window so
      these aren't great either)
    """

    for chan in chans:

        waveX = allChans[chan][0]
        waveY = allChans[chan][1]
        
        #make them start at the beginning

        #get them to overlay
#        argMax = np.argmax(waveY)
        argMax = np.argmax(tf.correlation(allChans["01BH"][1],allChans[chan][1]))
#        print argMax
        if chan!="01BH":
            waveY = np.roll(waveY,-argMax)
            
        f,fft = tf.genFFT(waveX,waveY)
        logMag = tf.calcLogMag(f,fft)
        
#        if any(x for x in badChans if chan in x):
#            continue

        lenWave = len(allChans[chan][0])

        if  (chan[2] == "T"):
            axIndex = 0
        elif (chan[2] == "M"):
            axIndex = 1
        elif (chan[2] == "B"):
            axIndex = 2

        if   (chan[3] == "V"):
                axV[axIndex].plot(waveX,waveY,label=chan[:2])
                axFFTV[axIndex].plot(f,logMag,label=chan[:2])

        elif (chan[3] == "H"):
#            if (chan[:2] == "13" and axIndex==2):
#                axH[axIndex].plot(waveX,waveY,label=chan[:2],color="red",lw=2)
#                axFFTH[axIndex].plot(f,logMag,label=chan[:2],color="red",lw=2)
#            else:
            axH[axIndex].plot(waveX,waveY,label=chan[:2])
            axFFTH[axIndex].plot(f,logMag,label=chan[:2])

    #labels and making it look pretty!
    
    axV[0].set_title("Vertical\nTop")
    axH[0].set_title("Horizontal\nTop")
    axV[1].set_title("Middle")
    axH[1].set_title("Middle")
    axV[2].set_title("Bottom")
    axH[2].set_title("Bottom")

    axFFTV[0].set_title("Vertical\nTop")
    axFFTH[0].set_title("Horizontal\nTop")
    axFFTV[1].set_title("Middle")
    axFFTH[1].set_title("Middle")
    axFFTV[2].set_title("Bottom")
    axFFTH[2].set_title("Bottom")

    axV[0].set_xlim([0,30])
    axH[0].set_xlim([0,30])

    axV[2].set_xlabel("Time (ns)")
    axH[2].set_xlabel("Time (ns)")
    axV[2].set_xticks(np.arange(0,30,1))
    axH[2].set_xticks(np.arange(0,30,1))

    if ant:
        axV[1].set_ylabel("Effective antenna height (m)")
        axH[1].set_ylabel("Effective antenna height (m)")

        axFFTH[1].set_ylabel("Effective antenna height $(log_{10}(m))$")
        axFFTV[1].set_ylabel("Effective antenna height $(log_{10}(m))$")

    elif sig:
        axV[1].set_ylabel("Gain (unitless)")
        axH[1].set_ylabel("Gain (unitless)")
        
        axFFTH[1].set_ylabel("Gain (dB)")
        axFFTV[1].set_ylabel("Gain (dB)")

    elif tf:
        axV[1].set_ylabel("Normalized instrument height (m)")
        axH[1].set_ylabel("Normalized instrument height (m)")
        
        axFFTH[1].set_ylabel("Normalized instument height  $(log_{10}(m))$")
        axFFTV[1].set_ylabel("Normalized instrument height $(log_{10}(m))$")


    axFFTV[0].set_xlim([0,2])
    axFFTH[0].set_xlim([0,2])

    axFFTV[2].set_xlabel("Frequency (GHz)")
    axFFTH[2].set_xlabel("Frequency (GHz)")

    

    for i in range(0,3):
        if ant:
            axFFTV[i].set_ylim([-30,5])
            axFFTH[i].set_ylim([-30,5])
        elif sig:
            axFFTV[i].set_ylim([30,65])
            axFFTH[i].set_ylim([30,65])
        else:
            axFFTV[i].set_ylim([15,50])
            axFFTH[i].set_ylim([15,50])


    axV[0].legend(frameon=True,ncol=2,bbox_to_anchor=(1.1,1.3))
    axH[0].legend(frameon=True,ncol=2,bbox_to_anchor=(1.1,1.3))

    axFFTV[0].legend(frameon=True,ncol=2,bbox_to_anchor=(1.1,1.3))
    axFFTH[0].legend(frameon=True,ncol=2,bbox_to_anchor=(1.1,1.3))

    figV.show()
    figH.show()
    figFFTV.show()
    figFFTH.show()

    if savePlots:
        figV.savefig("plotCompare_timeV.pdf")
        figH.savefig("plotCompare_timeH.pdf")
        figFFTV.savefig("plotCompare_fftV.pdf")
        figFFTH.savefig("plotCompare_fftH.pdf")

                
    return



def plotAntennaGain(allAnts,savePlots=False):

    #Produce a plot for the gain normalized to dBi (hence the 4pi)
    
    figGain,axGain = lab.subplots()
    figGain.suptitle("ANITA3 Antenna Gain")
    
    axGain.set_ylabel("Antenna Gain (dBi)")
    axGain.set_xlabel("Frequency (GHz)")
    axGain.set_ylim([-20,20])
    axGain.set_xlim([0,3])

    for chan in allAnts:
        x,y = allAnts[chan]
        f,fft = tf.genFFT(x,y)
        gainTF = tf.calcAntGain(f,fft)


        if chan == "01BH":
            axGain.plot(f,gainTF,label="H",color="red")
        if chan == "01BV":
            axGain.plot(f,gainTF,label="V",color="blue")
        if chan[-1] == "H":
            axGain.plot(f,gainTF,color="red")
        if chan[-1] == "V":
            axGain.plot(f,gainTF,color="blue")


    axGain.legend()
    figGain.show()

    if savePlots:
        figGain.savefig("plots/doPalAnt_AntGain"+chan+".png")
        



def plotOneAtATime(allChans):
    fig,ax = lab.subplots(2)
    fig.show()
    for chan in allChans:
        x,y = allChans[chan]
        ax[0].cla()
        ax[0].set_title(chan)
        ax[0].plot(x,y)
#        ax[0].set_ylim([-150,150])
        f,fft = tf.genFFT(x,y)
        ax[1].cla()
        ax[1].plot(f,tf.calcLogMag(f,fft))
#        ax[1].set_ylim([10,70])
#        ax[1].set_xlim([0,2])
        fig.canvas.draw()
        raw_input()


#==============================================================================
# Stupid test things
#==============================================================================

def phaseShifts(waveform):
    """
      The group delay and magnitude are what is effected by the signal chain
      and antennas, the absolute phase is not
      I'm just playing with creating more causal signals by uniformly shifting
      the phase of the entire complex domain
      It doesn't really work and isn't really the right thing to do, but is fun
      to mess with
    """

#    shifts = np.arange(0,np.pi*2,0.01)
    shifts = np.arange(0,np.pi/2,0.1)

    shiftedArrayMaxes = []
    waveformFFT = tf.fftw.rfft(waveform[1])
    cnt = 0
    fig2,ax2 = lab.subplots()
    for shift in shifts:
        print shift
        shiftedFFT = tf.fftPhaseShift(waveformFFT,shift)
        shiftedWave = tf.fftw.irfft(shiftedFFT)
        ax2.plot(waveform[0]+cnt*waveform[0][-1],np.roll(shiftedWave,512))
        shiftedArrayMaxes.append(np.max(shiftedWave))
        cnt += 1
        
    fig,ax = lab.subplots()
    ax.plot(shifts,shiftedArrayMaxes)
    fig.show()
    fig2.show()
    
    return shifts[np.argmax(shiftedArrayMaxes)]
        

def importAndPlotSurfs():
    """
      The surf waveforms are getting all cut off at the end of the window 
      (because of the ROOT correlation code)
      I gotta see how often it happens and ensure that it DOESN'T because it
      makes the TF wrong if the pulse occurs early in the window
    """

    chans = np.loadtxt("chanList.txt",dtype=str)

    fig,ax = lab.subplots(6,8)
    fig2,ax2 = lab.subplots(6,8)

    index=0
    for chan in chans:
        try:
            x,y,f,fft = importSurf(chan)
        except:
            print chan+" missing"

        if index<48:
            ax[index/8][index%8].plot(x,y)
            ax[index/8][index%8].set_title(chan)
            ax[index/8][index%8].set_xticklabels("")
            ax[index/8][index%8].set_yticklabels("")
        else:
            ax2[(index-48)/8][index%8].plot(x,y)
            ax2[(index-48)/8][index%8].set_title(chan)
            ax2[(index-48)/8][index%8].set_xticklabels("")
            ax2[(index-48)/8][index%8].set_yticklabels("")
            
        index += 1

    fig.show()
    fig2.show()

    return


def fourierZeroPadding():

    a,b,f,fft = getCables("A-B_PULSER-SCOPE.s2p")
    print len(fft)
    f2,fft2 = tf.fourierZeroPad(f,fft,len(fft)*3)

    fig,ax = lab.subplots(2)
    ax[0].plot(f,np.real(fft),color="blue")
    ax[0].plot(f,np.imag(fft),'--',color="blue")

    ax[0].plot(f2,np.real(fft2),color="red")
    ax[0].plot(f2,np.imag(fft2),'--',color="red")

    ax[1].plot(*tf.genTimeSeries(f,fft))
    ax[1].plot(*tf.genTimeSeries(f2,fft2))

    fig.show()

    return f, fft, f2, fft2


def fourierCausalPadding():
    f,cableGainLin,cablePhase = s2pParser(cablesBaseDir + "A-C_PULSER-TEST_40DB.s2p")
#    cablePhase -= cablePhase[1]
    fft = tf.gainAndPhaseToComplex(cableGainLin,cablePhase)
    
    ampliToExtrapolate = np.absolute(fft[-1])
    print ampliToExtrapolate
    f2,fft2 = tf.fourierPad(f,fft,len(fft)*9,ampliToExtrapolate)

    fig,ax = lab.subplots(3)
    ax[0].plot(f,np.real(fft),label="real")
    ax[0].plot(f,np.imag(fft),'.',label="imag")

    ax[1].plot(*tf.genTimeSeries(f,fft),ls='-.',label="noPad")
    ax[1].plot(*tf.genTimeSeries(f2,fft2*5),label="padding")
    ax[1].legend()

    ax[2].plot(f2,np.real(fft2),label="real")
    ax[2].plot(f2,np.imag(fft2),'.',label="imag")
    ax[2].legend()

    fig.show()

    return f, fft, f2, fft2

    
def calibrationPulseSignal(chan):

    #get scope data from waveforms (extracted from raw data from another script)
    scopeRawX,scopeRawY,scopeRawF,scopeRawFFT = importScope(chan)
    scopeX,scopeY = processWaveform(scopeRawX,scopeRawY,"calScope")
    scopeF,scopeFFT = tf.genFFT(scopeX,scopeY)
    #1.26 is the max coherance phase?
#    scopeFFT = tf.fftPhaseShift(scopeFFT,1.26)

    #Get the cable's (H(f) transfer function for pulser to scope)
    #The 17 (for the "center") is from tf.compPhaseShifts3(), which makes a nice film of where the 
    # phase center is
    print "Getting Pulser to Scope Cables..."
    cableScopeF,cableScopeFFT= getCables("A-B_PULSER-SCOPE.s2p",tZero=17,resample="interp")

    #deconvolve cable pulse * cable = scope -> pulse = scope/cable
    #finds just the pulser impulse
    pulseDeconvFFT = scopeFFT/cableScopeFFT

    return cableScopeF, pulseDeconvFFT


def weinerDeconv(sigInX, sigInY, chan):
    """
      A weiner deconvolution!
      sigIn: The signal you want to deconvolve (from SURF probably, in time domain?)
      chan, which loads:
      tfIn:  The transfer function of the system (from autoPlots, in time domain)
      noiseIn: The snr for each frequency (from surfNoise, in log mag)
      
      and uses those to generate snrIn, which is then used to do deconvolution
    """
    baseDir = "/Users/brotter/Science/ANITA/ANITA3/benCode/analysis/impulseResponse/integratedTF/"
    tfIn = np.loadtxt(baseDir+"autoPlots/A3ImpulseResponse/"+chan+".txt").T
    noiseInF,noiseLogMag = np.loadtxt(baseDir+"noiseCal/chan_"+chan+".txt").T

    tfF,tfFFT = tf.genFFT(*tfIn)
    tfMag = np.absolute(tfFFT)

    #noise from the surf is stored in dBm!!!! (not dBm/Hz)
    noiseF = noiseInF/1000. #in us / MHz which is stupid, so I'll change it (power already f dependant)
    
    sigX,sigY,sigF,sigFFT = importSurf(chan) #this has sigY in Volts (thus sigFFT is in that unit scale too)
    signalLogMag = tf.calcLogMag(sigF,sigFFT)+70
    signalLogMag[:16] = np.ones(16)*signalLogMag[0]
    #lets just say this pulser thing has no low frequency power

    snrMag = signalMag / noiseMag

    #weiner deconv!
    #G is the deconvolution multiplier I guess
    weiner = (tfMag**2)
    weiner /= (tfMag**2)+(1./snrMag)
    G = weiner/tfFFT

    sigInF,sigInFFT = tf.genFFT(sigInX,sigInY)

    sigOutFFT = sigInFFT * G
    sigOutX,sigOutY = tf.genTimeSeries(sigInF,sigOutFFT)
    
    fig,ax = lab.subplots(2)
    ax[0].plot(sigInX,sigInY)
    ax[1].plot(sigOutX,sigOutY)
    fig.show()

    print "inputPeakRatio",np.max(sigInY**2)/np.mean(sigInY**2)
    print "outputPeakRatio",np.max(sigOutY**2)/np.mean(sigOutY**2)

    return sigOutX,sigOutY,sigInFFT,sigOutFFT,G,snrMag,tfMag


def testWeiner():

    chan = "01BH"
    x,y,f,fft = importSurf(chan)
    outX,outY,a,b,c,d,e = weinerDeconv(x,y,chan)

    inOutScale = np.max(y) / np.max(outY)

#    lab.close("all")
    fig,ax = lab.subplots(2)
    ax[0].plot(np.log10(a),label="in")
    ax[0].plot(np.log10(b),label="out")
    ax[0].plot(np.log10(c),label="g")
    ax[0].plot(np.log10(d),label="snr")
    ax[0].plot(np.log10(e),label="tf")
    ax[0].legend()
    ax[1].plot(x,y,label="measured Pulse")
    ax[1].plot(outX,outY*inOutScale,label="deconvolved")
    ax[1].legend()

    fig.show()

    return outX, outY


def makeWaveform(length=513,dF=5./512):
    """
      Can I like make a pulse?
      
      the magnitude should be 1 in our range, and 0 out of it
      the phase should meet the requirement of:
      group delay increases by 30ns across our band (I guess) and is always positive

      Tg(w) = -dPhase / dw
      w = 2*pi*f -> f = w/(2*pi)
      Tg(f) = -dPhase*2*pi / df
      
      anyway this isn't how it actually works but is fun to play with
    """
    
    #okay so this should get 10GS/s sampling for 1024 samples
    #in units of GHz
    fArray = np.arange(0,length)*dF

    magArray = np.ones(length)

    Tgw = -1
    for i in range(0,len(fArray)):
        f = fArray[i]
        if (f > .15) and (f < .2):
            magArray[i] *= (-np.cos(((f-.15)*np.pi)/.05)+1)/2.
        elif (f > .2) and (f < 1.2):
            pass
        elif (f > 1.2) and (f < 1.25):
            magArray[i] *= (np.cos(((f-1.2)*np.pi)/.05)+1)/2.
        else:
            magArray[i] = 0


# some random guess as to phase structure
#    for i in range(0,len(fArray)):
#        f = fArray[i]
#        Tgw += ((3/102.4)*(1-magArray[i]))/(np.pi*2.)
#        phaseArray.append(phaseArray[-1]-Tgw)


    #constant second order quadratic phase increase
#    phaseArray = (25*(np.arange(0,len(fArray))/float(len(fArray))))**2

            
    phaseArray = (2*fArray)**(1-np.absolute(np.cumsum(np.gradient(magArray))))

    print phaseArray

        
    fft = tf.gainAndPhaseToComplex(magArray,phaseArray)

    outX,outY = tf.genTimeSeries(fArray,fft)

    fig,ax = lab.subplots(2)

    ax[0].set_title("Making a waveform from scratch")
    ax[0].plot(outX,np.roll(outY,300))
    ax[0].set_xlabel("time (ns)")
    ax[0].set_ylabel("voltage (arbitrary)")
    ax[1].plot(fArray,phaseArray,color="blue")
    ax[1].set_xlabel("frequency (GHz)")
    ax[1].set_ylabel("phase (unwrapped radians)",color="blue")
    ax11 = ax[1].twinx()
    ax11.set_ylabel("linear gain",color="red")
    ax11.plot(fArray,magArray,'.',color="red")
    
    fig.show()


    return outX,outY


def compPhase(allChans):

    fig,ax = lab.subplots(2,2)
    
    for chan in allChans.keys():
        axIndex=0
        if (chan[-1] == "H"): axIndex=1
        f,fft = tf.genFFT(allChans[chan][0],allChans[chan][1])
        if (chan=="05MH"):
            ax[1][axIndex].plot(f,tf.calcLogMag(f,fft),label=chan,color="red")
            ax[0][axIndex].plot(f,tf.calcPhase(tf.minimizeGroupDelayFromFFT(f,fft,highLim=256)),color="red")
        else:
            ax[1][axIndex].plot(f,tf.calcLogMag(f,fft),label=chan)
            ax[0][axIndex].plot(f,tf.calcPhase(tf.minimizeGroupDelayFromFFT(f,fft,highLim=256)))
            



    ax[0][0].legend()
    ax[1][0].legend()

    fig.show()






def plotSavedFiles(baseName="avgSurfWaveform",dir="waveforms_new"):
    #Just plots the saved waveforms against each other overlayed

    sns.set_palette(sns.color_palette("husl",16))

    fig,ax = lab.subplots(2,3,sharex=True)

    files = glob(dir+"/*"+baseName+"*")


    allChans = {}

    for file in files:
        name = file.split("/")[1].split("_")[1].split(".")[0]
        print name
        data = np.loadtxt(file).T
        allChans[name] = data
        if name[-1] == "H":
            if name[-2] == "T":
                ax[0][0].plot(data[0],data[1],label=name[:2])
            if name[-2] == "M":
                ax[0][1].plot(data[0],data[1],label=name[:2])
            if name[-2] == "B":
                ax[0][2].plot(data[0],data[1],label=name[:2])
        else:

            if name[-2] == "T":
                ax[1][0].plot(data[0],data[1],label=name[:2])
            if name[-2] == "M":
                ax[1][1].plot(data[0],data[1],label=name[:2])
            if name[-2] == "B":
                ax[1][2].plot(data[0],data[1],label=name[:2])

    ax[0][0].set_title("H Top")
    ax[0][1].set_title("H Middle")
    ax[0][2].set_title("H Bottom")
    ax[1][0].set_title("V Top")
    ax[1][1].set_title("V Middle")
    ax[1][2].set_title("V Bottom")


    ax[0][2].legend(bbox_to_anchor=(1.1, 1.05))

    fig.show()

    return allChans
