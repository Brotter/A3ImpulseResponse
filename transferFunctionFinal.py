import numpy as np
import pylab as lab

from scipy import signal
from scipy.interpolate import interp1d
from scipy.interpolate import Akima1DInterpolator

from glob import glob

import tfUtils as tf

import cmath #for complex sqrt

import time #for drawing with a pause

try:
    import seaborn as sns #just a pretty plotting package, it isn't important
except:
    print "You don't have seaborn but thats okay!"


"""
Ben Rotter
University of Hawaii at Manoa
ANITA3
February 2016 - through September 2016
"""

"""
So this is a ton of my code all glommed together into one giant library.  This should make it easy to generate the full transfer function from all the different waveforms, and for all the channels easily

It uses tfUtils (which I might turn into a auto-linked brotterlib at some point) which has a lot of FFT analysis and fun tools to mess with

It uses data from avgSigChainWaveform and avgAntWaveform, which imports, correlates and averages the waveforms using eventReaderRoot and libFftwWrapper, things that I don't want to recreate.


This is a selection of a much larger code base that does all the work, which is available upon request if you want to see a ton of debugging and testing functions, in addition to some classes that make plotting easier.

"""


#where the cable info is stored
cablesBaseDir = "/Volumes/ANITA3Data/antarctica14/S21ExternalRFChainForSingleChannelCallibration/"


#signal chain data (57dB seems to be a happy middle power)
waveformDir = "/Users/brotter/benCode/impulseResponse/integratedTF/waveforms_57dB/"

#directories with antenna pulse data
localDir = "/Volumes/BenANITA3Data/"

roofDir = "Rooftop_Seevey_Antenna/Antenna_Impulse_Testing/"
currLocalRoofDir=localDir+roofDir+"Attempt_2_02July2012/Hpol/"
remoteRoofDir = "/storageb/ANITA-3/"+roofDir
currRemoteRoofDir=remoteRoofDir+"Attempt_2_02July2012/Hpol/"

chamberDir="Anechoic_Chamber_New_Horn_Impulse/ANITA2Antennas/HPol_Pulse"
remoteChamberDir="/storage/Krauss_Backup/ANITA/"+chamberDir

chamberIdentDir="Anechoic_Chamber_New_Horn_Impulse/IdenticalAntennas/Ant2Blue_to_Ant0/HPol_Pulse/HPol_Pulse_New_Cable/"
chamberRefPulse = localDir+chamberIdentDir+"Impulse_NewRedCable_23May2012.dat"

chamberBoresightPulse = localDir+chamberIdentDir+"Impulse_p180Deg_t090Deg_22May2012_H_Two.dat"

###############################################################################################################
#Importing things


def importANITA1():
    a1 = np.loadtxt("ANITA1ImpulseResponse.dat")
    return a1.T[0]*1e9,a1.T[1]

###################################
def findPalestineAntennaFile(chan,inOrOut):
    
    antennaInfo = np.loadtxt("antennaInformation_cut.csv",dtype=str,delimiter=",")

    antennaNumber = int(antennaInfo.T[1][np.argwhere(antennaInfo.T[2]==chan[:3])[0][0]][1:])

    print chan+" maps to antenna number "+str(antennaNumber)
        
    dir = "/Volumes/ANITA3Data/palestine14/SeaveyAntennas/S21s/"
        
    fileName = dir+chan[-1].lower()+"pol_ezLinks/rxp"+str(antennaNumber).zfill(2)
    if inOrOut == "in":
        fileName = fileName +"_inputPulse.csv"
    if inOrOut == "out":
        fileName = fileName +"_coPol.csv"
 
    return fileName


##########################
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

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)


    return dataX,dataY,dataF,dataFFT

#############################
def importPalAntOut(chan):

    fileName = findPalestineAntennaFile(chan,"out")

    dataX,dataY = np.loadtxt(fileName,comments="\"",delimiter=",",usecols=(3,4)).T

    dataX *= 1e9 #s ->ns
    dataY *= 100 #there is something super wrong with the scaling, 20mV Vpp, 15uV resolution?  Weird...
    dataY *= -1 #also the polarity is flipped?  Look in log book to see if this is right...


    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)


    return dataX,dataY,dataF,dataFFT


#############
def importChamberAnts(fileName=localDir+roofDir+"HPulse_10degCCW_148p_090t_10dBatt.dat",numEvents=-1,channel=0):
    """
    I always need to write more parsers, always be parsing (ABP)
    """
    print "Importing locally -> "+fileName

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


######################
def importRoofAntIn():
    #Hawaii Rooftop Antenna Input Pulse (FID)
    fileName = "avgFid.txt"
    dataX,dataY = np.loadtxt(fileName).T

    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)
        
    return dataX,dataY,dataF,dataFFT
                
#######################
def importRoofAntOut():
    #Hawaii Rooftop Antenna Input Pulse (Ant)
    fileName = "avgAnt.txt"
    dataX,dataY = np.loadtxt(fileName).T

    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)
        
    return dataX,dataY,dataF,dataFFT


########################
def importSurf(chan):
    #LAB chip from the SURF
    fileName = waveformDir+chan+"_avgSurfWaveform.txt"
    dataX,dataY = np.loadtxt(fileName).T
    
    dataY /= 1000 #mv->V

    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)
        
    return dataX,dataY,dataF,dataFFT
    

#########################
def importScope(chan):
    #Calibration Scope
    fileName = waveformDir+chan+"_avgScopeWaveform.txt"
    dataX,dataY = np.loadtxt(fileName).T

    dataX *= 1e9 #s ->ns

    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)
        
    return dataX,dataY,dataF,dataFFT


#################################################
def getCables(fileName,dF=5000./512.,length=513):
    """
    New getCables():

    imports the data (correctly with the s2pParser fix), then resamples the gain and phase
    to meet the specified dF and length (which should stay constant because the time domain I choose defines
    them)
    phase -> just a constant group delay (and I don't care about absolute phase)
    gain -> interpolated inside the range you can, zero outside? (sure why not)

    Lets skip causality for now (that seems to be a bigger problem at the end anyway)

    """


    cableFreq,cableGainLin,cablePhase = s2pParser(cablesBaseDir + fileName)

    cableNewFreq = np.arange(0,length)*dF
    cableNewPhase = tf.regenerateCablePhase(cablePhase)
    cableNewGainLin = tf.regenerateCableLinMag(cableFreq,cableGainLin)

    cableFFT = tf.gainAndPhaseToComplex(cableNewGainLin,cableNewPhase)

    return cableNewFreq,cableFFT




#####################################
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




    #Resampling
    ##########################
    #this might do bad things...
    #the goal of it was to get the time binning and window length the same
    #for an irfft, but this just worstens the fractional phase offset error at first...
    #-> Yeah even for a "phase aligned" waveform it messes it up
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

########################
def s2pParser(fileName):

    freq = []
    gainLin = []
    phase = []
    with open(fileName) as inFile:
        for line in inFile:
            if (line[0] !='!' and line[0] != '#'):
                split = line.split()
                freq.append(float(split[0])/1e9)
                phase.append(float(split[4])*(np.pi/180.))
                gainLin.append(10**(float(split[3])/20.))
                
#    fig,ax = lab.subplots(2,sharex=True)
#    ax[0].plot(phase)
#    phase = tf.unwrapPhase(phase)
#    ax[1].plot(phase)
#    fig.show()
                

    return np.array(freq),np.array(gainLin),np.array(phase)


###############################################################################################################
#Processing

##########################################
def processWaveform(graphT,graphV,source):
    """
    Single functions to do lots of things!  Input what the "source" is and it will do the type
    of processing that you want :)  This is also a good way to make a sort of "table" to deterimine
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
#        kHanningWindow = True //this cuts off a reflection pulse that is only like 35ns late
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
        kHanningWindow = True
        kPostZeroMean  = False #created problems with zero division at f=0

    #          Palestine antenna testing output
    #Raw: record length of 8000 (7994?) | samping rate of 10GS/s
    #Now: same (only one waveform, so read right into python)
    elif source.lower() == "palout":
        kPreZeroMean   = True
        kBandPassFilter= True
        kZeroPad       = True
        kHanningWindow = True
        kPostZeroMean  = True


    print "---->Starting processWaveform() on "+source

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
        graphT,graphV = tf.zeroPadEqual(graphT,graphV,1024*3)
        graphMaxV = np.argmax(graphV)
        graphT,graphV = tf.hanningWindow(graphT,graphV,graphMaxV+240,512,slope=1)
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



###############################################################################################################
#Generating impulse responses


#getting the cables is really resource intensive, so I'd like to do it once and then be done with it
#so I guess I can save it within the module?  Initialize them as False and then fill them if they
# aren't an numpy ndarray (which a boolean isn't)
P2SF = False
P2SFFT = False
P2AF = False
P2AFFT = False

##################################
def doSigChainWithCables(chan,savePlots=False):
    #phase shifts are likely pointless!
    #they were: scopeFFT=1.26
    #           phaseShift=1.698


    
    #get scope data from waveforms (extracted from raw data from another script)
    scopeRawX,scopeRawY,scopeRawF,scopeRawFFT = importScope(chan)
    scopeX,scopeY = processWaveform(scopeRawX,scopeRawY,"calScope")
    scopeF,scopeFFT = tf.genFFT(scopeX,scopeY)


    if savePlots:
        fig,ax = lab.subplots(2,figsize=(11,8.5))
        ax[0].set_title("Raw Signal Chain Calibration Pulser")
        ax[0].plot(scopeRawX,scopeRawY,label="raw scope pulse")
        ax[0].plot(scopeX,scopeY,label="processed scope pulse")
        ax[0].set_xlabel("Time (ns)")
        ax[0].set_ylabel("Voltage (V)")
        ax[0].legend()
        ax[1].plot(scopeRawF,tf.calcLogMag(scopeRawF,scopeRawFFT),label="raw scope pulse")
        ax[1].plot(scopeF,tf.calcLogMag(scopeF,scopeFFT),label="processed scope pulse")
        ax[1].set_xlabel("Frequency (GHz)")
        ax[1].set_ylabel("Spectral Power (dBm/Hz)")
        ax[1].legend()
        fig.savefig("plots/doSigChainWithCables_A.png")
        


    #1.26 is the max coherance phase?
#    scopeFFT = tf.fftPhaseShift(scopeFFT,1.26)

    #Get the cable's (H(f) transfer function for pulser to scope)
    #The 17 (for the "center") is from tf.compPhaseShifts3(), which makes a nice film of where the 
    # phase center is
    global P2SF
    global P2SFFT
    if type(P2SF) != np.ndarray:
        print "Getting Pulser to Scope Cables..."
#        P2SF,P2SFFT= getCables("A-B_PULSER-SCOPE.s2p",tZero=17,resample="interp")
        P2SF,P2SFFT= getCables("A-B_PULSER-SCOPE.s2p")



    #deconvolve cable pulse * cable = scope -> pulse = scope/cable
    #finds just the pulser impulse
    pulseDeconvFFT = scopeFFT/P2SFFT


    #Get the cable's (H(f) transfer function for pulser to AMPA)
    # again, the 973 is from tf.compPhaseShifts3()
    global P2AF
    global P2AFFT
    if type(P2AF) != np.ndarray:
        print "Getting Pulser to Ampa Cables..."
#        P2AF,P2AFFT= getCables("A-C_PULSER-TEST_66DB.s2p",tZero=True,hanning=True,resample="fftFit")
        P2AF,P2AFFT= getCables("A-B_PULSER-TEST_0DB.s2p")

    
    #convolve it with that transfer function to get pulse at AMPA
    ampaInputFFT = P2AFFT*pulseDeconvFFT


    #get the surf (extracted from ROOT data from another script)
    surfRawX,surfRawY,surfRawF,surfRawFFT = importSurf(chan)
    surfX,surfY = processWaveform(surfRawX,surfRawY,"surf")

    surfF,surfFFT = tf.genFFT(surfX,surfY)

    return surfRawX,surfRawY,surfRawF,surfRawFFT

    #deconvolve signal chain transfer function out of that!
    tfFFT = surfFFT/ampaInputFFT

    #fix the infinities and nans
    tfFFT[np.isposinf(np.real(tfFFT))] = 0
    tfFFT[np.isneginf(np.real(tfFFT))] = 0
    tfFFT = np.nan_to_num(tfFFT)

    #zero out everything above 1.3GHz because that's our nyquist limit
#    for i in range(0,len(surfF)):
#        if surfF[i] > 1.3:
#            tfFFT[i] /= 1e6

    #change it back to time domain
    tfY = tf.fftw.irfft(tfFFT)
    
    #clean up the tail and make it start at the beginning
#    tfY = np.roll(tfY,30-np.argmax(tfY))
#    tfY = tf.hanningTail(tfY,200,100)
    

    return surfX,tfY,surfF,tfFFT 





##########################
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


###################
def doPalAnt(chan):
    #
    # Go through ALL the antennas measured in palestine 2014
    #



    inX,inY,inF,inFFT = importPalAntIn(chan)
    inX,inY = processWaveform(inX,inY,"palin")
    inF,inFFT = tf.genFFT(inX,inY)

    outX,outY,outF,outFFT = importPalAntOut(chan)
    outX,outY = processWaveform(outX,outY,"palout")
    outF,outFFT = tf.genFFT(outX,outY)

    antTFX,antTFY,antTFF,antTFFFT = computeTF(inF,inFFT,outF,outFFT)


    #one antenna should be the square root of the response
    #cmath can't operate on numpy matricies so I have to iterate by hand
    #this also doesn't work... maybe I can just do it in group delay and magnitude?
    #basically since square roots are poorly defined in complex space, you NEED to do it with
    # the sines and cosines, which I added to tfUtils! (it works and is correct now)
    antTFFFT = tf.sqrtOfFFT2(antTFFFT)


    #Make causal?  no probably not, that seems like a "last thing you do" sort of technique
#    antTFFFT = tf.makeCausalFFT(antTFFFT,np.argmax(tf.fftw.irfft(antTFFFT)))

    #I have to get the time domain again after that
    antTFY = tf.fftw.irfft(antTFFFT)

    
    #clean up the tail and make it start at the beginning
    print "doing hanning window"
    antTFY = np.roll(antTFY,20-np.argmax(antTFY))
    antTFY = tf.hanningTail(antTFY,370,30)
    antTFFFT = tf.fftw.rfft(antTFY)


    return antTFX,antTFY,antTFF,antTFFFT






########################
def doSigChainAndAntenna(chan):
    #get antenna
#    antX,antY,antF,antFFT = doRoofAntWithCables()
    antX,antY,antF,antFFT = doPalAnt(chan)
    
    #get sig chain
    sigChainX,sigChainY,sigChainF,sigChainFFT = doSigChainWithCables(chan)


    #convolve the two (full ANITA3 transfer function!)
    a3F = sigChainF
    a3FFT = sigChainFFT * antFFT

    #clean it up a bit
    #something was adding a ton of 1.25GHz noise
#    for i in range(0,len(a3FFT)):
#        if a3F[i] > 1.1:
#            a3FFT[i] /= 1e4

#    a3FFT = np.concatenate((a3FFT[:171],np.zeros(342))) #this isn't causal...
    a3X = antX
    a3Y = tf.fftw.irfft(a3FFT)


    

    #make it look nice and cut off the junk at the end
#    a3Y = np.roll(a3Y,40-np.argmax(a3Y))
#    a3Y = tf.hanningTail(a3Y,300,200)
#    a3X,a3Y = tf.zeroPadEqual(a3X,a3Y,1024)
    
    

    return a3X,a3Y


#####################################
def computeTF(inF,inFFT,outF,outFFT):
    """
    Generating the transfer function is always the same, so it gets its own class
    Comparable to "deconvwnr()" in matlab (nearly identical)
    
    It has a bunch of other things you can uncomment to use, but most end up being dumb
    """
        #then generate the transfer function!
        #In(f)*H(f) = Out(f) -> H(f) = Out(f)/In(f)
    tfFFT = np.divide(outFFT,inFFT)
    tfF = inF
    
    tfY = tf.fftw.irfft(tfFFT)        

#    print tfY

    tfX = np.arange(0,len(tfY))*(1./(2.*tfF[-1]))

        #We can also try "cleaning it up" if we want, these lettered things basically do that.

        #A) the group delay offset of this fft is going to be weird, it is the addition of the phase offset for both of them (I think?)
        #so I have to roll it back a bit or it will wrap around the window, probably by the mean of the group delay!

    tfGrpDlyMean = np.mean(tf.calcGroupDelay(tfF,tfFFT))
    dT = tfX[1]-tfX[0]
    print "tfGrpDlyMeanPt="+str(tfGrpDlyMean)

        #roll it backwards, 1/8th off window (the mean is at the middle I think)
    try:
        tfY = np.roll(tfY,-(int(tfGrpDlyMean/dT)-(len(tfY)*3)/8)) 
    except:
        print "Warning in computeTF: coulnd't roll because the group delay is totally borked"
        


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

########################
def doTheWholeShebang():
    """
    Generate the transfer function for all the channels!
    """
    chans = np.loadtxt("chanList.txt",dtype=str)
    allChans = {}
    for chan in chans:
        try:
            allChans[chan] = doSigChainAndAntenna(chan)
        except:
            print chan+" FAILED"
    
    writeAll(allChans)

    return allChans

#####################
def doAllSigChains():
    """
    just does all the signal chains without the antennas
    """
    chans = np.loadtxt("chanList.txt",dtype=str)
    allChans = {}
    for chan in chans:
        try:
            allChans[chan] = doSigChainWithCables(chan)
        except:
            print chan+" FAILED************************************************"

    return allChans

####################
def doAllAntennas():
    """
    just does all the palestine antennas
    """
    chans = np.loadtxt("chanList.txt",dtype=str)
    allChans = {}
    for chan in chans:
        try:
            allChans[chan] = doPalAnt(chan)
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


########################
def findSignalToNoise():
    """
    I need to find the signal to noise ratio of the signal to do a Weiner Deconvolution
    Basically, out of band our deconvolutions introduce noise (small/small = big possibly)
    so we need to take the SNR into account.

    The signal in our transfer function is actually the signal in the input pulse!
    It is super well averaged (so there is basically no noise)

    The noise in our transfer function then will be whatever the noise spectrum of the terminated
    SURF input is?  Basically doing an average of all the FFTs of a bunch of waveforms with no signal...

    """

    



###############################################################################################################
#Plotting and writing txt files



#######################
def writeAll(allChans):
    """
    Writes all the files to disk in a space seperated variable format so that I can post them

    """
    for chan in allChans.keys():
        try:
            file = open("autoPlots/"+chan+".txt","w")
            for i in range(0,len(allChans[chan][0])):
                file.write(str(allChans[chan][0][i])+" "+str(allChans[chan][1][i])+"\n")
            file.close()
        except:
            print chan+" FAILED"





###############################
def saveAllNicePlots(allChans):
    """
    Saves a bunch of plots that I made so I can find channels that didn't work correctly
    """

    lab.rcParams['figure.figsize'] = [16.0,11.0]
    fig,ax = lab.subplots(2)

    #get the ANITA1 transfer function
    a1X,a1Y = importANITA1()
    a1Y *= 1000 #V to mV (and a nudge upwards for some reason) (also is it phase flipped(or maybe I am)?)
    
    a1F,a1FFT = tf.genFFT(a1X,a1Y)


    for chan in allChans:
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
        ax[1].plot(a1F,tf.calcLogMag(a1F,a1FFT)-65,label="ANITA1",color="green") #line em up with nudge
        ax[1].legend()
        ax[1].set_xlabel("frequency (GHz)")
        ax[1].set_ylabel("gain (dB)")
        ax[1].set_xlim([0,2])
        ax[1].set_ylim([-10,60])
        ax[1].set_autoscale_on(False)
    
        fig.savefig("autoPlots/"+chan+".png")




######################
def plotCompare(allChans):
    """
    Plot a comparison of all the channels in a way that is sort of easy to see
    Also has a good way to exclude channels that you think are "bad" (so you can fix them later of course)

    """

    lab.close("all")

    figV,axV = lab.subplots(3,sharex=True)
    figH,axH = lab.subplots(3,sharex=True)

    figFFT,axFFT = lab.subplots()

    try:
        sns.set_palette("viridis", n_colors=16)
    except:
        pass

    chans = allChans.keys()
    chans.sort()

    #why are these bad?????
    #calibration runs don't seem to have pulses in them!
    #for 67dB runs
#    badChans = ["01MH","02TH","03MH","04MV","04TV",
#                "05MH","05TH","05TV","08MH","11TV",
#                "12BV","12TV","13MH","15TH","16TH"]
                
    #for 47dB runs (which have a tail that runs off the end of the window so these aren't great either)


    for chan in chans:

        waveX = allChans[chan][0]
        waveY = allChans[chan][1]
        
        #make them start at the beginning
        argMax = np.argmax(waveY)
        waveY = np.roll(waveY,-argMax+100)


        f,logMag = tf.genLogMag(waveX,waveY)
        axFFT.plot(f,logMag)
        
#        if any(x for x in badChans if chan in x):
#            continue

        argMax = np.argmax(allChans[chan][1])
        lenWave = len(allChans[chan][0])

        if  (chan[2] == "T"):
            axIndex = 0
        elif (chan[2] == "M"):
            axIndex = 1
        elif (chan[2] == "B"):
            axIndex = 2

        if   (chan[3] == "V"):
                axV[axIndex].plot(waveX,waveY,label=chan[:2])
        elif (chan[3] == "H"):
                axH[axIndex].plot(waveX,waveY,label=chan[:2])



    #labels and making it look pretty!
    axV[0].legend(bbox_to_anchor=(1.1, 1.05))
    axH[0].legend(bbox_to_anchor=(1.1, 1.05))
    
    axV[0].set_title("Vertical\nTop")
    axH[0].set_title("Horizontal\nTop")
    axV[1].set_title("Middle")
    axH[1].set_title("Middle")
    axV[2].set_title("Bottom")
    axH[2].set_title("Bottom")

    figV.show()
    figH.show()
    figFFT.show()
                
    return



###############################################################################################################
#Stupid test things

##########################
def phaseShifts(waveform):
    """
    The group delay and magnitude are what is effected by the signal chain and antennas, the absolute phase is not
    I'm just playing with creating more causal signals by uniformly shifting the phase of the entire complex domain
    It doesn't really work and isn't really the right thing to do, but is fun to mess with
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
        

#########################
def importAndPlotSurfs():
    """
    The surf waveforms are getting all cut off at the end of the window(because of the ROOT correlation code)
    I gotta see how often it happens and ensure that it DOESN'T because it makes the TF wrong if the pulse occurs early in the window
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


    return f,fft,f2,fft2


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


    return f,fft,f2,fft2

    
    


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
    global P2SF
    global P2SFFT
    if type(P2SF) != np.ndarray:
        print "Getting Pulser to Scope Cables..."
        P2SF,P2SFFT= getCables("A-B_PULSER-SCOPE.s2p",tZero=17,resample="interp")



    #deconvolve cable pulse * cable = scope -> pulse = scope/cable
    #finds just the pulser impulse
    pulseDeconvFFT = scopeFFT/P2SFFT

    return P2SF,pulseDeconvFFT


def weinerDeconv(sigInX,sigInY,chan):
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


    snrMag = 10**(signalLogMag - noiseLogMag)


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

    lab.close("all")
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

    return outX,outY


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

    phaseArray = [0]
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

    for i in range(0,len(fArray)):
        f = fArray[i]
        Tgw += ((3/102.4)*(1-magArray[i]))/(np.pi*2.)
        phaseArray.append(phaseArray[-1]-Tgw)
        
    return fArray,magArray,phaseArray

    fft = tf.gainAndPhaseToComplex(magArray,phaseArray)

    return fArray,fft

        
    
