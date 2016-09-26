import numpy as np
import pylab as lab
import scipy.signal as signal
from scipy.interpolate import interp1d

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


#signal chain data
waveformDir = "/Users/brotter/benCode/impulseResponse/integratedTF/waveforms/"

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
def findPalestineAntennaFile(chan):
    
    antennaInfo = np.loadtxt("antennaInformation_cut.csv",dtype=str,delimiter=",")

    antennaNumber = int(antennaInfo.T[1][np.argwhere(antennaInfo.T[2]==chan[:3])[0][0]][1:])

    print chan+" maps to antenna number "+str(antennaNumber)
        
    dir = "/Volumes/ANITA3Data/palestine14/SeaveyAntennas/S21s/"
        
    fileName = dir+"hpol_ezLinks/rxp"+str(antennaNumber).zfill(2)+"_coPol.csv"

    return fileName


##########################
def importPalAntIn(chan):
    #Palestine Antennas
        # there is only one waveform :( 
        # so I can't do any averaging and I just import the raw data I guess
        # Also this means I probably need to link the whole dir (since they are far away)

    fileName = findPalestineAntennaFile(chan)

    dataX,dataY = np.loadtxt(fileName,comments="\"",delimiter=",",usecols=(3,4)).T


    #make the time range start at zero
    dataX -= dataX[0]

    #also I always want to have the fourier spaces of these things
    dataF,dataFFT = tf.genFFT(dataX,dataY)


    return dataX,dataY,dataF,dataFFT

#############################
def importPalAntOut(chan):

    fileName = findPalestineAntennaFile(chan)

    dataX,dataY = np.loadtxt(fileName,comments="\"",delimiter=",",usecols=(3,4)).T

    dataX *= 1e9 #s ->ns
    if source.lower() == "palout":
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





#####################################
def getCables(fileName,phaseShift=0):
    cableFreq,cableGainLin,cablePhase = s2pParser(cablesBaseDir + fileName)


    #shift phase?
    cablePhase = cablePhase+phaseShift

    cableFFT = tf.gainAndPhaseToComplex(cableGainLin,cablePhase)
    cableF2,cableFFT = tf.complexZeroPadAndResample(cableFreq,cableFFT)




    return cableF2,cableFFT

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
    phase = tf.unwrapPhase(phase)
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


    #find max and do hanning window with SUPER SHART EDGES(cuts waveform to desired length too)
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

    
def genPalAntTF(chan):
        #This should get the aligned and averaged waveforms, then process them so they have the same record length and dTs
        #It does aliasing filtering and windowing as well, so everything should be pretty and wonderful to compare the two

        """
        You really only need to store everything as a class variable if you're going to plot it later
        or if you want to be able to access it later
        """


        print "Importing Data..."
        inRawX,inRawY,inRawF,inRawFFT = importPalAntIn
        outRawX,outRawY,outRawF,outRawFFT = importData(outSource)


        print "Processing Data..."
        inX,inY = processWaveform(inRawX,inRawY,inSource)
        inF,inFFT = tf.genFFT(inX,inY)

        outX,outY = processWaveform(outRawX,outRawY,outSource)
        outF,outFFT = tf.genFFT(outX,outY)


        #calculate group delay difference bewteen the two
        meanGD = calcMeanGD(outF,outFFT,inF,inFFT)

        #print out some sanity checks
        sanityCheck()

        #generate transfer function
        TFX,TFY,TFF,TFFFT = computeTF(inF,inFFT,outF,outFFT)


        return 


##################################
def doSigChainWithCables(chan):
    #phase shifts are likely pointless!
    #they were: scopeFFT=1.26
    #           phaseShift=1.698


    #get scope data from waveforms (extracted from raw data from another script)
    scopeRawX,scopeRawY,scopeRawF,scopeRawFFT = importScope(chan)
    scopeX,scopeY = processWaveform(scopeRawX,scopeRawY,"calScope")
    scopeF,scopeFFT = tf.genFFT(scopeX,scopeY)
    #1.26 is the max coherance phase?
#    scopeFFT = tf.fftPhaseShift(scopeFFT,1.26)

    #Get the cable's (H(f) transfer function for pulser to scope)
    #phaseShift was calculated using compPhaseShifts2(), it is when there is peak coherence
    P2SF,P2SFFT= getCables("A-B_PULSER-SCOPE.s2p")

    #deconvolve cable pulse * cable = scope -> pulse = scope/cable
    #finds just the pulser impulse
    pulseDeconvFFT = scopeFFT/P2SFFT


    #Get the cable's (H(f) transfer function for pulser to AMPA)
    P2AF,P2AFFT= getCables("A-C_PULSER-TEST_66DB.s2p")
    
    #convolve it with that transfer function to get pulse at AMPA
    ampaInputFFT = P2AFFT*pulseDeconvFFT

    #get the surf (extracted from ROOT data from another script)
    surfRawX,surfRawY,surfRawF,surfRawFFT = importSurf(chan)
    surfX,surfY = processWaveform(surfRawX,surfRawY,"surf")
    surfF,surfFFT = tf.genFFT(surfX,surfY)

    #deconvolve signal chain transfer function out of that!
    tfFFT = surfFFT/ampaInputFFT
    

    return surfF,tfFFT 

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



########################
def doSigChainAndAntenna(chan):
    #get antenna
    antX,antY,antF,antFFT = doRoofAntWithCables()

    #get sig chain
    sigChainF,sigChainFFT = doSigChainWithCables(chan)

    #convolve the two (full ANITA3 transfer function!)
    a3F = sigChainF
    a3FFT = sigChainFFT * antFFT
    #clean it up a bit
    a3FFT = np.concatenate((a3FFT[:171],np.zeros(342)))
    a3X = antX
    a3Y = tf.fftw.irfft(a3FFT)
    #get the pulse at the beginning
    yMax = np.argmax(a3Y)
    a3Y = np.roll(a3Y,50-yMax)
    #make it look nice and cut off the shit at the end
    a3X,a3Y = tf.hanningWindow(a3X,a3Y,np.argmax(a3Y)+100,totalWidth=500,slope=50)
    a3X,a3Y = tf.zeroPadEqual(a3X,a3Y,1024)

    

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
                file.write(str(allChans[chan][0][i])+" "+str(allChans[chan][1][i])+"n")
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
    

        #line them up
        a1Max = a1X[np.argmax(a1Y)]
        a3Max = a3X[np.argmax(a3Y)]
        maxDiff = a1Max-a3Max
        maxDiff += 0.3 #just a litte nudge

        #scale them
        a1MaxVal = np.max(a1Y)
        a3MaxVal = np.max(a3Y)
        scale = a1MaxVal/a3MaxVal
        

        #plot it
        ax[0].cla()
        ax[0].set_title(chan)
        ax[0].plot(a3X,a3Y*scale,'.-',lw=3,label="ANITA3",color="red")
        ax[0].plot(a1X-maxDiff,a1Y,'.',label="ANITA1",color="green")
        ax[0].plot(a3X-maxDiff,a1YInterp,'-',color="green")
        ax[0].legend()
        ax[0].set_xlabel("time (ns)")
        ax[0].set_ylabel("voltage (arbitrary)")
        ax[0].set_xlim([-20,60])
        
        ax[1].cla()
        ax[1].plot(a3F,tf.calcLogMag(a3F,a3FFT),label="ANITA3",color="red")
        ax[1].plot(a1F,tf.calcLogMag(a1F,a1FFT)+7.5,label="ANITA1",color="green") #line em up with nudge
        ax[1].legend()
        ax[1].set_xlabel("frequency (GHz)")
        ax[1].set_ylabel("gain (dB)")
        ax[1].set_xlim([0,2])
        ax[1].set_ylim([0,60])
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

        f,logMag = tf.genLogMag(allChans[chan][0],allChans[chan][1])
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
                axV[axIndex].plot((np.arange(0,lenWave)-argMax)*0.1,allChans[chan][1],label=chan[:2])
        elif (chan[3] == "H"):
                axH[axIndex].plot((np.arange(0,lenWave)-argMax)*0.1,allChans[chan][1],label=chan[:2])



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
