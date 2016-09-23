import numpy as np
import pylab as lab
import scipy.signal as signal
from scipy.interpolate import interp1d
from scipy.interpolate import Akima1DInterpolator

import copy

import pyfftw.interfaces.numpy_fft as fftw

import benMathFit as mf

#######################################################################
#
#  Ben Rotter - BenJRotter@gmail.com - University of Hawaii at Manoa
#
#  A collection of functions for doing fourier domain manipulation
#  Used primarily for generating ANITA3 transfer functions
#  Recreates many things in Ryan Nichol's FFTtools, plus a bunch of
#       useful things that I couldn't find in there
#
######################################################################





def genFFT(graphX,graphY):
    outFFT = fftw.rfft(graphY)
    outF = genFreqArray(graphX)
    
    return outF,outFFT

def genFreqArray(graphX):
    length = len(graphX)/2. + 1
    dF = 1./((graphX[1]-graphX[0])*len(graphX))

    fMax = length * dF

    return np.arange(0,fMax,dF)



#def genFreqArray(graphX):
#    dF = 1./((graphX[1]-graphX[0])*len(graphX))

#    N=len(graphX)
#    fMax = N/2 * dF

#    return np.arange(0,fMax+dF,dF)


def genLogMag(graphX,graphY):
    graphF,graphFFT = genFFT(graphX,graphY)
    logMag = calcLogMag(graphF,graphFFT)
    return graphF,logMag


def calcLogMag(graphF,graphFFT):
    """
    Returns Power Spectral Density (dBm/Hz)
    """

    dF = graphF[1]-graphF[0]
    #unfold spectrum:
    #1) discard negative frequencies (second half)
    #2) double all positive frequencies (but not DC (bin 0) and nyquist (bin length/2+1)
    #3) now you can do the log mag

    return 10.*np.log10((np.abs(graphFFT)**2/(50.*1000.*dF)))
    



def genGaussian(graphX,width=10,center=False):
    length = len(graphX)
    if (center==False or center > length or center < 0):
        center = length/2
        print "genGaussian: default length ("+str(center)+")"
        
    dT = graphX[1] - graphX[0]
    norm=1./(np.sqrt(2.*np.pi)*width)
    graphY = norm*np.exp(-(graphX-graphX[center])**2./(width)**2.)
    return graphY
    

def genDelta(graphX):
    outGraph = np.zeros(len(graphX))
    outGraph[len(graphX)/2] = 1

    return outGraph


def zeroPadBeginning(inXArray,inYArray,numPoints):
    tempXList = list(inXArray)
    tempYList = list(inYArray)

    dT = inXArray[1]-inXArray[0]
    for i in range(0,numPoints):
        tempXList.append(tempXList[-1]+dT)
        tempYList.append(np.mean(tempYList[:20]))

    return np.array(tempXList),np.array(tempYList)


def zeroPadEnd(inXArray,inYArray,numPoints):
    tempXList = list(inXArray)
    tempYList = list(inYArray)

    dT = inXArray[1]-inXArray[0]
    for i in range(0,numPoints):
        tempXList.append(tempXList[-1]+dT)
        tempYList.append(np.mean(tempYList[:20]))

    return np.array(tempXList),np.array(tempYList)


def zeroPadEqual(inXArray,inYArray,outPoints):
    """
    Zero pad the input arrays into lengths specified by outPoints
    Be aware that there are many values you can "zero" it with, but, if you're constantly zero meaning, zero might be best

    """

    switch = 1

    tempXList = list(inXArray)
    tempYList = list(inYArray)

    dT = inXArray[1]-inXArray[0]

    while len(tempXList) < outPoints:
        if switch==1:
            tempXList.append(tempXList[-1]+dT)
#            tempYList.append(np.mean(tempYList[-20:]))
            tempYList.append(0)
            switch = -1
        elif switch==-1:
            tempXList.insert(0,tempXList[0]-dT)
#            tempYList.insert(0,np.mean(tempYList[:20]))
            tempYList.insert(0,0)
            switch = 1
        else:
            print "lol how did this happen??"

    return np.array(tempXList),np.array(tempYList)

def hanningWindow(inXArray,inYArray,center,totalWidth=1000,slope=200):
    """
    Takes two input arrays (x and y, since they both will change length), a total width, and a center (and optional slope).
    Returns two arrays, a hanning windowed Y axis with a width of totalWidth centered at center,
    and the corresponding X axis with the same width.
    The slope of the hanning window is defined by slope, and defaults to 200.
    You can test it with just a flat DC line! :)
    """

#    print "Hanning: center="+str(center)+" totalWidth="+str(totalWidth)
    outXArray = []
    outYArray = []

    lenGraph = len(inYArray)
    
    #I'm doing this weird, so width is where the hanning "turns on" and totalWidth is the end length of the waveform
    width = totalWidth-(2*slope)
    
    if ( center-(width/2+slope) < 0 ):
        print "Warning in hanningWindow: center ("+str(center)+") too close to zero, moving to min ("+str(width/2)+")"
        center = totalWidth/2
        
    if ( center+(width/2+slope) > lenGraph):
        print "Warning in hanningWindow: center ("+str(center)+") too large, moving to max ("+str(lenGraph-width/2)+")"
        center = lenGraph-totalWidth/2
        
   
    
    value = 0.0;
    for pt in range(center-(width/2+slope),center+(width/2+slope)):
        if ( pt < (center-(width/2)-slope) ):
            outYArray.append(0)
            outXArray.append(inXArray[pt])

        elif ( ( pt >= (center-(width/2)-slope) ) and ( pt < center-(width/2) ) ):
            pt_corr = pt - (center-(width/2.)-slope) + 1
            value = ( 0.5 - 0.5*np.cos((np.pi*(-pt_corr/slope) ) ))
            outYArray.append(value*inYArray[pt])
            outXArray.append(inXArray[pt])

        elif ( ( pt >= center-(width/2) ) and ( pt < center+(width/2) ) ):
            outYArray.append(inYArray[pt])
            outXArray.append(inXArray[pt])
            continue

        elif ( ( pt >= center+(width/2) ) and ( pt < center+(width/2)+slope ) ):
            pt_corr = (center+(width/2.)+slope) - pt
            value = ( 0.5 - 0.5*np.cos((np.pi*(pt_corr/slope)) ) )
            outYArray.append(value*inYArray[pt])
            outXArray.append(inXArray[pt])

        elif ( pt >= center+(width/2)+slope ):
            outYArray.append(0)
            outXArray.append(inXArray[pt])

          
    return np.array(outXArray),np.array(outYArray)




def nyquistLimit(inputWaveVolts,fMax):
    
    #this is a real butterworth filter
    lpFilt = 1.25 #GHz (1.25 is nyquist of LAB3B, so this is a const)
    lpFilt /= fMax #scale it by nyquist
    b,a = signal.butter(5,lpFilt,'lowpass')

    filtered = signal.lfilter(b,a,inputWaveVolts,axis=0)

        
    return filtered

def lowPassFilter(inputX,inputY,lpFilter=1.25):
    fMax = 2./(inputX[1]-inputX[0])
    lpFilter /= fMax

    b,a = signal.butter(5,lpFilter,"lowpass")

    filtered = signal.lfilter(b,a,inputY)
    
    return filtered


def calcPhase(inputFFT):

    phase = np.angle(inputFFT)
    phase = unwrapPhase(phase)

    return phase


def genPhase(inputX,inputY):
    inputF,inputFFT = genFFT(inputX,inputY)
    grpDly = calcPhase(inputFFT)

    return inputF,grpDly



def unwrapPhase(phase):

    for i in range(1,len(phase)):
        
        if (phase[i] - phase[i-1] > np.pi):
            for j in range(i,len(phase)):
                phase[j] -= 2 * np.pi
        elif (phase[i] - phase[i-1] < -np.pi):
            for j in range(i,len(phase)):
                phase[j] += 2 * np.pi

    return phase

def fftPhaseShift(fft,shift):
    gainLin,phase = complexToGainAndPhase(fft)
    phase = phase+shift

    fft = gainAndPhaseToComplex(gainLin,phase)
    return fft


def gainAndPhaseToComplex(gainLin,phase):
    real = []
    imag = []
    for i in range(0,len(gainLin)):
        real.append(gainLin[i]*np.sin(phase[1]))
        imag.append(gainLin[i]*np.cos(phase[1]))
  
    return np.array(real)+np.array(imag)*1j


def complexToGainAndPhase(fft):
    gainLin = []
    phase = []
    for i in range(0,len(fft)):
        gainLin.append(np.absolute(fft[i]))
        phase.append(np.angle(fft[i]))
  
    phase = unwrapPhase(phase)

    return np.array(gainLin),np.array(phase)


def calcGroupDelay(inputF,inputFFT):

    phase = calcPhase(inputFFT)
    
    dF = inputF[1]-inputF[0]
#    dF = 1
    GrpDly = np.abs(np.diff(phase)/(2*np.pi*dF))

    return GrpDly


def genGroupDelay(inputX,inputY):
    inputF,inputFFT = genFFT(inputX,inputY)
    grpDly = calcGroupDelay(inputF,inputFFT)

    return inputF[:-1],grpDly


def bandPassFilter(inputWaveX,inputWaveY,lpFilter=1.4,hpFilter=0.175):
    """
    l/hpFilter in GHz (1.25 is nyquist of LAB3B, so this is a const)
    fMax = 2./(inputWaveX[1]-inputWaveX[0])
    """
    fMax = 1/(2*(inputWaveX[1]-inputWaveX[0]))

    lpFilter /= fMax
    hpFilter /= fMax

    b,a = signal.butter(3,[hpFilter,lpFilter],"bandpass")

    filtered = signal.lfilter(b,a,inputWaveY,axis=0)

    return filtered


def highPass(inputWaveX,inputWaveY,hpFilter=0.150):
    """
    lpFilter in GHz (1.25 is nyquist of LAB3B, so this is a const)
    """

    fMax = 1/(2*(inputWaveX[1]-inputWaveX[0]))
    
    #this is a real butterworth filter
    hpFilter /= fMax #scale it by nyquist
    b,a = signal.butter(5,hpFilter,'highpass')

    filtered = signal.lfilter(b,a,inputWaveY,axis=0)


def lowPass(inputWaveX,inputWaveY,lpFilter=1.50):
    """
    hpFilter in GHz (1.25 is nyquist of LAB3B, so this is a const)
    """

    fMax = 1/(2*(inputWaveX[1]-inputWaveX[0]))
    
    #this is a real butterworth filter
    lpFilter /= fMax #scale it by nyquist
    print "lowpass filter value: "+str(lpFilter)
    b,a = signal.butter(5,lpFilter,'lowpass')

    filtered = signal.lfilter(b,a,inputWaveY,axis=0)

        
    return filtered


def printTimeDomainToFile(x,y,file):
    file = open(file,"w")
    
    if ( len(x) != len(y) ):
        print "Warning: can't print, different lengths (+"+str(len(x))+" != "+str(len(y))+")"
        return -1

    for pt in zip(x,y):
        file.write(str(pt[0])+" "+str(pt[1])+"\n")

    file.close()

    return



def convolve(grA,grB):
    """
    I should write my own deconvolution algorithm if it is so easy

    Takes grA, tacks on len(grB) to each end, and then sweeps grB across it,
    integrating the multiplied value all the while
    returns whatever that gives you (so it will have a length of len(A)+2*len(B)
    grA and grB don't need to have the same length
    """

    oldLength = len(grA)
    newLength = len(grA)+2*len(grB)
    #zero pad wants an x array too, don't need it though
    fart,grA = zeroPadEqual(np.arange(0,oldLength),grA,newLength)

    grOut = []

    #now we compare grA vs grB
    stop =0
    pt = 0
    while (stop != 1):
        try:
            value = np.sum(grA[pt:pt+len(grB)]*grB[::-1])
            grOut.append(value)
            pt += 1
        except:
            stop =1

    return np.array(grOut)

    

def laplace(graphX,graphY):
    """
    for mapping "real valued functions to real valued functions"
    A particular real form of the Fourier transform
    F(s) = integral(inf)(0){exp(-s*t)*f(t)dt}
    
    I want to see how this reacts to our real signals
    """

    t0 = graphX[0]
    dT = graphX[1]-graphX[0]
    dF = 1/(dT*len(graphX))
    dW = dF*2*np.pi

    outGraphX = []
    outGraphY = []

    for s in np.arange(0,len(graphX))*dW*1j:
        outGraphX.append(s/(1j*2*np.pi))
        value = 0
        for i in range(len(graphX)):
            t = graphX[i]-t0
            ft = graphY[i]
#            print str(s)+" "+str(t)+" "+str(ft)+" "+str(value)
            value += np.exp(-s*t)*ft
        outGraphY.append(value)
        

    return outGraphX,outGraphY



def correlation(graphA,graphB):
    """
    Takes in two EVENLY SAMPLED graphs, then does the cross correlation in fourier space

    NOTE: I need to roll this a bit! 
    Note to Note: Why?  to draw it away from zero?  Maybe no roll :(
    """

#    fart,graphA = zeroPadEqual(np.arange(0,len(graphA)),graphA,len(graphA)*3)
#    fart,graphB = zeroPadEqual(np.arange(0,len(graphB)),graphB,len(graphB)*3)
    

    fftA = fftw.rfft(np.array(graphA))
    fftB = fftw.rfft(np.array(graphB))

    xCorr = np.conj(fftA)*fftB

    iXCorr = fftw.irfft(xCorr)
    
    return np.roll(iXCorr,len(iXCorr)/2)
#    return iXCorr
   



def findPeakCorr(dataA,dataB,xMin=0,xMax=-1):
    """
    A wrapper for correlation and fitAndPinpoint basically

    inputs:
    dataA - the first Y data to correlate
    dataB - the second Y data to correlate against dataA
    xMin/xMax - in case you want to select a range in which the correlation is valid


    returns:
    max = fit maximum peak of the correlation (within the range)
    rSq = how good the fit was
    peak = just the raw correlation peak
    """

    xCorrY = correlation(dataA,dataB)[xMin:xMax]

    params,rSq,peak = fitAndPinpoint(xCorrY)
    max = params[1]
#    if max > len(xCorrY)/2.:
#        max -= len(xCorrY)
        
    return max,rSq,peak


def findPeakCorrs(data,scopeA=0,scopeB=1,chanA=0,chanB=0,xMin=0,xMax=-1,roll=0):
    """
    will eventually be the same as findPeakCorr, except returns it in an array..

    this is also oddly specific to the transfer function stuff I was doing so I 
    shouldn't touch it too much until I split it off
    """


    maxes = []
    rSqs = []
    peaks = []
    
    dT = data[0][scopeA][chanA][0][1] - data[0][scopeA][chanA][0][0]

    horReco = len(data[0][scopeA][chanA][0])
    if xMax==-1: 
        xMax=horReco

    for eventNum in range(0,len(data)):
        xArray = np.roll(np.arange(0,horReco)*dT,roll)
        dataA = data[eventNum][scopeA][chanA][1]
        dataB = data[eventNum][scopeB][chanB][1]
        xCorr = np.roll(tf.correlation(dataA,dataB),roll)[xMin:xMax]
        params,rSq,peak = fitAndPinpoint(xArray,xCorr)
        max = params[1] 
        rSqs.append(rSq)
        peaks.append(peak)
        if max > horReco/2.:
            max -= horReco
        
        #print str(eventNum)+" "+str(max)
        maxes.append(max)


    return np.array(maxes),np.array(rSqs),np.array(peaks)



def fitAndPinpoint(xCorrY,windowSize=2):
    """
    Takes a cross correlation graph, fits a gaussian to it, and returns the offset
    fits to (from benMathFit):
    y = D*exp(((x-B)**2)/A)+C
    
    Inputs:
    xCorrY=the correlation graph you want to find the maximum of
    window=the size of the guassian fit, defaults to +-10

    I want to do this in POINTS, not whatever the x axis is.  Other things can deal
    with the stupid X axis

    returns:
    B (so it doesn't really matter that all those other things are stupid)
    Actually lets return everything:
    params (p[1] is the peak)
    rSq of fit
    peak of correlation (to see how much improvement this makes)

    So basically this is just a wrapper for a wrapper to the fit function lol whoops

    """
    

    #find the absolute peak of the correlation
    peak = np.argmax(xCorrY)
    #remember that peak = len(xCorrY)/2 is zero!
    peakValue = peak-len(xCorrY)/2
    #print "peak="+str(peak)
    
    #We need a guess!
    guessA = 2 #sure why not, it might be like 2 points wide
    guessB = 0 #should be centered at zero (because of the windowing)
    guessC = 0 #should be zero meaned
    guessD = np.max(xCorrY) #amplitude should be the peak

    #set the window around the peak 
    #(I want there to be an odd number of points so it is even on both sides of zero)
    #You have to remember this offset!  Otherwise you end up not knowing the position.
    windowY = xCorrY[peak-windowSize:peak+windowSize+1]
    windowX = np.arange(-windowSize,windowSize+1)


    #actually do the fit
    params,rSq = mf.fitGaussian(windowX,windowY,[guessA,guessB,guessC,guessD])


    #plot it for a sanity check
#    lab.plot(xCorrY)
    upsampleX = np.arange(-windowSize*10,windowSize*10+0.01,0.01)
#    print "max="+str(params[1])+" peak="+str(peakValue)


    #Now the params you are returning are
    
    return params,rSq,peakValue

def shiftWaveform(yData,offset):
    """
    shifts the yData over by a set offset <- -|+ ->
    zero pads new points that are created by the shift and dont overlay on the new graph

    """

    length = len(yData)

    if (offset>0):
        yData = yData[:length-offset]
        for i in range(0,offset):
            yData = np.insert(yData,0,0)
        
    elif (offset<0):
        yData = yData[-offset:]
        for i in range(0,-offset):
            yData = np.append(yData,0)


#    print length," ",len(yData)
    return yData

def resampleAtOffset(yData,offset):
    """
    does a cubic (or akima spline) of the yData and returns the graph as if it were
    resampled at whatever "offset" is.

    input:
    yData - input data array
    offset - point offset of grid, should be in fractions of a point
    """


    
    origXGrid = np.arange(0,len(yData))
    yInterp = Akima1DInterpolator(origXGrid,yData)

    newXGrid = origXGrid+offset
    yNew = np.nan_to_num(yInterp(newXGrid))

    #you can't predict with a spline, so you end up losing a point
    #which point you lose is dependant of the sign of the offset!

    #so for now, lets just set the point you are "missing" to zero (we are averaging)
    #otherwise it returns "nan" which screws things up
    #python has a thing that does that!! (nan_to_num)

    return yNew



def correlateAndAverage(events,makePlots=False):
    """
    Takes a bunch of events that you've already imported a bunch of and correlates and averages them
    """

    avgEvent = [[],[]]
    currEvent = [[],[]]
    numEvents = len(events)
    firstFlag = True
    eventsCopy = copy.deepcopy(events)
    for currEvent in eventsCopy:
        if firstFlag == True:
            #print "lets DO THIS THING"
            #print "copying over first event"
            avgEvent[0] = currEvent[0][:]
            avgEvent[1] = currEvent[1][:]/numEvents
            firstFlag = False
        else:
            max,rSq,peak = findPeakCorr(currEvent[1],avgEvent[1])
            shiftedEvent = shiftWaveform(currEvent[1],peak)
            resampledEvent = resampleAtOffset(shiftedEvent,max)
            avgEvent[1] += resampledEvent/numEvents
            
            #test plotting
            if makePlots==True:
                if len(lab.get_fignums())==0:
                    fig = lab.figure()
                    axA = fig.add_subplot(211)
                    axB = fig.add_subplot(212,sharex=axA)
                    fig.show()
                else:
                    fig = lab.figure(lab.get_fignums()[0])
                    axA = fig.get_axes()[0]
                    axB = fig.get_axes()[1]
                    axB.cla()
                    axA.cla()

                axA.plot(currEvent[1],label="curr")
                axA.plot(resampledEvent,label="resamp")
                axA.legend()
                axB.plot(avgEvent[1],label="avg")
                axB.legend()
                fig.canvas.draw()



        currEvent = [[],[]]
    
    return avgEvent[0],avgEvent[1]




def CheckParseval(dT,):
    """
    The result of a fourier transform is unitary:
    "The sum (or integral) of the square of the function is equal to the sum (or integral) of the square of it's transform"
    This keeps the total energy the same in both
    """

    return

from scipy import interpolate

def interp(xNew,xOld,yOld,right=0):
    """
    remaps the scipy interpolation scheme (which has two lines)
    to the numpy interpolation syntax (which is one line!)
    "right" in this case maps to both right and left because I am lazy!
    """

    interpObj = interpolate.Akima1DInterpolator(xOld, yOld)
    
    interpWave = interpObj(xNew)

    #the akima returns nan out of it's range
    interpWave[np.isnan(interpWave)] = right
    
    return interpWave


def complexZeroPadAndResample(inputF,inputFFT,sampleLength=1024,sampleResolu=0.1):
    """
    This is the important function for matching so I'll explain it:
    So far I've done everything to 1024 samples at 0.1ns in time domain

    ->>  That translates to 513 samples (len/2+1) at 0.009765625GHz (1/Ttot)
         Fmax = 5GHz (zero pad to reach this value)
    """

    finaldF = 1./(sampleLength*sampleResolu);

    outputF = np.arange(0,(sampleLength/2.)+1)*finaldF


    #numpy 1.11 can't do complex interp, 1.12 can though! (requires python 2.7)
    #so I have to do this by hand (wonderful)
    inputReal = np.real(inputFFT)
    inputImag = np.imag(inputFFT)
    outputReal = interp(outputF,inputF,inputReal,right=1)
    outputImag = interp(outputF,inputF,inputImag,right=1)

    outputFFT = outputReal+outputImag*1j

#    print outputReal


#    fig,ax = lab.subplots()
#    ax.plot(outputF,outputReal,'-')
 #   ax.plot(inputF,inputReal,'.')
 #   fig.show()

    
    return outputF,outputFFT


def fourierZero(FFT):
    newFFT = []
    for i in range(0,len(FFT)):
        if i < 10:
            newFFT.append(0)
        elif i > 140:
            newFFT.append(0)
        else:
            newFFT.append(FFT[i])

    return np.array(newFFT)
