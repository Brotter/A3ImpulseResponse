"""
  Ben Rotter - BenJRotter@gmail.com - University of Hawaii at Manoa

  A collection of functions for doing fourier domain manipulation
  Used primarily for generating ANITA3 transfer functions
  Recreates many things in Ryan Nichol's FFTtools, plus a bunch of
  useful things that I couldn't find in there
  
  Edited by John Russell.
"""

import numpy as np
import pylab as lab
import scipy.signal as signal
from scipy import interpolate
from scipy.interpolate import interp1d, Akima1DInterpolator

import copy

import seaborn as sns

try:
    import pyfftw.interfaces.numpy_fft as fftw
except:
    print "You don't have the pyfftw package, so I'll use numpy's FFT package FFTPACK"
    import numpy.fft as fftw

import cmath #for doing complex square root

import benMathFit as mf


def genFFT(graphX,graphY):
    outFFT = fftw.rfft(graphY)
    outF = genFreqArray(graphX)
    
    return outF,outFFT


def genFreqArray(graphX):
    length = len(graphX)/2. + 1
    dF = 1./((graphX[1]-graphX[0])*len(graphX))

    fMax = length * dF

    return np.arange(0,fMax,dF)


def genTimeSeries(graphF,graphFFT):
    outY = fftw.irfft(graphFFT)
    outX = genTimeArray(graphF)

    return outX,outY


def genTimeArray(graphF):
    lengthT = (len(graphF) - 1) * 2.
    
    dT = 1./((graphF[1]-graphF[0])*lengthT)
    
    tMax = lengthT * dT

    print dT,tMax

    return np.arange(0,tMax,dT)


def genLogMag(graphX,graphY):
    graphF,graphFFT = genFFT(graphX,graphY)
    logMag = calcLogMag(graphF,graphFFT)
    return graphF,logMag


def calcLogMag(graphF,graphFFT):
    """
    Returns Power Spectral Density (dBm/Hz)
    dF is likely in GHz (because thats how I roll) so you need to scale that for this to be dBm/Hz
    This also assumes you're giving it voltage in VOLTS

    http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf
    time-integral squared amplitude
    
    I somehow messed this up?
    """

    dF = (graphF[1]-graphF[0]) * 1e9  # Resolution bandwidth (RBW), converted from GHz to Hz.
    return 10. * np.log10(np.abs(graphFFT)**2 / 5e-2 / dF)
    #unfold spectrum (This is done for you in fftw.rfft()):
    #1) discard negative frequencies (second half)
    #2) double all positive frequencies (but not DC (bin 0) and nyquist (bin length/2+1)
    #3) now you can do the log mag

    #do the calculation
    #Power = V**2 / R
    #Watts = V**2 / ohms
    #Transfers in fourier space so don't worry
    #dBm = 10.*log10(Power(in mW)/1mW)
    #"spectral power" -> divide it by frequency
#    dBm = 10.*np.log10(power)

#    return dBm
    

def calcLinMag(graphF,graphSpecMag):
    """
    It might be nice to have a way to reverse this (and return power in watts)
    """
    dF = graphF[1]-graphF[0]

    graphLogMag = graphSpecMag*dF
    power = 10**(graphLogMag / 10.)
    
    return power


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
    Be aware that there are many values you can "zero" it with, but, if you're
    constantly zero meaning, zero might be best
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


def hanningTail(inYArray,start,slope):
    """
    If you roll the waveform to the front of the window, you only want to hanning filter the end of it, so this does that.

    Comparable to doing a hanning window with a step function on the left side
    """
    
    outYArray = []
    for pt in range(0,len(inYArray)):
        if pt<=start:
            outYArray.append(inYArray[pt])
        if pt>start and pt<=start+slope: 
            pt_corr = pt - start
            value = 0.5 + 0.5*np.cos((np.pi*(-pt_corr/float(slope))))
            outYArray.append(inYArray[pt]*value)
        if pt>(start+slope):
            outYArray.append(0)

    return np.array(outYArray)


def hanningWindow(inXArray,inYArray,center,totalWidth=1000,slope=200):
    """
    Takes two input arrays (x and y, since they both will change length), a
    total width, and a center (and optional slope). Returns two arrays, a
    hanning windowed Y axis with a width of totalWidth centered at center,
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

def fftPhaseShift(fft,shift,gdShift=0):
    gainLin,phase = complexToGainAndPhase(fft)
    phase = phase+shift+np.arange(len(fft))*gdShift

    fft = gainAndPhaseToComplex(gainLin,phase)
    return fft


def timePhaseShift(Y,shift):
    fft = fftw.rfft(Y)
    gainLin,phase = complexToGainAndPhase(fft)
    phase += shift
#    phase[phase>=np.pi] -= np.pi
#    phase[phase<=-np.pi] += np.pi
    fft = gainAndPhaseToComplex(gainLin,phase)
    outY = fftw.irfft(fft)
    return outY


def accumulatePhase(phase):
    """
    If you're going to resample the phase, you need it to take the wraps into
    account! The easiest way to do this is to "accumulate" the phase.
    Basically, when it wraps upwards add +2pi to all all subsequent values,
    and vice versa

    This takes in the phase array and returns that accumulated array
    """

    phaseOut = copy.deepcopy(phase)

    for pt in range(1,len(phase)):
        if ( (phase[pt-1] - phase[pt]) > np.pi ):
            for pt2 in range(pt,len(phase)):
                phaseOut[pt2] += np.pi*2
        elif ( (phase[pt] - phase[pt-1]) > np.pi ):
            for pt2 in range(pt,len(phase)):
                phaseOut[pt2] -= np.pi*2

    return phaseOut


def regenerateCablePhase(phase,dF=5000./512.,length=513):
    """
    A cable has a very flat group delay (Tg(w)), so we should, instead of trying to resample it,
    just regenerate the phase entirely of a cable using the measured phase
    
    Also, since I'm doing this for the transfer function most likely, I'm going to set the defaults
    that everything else has:
    time domain: 10GS/S and n=1024
    freq domain: dF = 5000/512 and length=513
    """

    aPhase = accumulatePhase(phase)
    groupDelay = np.diff(aPhase)
    #only use the first half for determining the mean, since it gets noisier at the ends 
    groupDelayMean = np.mean(groupDelay[:len(groupDelay)/2])
    
    #now, the DC is going to be zero (just because) and everything else is going to DECREASE from there
    phaseOut = np.arange(length)*(-groupDelayMean)*dF

    return phaseOut


def regenerateCableLinMag(f,linMag,dF=5000./512.,length=513):
    """
    Just like the phase, I need to be able to regenerate the cable linear magnitude.

    This one gets the frequency too because I don't know where the boundries of the measurement are
    offhand and want to be able to determine them
    """

    linMagInterp = Akima1DInterpolator(f*1000.,linMag)

    outF = np.arange(length)*dF
    newLinMag = linMagInterp(outF)

    #DC should be zero
    newLinMag[0] = 0

    #out of range it should be zero (nan_to_num does that hooray!)
    newLinMag = np.nan_to_num(newLinMag)

    return newLinMag


def complexToGainAndPhase(fft):
    """
    the gain and phase are polar representations of the real/imag values on the
    complex plane

    np.angle = "the counterclockwise angle from the positive real axis on the
    complex plane"

    That seems like a fine convention to stick to, I'm sure other people use it
    too
    """

    gainLin = []
    phase = []
    for i in range(0,len(fft)):
        gainLin.append(np.absolute(fft[i]))
        phase.append(np.angle(fft[i]))
  
    #can't do this!  It distorts which phasor quadrant you are in
#    phase = unwrapPhase(phase)

    return np.array(gainLin),np.array(phase)


def gainAndPhaseToComplex(gainLin,phase):
    """
    Just need to follow the conventions from complexToGainAndPhase() which is
    defined by np.angle
    """

    real = []
    imag = []
    for i in range(0,len(gainLin)):
        imag.append(gainLin[i]*np.sin(phase[i]))
        real.append(gainLin[i]*np.cos(phase[i]))
  
    return np.array(real)+np.array(imag)*1j


def sqrtOfFFT2(fft):
    """
    This is the actual way to do it!  The phase is divided in two and the
    magnitude is rooted.  
    I THINK this preserves the signal (though it does time invert the causality
    so you have to reverse the inverse fft)
    """

    gain,phase = complexToGainAndPhase(fft)
    gain = np.sqrt(gain)
    gphase = phase/2.
    fftOut = gainAndPhaseToComplex(gain,phase)

    return fftOut


def sqrtOfFFT1(fft):
    """
    The dumbest and first approach at taking the square root of a complex
    transfer function
    Since sqrt isn't well defined for a complex number, this doesn't preserve
    the time domain signal
    """
    out = []
    for i in range(0,len(fft)):
        out.append(cmath.sqrt(fft[i]))

    return np.array(out)


def calcGroupDelay(inputFFT,inputF=-1,dF=-1):

    phase = calcPhase(inputFFT)
    
    if inputF != -1 and dF != -1:
        print "WARNING IN tfUtils::calcGroupDelay - I'm going to use inputF to generate dF even though you also gave me dF"

    if  inputF != -1:
        dF = inputF[1]-inputF[0]
    if dF == -1:
        dF = 0.1 #sure why not

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


def lowPassFilter(inputX,inputY,lpFilter=1.25):
    fMax = 2./(inputX[1]-inputX[0])
    lpFilter /= fMax

    b,a = signal.butter(5,lpFilter,"lowpass")

    filtered = signal.lfilter(b,a,inputY)
    
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
    fart,grA = zeroPadEqual(np.arange(oldLength), grA, newLength)

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

    for s in np.arange(len(graphX))*dW*1j:
        outGraphX.append(s/(1j*2*np.pi))
        value = 0
        for i in range(len(graphX)):
            t = graphX[i]-t0
            ft = graphY[i]
#            print str(s)+" "+str(t)+" "+str(ft)+" "+str(value)
            value += np.exp(-s*t)*ft
        outGraphY.append(value)
        
    return outGraphX, outGraphY


def correlation(graphA,graphB):
    """
    Takes in two EVENLY SAMPLED graphs, then does the cross correlation in
    fourier space

    NOTE: I need to roll this a bit! 
    Note to Note: Why?  to draw it away from zero?  Maybe no roll :(
    """

#    fart,graphA = zeroPadEqual(np.arange(len(graphA)),graphA,len(graphA)*3)
#    fart,graphB = zeroPadEqual(np.arange(len(graphB)),graphB,len(graphB)*3)

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
    xMin/xMax - in case you want to select a range in which the correlation is
    valid

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
        xArray = np.roll(np.arange(horReco)*dT,roll)
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
    Takes a cross correlation graph, fits a gaussian to it, and returns the
    offset fits to (from benMathFit):
    y = D*exp(((x-B)**2)/A)+C
    
    Inputs:
    xCorrY=the correlation graph you want to find the maximum of
    window=the size of the guassian fit, defaults to +-10

    I want to do this in POINTS, not whatever the x axis is.  Other things can
    deal with the stupid X axis

    returns:
    B (so it doesn't really matter that all those other things are stupid)
    Actually lets return everything:
    params (p[1] is the peak)
    rSq of fit
    peak of correlation (to see how much improvement this makes)

    So basically this is just a wrapper for a wrapper to the fit function lol
    whoops
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
    upsampleX = np.arange(-windowSize * 10, windowSize * 10 + 0.01, 0.01)
#    print "max="+str(params[1])+" peak="+str(peakValue)

    #Now the params you are returning are
    
    return params,rSq,peakValue


def shiftWaveform(yData,offset):
    """
    shifts the yData over by a set offset <- -|+ ->
    zero pads new points that are created by the shift and dont overlay on the
    new graph
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
    does a cubic (or akima spline) of the yData and returns the graph as if it
    were resampled at whatever "offset" is.

    input:
    yData - input data array
    offset - point offset of grid, should be in fractions of a point
    """
    
    origXGrid = np.arange(len(yData))
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
    Takes a bunch of events that you've already imported a bunch of and
    correlates and averages them
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
    "The sum (or integral) of the square of the function is equal to the sum
    (or integral) of the square of it's transform"
    This keeps the total energy the same in both
    """

    return


def interp(xNew,xOld,yOld,right=0):
    """
    remaps the scipy interpolation scheme (which has two lines)
    to the numpy interpolation syntax (which is one line!)
    "right" in this case maps to both right and left because I am lazy!
    """

    interpObj = Akima1DInterpolator(xOld, yOld)
    
    interpWave = interpObj(xNew)

    #the akima returns nan out of it's range
    interpWave[np.isnan(interpWave)] = right
    
    return interpWave


def fitAndRegenFFT(inputF,inputFFT,sampleLength=1024,sampleResolu=0.1):
    """
    So far I can't find a good way to resample the cables, so I'm going to fit
    the mangitude and then regenerate a causal waveform from that

    This ends up doing weird things, so lets fit the last 1/4th, then hanning
    window the end
    """

    length = len(inputF)
    dF = inputF[1]-inputF[0]

    extrapF = np.arange(inputF[-1]+dF,5,dF)
    
    inputF   = np.concatenate((inputF,extrapF))
    inputFFT = np.concatenate((inputFFT,np.zeros(len(extrapF))))

    magnitude = np.absolute(inputFFT)

    pFit,rSq = mf.fitPoly3(inputF,magnitude,[0,0,0,0])
    
    timeSeries = np.arange(sampleLength)*sampleResolu
    freqSeries = genFreqArray(timeSeries)

    newMag = mf.lambdaPoly3(pFit,freqSeries)

    newMag = hanningTail(newMag,len(newMag)/2,20)

    newMag[newMag<0] = 0

    fig,ax = lab.subplots()
    ax.plot(inputF,magnitude,'.')
    ax.plot(freqSeries,newMag)
    fig.show()

    newReal = np.sqrt(newMag)
    newReal[::2] *= -1
    newImag = np.sqrt(newMag)
    newReal[1::2] *= -1

    return freqSeries,newReal+newImag*1j

    
def complexZeroPadAndResample(inputF,inputFFT,sampleLength=1024,sampleResolu=0.1):
    """
    This is an important function for matching so I'll explain it:
    So far I've done everything to 1024 samples at 0.1ns in time domain

    ->>  That translates to 513 samples (len/2+1) at 0.009765625GHz (1/Ttot)
         Fmax = 5GHz (zero pad to reach this value)

    This BREAKS THINGS!! :( It makes a second weird pulse occur?  Maybe I
    should just do it in time domain...
    It is because setting values out of the range of the input pulse to 1 is
    adding power in weird places
    So I changed the "right=1" to "right=0" and that fixes things somewhat
    """

    finaldF = 1./(sampleLength*sampleResolu);

    outputF = np.arange(sampleLength/2.+1)*finaldF

    #numpy 1.11 can't do complex interp, 1.12 can though! (requires python 2.7)
    #so I have to do this by hand (wonderful)
    inputReal = np.real(inputFFT)
    inputImag = np.imag(inputFFT)
    outputReal = interp(outputF,inputF,inputReal,right=0)
    outputImag = interp(outputF,inputF,inputImag,right=0)

    outputFFT = outputReal+outputImag*1j

#    print outputReal

#    fig,ax = lab.subplots()
#    ax.plot(outputF,outputReal,'-',label="output")
#    ax.plot(inputF,inputReal,'.',label="input")
#    ax.legend()
#    fig.show()
    
    return outputF,outputFFT


def fourierZeroPad(f,fft,numPads):
    #this doesn't work and does bad things
    print "fourierPad doesn't work and does bad things"
    fftOut = np.concatenate((fft,np.zeros(numPads)))
    fOut = np.concatenate((f,np.arange(numPads)*(f[1]-f[0])+f[-1]))

    return fOut,fftOut


def fourierPad(f,fft,numPads,value):
    #this doesn't work and does bad things
    print "fourierPad doesn't work and does bad things"

#    realPad = np.arange(numPads,0,-1)*(value/numPads)
    realPad = np.ones(numPads)*value
#    realPad[::2] *= -1
    imagPad = np.imag(signal.hilbert(realPad))

    real = np.concatenate((np.real(fft),realPad))
    imag = np.concatenate((np.imag(fft),imagPad))

#    fftOut = np.concatenate((fft,(np.ones(numPads)-np.ones(numPads)*1j)*value))
    fOut = np.concatenate((f,np.arange(numPads)*(f[1]-f[0])+f[-1]))

    return fOut,real+imag*1j


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


def fourierExtrapolate(fftIn,numPads,fitStart=-1):
    if fitStart==-1:
        fft = fftIn
    else:
        fft = fftIn[:fitStart]
    p = mf.fitPoly3(np.arange(len(fft)),np.absolute(fft),[0,0,0,0])
    print p
    extrap = mf.lambdaPoly3(p[0],np.arange(len(fftIn)+numPads))
    extrap = hanningTail(extrap,len(fft),len(fft)+50)

    #I put work into this so I want to keep it but it doesn't work!
#    decay = mf.lambdaDecay([1.,100.],np.arange(numPads,dtype=float))
#    decay = np.concatenate((np.ones(len(fft)),decay))
#    print len(fft),len(decay),len(extrap)
#    extrap *= decay

    fig,ax = lab.subplots()
    ax.plot(np.absolute(fft))
    ax.plot(extrap)
    fig.show()

    return extrap


#==============================================================================
# Using phase shifting to make a causal signal
#==============================================================================

def makeCausalFFT(fft,tZero):
    """
    http://cdn.teledynelecroy.com/files/whitepapers/14-wp6paper_doshi.pdf
    
    Basically when taking the S21 of a cable or signal (which we are doing,
    just a forward gain complex phasor,), Fractional sample offsets distort the
    signal from causal to an acausal sinc function in the time domain.  To
    alieviate this, it is allowed to shift the total phase of the S21 in order
    to achieve maximum causality, which in this case would be a high pre-T0 to
    post-T0 power ratio.

    I wrote a section about this in the impulse response white paper citing the
    above linked paper

    I think the max causality has to be in between -pi and pi... though
    probably smaller
    HOW IS THE CENTER DEFINED???? (17 is just for the cables which I found by
    looking, but it has to be calculatable somehow).  You have to provide this
    because I'm strapped for time

    Also at some point I should fit the "causalityRatio" thing and actually
    find the perfect point, but that is polishing the cannon for sure.

    Also note that this can 180 degree phase shift things (flip the polarity)
    so BE CAREFUL
    """

    shifts = np.arange(-np.pi,np.pi,0.01)
    
    causalityRatio = []

    for i in range(0,len(shifts)):
        shiftedFFT = fftPhaseShift(fft,shifts[i])
        shifted = fftw.irfft(shiftedFFT)
        causalityRatio.append(np.sum(shifted[:tZero]**2)/np.sum(shifted[tZero:]**2))

    #minimum is when things are "most causal" (second half biggest than first half)
    maxCausal = np.argmin(causalityRatio)
    
    print "Maximum causal waveform found at: ",shifts[maxCausal]
    shiftedFFT = fftPhaseShift(fft,shifts[maxCausal])

    return shiftedFFT


def makeCausalTime(y,tZero):
    """
    If you have the time series, this makes it easier (but slightly slower!
    Two more FFTS!)
    
    you have to provide the pulse (in bins)

    """

    fft = fftw.rfft(y)
    shifted = makeCausalFFT(fft,tZero)
    yOut = fftw.irfft(shifted)

    return yOut


def compPhaseShiftSimple(cableName):
    shifts = np.arange(-np.pi*2,np.pi*2,0.1)
    shiftedArrayMaxes = []
    for shift in shifts:
        a,b = getCables(cableName,phaseShift=shift)
        fft = tf.fftw.irfft(b)
        shiftedArrayMaxes.append(np.max(fft))

    fig,ax = lab.subplots()
    ax.plot(shifts,shiftedArrayMaxes)
    fig.show()
    
    return shiftedArrayMaxes


def compPhaseShifts(cableName):
    shifts = np.arange(0,np.pi*2,0.001)
    shiftedArrayMaxes = []
    for shift in shifts:
        a,b = getCables(cableName,phaseShift=shift)
        fft = tf.fftw.irfft(b)
        shiftedArrayMaxes.append(np.max(fft))

    fig,ax = lab.subplots()
    ax.plot(shifts,shiftedArrayMaxes)
    fig.show()
    
    return shifts[np.argmax(shiftedArrayMaxes)]
        

def compPhaseShifts2(cableName="A-C_PULSER-TEST_66DB.s2p"):
    max = compPhaseShifts(cableName)
    print max
    shiftedArrayMaxes = []
    shifts = np.arange(max-1,max+1,0.25)
    fig,ax = lab.subplots()

    for shiftNum in range(0,len(shifts)):
        a,b = getCables(cableName,phaseShift=shifts[shiftNum])
        fft = tf.fftw.irfft(b)
        ax.plot(np.arange(len(fft))+shiftNum*len(fft),np.roll(fft,400))

    fig.show()
    return


def compPhaseShifts3(y=[],center=[],save=False):

    if y == []:
        y = np.zeros(1024)
        y[50] = 1
        center = 50

    if center == []:
        center = len(y)/2

    fft = fftw.rfft(y)

    maxes  = []
    mins   = []
    maxLoc = []
    minLoc = []
    causalityRatio = []
    shifts = np.arange(-2*np.pi,2*np.pi,0.1)
#    shifts = np.arange(-np.pi/4.,np.pi/4.,0.01)
    fig,ax = lab.subplots(3,1,figsize=(15,7))
    fig.show()
    for i in range(0,len(shifts)):
        shiftedFFT = fftPhaseShift(fft,shifts[i])
        shifted = fftw.irfft(shiftedFFT)

        #time domain pulse
        ax[0].cla()
        ax[0].set_title("time domain");
        ax[0].plot(shifted,lw = 3)
        ax[0].plot(range(0,center),shifted[:center],'.',color="red")
        ax[0].plot(range(center,len(shifted)),shifted[center:],'.',color="blue")
#        ax[0].set_xlim([0,200])
        ax[0].set_ylim([-0.5,0.5])

        #phase
        ax[1].cla()
        ax[1].set_title("Phase")
        ax[1].plot(np.angle(shiftedFFT[1:-1]))
        ax[1].set_ylim([-np.pi,np.pi])

        #group delay
        ax[2].cla()
        ax[2].set_title("Group Delay")
        ax[2].plot(calcGroupDelay(shiftedFFT))

#        ax[1][0].plot(10.*np.log10((np.abs(shiftedFFT[1:-1])**2)))

        #fourier components
#        ax[2][0].cla()
#        ax[2][0].set_title("Group Delay")
#        ax[2][0].plot(np.real(shiftedFFT),label="real")
#        ax[2][0].plot(np.imag(shiftedFFT),label="imag")
#        ax[2][0].legend()

        #time domain "zeropoint" peak ratio
#        maxes.append(np.max(shifted))
#        mins.append(np.min(shifted))
#        maxLoc.append(np.argmax(shifted))
#        minLoc.append(np.argmin(shifted))
#        causalityRatio.append(np.sum(shifted[:center]**2)/np.sum(shifted[center:]**2))
#        ax[0][1].set_title("Time Domain Zeropoint Peak Ratio");
#        ax[0][1].plot(shifts[i],maxes[-1],'.',color="blue")
#        ax[0][1].plot(shifts[i],mins[-1],'.',color="red")
#
#        #location of peak
#        ax[1][1].set_title("peak location")
#        ax[1][1].plot(shifts[i],maxLoc[-1],'.',color="blue")
#        ax[1][1].plot(shifts[i],minLoc[-1],'.',color="red")
#
#        #ratio of pre-peak to post-peak power
#        ax[2][1].set_title("Ratio of pre-peak to post-peak power");
#        ax[2][1].plot(shifts[i],causalityRatio[-1],'.',color="purple")

        #draw it
        fig.canvas.draw()
        #save it for moviemakin'
        if save==True:
            fig.savefig("phaseShiftMovie/"+str(i).zfill(3)+".png")

#    maxCausal = np.argmax(causalityRatio)
    
#    shiftedFFT = fftPhaseShift(fft,shifts[maxCausal])
#    shifted = fftw.irfft(shiftedFFT)

    return causalityRatio, shifted


def compPhaseShifts4(y):

    sns.set_palette(sns.color_palette("husl",20))

    fig,ax = lab.subplots()
    
    fft = fftw.rfft(y)
    shifts = np.arange(-np.pi,np.pi,0.01)

    for tZero in range(310,330):
        print tZero
        causalityRatio = []

        for i in range(0,len(shifts)):
            shiftedFFT = fftPhaseShift(fft,shifts[i])
            shifted = fftw.irfft(shiftedFFT)
            causalityRatio.append(np.sum(shifted[:tZero]**2)/np.sum(shifted[tZero:]**2))

        ax.plot(shifts,causalityRatio,label=tZero)

    ax.legend()
    fig.show()

    return