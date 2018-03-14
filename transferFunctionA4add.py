"""
Script written by John Russell which includes functions specific to ANITA-4 over
ANITA-3. Takes some functions from Ben Rotter's "transferFunctionFinal.py" and
modifies them for ANITA-4 usage.
"""
import os as os # import this to use environmental variables for portability (define $PALESTINE_DATA where all the CSBF16 stuff is)
import numpy as np
import pylab as lab
import tfUtils as tfu
from scipy import interpolate
import transferFunctionFinal as tff  #  This is to reference what Ben has already done.


#==============================================================================
# Directories which are potentially referred to in this script.
#==============================================================================


# Directory where data saved in GitHub reposity for A3ImpulseResponse.
A4Dir = "./"

# Directories with antenna pulse data.
localDir = "../calibrationLinks/"

# Signal chain data (57dB seems to be a happy middle power).
waveformDir = localDir + "waveforms/"


#==============================================================================
# Importing things
#==============================================================================


def findPalestineAntennaFile(chan, inOrOut):
    antennaInfo = np.loadtxt(A4Dir + "ANITA4-Antenna-assignments.csv",
                             dtype = bytes, delimiter = ",").astype(str)

    antennaSN = str(antennaInfo.T[3][np.argwhere(antennaInfo.T[2]==chan[:3])[0][0]])
    pol = chan[3]

    antdir = '{0}'.format(os.environ['PALESTINE_DATA'])

    if inOrOut == "S21":
        if pol == "H":
            fileName = antdir + "/BS/CSBF16_SN" + antennaSN + ".H.C.csv"
        if pol == "V":
            fileName = antdir + "/BS/CSBF16_SN" + antennaSN + ".V.C.csv"

    if inOrOut == "in":
        fileName = antdir + "/OFF-BS-PLSR/CALIBRATION/06_22-INPUT_PULSE_20dB_DOWN-waveform.csv"

    if inOrOut == "out":
        if (antennaSN == "216501" or antennaSN == "216507") and pol == "H":
            fileName = antdir + "/OFF-BS-PLSR/SN" + antennaSN + "/H-TRANS/waveform/06_28-el_0-az_0-" + pol + "-C-waveform.csv"
        if (antennaSN == "208552" or antennaSN == "216513") and pol == "H":
            fileName = antdir + "/OFF-BS-PLSR/SN" + antennaSN + "/H-TRANS/waveform/06_29-el_0-az_0-" + pol + "-C-waveform.csv"
        if antennaSN == "216507" and pol == "V":
            fileName = antdir + "/OFF-BS-PLSR/SN" + antennaSN + "/V-TRANS/waveform/06_23-el_0-az_0-" + pol + "-C-waveform.csv"
        if antennaSN != "216507" and pol == "V":
            fileName = antdir + "/OFF-BS-PLSR/SN" + antennaSN + "/V-TRANS/waveform/06_30-el_0-az_0-" + pol + "-C-waveform.csv"
    return fileName


def importPalAntIn(chan):
    """
      Palestine Antennas
        there is only one waveform :( 
        so I can't do any averaging and I just import the raw data I guess
        Also this means I probably need to link the whole dir (since they are far away)
    """

    fileName = findPalestineAntennaFile(chan,"in")

    dataX,dataY = np.loadtxt(fileName, delimiter = ",", skiprows = 6, usecols=(3,4)).T
    dataX *= 1e9
    #make the time range start at zero
    dataX -= dataX[0]

    dataY = np.roll(dataY,30-np.argmax(dataY))
    #also I always want to have the fourier spaces of these things
    dataF, dataFFT = tfu.genFFT(dataX,dataY)

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
    dataF,dataFFT = tfu.genFFT(dataX,dataY)

    return dataX, dataY, dataF, dataFFT

def importPalAntS21(chan):
    """
      Palestine Antennas on BS data
    """

    fileName = findPalestineAntennaFile(chan,"S21")

    freq = []
    mag = []
    n = 0
    with open(fileName) as inFile:
        for line in inFile:
            if(line[0] == 'E' and n > 0):
                break
            if(line[0] != '!' and line[0] != 'B'):
                split = line.split(',')
                freq.append(float(split[0])/1e9)
                mag.append(float(split[1].strip()))
                n+=1

    return np.array(freq), np.array(mag)

def importAverageBRotter():
    """
    Ben Rotter's average antenna
    """

    fname = '{0}/share/AnitaAnalysisFramework/responses/SingleBRotter/all.imp'.format(os.environ['ANITA_UTIL_INSTALL_DIR'])

    t,v = np.loadtxt(fname).T
    x = np.array(t)
    y = np.array(v)
    f,fft = tfu.genFFT(x,y)
    return x,y,f,fft 

#==============================================================================
# Processing.
#==============================================================================

"""
Used by makePalAntAverageTransferFunction, unimportant
"""
def addToArray(array, fname):
    tempIn = np.genfromtxt(fname, delimiter=" ")
    array += tempIn.T[1]
    return array

"""
Used by makePalAntAverageWaveform, unimportant
"""
def addToArray2(array, chan):
    inX,inY=importPalAntS21(chan)
    array += inY
    return array

"""
  Develop an average transfer function from the transfer functions of antennas measured in palestine 2016.
  Normalized to V/m (saved in V - ns fomat)
"""
def makePalAntAverageTransferFunction(pol, showPlots=False, savePlots=False, writeFiles=False):
    
    tempIn = np.genfromtxt("transferFunctions/palAntA4_02T" + pol + ".txt", delimiter=" ")
    inX = tempIn.T[0]
    inY = tempIn.T[1]

    inY = addToArray(inY, "transferFunctions/palAntA4_07T" + pol + ".txt")
    inY = addToArray(inY, "transferFunctions/palAntA4_12M" + pol + ".txt")
    inY = addToArray(inY, "transferFunctions/palAntA4_13M" + pol + ".txt")
    inY /= 4
    inY -= np.mean(inY)
    
    inF,inFFT = tfu.genFFT(inX,inY)
    
    if savePlots or showPlots:
        fig0,ax0 = lab.subplots(3,1)

        ax0[0].plot(inX,inY,label="average wf",color="red")
        ax0[0].set_xlim([0,25])
        ax0[0].set_xlabel("Time (ns)")
        ax0[0].set_ylabel("Voltage (V)")
        ax0[0].legend(loc="lower right")

        ax0[1].plot(inF,tfu.calcSpecMag(inF,inFFT),label="average transfer fn",color="red")
        ax0[1].set_ylabel("Spectral Magnitude \n [dBm/Hz]")
        ax0[1].legend()

        ax0[2].plot(inF,tfu.calcGroupDelay(inFFT,inputF=inF),label="average transfer fn",color="blue")
        ax0[2].set_ylabel("Group Delay (ns)")
        ax0[2].legend()

        if showPlots:
            fig0.show()
            fig0.tight_layout()
            fig0.canvas.draw()
        if savePlots:
            fig0.savefig("plots/doPalAnt_averageTF.png")

    if writeFiles:
        tff.writeOne([inX,inY], pol, "palAntA4Average")

    return inX, inY

"""
  Create an average waveform function from the waveforms of antennas measured in palestine 2016.
"""
def makePalAntAverageS21(pol, writeFiles=False):
    
    inX,inY=importPalAntS21("02T" + pol)

    inY = addToArray2(inY, "07T" + pol)
    inY = addToArray2(inY, "12M" + pol)
    inY = addToArray2(inY, "13M" + pol)
    inY /= 4
    
    if writeFiles:
        tff.writeOne([inX,inY],pol,"palAntA4AverageS21")

    return inX, inY


"""
  Develop a transfer function for an antenna measured in palestine 2016, mapped to channel.
  Normalized to V/m (saved in V - ns fomat)
"""
def doPalAnt(chan,showPlots=False,savePlots=False,writeFiles=False):

    if showPlots or savePlots:
        lab.close("all")

    #input pulse
    inRawX,inRawY,inRawF,inRawFFT = importPalAntIn(chan)

    inX,inY = tff.processWaveform(inRawX,inRawY,"palin")
    inY = np.roll(inY,30-np.argmax(inY))
    inX = inX[:1024]
    inY = inY[:1024]

    inF,inFFT = tfu.genFFT(inX,inY)
    # Increase by 20dB (from coupler)
    inFFT *= 10
#    #increase by 10dB (from coupler)
##    inFFT *= 10**(10./20)
    inX,inY = tfu.genTimeSeries(inF,inFFT)

#    inFFT = tfu.minimizeGroupDelayFromFFT(inF,inFFT)

    #output pulse
    outRawX,outRawY,outRawF,outRawFFT = importPalAntOut(chan)
    outX,outY = tff.processWaveform(outRawX,outRawY,"palout")
    outY = np.roll(outY,50-np.argmax(outY))
#    outY = tfu.hanningTail(outY,130,30)
    outX = outX[:1024]
    outY = outY[:1024]
#    outY = tfu.highPass(outX,outY)
  
    #Some are upside down.  this attempts to constrain polarity
    if np.argmin(outY) < np.argmax(outY):
        print ("Flipped!")
        outY *= -1
        outY = np.roll(outY,50-np.argmax(outY))

    outF,outFFT = tfu.genFFT(outX,outY)

#    outFFT = tfu.minimizeGroupDelayFromFFT(outF,outFFT)

    #cables (i cant find the cables files, not sure they exist at all)
#    cablesF,cablesFFT = palAntCables(showPlots=showPlots)
#    outFFT /= cablesFFT
#    outX,outY = tfu.genTimeSeries(outF, outFFT)

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

        ax0[1][0].plot(inRawF,tfu.calcSpecMag(inRawF,inRawFFT),label="Raw Input Pulse",color="red")
        ax0[1][0].plot(inF,tfu.calcSpecMag(inF,inFFT),label="Processed Pulse",color="blue")
        ax0[1][0].set_ylabel("Spectral Magnitude \n [dBm/Hz]")
        ax0[1][0].legend()

        ax0[2][0].plot(inRawF,tfu.calcGroupDelay(inRawFFT,inputF=inRawF),label="Raw Input Pulse",color="red")
        ax0[2][0].plot(inF,tfu.calcGroupDelay(inFFT,inputF=inF),label="Processed Pulse",color="blue")
        ax0[2][0].set_ylabel("Group Delay (ns)")
        ax0[2][0].legend()

        ax0[1][1].plot(outRawF,tfu.calcSpecMag(outRawF,outRawFFT),label="Raw Output Pulse",color="red")
        ax0[1][1].plot(outF,tfu.calcSpecMag(outF,outFFT),label="Processed Pulse",color="blue")
        ax0[1][1].set_ylabel("Gain (dB)")
        ax0[1][1].set_xlabel("Frequency (GHz)")
        ax0[1][1].legend()

        ax0[2][1].plot(outRawF,tfu.calcGroupDelay(outRawFFT,inputF=outRawF),label="Raw Output Pulse",color="red")
        ax0[2][1].plot(outF,tfu.calcGroupDelay(outFFT,inputF=outF),label="Processed Pulse",color="blue")
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

    inF[0] = 1e-5 #to avoid division by zero errors
    antTFFFT = np.divide(outFFT,inFFT)*(0.3/(1j*inF))
    antTFFFT[0] = 0
    inF[0] = 0

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
    # According to notes, ~5.203m face to face flight distance
    # Antennas have a depth of 22" = 0.5588m, so r(f) should have two points, (0.180,5.203) and (1.2,5.203 + 2*0.5588)
    # makes it in units of "meters squared"
    #I'm not convinced ben did this quite right so i changed it
    distAntToAnt = 5.203
    antDepth = 0.5588
    maxFlightDist = distAntToAnt + 2 * antDepth
    maxAntFreq = 1.3
    minAntFreq = .18
    dist = 2*antDepth*(1 - (maxAntFreq - antTFF)/(maxAntFreq - minAntFreq)) + distAntToAnt
    for i in range(0,len(dist)):
        if dist[i] > maxFlightDist:
            dist[i] = maxFlightDist
    
    antTFFFT *= dist

    print ("length antTFFFT:",len(antTFFFT))

    if showPlots and 1:
        figDist,axDist = lab.subplots()
        axDist.plot(antTFF,dist)
        axDist.set_xlabel("Frequency (GHz)")
        axDist.set_ylabel("Distance (m)")
        figDist.show()
    
    #take the square root to get the normalized complex height
    # units of "meters again!"
    antTFFFT = tfu.sqrtOfFFT2(antTFF,antTFFFT)

    #causality conditions aren't good for measured signals (artifically distort phase)!
#    antTFFFT = tfu.makeCausalFFT(antTFFFT,np.argmax(tfu.fftw.irfft(antTFFFT)))

    #I have to get the time domain again after that
    # units of "meters / nanosecond"
    antTFX,antTFY = tfu.genTimeSeries(antTFF,antTFFFT)

    #remove all the stupid power above 1.3GHz since ANITA can't measure it (probably should be lower)
    #butterworth filter @ 1.3
#    antTFY = tfu.highPass(antTFX,antTFY)
    antTFY = tfu.nyquistLimit(antTFY,5)
   
    #make sure xfer function is also all the same polarity
    if np.max(antTFY) < abs(np.min(antTFY)):
        print ("Flipped!")
        antTFY *= -1

    antTFY = tfu.hanningTail(antTFY, np.argmax(antTFY)+40, 100)
    antTFY = tfu.hanningTailNeg(antTFY, np.argmax(antTFY)-20, 20)
    antTFY -= np.mean(antTFY)
    
    antTFF,antTFFFT = tfu.genFFT(antTFX,antTFY)

    #Produce a plot for the gain normalized to dBi (hence the 4pi)
    if showPlots or savePlots:
        gainTF = ((4*np.pi*antTFF**2)/(0.3**2))*np.abs(antTFFFT)**2
        gainTF[0] = 1e-15 #to avoid division by zero error
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
#    antTFY = tfu.hanningTail(antTFY,20,20)
#    antTFFFT = tfu.fftw.rfft(antTFY)

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

        ax[1][0].plot(inF,tfu.calcSpecMag(inF,inFFT),label="input Pulse",color="red")
        ax[1][1].plot(outF,tfu.calcSpecMag(outF,outFFT),label="output Pulse",color="blue")
        ax[1][2].plot(antTFF,tfu.calcLogMag(antTFF,antTFFFT),label="transfer Function",color="black")
        ax[1][0].legend()
        ax[1][1].legend()
        ax[1][2].legend()
       
        ax[2][0].plot(inF,tfu.calcGroupDelay(inFFT,inputF=inF),label="input Pulse",color="red")
        ax[2][1].plot(outF,tfu.calcGroupDelay(outFFT,inputF=inF),label="output Pulse",color="blue")
        ax[2][2].plot(antTFF,tfu.calcGroupDelay(antTFFFT,inputF=antTFF),label="transfer Function",color="black")
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
        tff.writeOne([antTFX,antTFY],chan,"palAntA4")
        
    return antTFX,antTFY,antTFF,antTFFFT

def palAntCables(showPlots=False):
    """
        

    """

    #input pulse (with coupler like always)
    inputX,inputY = tff.palAntImporter('OFF-BS-PLSR/CALIBRATION/06_22-INPUT_PULSE_20dB_DOWN-waveform.csv')
    #inputY *= -1 #it is flipped?
    inputY = np.roll(inputY,25-np.argmax(inputY))
    inputX = np.arange(0,1024)*np.diff(inputX)[0]
    inputY = inputY[:1024]
    zeroMean = np.mean(inputY)
    inputY -= zeroMean
    print( "input: len=",len(inputX), " dT=",np.diff(inputX)[0])
    inputF,inputFFT = tfu.genFFT(inputX,inputY)

    #output pulse (after going through ALL cables)
    outputX,outputY = tff.palAntImporter('OFF-BS-PLSR/CALIBRATION/06_22-BARREL_CAL-V-C-waveform.csv')
    outputY = np.roll(outputY,25-np.argmax(outputY))
    #outputX = outputX[:1024]
    outputX = np.arange(0,1024)*np.diff(outputX)[0]
    outputY = outputY[:1024]
    zeroMean = np.mean(outputY)
    outputY -= zeroMean
    print( "output: len=",len(outputX), " dT=",np.diff(outputX)[0])
    outputF,outputFFT = tfu.genFFT(outputX,outputY)

    #output pulse had a 20dB attenuator on it that isn't in the full setup
#    atten = 20.
#    outputFFT *= 10**(atten/20)

    outputX,outputY = tfu.genTimeSeries(outputF,outputFFT)
    
    tfF = outputF
    tfFFT = outputFFT/inputFFT

    tfX,tfY = tfu.genTimeSeries(tfF,tfFFT)

    if (showPlots):
        fig,ax = lab.subplots(2,2)
        
        ax[0][0].plot(inputX,inputY,label="direct from coupler")
        ax[0][0].plot(outputX,outputY,label="through all cables")
        ax[0][0].legend()
        
        inputLogMag = tfu.calcLogMag(inputF,inputFFT)
        outputLogMag = tfu.calcLogMag(outputF,outputFFT)
        ax[0][1].plot(inputF,inputLogMag,label="direct from coupler")
        ax[0][1].plot(outputF,outputLogMag,label="through all cables")
        ax[0][1].legend()
                
#        ax[1][0].plot(tfX,tfY,label="transFunc")
#        ax[1][0].legend()
        
        tfLogMag = tfu.calcLogMag(tfF,tfFFT)
        ax[1][1].plot(tfF,tfLogMag,label="trans func")
        ax[1][1].legend()
        
        fig.show()

    return tfF,tfFFT

"""
This is for generating the antenna responses for antennas we didn't record phase info for.  It keeps the phase information from the average transfer function and scales the gain by the antenna gain we measured 
"""

def doPalAntWithAveragePhase(chan, showPlots=False, savePlots=False, writeFiles=False):
    pol = chan[3]
    # load up average transfer function
    tempIn = np.genfromtxt("transferFunctions/palAntA4Average_" + pol + ".txt", delimiter=" ")
    tfX = tempIn.T[0]
    tfY = tempIn.T[1]
    tfF,tfFFT = tfu.genFFT(tfX,tfY)
    
    # load up average waveform
    tempInWF = np.genfromtxt("transferFunctions/palAntA4AverageS21_" + pol + ".txt", delimiter=" ")
    aveF = tempInWF.T[0]
    aveFFTLog = tempInWF.T[1]
    aveFFT = 10.**(aveFFTLog/20.)

    #load up S21 of the channel you want
    antF,antFFTLog=importPalAntS21(chan)
    antFFT = 10.**(antFFTLog/20.)
    antFFT /= aveFFT
    fnew = interpolate.interp1d(antF, antFFT, "cubic")
    
    for i in range(len(tfF)):
        try:
            m = np.absolute(tfFFT[i]) * fnew(tfF[i])
        except ValueError:
            m = np.absolute(tfFFT[i])
        p = np.angle(tfFFT[i])
        tfFFT[i] = m*np.cos(p) + m*np.sin(p)*1j

    tfX,tfY = tfu.genTimeSeries(tfF, tfFFT)
    tfY = tfu.hanningTail(tfY, np.argmax(tfY)+400, 100)
    tfY = tfu.hanningTailNeg(tfY, np.argmax(tfY)-20, 20)
    tfY -= np.mean(tfY)

    if writeFiles:
        tff.writeOne([tfX,tfY],chan,"palAntA4")
        
    return tfX,tfY,tfF,tfFFT

"""
  Evaluating the transfer function chain for a TUFF. Assumes input frequency array
  (f) is in Hz. By initially not setting the pF capacitance for each notch's
  variable capacitor (varCap1, varCap2, varCap3), the default behavior is for the
  notches to be off. Also default behavior is for phase to be in radians.
"""
def doTUFF(f, varCap1 = None, varCap2 = None, varCap3 = None, deg = False):
    R = 6  #  Parasitic notch resistance in Ohms.
    C = 0.6  #  Parasitic notch capacitance in pF.
    L = 56e-9  #  Notch inductance in H.
    
    # Calculate total admittance (Y = 1 / Z, Z being impedance) from parallel notch components.
    Ynotches = np.zeros_like(f, 'complex')  #  Start of total notch admittance.
    def Ynotch(capVal):  #  Input admittance function for notch filters.
        Znotch = R + 1j * 2 * np.pi * f * L + (1j * 2 * np.pi * f * (C + capVal) * 1e-12)**-1
        return Znotch**-1
    #  Starting here, we apply additional admittances when notches switched on.
    if varCap1 is not None: Ynotches += Ynotch(1.8 + varCap1)
    if varCap2 is not None: Ynotches += Ynotch(12 * varCap2 / (12 + varCap2))
    if varCap3 is not None: Ynotches += Ynotch(1.5 * varCap3 / (1.5 + varCap3))
    
    #  Calculate complex gain from TUFF passive components.
    GTUFF = (2 + 50 * Ynotches)**-1
    
    #  Find transfer function in dB and its phase.
    TdBTUFF = 40 + 20 * np.log10(np.abs(GTUFF))  #  Accounting for 40 dB amplification.
    phaseTUFF = np.angle(GTUFF, deg)  #  Phase in radians by default.
    
    return TdBTUFF, phaseTUFF


"""
  Essentially the same as doTUFF() above, except it replaces the variable
  capacitance input with resonant frequencies in Hz for driven series RLC
  circuits. To switch off a given notch, set the notch's corresponding resonant
  frequency equal to 'None' (ex. doTUFFFreq(f, resFreq1 = None)).
"""
def doTUFFFreq(f, resFreq1 = 260e6, resFreq2 = 375e6, resFreq3 = 460e6, deg = False):
    R = 6  #  Parasitic notch resistance in Ohms.
    L = 56e-9  #  Notch inductance in H.
    
    # Calculate total admittance (Y = 1 / Z, Z being impedance) from parallel notch components.
    Ynotches = np.zeros_like(f, 'complex')  #  Start of total notch admittance.
    def Ynotch(freqVal):  #  Input admittance function for notch filters.
        Znotch = R + 1j * 2 * np.pi * f * L
        Znotch += (0.5 * R**2 + (2 * np.pi * freqVal * L)**2) / (1j * 2 * np.pi * f * L)
        return Znotch**-1
    #  Starting here, we apply additional admittances when notches switched on.
    if resFreq1 is not None: Ynotches += Ynotch(resFreq1)
    if resFreq2 is not None: Ynotches += Ynotch(resFreq2)
    if resFreq3 is not None: Ynotches += Ynotch(resFreq3)
    
    #  Calculate complex gain from TUFF passive components.
    GTUFF = (2 + 50 * Ynotches)**-1
    
    #  Find transfer function in dB and its phase.
    TdBTUFF = 40 + 20 * np.log10(np.abs(GTUFF))  #  Accounting for 40 dB amplification.
    phaseTUFF = np.angle(GTUFF, deg)  #  Phase in radians by default.
    
    return TdBTUFF, phaseTUFF

def doAllAntennas(showPlots = False, savePlots = False, writeFiles = False):

    finishedChannels = []
    chans = np.loadtxt("chanList.txt", dtype = bytes).astype(str)
    channelsWithPhaseInfo = ["02TH","07TH","12MH","13MH","02TV","07TV","12MV","13MV"]
    for chan in channelsWithPhaseInfo:
        doPalAnt(chan, showPlots, savePlots, writeFiles)
        finishedChannels.append(chan)

    makePalAntAverageTransferFunction("H", showPlots, savePlots, writeFiles)
    makePalAntAverageTransferFunction("V", showPlots, savePlots, writeFiles)
    makePalAntAverageS21("H", writeFiles)
    makePalAntAverageS21("V", writeFiles)
    
    for chan in chans:
        if chan not in channelsWithPhaseInfo:
            print(chan)
            doPalAntWithAveragePhase(chan, showPlots, savePlots, writeFiles)
            finishedChannels.append(chan)

    return finishedChannels

















