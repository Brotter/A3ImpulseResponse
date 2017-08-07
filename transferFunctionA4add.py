"""
Script written by John Russell which includes functions specific to ANITA-4 over
ANITA-3. Takes some functions from Ben Rotter's "transferFunctionFinal.py" and
modifies them for ANTIA-4 usage.
"""

import numpy as np
import pylab as lab
import tfUtils as tfu
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
                             dtype = str, delimiter = ",")

    antennaSN = str(antennaInfo.T[3][np.argwhere(antennaInfo.T[2]==chan[:3])[0][0]])

    dir = localDir + "CSBF16/OFF-BS-PLSR/"
    if inOrOut == "in":
        fileName = dir + "CALIBRATION/06_22-INPUT_PULSE_20dB_DOWN-waveform.csv"
    if inOrOut == "out":
        fileName = dir + "SN" + antennaSN + "/V-TRANS/waveform/06_23-el_0-az_0-V-C-waveform.csv"
 
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


#==============================================================================
# Processing.
#==============================================================================


"""
  Develop a transfer function for an antenna measured in palestine 2016, mapped to channel.
  Normalized to V/m (saved in V - ns fomat)
"""
def doPalAnt(chan,savePlots=False,showPlots=False,writeFiles=False):

    if showPlots or savePlots:
        lab.close("all")

    #input pulse
    inRawX,inRawY,inRawF,inRawFFT = importPalAntIn(chan)

    inX,inY = tff.processWaveform(inRawX,inRawY,"palin")
    inY = np.roll(inY,30-np.argmax(inY))

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
##    outX = outX[:1024]
##    outY = outY[:1024]
#    outY = tfu.highPass(outX,outY)
    outF,outFFT = tfu.genFFT(outX,outY)

#    outFFT = tfu.minimizeGroupDelayFromFFT(outF,outFFT)

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
    antTFFFT = tfu.sqrtOfFFT2(antTFF,antTFFFT)

    #causality conditions aren't good for measured signals (artifically distort phase)!
#    antTFFFT = tfu.makeCausalFFT(antTFFFT,np.argmax(tfu.fftw.irfft(antTFFFT)))

    #I have to get the time domain again after that
    # units of "meters / nanosecond"
    antTFX,antTFY = tfu.genTimeSeries(antTFF,antTFFFT)

    #remove all the stupid power above 1.3GHz since ANITA can't measure it (probably should be lower)
    #butterworth filter @ 1.3
    antTFY = tfu.highPass(antTFX,antTFY)
    antTFY = tfu.nyquistLimit(antTFY,5)

    #10BH is upside down?
    #A bunch are actually!  Nothing here constrains absolute polarity
    if np.max(antTFY) < -np.min(antTFY):
        print "Flipped!"
        antTFY *= -1
    
    antTFF,antTFFFT = tfu.genFFT(antTFX,antTFY)

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
            tff.writeOne([antTFX,antTFY],chan,"palAnt")
        
        return antTFX,antTFY,antTFF,antTFFFT


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
  circuits.
"""
def doTUFFFreq(f, resFreq1 = 260e6, resFreq2 = 375e6, resFreq3 = 460e6, deg = False):
    R = 6  #  Parasitic notch resistance in Ohms.
    L = 56e-9  #  Notch inductance in H.
    
    # Calculate total admittance (Y = 1 / Z, Z being impedance) from parallel notch components.
    Ynotches = np.zeros_like(f, 'complex')  #  Start of total notch admittance.
    def Ynotch(freqVal):  #  Input admittance function for notch filters.
        Znotch = R + 1j * 2 * np.pi * f * L
        Znotch += ((2 * np.pi * freqVal * L)**2 + 0.5 * R**2) / (1j * 2 * np.pi * f * L)
        return Znotch**-1
    #  Starting here, we apply additional admittances when notches switched on.
    if resFreq1: Ynotches += Ynotch(resFreq1)
    if resFreq2: Ynotches += Ynotch(resFreq2)
    if resFreq3: Ynotches += Ynotch(resFreq3)
    
    #  Calculate complex gain from TUFF passive components.
    GTUFF = (2 + 50 * Ynotches)**-1
    
    #  Find transfer function in dB and its phase.
    TdBTUFF = 40 + 20 * np.log10(np.abs(GTUFF))  #  Accounting for 40 dB amplification.
    phaseTUFF = np.angle(GTUFF, deg)  #  Phase in radians by default.
    
    return TdBTUFF, phaseTUFF


