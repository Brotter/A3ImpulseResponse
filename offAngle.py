#/usr/bin/python


"""
My transferFunction.py code (in integratedTF) is sort of a big mess, so this code is specifically for 
looking at the differences between the angles.  It probably relies on that other code
though, since it is nice code.


This one does the stuff in the chamber too!  I don't know what I am doing in terms of organization, I just have such a big mess of code. I should clean it up.

Update Nov 1st 2016:  This code is a disaster.  Some of it works, some of it doesn't, it's been awhile since I've messed with it.

"""


import numpy as np
import pylab as lab
import glob
import copy #im suprised how little I end up using this...

import tfUtils as tf
import transferFunctionFinal as tFunc #transferFunction.py is old, transferFunctionFinal is new

from scipy import signal

import seaborn as sns

import remoteConnect as rc



#for saving multiple plots to a pdf
from matplotlib.backends.backend_pdf import PdfPages


#directories with antenna pulse data
localDir = "/Volumes/ANITA3Data/"

roofDir = "Rooftop_Seevey_Antenna/Antenna_Impulse_Testing/"
currLocalRoofDir=localDir+roofDir+"Attempt_2_02July2012/Hpol/"
remoteRoofDir = "/storageb/ANITA-3/"+roofDir
currRemoteRoofDir=remoteRoofDir+"Attempt_2_02July2012/Hpol/"

chamberDir="Anechoic_Chamber_New_Horn_Impulse/ANITA2Antennas/HPol_Pulse"
remoteChamberDir="/storage/Krauss_Backup/ANITA/"+chamberDir

chamberIdentDir="Anechoic_Chamber_New_Horn_Impulse/IdenticalAntennas/Ant2Blue_to_Ant0/HPol_Pulse/HPol_Pulse_New_Cable/"
chamberRefPulse = localDir+chamberIdentDir+"Impulse_NewRedCable_23May2012.dat"

chamberBoresightPulse = localDir+chamberIdentDir+"Impulse_p180Deg_t090Deg_22May2012_H_Two.dat"

def main():
    """
    does whatever silly thing I am trying on the command line
    """

    files = localFileList_ANITA2chamberHPol()

    fig = lab.figure()
    axes = []

    colors=["red","blue"]

    numFiles = len(files)
    index=1
    logMags = {}
    for file in files:
        if index!=1:
            axes.append(fig.add_subplot(numFiles,1,index,sharex=axes[0]))
        else:
            axes.append(fig.add_subplot(numFiles,1,index))
        
        for chan in range(0,2):
            events = importAll(file,numEvents=100,channel=chan)
            avgEvent = tf.correlateAndAverage(events)
            peak = np.argmax(np.abs(avgEvent[1])) + 500
            procAvgEvent = processAllEvents([avgEvent],totalWidth=2000,peak=peak)[0]        
            logMag = tf.genLogMag(procAvgEvent[0],procAvgEvent[1]*1000)
            grpDly = tf.genGroupDelay(procAvgEvent[0]*1e9,procAvgEvent[1]*1000)

            offPhi = file.split("/")[-1].split("_")[1][1:4]
            offTheta = file.split("/")[-1].split("_")[2][1:4]
            print "Phi="+str(offPhi)+" Theta="+str(offTheta)

            
            label = str(offPhi)+","+str(offTheta)
            axes[-1].plot(grpDly[0][1:],grpDly[1][1:],label=label,color=colors[chan])
        index += 1
        

    for ax in axes:
        ax.legend()
    fig.show()



    return


def ANITA2_HchamberAnglePlot():


    files = localFileList_ANITA2chamberHPol()

    #Where on the grid the antenna was pointing for the chamber data
    graphPositions = (2,4,5,6,8)

    fig = lab.figure()
    axes = []

    colors=["red","blue"]

    index=0
    logMags = {}
    for file in files:
        if index!=0:
            axes.append(fig.add_subplot(3,3,graphPositions[index],sharex=axes[0]))
        else:
            axes.append(fig.add_subplot(3,3,graphPositions[index]))
        
        for chan in range(0,2):
            events = importAll(file,numEvents=100,channel=chan)
            avgEvent = tf.correlateAndAverage(events)
            peak = np.argmax(np.abs(avgEvent[1])) + 500
            procAvgEvent = processAllEvents([avgEvent],totalWidth=2000,peak=peak)[0]        
            logMag = tf.genLogMag(procAvgEvent[0],procAvgEvent[1]*1000)
            grpDly = tf.genGroupDelay(procAvgEvent[0]*1e9,procAvgEvent[1]*1000)

            offPhi = file.split("/")[-1].split("_")[1][1:4]
            offTheta = file.split("/")[-1].split("_")[2][1:4]
            print "Phi="+str(offPhi)+" Theta="+str(offTheta)

            
            label = str(offPhi)+","+str(offTheta)
#            axes[-1].plot(grpDly[0][1:],grpDly[1][1:],label=label,color=colors[chan])
            axes[-1].plot(logMag[0][1:200],logMag[1][1:200],label=label,color=colors[chan])
        index += 1
        

    for ax in axes:
        ax.legend()
    fig.show()



    return

    

def ANITA3_chamberResponse(files=-1,rotate="phi",channel=1,eventAvgs=50):
    #imports and generates some arrays for angle, frequency, and spectral magnitude of
    # correlated and averaged waveforms in a series of chamber data

    #the defaults are for HPol rotating in phi, but you can do others.
    # channel==0 is Vpol and channel==1 is Hpol (I think)


    if files==-1:
        files = localFileList_ANITA3chamberHPol()

    freqs = []
    allFreqs = []
    mags = []
    allMags = []
    angles = []

    #params
    numFreqs = 25
    stop = 100 #~1GHz
    start = 2
    step = stop/float(numFreqs)


    sns.set_palette(sns.color_palette("husl",numFreqs))

    inPulse =  antChamberResponse(chamberRefPulse,channel=0)

    for file in files:
        events = importAll(file,channel=channel)
        print len(events)
        fftAvgs = len(events)/eventAvgs
        for fftAvgNum in range(0,fftAvgs):
            avgEvent = tf.correlateAndAverage(events[eventAvgs*fftAvgNum:eventAvgs*(fftAvgNum+1)])
            peak = np.argmax(np.abs(avgEvent[1])) + 500
            print "I'M CRAB V.v.V"
            procAvgEvent = processAllEvents([avgEvent],totalWidth=2000,peak=peak)[0]
            transFunc = HpolTransferFunctionAnalysis(pulse=inPulse,respo=procAvgEvent)
            #friis equation transforms to S21!  (Gr*Gt)/(4*pi*R/lambda)**2
            R=5.9
            f,fft = tf.genFFT(transFunc[0],transFunc[1])
            friis = (fft**2)/(4.*np.pi*R*f*3e8)**2
            logMag = tf.calcLogMag(f,friis*1000)
            if freqs == []:
                freqs = np.ceil(f[start:stop:step]/1e6)
                allFreqs = np.ceil(f[start:stop]/1e6)
            if fftAvgNum == 0:
                logMagAvg = logMag/fftAvgs
            else:
                logMagAvg += logMag/fftAvgs

            

        mags.append(logMagAvg[start:stop:step])
        allMags.append(logMagAvg[start:stop])

        if rotate=="phi":
            angles.append(file.split("/")[-1].split("_")[1][1:4])
        if rotate=="theta":
            angles.append(file.split("/")[-1].split("_")[2][1:4])

        print(len(freqs))

    lab.close("all")
    fig = lab.figure()
    ax = fig.add_subplot(111)
    for freqNum in range(0,numFreqs):
        curve = []
        for angleNum in range(0,len(angles)):
            curve.append(mags[angleNum][freqNum])
        curve = np.array(curve)-np.max(curve)
        ax.plot(angles,curve,label=freqs[freqNum])
        

                  
    ax.legend()
    fig.show()


    return angles,allFreqs,allMags

def antChamberResponse(file,channel=1,eventAvgs=50):
    events = importAll(file,channel=channel)
    print len(events)
    fftAvgs = len(events)/eventAvgs
    for fftAvgNum in range(0,fftAvgs):
        avgEvent = tf.correlateAndAverage(events[eventAvgs*fftAvgNum:eventAvgs*(fftAvgNum+1)])
        peak = np.argmax(np.abs(avgEvent[1])) + 500
        procAvgEvent = processAllEvents([avgEvent],totalWidth=2000,peak=peak)[0]            

    return procAvgEvent


def HpolTransferFunctionAnalysis(pulse=False,respo=False,makeFigs=False):
    #import stuff
    if type(pulse)==bool and type(respo)==bool:
        print "type(pulse)"
        pulse = antChamberResponse(chamberRefPulse,channel=0)
        respo = antChamberResponse("/Volumes/BenANITA3Data/Anechoic_Chamber_New_Horn_Impulse/IdenticalAntennas/Ant2Blue_to_Ant0/HPol_Pulse/Impulse_p180Deg_t090Deg_16May2012_H.dat")

    #quickly making them have the same time steps
    pulseArgMax = np.argmax(pulse[1])
    respoArgMax = np.argmax(respo[1])
    offset = pulseArgMax - respoArgMax
    dT = (pulse[0][1] - pulse[0][0])
    pulseNewT = (np.arange(0,len(pulse[0]))-pulseArgMax)*dT
    respoNewT = (np.arange(0,len(respo[0]))-respoArgMax)*dT
    
    #make a quick plot of that
    if makeFigs==True:
        figInWaves = lab.figure(figsize=(16,8))
        axInWaves1 = figInWaves.add_subplot(111)
        axInWaves1.plot(pulseNewT,pulse[1],color="red",lw=1,alpha=0.6)
        axInWaves1.set_xlabel("Time (ns)")
        axInWaves1.set_ylabel("Input Pulse (V)",color="red")
        axInWaves1.set_title("HPol chamber input and output pulses")
        
        axInWaves2 = axInWaves1.twinx()
        axInWaves2.plot(respoNewT,respo[1],color="blue",lw=1,alpha=0.6)
        axInWaves2.set_ylabel("Recieved Pulse (V)",color="blue")
        axInWaves2.grid(b=True,which="major",color="white",linestyle="--")
        
        figInWaves.show()
        figInWaves.savefig("HpolBoresightWaveforms.png")

    #take the fourier transform
    pulseF,pulseFFT = tf.genFFT(pulseNewT,pulse)

    respoF,respoFFT = tf.genFFT(respoNewT,respo)


    #plot the log magnitude of that
    pulseLogMag = tf.calcLogMag(pulseF,pulseFFT)
    respoLogMag = tf.calcLogMag(respoF,respoFFT)

    if makeFigs == True:
        figLogMag = lab.figure(figsize=(16,8))
        axLogMag1 = figLogMag.add_subplot(111)
        axLogMag1.plot(pulseF,pulseLogMag,color="red",lw=1,alpha=0.6)
        axLogMag1.set_xlabel("Frequency (GHz)")
        axLogMag1.set_ylabel("Input Pulse (dBm)",color="red")
        axLogMag1.set_title("HPol chamber input and output spectral magnitude")
        
        axLogMag2 = axLogMag1.twinx()
        axLogMag2.plot(respoF,respoLogMag,color="blue",lw=1,alpha=0.6)
        axLogMag2.set_ylabel("Recieved Pulse (dBm)",color="blue")
        axLogMag2.grid(b=True,which="major",color="white",linestyle="--")
        
        figLogMag.show()
        figLogMag.savefig("HpolBoresightMagnitudes.png")

    #generate a transfer function
    transFFT = respoFFT/pulseFFT

    #Also cut off the top frequencies (since there isn't any power there)
    cutFraction = 4
    trans = tf.fftw.irfft(transFFT[:len(respo[0])/(cutFraction*2)+1])
    transNewT = respoNewT[::cutFraction]

    #it is phase aligned at zero and I want it to start in like 10a
    trans = np.roll(trans,len(respo[0])/10)

    if makeFigs == True:
        figTrans = lab.figure(figsize=(16,8))
        axTrans = figTrans.add_subplot(111)

        axTrans.plot(respoNewT[::4],trans)
        axTrans.set_xlabel("Time (ns)")
        axTrans.set_ylabel("Transfer Function (V)",color="black")
        axTrans.set_title("HPol Chamber Boresight Transfer Function")
        
        figTrans.show()
        figTrans.savefig("HpolBoresightTransferFunction.png")

    return transNewT,trans


def normalizeHeatmap(mags):
    mags = np.array(mags)

    normMags = []
    for i in range(0,len(mags.T)):
        column = np.array(mags.T[i]) - np.max(mags.T[i])
        normMags.append(column)

    normMags = np.array(normMags)
    return normMags.T


def plotHeatmapAndContour(angles,freqs,mags,title="test",center=180,invert=True):
    print center
    fig = plotHeatmap(angles,freqs,mags,title=title,center=center,invert=invert)
    plotContour(angles,freqs,mags,fig=fig,center=center)
    return fig


def plotHeatmap(angles,freqs,mags,title="test",center=180,invert=True):

    sns.set_style("whitegrid")

    #because I don't do this somewhere else
    angles = np.array(angles,dtype=np.float)
    freqs = np.array(freqs,dtype=np.float)
    mags = np.array(mags,dtype=np.float)

    #the ytick labels should be around zero probably
    #also note: angles[0] is negative and angles[-1] is positive
    #also, for some reason, are inverted from mags?
    #probably has to do with how imshow draws...
    angles = np.array(angles,dtype=np.int) - center


    #plot it!
    #note: need to turn off colorbar at first, since to access it's params you need to create it

    fig = lab.figure(figsize=(16,8))

    print freqs[0],freqs[-1],angles[0],angles[-1]
    dAngle = angles[1]-angles[0]
    dFreq = freqs[1]-freqs[0]
    #imshow makes 2d histograms :O  It just is sort of finicky

    #invert means that the angles go from positive to negative, so you have to flip everything
    # since imshow ALWAYS goes from top to bottom
    if invert==True:
        ax=lab.imshow(mags,cmap=lab.cm.gist_heat,aspect="auto",
                      extent=np.array((freqs[0],freqs[-1]+dFreq,
                                       angles[-1],angles[0]),dtype=np.int),
                      interpolation="none")
    elif invert==False:
        ax=lab.imshow(mags,cmap=lab.cm.gist_heat,aspect="auto",
                      extent=np.array((freqs[0],freqs[-1]+dFreq,
                                       0,len(angles)),dtype=np.int),
                      interpolation="none")
        
    ax.axes.set_title(title)

    
    
    if invert==True:
        ax.axes.set_yticks(np.linspace(angles[-1]-dAngle/2.,angles[0]+dAngle/2.,len(angles)))
        angles = np.array(angles,dtype=np.str)
        ax.axes.set_yticklabels(angles[::-1])

    elif invert==False:
        ax.axes.set_yticks(np.linspace(0+0.5,len(angles)-0.5,len(angles)))
        angles = np.array(angles,dtype=np.str)
        ax.axes.set_yticklabels(angles)

    print dFreq
    ax.axes.set_xticks(np.arange(freqs[0]+dFreq/2.,freqs[-1]+dFreq/2.,10*dFreq))
    ax.axes.set_xticklabels(freqs[::10])
    
    #label it!
    lab.xlabel("Frequency (MHz)")
    lab.ylabel("E-field rotation away from boresight")
    cbar =fig.colorbar(ax)
    cbar.set_label("Spectral power below maximum (dB)")

    fig.show()

    return fig
    


def plotHeatmapSeaborn(angles,freqs,mags):
    #Seaborn has a heatmap that inputs datastructs (whatever those are)
    # THIS IS DUMB because it doesn't map well to the lab.contour() plot, you need imshow for that


    #annoying way to set custom xticks in a heatmap >.<
    xticklabels = []
    for freq in freqs:
        if freq%100 != 0:
            xticklabels.append("")
        else:
            xticklabels.append(freq)

    #the ytick labels should be around zero probably
    angles = np.array(angles,dtype=np.int) - 180
    angles = np.array(angles,dtype=np.str)


    #plot it!
    #note: need to turn off colorbar at first, since to access it's params you need to create it

    fig = lab.figure()
    ax = sns.heatmap(-mags[::-1],yticklabels=angles[::-1],xticklabels=xticklabels,cbar=False)
    
    #label it!
    lab.xlabel("Frequency (MHz)")
    lab.ylabel("Phi rotation away from boresight")
    cbar = ax.figure.colorbar(ax.collections[0])
    cbar.set_label("Spectral power below maximum (dB)")

    fig.show()

    return fig


def plotContour(angles,freqs,mags,fig=-1,center=180):

    if fig==-1:
        newFig = True
        print "making new figure"
        fig = lab.figure()
        ax = fig.add_subplot(111)
    else:
        newFig = False
        ax = fig.get_axes()[0]


    #the ytick labels should be around zero probably
    angles = np.array(angles,dtype=np.int) - center
    angles = np.array(angles,dtype=np.str)

    dFreq = freqs[1]-freqs[0]
    contours = (-20,-10,-6,-3,-2,-1,-0.1)
    print "fart"
    colors = ("lightgrey","silver","grey","red","dimgrey","black","green")
    cs = ax.contour(freqs+dFreq/2.,angles,mags,contours,colors=colors)
    lab.clabel(cs, inline=1, fontsize=10)

    if newFig==False:
        fig.canvas.draw()
    else:
        fig.show()

##
#These fuctions are for generating the file lists of what I want to look at,
# that is why they have a weird syntax (the fileList_whatever).
# normally I hate underscores but I'll do this for now
##
def remoteFileList_chamberHPol(host="charm.phys.hawaii.edu",port=2505):
    baseDir = "/storage/Krauss_Backup/ANITA/Anechoic_Chamber_New_Horn_Impulse/ANITA2Antennas/HPol_Pulse/"
    
    remote = rc.remoteConnection(host,port=port)
    files = remote.sftpClient.listdir(baseDir)

    returnFiles = []
    for fileName in files:
        returnFiles.append(baseDir+fileName)
    
    if len(returnFiles) == 1:
        returnFiles = [returnFiles]

    return returnFiles

def localFileList_ANITA2chamberHPol():
    baseDir = "/Volumes/BenANITA3Data/Anechoic_Chamber_New_Horn_Impulse/ANITA2Antennas/HPol_Pulse/"

    returnFiles = (baseDir+"Impulse_p180Deg_t120Deg_15June2012.dat",
                   baseDir+"Impulse_p150Deg_t090Deg_15June2012.dat",
                   baseDir+"Impulse_p180Deg_t090Deg_15June2012.dat",
                   baseDir+"Impulse_p210Deg_t090Deg_15June2012.dat",
                   baseDir+"Impulse_p180Deg_t050Deg_15June2012.dat")

#    returnFiles = glob.glob(baseDir+"*.dat")

    return returnFiles


def localFileList_ANITA3chamberHPol():
    baseDir = "/Volumes/BenANITA3Data/Anechoic_Chamber_New_Horn_Impulse/IdenticalAntennas/Ant2Blue_to_Ant0/HPol_Pulse/"

    returnFiles = glob.glob(baseDir+"*.dat")

    #duplicate 180 degree files
    returnFiles.pop(5)
    returnFiles.pop(5)

    for file in returnFiles:
        print file

    return returnFiles


def localFileList_ANITA3chamberVPol():
    baseDir = "/Volumes/BenANITA3Data/Anechoic_Chamber_New_Horn_Impulse/IdenticalAntennas/Ant2Blue_to_Ant0/VPol_Pulse/"

    returnFiles = glob.glob(baseDir+"*p180Deg_*.dat")


    for file in returnFiles:
        print file

    return returnFiles


def localFileList_ANITA3chamberHPolFlipped():
    baseDir = "/Volumes/BenANITA3Data/Anechoic_Chamber_New_Horn_Impulse/IdenticalAntennas/Ant2Blue_to_Ant0/HPol_Pulse/"

    returnFiles = glob.glob(baseDir+"*.dat")

    #duplicate 180 degree files
    returnFiles.pop(5)
    returnFiles.pop(5)

    for file in returnFiles:
        print file

    return returnFiles
    

def ANITA3_rooftopResponse(files=-1,rotate="phi",channel=0,eventAvgs=50):
    #imports and generates some arrays for angle, frequency, and spectral magnitude of
    # correlated and averaged waveforms in a series of chamber data

    #the defaults are for HPol rotating in phi, but you can do others.
    # channel==0 is Vpol and channel==1 is Hpol (I think)


    if files==-1:
        files = localFileList_ANITA3rooftopHPol()

    print "number of files:"+str(len(files))

    freqs = []
    allFreqs = []
    mags = []
    allMags = []
    angles = []

    #params
    numFreqs = 25
    stop = 100 #~1GHz
    start = 2
    step = stop/float(numFreqs)


    sns.set_palette(sns.color_palette("husl",numFreqs))

    inPulse =  antChamberResponse(chamberRefPulse,channel=0)

    for file in files:
        events = importAll(file,channel=channel)
        print len(events)
        fftAvgs = len(events)/eventAvgs
        for fftAvgNum in range(0,fftAvgs):
            avgEvent = tf.correlateAndAverage(events[eventAvgs*fftAvgNum:eventAvgs*(fftAvgNum+1)])
            peak = np.argmax(np.abs(avgEvent[1])) + 500
            procAvgEvent = processAllEvents([avgEvent],totalWidth=2000,peak=peak).T
            transFunc = HpolTransferFunctionAnalysis(pulse=inPulse,respo=procAvgEvent)
            logMag = tf.genLogMag(transFunc[0],transFunc[1]*1000)
            if freqs == []:
                freqs = np.ceil(logMag[0][start:stop:step]/1e6)
                print freqs
                allFreqs = np.ceil(logMag[0][start:stop]/1e6)
            if fftAvgNum == 0:
                logMagAvg = logMag[1]/fftAvgs
            else:
                logMagAvg += logMag[1]/fftAvgs

        mags.append(logMagAvg[start:stop:step])
        allMags.append(logMagAvg[start:stop])

        try:
            angle = file.split("/")[-1].split("_")[1].split("deg")[0]
            direction = file.split("/")[-1].split("_")[1].split("deg")[1]
            if direction == "CW":
                angles.append(float(angle))
            elif direction == "CCW":
                angles.append(-1*float(angle))
        except:
            angles.append("0")

        print angles[-1]

    lab.close("all")
    fig = lab.figure()
    ax = fig.add_subplot(111)
    for freqNum in range(0,numFreqs):
        curve = []
        for angleNum in range(0,len(angles)):
            curve.append(mags[angleNum][freqNum])
        curve = np.array(curve)-np.max(curve)
        ax.plot(angles,curve,label=freqs[freqNum])
        

                  
    ax.legend()
    fig.show()


    return angles,allFreqs,allMags

def localFileList_ANITA3rooftopHPol():
    baseDir = localDir+roofDir+"Attempt_2_02July2012/Hpol/"

    returnFiles = glob.glob(baseDir+"*.dat")

    return returnFiles


def importAll(fileName=localDir+roofDir+"HPulse_10degCCW_148p_090t_10dBatt.dat",numEvents=-1,channel=0):
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


def remoteImportAll(fileName=currRemoteRoofDir+"HPulse_10degCCW_148p_090t_10dBatt.dat",
                    host="charm.phys.hawaii.edu",port=2505):
    """
    I always need to write more parsers, always be parsing (ABP)
    """

    print "Remotely importing from "+host+":"+str(port)+" -> "+fileName

    events = []
    event = None

    remote = rc.remoteConnection(host,port=port)
    file = remote.importRemoteFile(fileName)

    for line in file:
        if line[0] == '#':
            if event != None:
                events.append(np.array(event))
            event = [[],[]]
        else:
            event[0].append(float(line.split()[0]))
            event[1].append(float(line.split()[1]))
    
    events.append(event)

    remote.close()

    return events

def importAndMakeDoubleAveragedLogMags(fileList=-1):
    """
    Lets just average together a few in the time domain, then a few in the freq domain
    """
    eventAvgNum = 16

    if fileList==-1:
        baseDir=localDir+"/Attempt_3_23July2012/NewSeavey/Hpol/*degCW*"
        fileList = glob.glob(baseDir+"*.dat")
    else:
        try:
            fileList[0]
        except:
            fileList = [fileList]


    avgFFTs = {}
    for fileName in fileList:
        shortName = fileName.split("/")[-1]
        print shortName
        allEvents = importAll(fileName)

        avgFFT = ""
        numAvgs = len(allEvents)/eventAvgNum
        print "Number of Big Averages: "+str(numAvgs)
        print "Using "+str(numAvgs*eventAvgNum)
        for i in range(0,numAvgs):
            avgEvent = tf.correlateAndAverage(allEvents[i*eventAvgNum:(i+1)*eventAvgNum])
            if avgFFT == "":
                avgFFT = np.array(tf.genLogMag(avgEvent[0]*1e9,avgEvent[1]))
                print avgFFT[1]
            else:
                avgFFT[1] += tf.genLogMag(avgEvent[0],avgEvent[1])[1]
        avgFFTs[shortName] = (avgFFT[0],avgFFT[1]/numAvgs)


    return avgFFTs

def importAndMakeDoubleAveragedData(fileList=-1,eventAvgNum=49):
    """
    Lets just average together a few in the time domain, then a few in the freq domain
    """

    if fileList==-1:
        baseDir=localDir+"/Attempt_3_23July2012/NewSeavey/Hpol/"
        fileList = glob.glob(baseDir+"*.dat")
    else:
        if not fileList == list:
            fileList = [fileList]

    print fileList[0]

    avgDelays = {}
    avgLogMags = {}
    for fileName in fileList:
        shortName = fileName.split("/")[-1]
        if (shortName == "HPulse_60degCW_043p_090t_10dBatt_23July2012.dat" or \
#            shortName == "HPulse_60degCW_043p_090t_ANITAhornTrigger_10dBatt_23July2012.dat" or \
            shortName == "HPulse_90degCW_073p_090t_10dBatt_23July2012.dat" or \
#            shortName == "HPulse_90degCW_073p_090t_ANITAhornTrigger_10dBatt_23July2012.dat" or \
            shortName == "HPulse_BoresightToBoresight_InversePorts_10dBatt_23July2012.dat"):
            continue

        print shortName
        allEvents = importAll(fileName)

        allEvents = processAllEvents(allEvents)

        averagedEvents = averageAndPhaseAlign(allEvents)
        
        groupDelays = getDelays(averagedEvents)
        logMags = getLogMags(averagedEvents)

        avgDelays[shortName] = averageEvents(groupDelays)
        avgLogMags[shortName] = averageEvents(logMags)

    
    return avgDelays,avgLogMags

def averageEvents(events):
    """
    Not to be confused with correlateAndAverage!
    """

    events = np.array(events)

    numEvents=float(len(events))

    for event in events:
        try:
            avgEvent += event[1]/numEvents
        except:
            avgEvent = event[1]/numEvents

    return events[0][0],avgEvent
        

def averageAndPhaseAlign(allEvents,eventAvgNum=24):
    """
    Prepares for a "double average" of things.  Basically just does a boxcar average of events
    1) Does a time domain average (n=eventAvgNum events) (removes incoherant noise)
    2) Phase aligns all those waveforms against the first
    3) Returns a smaller number of averaged and phase aligned waveforms
    """

    numAvgs = len(allEvents)/eventAvgNum
    outputEvents = []
    print "Number of Big Averages: "+str(numAvgs)
    print "Using "+str(numAvgs*eventAvgNum)+" of "+str(len(allEvents))+" events"
    for i in range(0,numAvgs):
        avgEvent = tf.correlateAndAverage(allEvents[i*eventAvgNum:(i+1)*eventAvgNum])
        try:
            max,rsq,peak = tf.findPeakCorr(oldEvent,avgEvent)
            avgEvent = tf.resampleAtOffset(avgEvent,max)
        except:
            oldEvent = avgEvent
        outputEvents.append(avgEvent)

    return outputEvents

def doubleAverageDelay(allEvents,eventAvgNum=24):
    """
    Does a "double average" of the group delay, basically:
    1) Does a time domain average (n=eventAvgNum events) (removes incoherant noise)
    2) Transforms each of those averages into a delay
    3) Averages those delays together and returns THAT average
    """

    numAvgs = len(allEvents)/eventAvgNum
    print "Number of Big Averages: "+str(numAvgs)
    print "Using "+str(numAvgs*eventAvgNum)+" of "+str(len(allEvents))+" events"
    for i in range(0,numAvgs):
        avgEvent = tf.correlateAndAverage(allEvents[i*eventAvgNum:(i+1)*eventAvgNum])
        try:
            max,rsq,peak = tf.findPeakCorr(oldEvent,avgEvent)
            avgEvent = tf.resampleAtOffset(avgEvent,max)
        except:
            oldEvent = avgEvent
            
        try:
            avgFFT[1] += tf.genGroupDelay(avgEvent[0]*1e9,avgEvent[1])[1]
        except:
            avgFFT = np.array(tf.genGroupDelay(avgEvent[0]*1e9,avgEvent[1]))

    return (avgFFT[0],avgFFT[1]/numAvgs)



def processAllEvents(allEvents,totalWidth=500,slope=100,peak="auto"):
    """
    I need to work on the events a LITTLE BIT, not much, just a little
    1) Do hanning window around pulse (offset 
    2) Zero mean? (do this first!)
    """
    outAllEvents = []
    for event in allEvents:
        zeroMean = event[1] - np.mean(event[1])
        if peak=="auto":
            peak = np.argmax(np.abs(event[1])) + 100 #want to catch more tail than front
        outX,outY = tf.hanningWindow(event[0],zeroMean,
                                     peak,totalWidth=totalWidth,slope=slope)
        outAllEvents.append((outX,outY))

    return np.array(outAllEvents)

def plotDelays(delays,compare=0,boxcar=1):
    
    fig = lab.figure()

    keyList = delays.keys()
    keyList = sorted(keyList)[::-1]
    compKey = keyList.pop(compare)
    print "Comparison to:"+str(compKey)
    keyList = keyList[3:]

    sns.set_palette(sns.color_palette("husl",len(keyList)))

    x = delays[delays.keys()[0]][0]

    axCW = fig.add_subplot(211)
    axCCW = fig.add_subplot(212,sharex=axCW)
    axCW.set_xlabel("Frequency (GHz)")
    axCCW.set_xlabel("Frequency (GHz)")
    axCW.set_ylabel("Delay (ns)")
    axCCW.set_ylabel("Delay (ns)")    

    delayX = delays[keyList[compare]][0]

    for i in range(3,len(keyList),2):
        if i==compare:
            continue
        key = keyList[i]
        print key.split("_")[1]
        delayY = delays[keyList[compare]][1]-delays[key][1]

        boxX,boxY = boxcarMean(delayX,delayY,boxcar=boxcar)

        axCW.plot(boxX,boxY,label=key.split("_")[1])

    for i in range(4,len(keyList),2):
        if i==compare:
            continue
        key = keyList[i]
        print key.split("_")[1]
        delayY = delays[keyList[compare]][1]-delays[key][1]

        boxX,boxY = boxcarMean(delayX,delayY,boxcar=boxcar)

        axCCW.plot(boxX,boxY,label=key.split("_")[1])

    axCW.legend()
    axCCW.legend()
    
    fig.show()

    return


def boxcarMean(inputX,inputY,boxcar=25):
    numBoxcars = len(inputY)/boxcar
    inputY = np.array(inputY)
    inputY.resize(numBoxcars,boxcar)
    boxcarY = np.sum(inputY,axis=1)/float(boxcar)


    return inputX[::boxcar],boxcarY

def makeFigAndAxes(axisDims=(1,1)):
    fig = lab.figure()
    axes=[]
    numAxes=axisDims[0]*axisDims[1]
    for i in range(0,numAxes):
        axes.append(fig.add_subplot(axisDims[0],axisDims[1],i+1))

    return fig,axes
        
    

def filterGrp(grpDlyArray,window=11,order=3):

    filtered = {}
    for key in grpDlyArray.keys():
        filtered[key] = grpDlyArray[key][0],signal.savgol_filter(grpDlyArray[key][1],window,order)


    return filtered

def plotMain(avgEvents):

    sns.set_palette(sns.color_palette("hls",len(avgEvents)))

    fig = lab.figure()
    ax = fig.add_subplot(111)

    keys = avgEvents.keys()
    keys.sort()
    for shortName in keys:
        ax.plot(avgEvents[shortName][0],avgEvents[shortName][1],label=shortName,linewidth=5)

    ax.legend()
    fig.show()

    return
    

def getDictLogMags(avgEvents):
    
    keys = avgEvents.keys()
    keys.sort()

    logMags = {}
    for shortName in keys:
        currPlot = avgEvents[shortName]
        logMags[shortName] = tf.genLogMag(np.array(currPlot[0]),np.array(currPlot[1]))

    return logMags


def getLogMags(avgEvents):
    
    logMags = []
    for event in avgEvents:
        logMags.append(tf.genLogMag(np.array(event[0]),np.array(event[1])))

    return logMags


def getDictDelays(avgEvents):
    
    keys = avgEvents.keys()
    keys.sort()

    grpDlys = {}
    for shortName in keys:
        currPlot = avgEvents[shortName]
        grpDlys[shortName] = tf.genGroupDelay(np.array(currPlot[0]),np.array(currPlot[1]))

    return grpDlys

def getDelays(avgEvents):
    
    grpDelays = []
    for event in avgEvents:
        grpDelays.append(tf.genGroupDelay(np.array(event[0]),np.array(event[1])))

    return grpDelays


def getPhase(avgEvents):
    
    keys = avgEvents.keys()
    keys.sort()

    phases = {}
    for shortName in keys:
        currPlot = avgEvents[shortName]
        phases[shortName] = tf.genPhase(np.array(currPlot[0]),np.array(currPlot[1]))

    return phases





def importCorrelateAndAverage(fileName=localDir+roofDir+"HPulse_10degCCW_148p_090t_10dBatt.dat"):
    """
    How about instead of importing a giant chunk I just go one at a time! (sounds like fun)

    Basically only useful if you have some crazy long file

    DEPRECIATED!!!!!!!!!

    """

    
    eventNum=0
    
    avgEvent = ""
    currEvent = [[],[]]
    file = open(fileName)
    for line in file:
        if line[0] != '#':
            currEvent[0].append(float(line.split()[0]))
            currEvent[1].append(float(line.split()[1]))            
        else:
            if eventNum==0:
                print "importCorrelateAndAverage starting!"
            elif eventNum==1:
                #print "copying over first event"
                avgEvent = currEvent[:]
            else:
                max,rSq,peak = tf.findPeakCorr(currEvent[1],avgEvent[1])
                resampledEvent = np.array(tf.resampleAtOffset(currEvent[1],max)) 
                avgEvent[1] += resampledEvent
                #print str(eventNum)+" "+str(max)
            eventNum += 1
            currEvent = [[],[]]

            
    return np.array(avgEvent)/(eventNum+1)


def multipage(filename, figs=None):
    pp = PdfPages(filename)
    if figs is None:
        figs = [lab.figure(n) for n in lab.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()
