import numpy as np
import pylab as lab

import glob

import seaborn as sns

"""

Seavey gave us some data on the antennas long ago, I'm going to plot it nicely to compare


"""



dataDir = "/Users/brotter/Science/ANITA/ANITA3/ANITA_NewHorn_SeaveyMeasurements/"



def loadVEFiles():
    fileList = glob.glob(dataDir+"VportEplane*.csv")
    data = {}

    for file in fileList:
        name = file.split("_")[-1].split(".")[0][:-3]
        data[name] = np.loadtxt(file,skiprows=1,delimiter=",").T

    return data


def loadHEFiles():
    fileList = glob.glob(dataDir+"HportEplane*.csv")
    data = {}

    for file in fileList:
        name = file.split("_")[-1].split(".")[0][:-3]
        data[name] = np.loadtxt(file,skiprows=1,delimiter=",").T

    return data




def plotSingleData(data,title=""):

    fig,ax = lab.subplots()

    fig.suptitle(title)
    for key in data:
        x,y = data[key]
        ax.plot(x,y,label=key)

    ax.legend()

    fig.show()



def plotAll(savePlot=False):

    sns.set_palette(sns.color_palette("husl",5)[::-1])

    fig,ax = lab.subplots(1,2,figsize=[11,8.5])
    fig.suptitle("ANITA3/4 Seavey Factory Antenna Data")


    ax[0].set_title("Vertical")
    vData = loadVEFiles()
    keys = vData.keys()
    keys.sort(key=int)
    for key in keys:
        x,y = vData[key]
        ax[0].plot(x,y,label=key+" MHz")


    ax[1].set_title("Horizontal")
    hData = loadHEFiles()
    keys = hData.keys()
    keys.sort(key=int)
    for key in keys:
        x,y = hData[key]
        ax[1].plot(x,y,label=key+" MHz")


    ax[0].set_xlim([-60,60])
    ax[1].set_xlim([-60,60])

    ax[0].set_ylim([-6,0.2])
    ax[1].set_ylim([-6,0.2])

    ax[0].set_xlabel("Degrees From Boresight")
    ax[1].set_xlabel("Degrees From Boresight")

    ax[0].set_ylabel("Gain from Maximum (dB)")
    ax[1].set_ylabel("Gain from Maximum (dB)")

    ax[1].legend(frameon=True,bbox_to_anchor=(1.2,1))

    if savePlot:
        fig.savefig("seaveyDataPlotter.pdf")



    fig.show()

    
