#!/usr/bin/python

"""
Fits a variety of stuff to a spaceish separated input file
-brotter March 2013
"""
#import useful stuff
import sys
import os
import pylab as lab
import numpy as np
import scipy
from scipy.optimize import leastsq
import math as m

"""
#open file
data_file = open(sys.argv[1])
data_lines = data_file.readlines()

#initialize the empty data arrays
data_time=[]
data_voltage_pwrdtr=[]
data_voltage_input=[]
"""

#define function we want to fit to
#pe md as
lambdaLin      = lambda p, x: p[0]*(x-p[1])+p[2]
lambdaExpo     = lambda p, x: p[0]*(1-np.exp((x-p[1])*-p[2]))
lambdaRayleigh = lambda p, x: ((x-p[1])/p[0]**2)*np.exp(-((x-p[1])**2)/(2*p[0]**2))
lambdaMaxBoltz = lambda p, x: ((x-p[1])**2)*np.exp(-((x-p[1])**2)/(2*p[0]**2))/(p[0]**3)
lambdaGaussian = lambda p, x: p[3]*np.exp(((x-p[1])/p[0])**2)+p[2]
lambdaNormal   = lambda p, x: (1/(p[1]*np.sqrt(2*m.pi)))*np.exp(-1*((x-p[0])**2/(2*p[1]**2)))+p[2]
lambdaSin      = lambda p, x: p[1]*sin((2*m.pi*x)/p[0]+p[3])+p[2]
lambdaPoly2    = lambda p, x: p[2]*(x**2)+p[1]*x+p[0]
lambdaPoly3    = lambda p, x: p[3]*(x**3)+p[2]*(x**2)+p[1]*x+p[0]
lambdaPoly4    = lambda p, x: p[4]*(x**4)+p[3]*(x**3)+p[2]*(x**2)+p[1]*x+p[0]
lambdaDecay    = lambda p, x: p[0]*np.exp(-x/p[1])


def fitLin(aDataX,aDataY,pGuess):
    return fitFunction(aDataX,aDataY,lambdaLin,pGuess)

def fitExpo(aDataX,aDataY,pGuess):
    return fitFunction(aDataX,aDataY,lambdaExpo,pGuess)

def fitRayleigh(aDataX,aDataY,pGuess):
    return fitFunction(aDataX,aDataY,lambdaRayleigh,pGuess)

def fitMaxBoltz(aDataX,aDataY,pGuess):
    return fitFunction(aDataX,aDataY,lambdaMaxBoltz,pGuess)

def fitGaussian(aDataX,aDataY,pGuess): #4
    return fitFunction(aDataX,aDataY,lambdaGaussian,pGuess)

def fitNormal(aDataX,aDataY,pGuess): #3
    return fitFunction(aDataX,aDataY,lambdaNormal,pGuess)

def fitSin(aDataX,aDataY,pGuess):
    return fitFunction(np.array(aDataX),np.array(aDataY),lambdaSin,pGuess)

def fitPoly2(aDataX,aDataY,pGuess):
    return fitFunction(np.array(aDataX),np.array(aDataY),lambdaPoly2,pGuess)

def fitPoly3(aDataX,aDataY,pGuess):
    return fitFunction(np.array(aDataX),np.array(aDataY),lambdaPoly3,pGuess)

def fitPoly4(aDataX,aDataY,pGuess):
    return fitFunction(np.array(aDataX),np.array(aDataY),lambdaPoly4,pGuess)

def fitDecay(aDataX,aDataY,pGuess):
    return fitFunction(np.array(aDataX),np.array(aDataY),lambdaDecay,pGuess)



def fitFunction(aDataX,aDataY,fitFunc,pGuess):

    """
    Generalized fitting function for any function, returns fit parameters and rSq value
    """

    #and fuction to see difference between fitted function and data(y)
    errfunc = lambda p, x, y: (y-fitFunc(p,x))

    #optimize (fit)
    pFit, cov_x, infodict, mesg, ier = scipy.optimize.leastsq(errfunc,pGuess,args=(aDataX,aDataY),full_output=True)

    #calc goodness
    rSq = rSquared(aDataX,aDataY,infodict)

    return pFit,rSq

def rSquared(aDataX,aDataY,infodict):
    """
    Calculating rSquared for the fit function
    from http://stackoverflow.com/questions/7588371/scipy-leastsq-goodness-of-fit-estimator
    """

    ss_err = (infodict['fvec']**2).sum()
    ss_tot = ((aDataY-aDataY.mean())**2).sum()
    
    rSquared = 1-(ss_err/ss_tot)
    
    return rSquared




"""
#show results to world
print success,p1_fit,success2,p1_fit2
plot(func2plotx,func2ploty,'ro',label="AD8362 Vout")
plot(func2plotx,fitfunc(p1_fit,func2plotx), label=str(round(p1_fit[0],5))+"*(1-exp((x+"+str(round(p1_fit[1],5))+")*-"+str(round(p1_fit[2],5))+"))")
#     func2plotx,fitfunc2(p1_fit2,func2plotx))
grid(True)
legend(loc=4)
title("AD8362 Vout response from ~-90dBm to -26.75dBm \n w/ .1uF CLPF cap")
xlabel("Time (s)")
ylabel("Voltage (V)")
show()
"""



def findHistCenters(histoBins):
    """
    Finds the center points of a array of equally spaced bins
    """
    length = len(histoBins)
    separation = histoBins[1]-histoBins[0]
    beginCenter = histoBins[0]+(separation/2)

    
    histoCenters = []
    for i in range(0,length-1):
        histoCenters.append(beginCenter+separation*i)

    return histoCenters

def turnOnTime(data):
#find "turn on time" (when data starts to rise)
    for n in range(0,len(data_time)):
        if data_voltage_pwrdtr[n]>=0.2:
            turnon_point=n
            print data_voltage_pwrdtr[n],data_time[n],n
            break
    turnoff_point=turnon_point+100
    return turnon_point


def selectExpoZone(dataTime,dataVolt):
    #select part of plot we want to fit and put into numpy arrays
    adataTime = array(dataTime[turnon_point:turnoff_point])
    adataVolt = array(dataVoltage_pwrdtr[turnon_point:turnoff_point])
    return adataTime,adataVolt





def weightedLeastSqrs(x,y,w):

    #from https://www.che.udel.edu/pdf/FittingData.pdf

    sumX = np.sum(x*w**2)
    sumX2= np.sum(x**2*w**2)
    sumY = np.sum(y*w**2)
    sumXY = np.sum((x*y)*w**2)
    sumW = np.sum(w**2)
    

    slope = ( (sumX*sumY)-(sumXY*sumW) ) / ( sumX**2 - (sumX2)*(sumW) )
    int = ( sumXY - slope*sumX2 ) / sumX

    return slope,int
