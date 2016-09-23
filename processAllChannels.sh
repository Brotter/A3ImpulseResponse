#!/bin/bash

##################################################
#
# Ben Rotter - September 2016 - University of Hawaii at Manoa
# BenJRotter@gmail.com
#
#  Just loops over the code that finds and generates averaged waveforms for the impulse response calibration
#
###################################################



for ant in `cat chanList.txt`; do
    ./avgSigChainWaveform ${ant} >> allChannels.log 2>&1
done
    