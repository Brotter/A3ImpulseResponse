#!/bin/bash

##################################################
#
# Ben Rotter - September 2016 - University of Hawaii at Manoa
# BenJRotter@gmail.com
#
#  Just loops over the code that finds and generates averaged waveforms for the impulse response calibration
#
###################################################

suffix=${1}

for ant in `cat chanList.txt`; do
    ./avgSigChainWaveform ${ant} ${suffix}
done
    