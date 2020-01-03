#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Alexis Perez-Bellido v1.0

version = 1.0

import numpy as np
import matplotlib.pyplot as plt
import os, sys
import exp_func as exp
import stimuli as st

from psychopy import visual, logging, core, event,  gui, data
from psychopy.tools.filetools import fromFile, toFile # wrappers to save pickles
from random import random
from numpy import sin, pi
from numpy.random import vonmises
from scipy import signal, stats
from psychopy.preferences import prefs

# Some general presets
#event.globalKeys.clear() # implementing a global event to quit the program at any time by pressing ctrl+q
#event.globalKeys.add(key='q', modifiers=['ctrl'], func=core.quit)



prefs.hardware['audioLib']=['pyo']
from psychopy.sound import Sound

from psychopy.sound import Sound

Sound

med.play()
prefs.general['audioLib'] = ['pygame']

low = Sound(600, sampleRate=44100, secs=0.1, stereo=True ,loops=300,  blockSize=128, preBuffer=-1)
med = Sound(800, sampleRate=44100, secs=0.1 )
high = Sound(1000, sampleRate=44100, secs=0.1, stereo=True )

med = Sound(800, sampleRate=44100, secs=0.1, stereo=True )


low = Sound(600, sampleRate=44100, secs=0.1, stereo=True ,loops=0,hamming=True) 
high = Sound(1000, sampleRate=44100, secs=0.1, stereo=True,loops=0,hamming=True )

low.play()

high.play()

med = Sound(800, sampleRate=44100, secs=0.1, stereo=True ,loops=0,hamming=True)
med.play() 

high.play()

standard # you can verify that the audiolib has changed

import 
psychopy.sound.backend_sounddevice.SoundDeviceSound(value='C', secs=0.5, octave=4, stereo=-1, volume=1.0, loops=0, sampleRate=None, blockSize=128, preBuffer=-1, hamming=True, startTime=0, stopTime=-1, name='', autoLog=True)