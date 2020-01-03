#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 15:40:38 2019
Lets code a functional localizer task
@author: alex
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Alexis Perez-Bellido

# exec(open('staircase_exp.py').read())
version = 1.1

#cd('/home/node2/Experiments/PreDCN-master/predcision')
import numpy as np

import matplotlib.pyplot as plt
import os, sys
import exp_func as exp
import stimuli as st
import instructions as instr


from psychopy import visual, logging, core, event,  gui, data
from psychopy.tools.filetools import fromFile, toFile # wrappers to save pickles
from random import random
from numpy import sin, pi
from numpy.random import vonmises
from scipy import signal, stats
from psychopy.preferences import prefs
from pandas import DataFrame
# Some general presets
event.globalKeys.clear() # implementing a global event to quit the program at any time by pressing ctrl+q
event.globalKeys.add(key='q', modifiers=['ctrl'], func=core.quit)

prefs.hardware['audioLib']=['pyo'] # use Pyo audiolib for good temporal resolution
from psychopy.sound import Sound # This should be placed after changing the audio library

# monitors.getAllMonitors()
#from general_settings import * # reads the variables in one script

# Collect subject data
expInfo, resultspath = exp.stcexp_subject_info(version)
subj_id = expInfo['subjInfo']['observer']

# Loading monitor definitions
monitores = st.monitor_def() 

mon, expInfo['monitor'] = exp.define_monitor(monitores[1]) # select the correct monitor

# Creating a new experimental window
monitor_features = {}
monitor_features['monitor'] = mon
monitor_features['units'] = 'deg' # units to define your stimuli
monitor_features['screen_id'] = 0 # when using a extended display 
monitor_features['full']  = True
monitor_features['Hz'] = 'auto' #60 # this can be set to "auto" to estimate the refreshing rate of the monitor, although it can fail often

   
win, monitor_features = exp.create_window(monitor_features)
ifi = monitor_features['ifi']

expInfo['monitor_features'] = monitor_features


# Experiment timings and characteristics 
stim = st.stim_config(ifi) #Loading stim characteristics

expInfo['stim'] = stim # saving
      
basic_stim = st.draw_basic(win,expInfo['stim'])

# Create the sounds that I am going to use
#lowf_s = Sound(600, sampleRate=44100, secs=0.1, stereo=True ,loops=0,hamming=True)
medf_s = Sound(800, sampleRate=44100, secs=0.1, stereo=True ,loops=0,hamming=True)
#highf_s = Sound(1000, sampleRate=44100, secs=0.1, stereo=True ,loops=0,hamming=True)


# Experimental condition preparation
# response mapping (shuffled in each experiment )
expInfo['resp_maps']  = np.array([0, 45]) # 0 -> cardinal; 45 -> diagonal
np.random.shuffle(expInfo['resp_maps']) # first response option (cardignal or diagonal) will be placed at left, and second at right 
 
# Show some instructions 
instr.localizer_instructions(win, basic_stim, ifi)
instr.block_start(win)

# Create an experiment clock
Clock = core.Clock()
ExpClockStart = Clock.getTime() # global experiment time

#trials = data.TrialHandler(stimList, ntrials_per_cond, method='random') # when calling next trial it continues to the next randomly ordered trial

black_resps = np.array([[-1, -1, -1],[-1, -1, -1]]) # default color resp options


loc_exp  = {}
loc_exp['nblocks']     = 6 # 6 possible blocks (we will only run 2 in theory)
loc_exp['trial_reps']  = 9 # 1 # 1 bloque de 1 repeticion son 90 s aprox. 5
loc_exp['Exp_blocks']  = [None] * loc_exp['nblocks'] # assigning memory for storing block data


# Experiment design
stimList = []
for props in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]: # if props = 1, we add oddball grating
    stimList.append({'prob': props}) #, 'n_reps': n_reps


SF_difficulty = 0.5 # it can vary 1 visual degree up or down
BlockClockStart = Clock.getTime() # block experiment time
block = {} # dummy variable to save block data
trialClocktimes = np.array([]) # saving whole times here
correct_seq = np.array([]) # saving seq. of correct responses per trial sequence

thr_trials_var = [['subj','cond', 'n_odd', 'resp', 'correct', 'RT']] # saving conditions here
thr_trials_ori = [['o1','o2','o3','o4','o5','o6','o7','o8']]


trials = data.TrialHandler(stimList, loc_exp['trial_reps'], method='random') # when calling next trial it continues to the next randomly ordered trial

# the 8 main orientations that i will use for the localizer. Now they should be randmomly permuted
orient_vector = np.array([0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5]) 
orientations = np.tile(orient_vector, [trials.nTotal,1])

for i in range(orientations.shape[0]): # shuffling order orientations per trial
     np.random.shuffle(orientations[i,:])

orientations = np.around(np.deg2rad(orientations), decimals=3) # to save some space just keep 3 decimals in memory
    
# collapse the circular angles on one side of the circle to make the analyses easier (see that 45 is equal to 225
# degrees when drawing grating orientations)
orientations[orientations < 0] += np.deg2rad(180)
orientations = np.around(orientations, decimals = 3)
decision_var = signal.sawtooth(4 * ((orientations)), 0.5)  # DU decision variable
decision_var_T = np.mean(decision_var, 1)
decision_var_std = np.std(decision_var, 1)
orientation_var_std = stats.circstd((orientations), axis=1)

trialClocktimes = np.array([]) # saving whole times here

#win.getMovieFrame()  

for thisTrial in trials: # iterate over trials
     
    trial_times = np.array([]) # logging trial event time stamps
    basic_stim['fixation_point'].color = [1, 1, 1]
    ExpClockTrial = Clock.getTime()
    stim['ITI_frames'] = round(np.random.uniform(0.9, 1.1)*1000 /ifi) # only for the localizer(inter-trial interval)
  
    if thisTrial['prob'] < 4: # only in these trials the n_odd can match one of the 8 stimulus
         n_odd = np.random.randint(8)
    else:
         n_odd = 30

    # Initialize some default paratemers for this trial
    thisResp=[]
    rt_deci = 0 # np.array(0) # empty resp time
    resp = None
    col_resp = [1, 1, 1]
    trialClockStart = Clock.getTime()
    t_orient = orientations[ trials.thisN, ]
    for i_si in range(stim['ISI1_frames']): # first period before the first beep
        st.fixation(win, basic_stim)
        t = win.flip()
        if (i_si ==0): trial_times = np.append(trial_times, t)

    medf_s.play() # reproduce 1st auditory cue

    for i_si in range(stim['ISI2_frames']): # second period before the second beep
        st.fixation(win, basic_stim)
       # st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
        st.draw_contour(win,basic_stim)
        t = win.flip()
        if (i_si ==0): trial_times = np.append(trial_times, t)

    medf_s.play() # reproduce 2nd auditory cue (here it is always the same)
    
    for i_si in range(stim['ISI3_frames']): # second period before the contour
        st.fixation(win, basic_stim)
        #st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
        st.draw_contour(win,basic_stim)
        t = win.flip()
        if (i_si ==0): trial_times = np.append(trial_times, t)
        

    for i_si in range(stim['ISI4_frames']): # second period before the trial starts (response mapping cues)
        st.fixation(win, basic_stim)
        #st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)  
        st.draw_contour(win,basic_stim)
        t = win.flip()
        if (i_si ==0): trial_times = np.append(trial_times, t)

    # Draw the stimuli
    for istim in range(stim['nstim']+2): # 2 stim for the masks sandwiching
        event.clearEvents()
        if (istim == 0) | (istim == stim['nstim']+1):
           for frame in range(stim['stim_frames']):  # drawing stim frames
                if frame < stim['stim_frames'] - 1:  # the last frame should be empty
                    st.draw_mask(win, basic_stim)
                    st.fixation(win, basic_stim)
                   # st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                    st.draw_contour(win,basic_stim)
                    t = win.flip()
                    if (frame ==0): trial_times = np.append(trial_times, t)
                    respClockStart = Clock.getTime()
        else:            
            basic_stim['grating'].ori =  np.rad2deg(t_orient[istim-1]) # change orientation for each stim
            basic_stim['grating'].phase = np.random.rand()
            basic_stim['grating'].sf = stim['SF'] # reset SF (just in case previous stimulus was an oddball)
            
            if (istim == (n_odd+1)): # select stimulus containing oddball
                 basic_stim['grating'].sf = stim['SF'] - SF_difficulty # np.random.choice((-1, 1))*
                 
            for frame in range(stim['stim_frames']): # drawing stim frames
                if frame < stim['stim_frames'] - 1:  # the last frame should be empty
                    basic_stim['grating'].draw()
                    st.fixation(win, basic_stim)
                   # st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                    st.draw_contour(win,basic_stim)
                    t = win.flip()
                    if (frame ==0): trial_times = np.append(trial_times, t)
                        
            #st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)       
            st.draw_contour(win,basic_stim)
            st.fixation(win, basic_stim)
            t = win.flip()
                                      
    st.draw_mask(win, basic_stim)
    st.fixation(win, basic_stim)
    st.draw_contour(win,basic_stim)
    #st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)     
    t = win.flip()
          
    # Get responses stim['ITI_frames'] 
    for i_si in range(stim['ITI_frames']): # 
            #thisResp = exp.getResponse(win, ["space"], Clock) # you have to pass a clock function
            if thisResp == []:
                thisResp = event.getKeys(keyList=["space"], timeStamped=Clock)
                # print(thisResp)
            st.draw_mask(win, basic_stim)
            st.fixation(win, basic_stim)
            st.draw_contour(win,basic_stim)
            win.flip()  
   
    if (thisResp != []):
        rt_deci = thisResp[0][1] -  respClockStart
        rt_deci = np.around(rt_deci, decimals = 3)
        resp =  thisResp[0][0]
        print(rt_deci) # you can make a function of this to assign the correctness

        
    if (n_odd < 8 and thisResp != []) or ( n_odd > 7 and thisResp == []):
        correct = 1  # correct
        print("correct")
    else:
        print("incorrect")
        correct = -1  # incorrect
     
    for i_si in range(stim['feedback_frames']): # second period before the trial starts (response mapping cues)
        st.draw_mask(win, basic_stim)
        st.fixation(win, basic_stim)
        if correct == 1:
             basic_stim['fixation_point'].color = [-1, 1, -1] # colors are expressed as deviation from grey red -> [1,-1,-1]
        else:
             basic_stim['fixation_point'].color = [1, -1, -1]
            
        basic_stim['fixation_point'].draw()
        #st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
        st.draw_contour(win,basic_stim)
        t = win.flip()
        if (i_si ==0): trial_times = np.append(trial_times, t)
        #print('Overall, %i frames were dropped.' % win.nDroppedFrames)
    
    trial_times = np.array(trial_times ) 
    trialClocktimes = np.vstack([trialClocktimes, trial_times]) if trialClocktimes.size else trial_times  
    
    # Storing data in variables
    
    thr_trials_var.append([subj_id, thisTrial['prob'], n_odd, resp, correct, rt_deci])# saving conditions here
    thr_trials_ori.append(t_orient.tolist())
    #trialClocktimes.append(trial_times)
# Get datafiles in pandas format and attack to main Exp.variable



headers =  thr_trials_ori.pop(0)
block['trial_orientations'] = DataFrame(thr_trials_ori , columns=headers)

headers =  thr_trials_var.pop(0)
block['data'] = DataFrame(thr_trials_var , columns=headers)
block['block_duration'] =  Clock.getTime() -  BlockClockStart


# Saving block data
thisBlock = expInfo['locInfo']['block']

loc_exp['Exp_blocks'][thisBlock] =  block  # saving block data
expInfo['loc_exp'] =  loc_exp


print('Overall, %i frames were dropped.' % win.nDroppedFrames)
loc_frame_log = {}
loc_frame_log['droppedframes'] = win.nDroppedFrames
loc_frame_log['timings'] = np.diff(trialClocktimes,axis = 1)
expInfo['loc_frame_log'] =  loc_frame_log 

# closing everything
instr.end_experiment(win)
win.close()
# data saved
toFile(resultspath, expInfo) #saving file to disk

core.quit
