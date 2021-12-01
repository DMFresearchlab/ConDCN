#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Alexis Perez-Bellido
# Code for the confirmation bias experiment (using predcision code)


version = 1.4

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


#prefs.hardware['audioLib']=['pyo'] # use Pyo audiolib for good temporal resolution
#from psychopy.sound import Sound # This should be placed after changing the audio library
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
monitor_features['Hz'] = 144 #'auto' #60 # this can be set to "auto" to estimate the refreshing rate of the monitor, although it can fail often

   
win, monitor_features = exp.create_window(monitor_features)
ifi = monitor_features['ifi']

expInfo['monitor_features'] = monitor_features


# Experiment timings and characteristics 
stim = st.stim_config(ifi) #Loading stim characteristics

expInfo['stim'] = stim # saving
      
basic_stim = st.draw_basic(win,expInfo['stim'])

# Create the sounds that I am going to use
#lowf_s = Sound(600, sampleRate=44100, secs=0.1, stereo=True ,loops=0,hamming=True)
#medf_s = Sound(800, sampleRate=44100, secs=0.1, stereo=True ,loops=0,hamming=True)
#highf_s = Sound(1000, sampleRate=44100, secs=0.1, stereo=True ,loops=0,hamming=True)


# Experimental condition preparation
# response mapping (shuffled in each experiment )
expInfo['resp_maps']  = np.array([0, 45]) # 0 -> cardinal; 45 -> diagonal


# set path of stim file. Load the oriented trials dataset
stimfile = os.path.join(os.getcwd(),'stim_matrix') # os.sep

if os.path.isfile(stimfile + ".npy"): # check if file exist, otherwise, create a new one.
    orientations = np.load(stimfile + ".npy")
else: # This is the procedure to generate the trials matrix. I RECOMEND TO RUN THIS MANUALLY TO VISUALIZE THE RESULTS ONLINE
    from create_stimuli import stim_creation
    orientations = stim_creation(stim['nstim'], stimfile)

# collapse the circular angles on one side of the circle to make the analyses easier (see that 45 is equal to 225
# degrees when drawing grating orientations)
orientations[orientations < 0] += np.deg2rad(180)
orientations = np.around(orientations, decimals = 3)
decision_var = signal.sawtooth(4 * ((orientations)), 0.5)  # DU decision variable
decision_var_T = np.mean(decision_var, 1)
decision_var_std = np.std(decision_var, 1)
orientation_var_std = stats.circstd((orientations), axis=1)


# Show some instructions 
instr.init_instructions(win, basic_stim, ifi)
instr.block_loc_start(win)

# Create an experiment clock
Clock = core.Clock()
ExpClockStart = Clock.getTime() # global experiment time

#trials = data.TrialHandler(stimList, ntrials_per_cond, method='random') # when calling next trial it continues to the next randomly ordered trial

guess = expInfo['subjInfo']['guess']
black_resps = np.array([[-1, -1, -1],[-1, -1, -1]]) # default color resp options

stepsize = [0.3/3] # according to paper is SDT

staircase = data.StairHandler(startVal = guess,
                          stepType = 'lin', stepSizes=np.array(stepsize), # this determines the number of reversals and therefore the number of trials
                          nUp=1, nDown=2,  # will home in on the 80% threshold
                          minVal = 0, maxVal = 0.6,
                          nTrials= 20)


thr_trials_var = [['subj','cond', 'DV', 'resp', 'r_map', 'correct', 'RT']] # saving conditions here
thr_trials_ori = [['o1','o2','o3','o4','o5','o6']]
trialClocktimes = np.array([]) # saving whole times here

    
for thisIncrement in staircase:  # will continue the staircase until it terminates!
    trial_times = np.array([]) # logging trial event time stamps
    cond = np.random.choice([-1, 1])  # will be cardinal or diagonal
    basic_stim['fixation_point'].color = [1, 1, 1]
    ExpClockTrial = Clock.getTime()
    # decision variable in this trial
    x =  cond*(thisIncrement)
    print(x)
   # x = -0.5
    #sel_trials = decision_var_T[np.where((decision_var_T > x-0.001) & (decision_var_T < x+0.001))] # to return an array and not tuple
    sel_trials = np.where((decision_var_T > x-0.025) & (decision_var_T < x+0.025))
    trial_sel_idx = np.random.choice(sel_trials[0]) # selecting orientation vector for this trial.
    t_orient = orientations[trial_sel_idx,]
    decision_avg = decision_var_T[trial_sel_idx]

    # Initialize some default paratemers for this trial
    thisResp=None
    col_resp = [1, 1, 1]
    trialClockStart = Clock.getTime()
  
    np.random.shuffle(expInfo['resp_maps']) # Shuffling resp_map -> first response option (cardignal or diagonal) will be placed at left, and second at right
    
    stim['ISI1_frames'] = round(np.random.randint(400,700)/ifi)  # fixation alone
    for i_si in range(stim['ISI1_frames']): # first period before the first beep
        st.fixation(win, basic_stim)
        t = win.flip()
        if (i_si ==0): trial_times = np.append(trial_times, t)

    # medf_s.play() # reproduce 1st auditory cue

    for i_si in range(stim['ISI2_frames']): # second period before the second beep
        st.fixation(win, basic_stim)
       # st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
        st.draw_contour(win,basic_stim)
        t = win.flip()
        if (i_si ==0): trial_times = np.append(trial_times, t)

    # medf_s.play() # reproduce 2nd auditory cue (here it is always the same)

    # Draw the stimuli
    for istim in range(stim['nstim']+2): # 2 stim for the masks sandwiching
        event.clearEvents()
        if (istim == 0) | (istim == stim['nstim']+1):
           for frame in range(stim['stim_frames']):  # drawing stim frames
                if frame < stim['stim_frames'] - 1:  # the last frame should be empty
                    st.draw_mask(win, basic_stim)
                    st.fixation(win, basic_stim)
                 #   st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                    st.draw_contour(win,basic_stim)
                    t = win.flip()
                    if (frame ==0): trial_times = np.append(trial_times, t)
                    respClockStart = Clock.getTime()
        else:
            basic_stim['grating'].ori =  np.rad2deg(t_orient[istim-1]) # change orientation for each stim
            basic_stim['grating'].phase = np.random.rand()

            for frame in range(stim['stim_frames']): # drawing stim frames
                if frame < stim['stim_frames'] - 1:  # the last frame should be empty
                    basic_stim['grating'].draw()
                    st.fixation(win, basic_stim)
                  #  st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                    st.draw_contour(win,basic_stim)
                    t = win.flip()
                    if (frame ==0): trial_times = np.append(trial_times, t)
          #  st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)       
            st.draw_contour(win,basic_stim)
            st.fixation(win, basic_stim)
            t = win.flip()
    
    for i_si in range(stim['ISI3_frames']): # second period before the second beep
                st.draw_mask(win, basic_stim)
                st.fixation(win, basic_stim)
               # st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                st.draw_contour(win,basic_stim)
                t = win.flip()
                trial_times = np.append(trial_times, t)
                
    # Mao response         
    if expInfo['resp_maps'][0] == 0:
        resp_pos = 1
    else:
        resp_pos = -1
        
         
    ticks = [-1,1] #['1','2','3','4','5']
    labels = ['Cardinal','','Diagonal'] # labels = labels,
    deci_slider = visual.Slider(win,ticks= ticks,  pos = [0,0],  
                              size = (6,1), style = 'radio', granularity = 1)
    

    while not deci_slider.rating:
        st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps* -1)     
        deci_slider.draw()  
        #visual.Line(win=win,lineWidth = 20,units="deg",lineColor=[1, 1, 1],start = [0, -1],end = [0, 1]).draw()
        st.fixation(win, basic_stim)
        win.flip()
       
    resp_ang = deci_slider.rating * resp_pos
    rt_deci = deci_slider.getRT()
      
    
    win.mouseVisible = False 
     
    if (x > 0 and resp_ang > 0) or (x < 0 and  resp_ang < 0):
        correct = 1  # correct
        print("correct")
        steps = staircase.stepSizes # change stepsize as suggested in Miguel Angel Perez paper
        steps = np.array(steps)
        staircase.stepSizes = 0.871*steps
    else:
        print("incorrect")
        correct = -1  # incorrect
    # trial_perform = np.append(trial_perform, correct) 
     
    staircase.addResponse(correct) # adding information to staircase

    for i_si in range(stim['feedback_quest_frames']): # second period before the trial starts (response mapping cues)
        #st.draw_mask(win, basic_stim)
        st.fixation(win, basic_stim)
        if correct == 1:
             basic_stim['fixation_point'].color = [-1, 1, -1] # colors are expressed as deviation from grey red -> [1,-1,-1]
        else:
             basic_stim['fixation_point'].color = [1, -1, -1]
            
        basic_stim['fixation_point'].draw()
        #st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
        #st.draw_contour(win,basic_stim)
        t = win.flip()
        if (i_si ==0): trial_times = np.append(trial_times, t)
        #print('Overall, %i frames were dropped.' % win.nDroppedFrames)
    
    trial_times = np.array(trial_times ) 
    trialClocktimes = np.vstack([trialClocktimes, trial_times]) if trialClocktimes.size else trial_times  
    
    # Storing data in variables
    
    thr_trials_var.append([subj_id, cond, x, resp_ang, expInfo['resp_maps'][0], correct, rt_deci])# saving conditions here
    thr_trials_ori.append(t_orient.tolist())
    #trialClocktimes.append(trial_times)
# Get datafiles in pandas format and attack to main Exp.variable

approxThreshold = np.average(staircase.reversalIntensities[-5:])

expInfo['subjInfo']['guess'] = approxThreshold  # save threshold for next experiment

headers =  thr_trials_ori.pop(0)

stc_exp = {}
stc_exp['staircase'] = staircase

stc_exp['trial_orientations'] = DataFrame(thr_trials_ori , columns=headers)

headers =  thr_trials_var.pop(0)
stc_exp['data'] = DataFrame(thr_trials_var , columns=headers)

expInfo['staircase_exp'] =  stc_exp



plt.plot(np.arange(len(staircase.intensities)),staircase.intensities,linewidth=1, color = 'black')
#plt.ylim((0,0.6))
plt.axhline(y=approxThreshold, linewidth=1, color='r')
plt.xlabel('Decision variable')
plt.xlabel('N trial')

print('Overall, %i frames were dropped.' % win.nDroppedFrames)
stc_frame_log = {}
stc_frame_log['droppedframes'] = win.nDroppedFrames
stc_frame_log['timings'] = np.diff(trialClocktimes,axis = 1)
#stc_frame_log['timings'] = (trialClocktimes)
expInfo['stc_frame_log'] =  stc_frame_log 

# closing everything
instr.end_experiment(win)
win.close()
# data saved
toFile(resultspath, expInfo) #saving file to disk

print(staircase.reversalIntensities)
print('The approximate threshold is %f' % approxThreshold)

core.quit

