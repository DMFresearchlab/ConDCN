#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Alexis Perez-Bellido
# Code for the confirmation bias experiment (using predcision code)


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


sst = False # Using parallel port to send triggers

if sst: 
    p_port = serial.Serial('COM3', 115200, timeout = 0) # this is for windows
    p_port.write(b'00')
    core.wait(0.2)
    p_port.write(b'RR')

# Trigger keys

tg_mblock =  '07'
tg_mtrial = '01'
tg_stim = '13'
tg_resp = '05'
tg_zero = '00'
#prefs.hardware['audioLib']=['pyo'] # use Pyo audiolib for good temporal resolution
#from psychopy.sound import Sound # This should be placed after changing the audio library
# monitors.getAllMonitors()
#from general_settings import * # reads the variables in one script

# Collect subject data
expInfo, resultspath = exp.mainexp_subject_info(version)
subj_id = expInfo['subjInfo']['observer']

# Loading monitor definitions
monitores = st.monitor_def() 

mon, expInfo['monitor'] = exp.define_monitor(monitores[4]) # select the correct monitor

# Creating a new experimental window
monitor_features = {}
monitor_features['monitor'] = mon
monitor_features['units'] = 'deg' # units to define your stimuli
monitor_features['screen_id'] = 0 # when using a extended display 
monitor_features['full'] = True
monitor_features['Hz'] = 'auto' #60 # this can be set to "auto" to estimate the refreshing rate of the monitor, although it can fail often

   
win, monitor_features = exp.create_window(monitor_features)
ifi = monitor_features['ifi']

# Experiment timings and characteristics 
stim = st.stim_config(ifi) #Loading stim characteristics


      
basic_stim = st.draw_basic(win,stim)

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
instr.main_instructions(win)

# Create an experiment clock
Clock = core.Clock()
ExpClockStart = Clock.getTime() # global experiment time

#trials = data.TrialHandler(stimList, ntrials_per_cond, method='random') # when calling next trial it continues to the next randomly ordered trial

guess = expInfo['subjInfo']['guess']
black_resps = np.array([[-1, -1, -1],[-1, -1, -1]]) # default color resp options

stepsize =  [0.3/3] # according to paper is SDT - other options [0.15, 0.1, 0.1, 0.05, 0.05] #

staircase = data.StairHandler(startVal = guess,
                          stepType = 'lin', stepSizes=stepsize, # this determines the number of reversals and therefore the number of trials
                          nUp=1, nDown=2,  # will home in on the 80% threshold
                          minVal = 0.01, maxVal = 0.6,
                          nTrials= 1000) # 40

main_exp  = {}
main_exp['nblocks']     = 4 # 4 # totaltime = 90 * 6 * 5
main_exp['Exp_blocks']  = [None] * main_exp['nblocks'] # assigning memory for storing block data
main_exp['trial_reps']  = 40 
   



for thisBlock in range(main_exp['nblocks']): # iterate over blocks
    instr.block_start(win)
    # lets trigger the beggining of the experiment
    if sst: win.callOnFlip(p_port.write, tg_mblock.encode())
    win.flip()
    if sst: win.callOnFlip(p_port.write, tg_zero.encode())
    win.flip()
    
    BlockClockStart = Clock.getTime() # block experiment time
    block = {} # dummy variable to save block data
    correct_seq = np.array([]) # saving seq. of correct responses per trial sequence
    trial_rep = 0 # update staricase after each iteraction
    thr_trials_var = [['subj','nblock', 'ntrial', 'nrep', 'cond', 'DV', 'resp', 'r_map', 'correct', 'RT']] # saving conditions here
    thr_trials_ori = [['o1','o2','o3','o4','o5','o6']]
    trialClocktimes = np.array([]) # saving whole times here
    
    for iTrial in range(main_exp['trial_reps']):  # will continue the staircase until it terminates!
        thisIncrement = np.around(next(staircase),decimals = 5)
        cond = np.random.choice([-1, 1])  # will be cardinal or diagonal
        ExpClockTrial = Clock.getTime()
        # decision variable in this trial
        x =  cond*(thisIncrement)
        print(x)
       # x = -0.5
        #sel_trials = decision_var_T[np.where((decision_var_T > x-0.001) & (decision_var_T < x+0.001))] # to return an array and not tuple
        sel_trials = np.where((decision_var_T > x-0.025) & (decision_var_T < x+0.025))
        trial_sel_idx = np.random.choice(sel_trials[0]) # selecting orientation vector for this trial.
        t_orient = orientations[trial_sel_idx]
        decision_avg = decision_var_T[trial_sel_idx] 
        
        fixation_color = [1, 1, 1]
        trial_perform = np.array([])
        for i_rep in range(stim['nreps']):
            # Initialize some default paratemers for this trial
            thisResp=None
            basic_stim['fixation_point'].color = fixation_color
            trial_times = np.array([]) # logging trial event time stamps          
            col_resp = [1, 1, 1]
            trialClockStart = Clock.getTime()        
            np.random.shuffle(expInfo['resp_maps']) # Shuffling resp_map -> first response option (cardignal or diagonal) will be placed at left, and second at right
            
            stim['ISI1_frames'] = round(np.random.randint(800,900)/ifi)  # fixation alone
            for i_si in range(stim['ISI1_frames']): # first period before the first beep
                st.fixation(win, basic_stim)
                basic_stim['feedback1'].color = [0.5, 0.5, 0.5]
                basic_stim['feedback1'].autoDraw = True
                basic_stim['feedback1'].draw()  
                
                if i_rep == 1:
                     basic_stim['feedback2'].color = [0.5, 0.5, 0.5]
                     basic_stim['feedback2'].autoDraw = True
                     basic_stim['feedback2'].draw()
                if i_rep == 2:
                     basic_stim['feedback3'].color = [0.5, 0.5, 0.5]
                     basic_stim['feedback3'].autoDraw = True
                     basic_stim['feedback3'].draw()
                     
                t = win.flip()
                
                if (i_si ==0):
                    basic_stim['white_test'].draw()  
                    trial_times = np.append(trial_times, t)
        
            # medf_s.play() # reproduce 1st auditory cue
        
            for i_si in range(stim['ISI2_frames']): # second period before the second beep
                st.fixation(win, basic_stim)
               # st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                st.draw_contour(win,basic_stim)
                if i_si == 0:
                    basic_stim['white_test'].draw()  
                    if sst: win.callOnFlip(p_port.write, tg_mtrial.encode()) # send trigger
                if i_si == 1:     
                    if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
                t = win.flip()
                if (i_si ==0): trial_times = np.append(trial_times, t)
        
            # medf_s.play() # reproduce 2nd auditory cue (here it is always the same)
        
            # Draw the stimuli
            for istim in range(stim['nstim']+2): # 2 stim for the masks sandwiching
                event.clearEvents()
                if (istim == 0) | (istim == stim['nstim']+1): # if first or last mask
                   for frame in range(stim['stim_frames']):  # drawing stim frames
                        if frame == 0:
                            if sst: win.callOnFlip(p_port.write, tg_stim.encode()) # send trigger
                        if frame == 1:     
                            if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
                        if (frame == stim['stim_frames']) and (istim == 0):  # the last frame should be emptyin the first mask
                           # -1 and istim == 0
                            st.fixation(win, basic_stim)
                            st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                            st.draw_contour(win,basic_stim)
                            t = win.flip()                    
                        else: # flip empty frame
                            st.draw_mask(win, basic_stim)
                            st.fixation(win, basic_stim)
                            st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                            st.draw_contour(win,basic_stim)
                            if frame == 0:
                                basic_stim['white_test'].draw()  
                            t = win.flip()
                        if (frame ==0): 
                            
                            trial_times = np.append(trial_times, t)
                else:
                    basic_stim['grating'].ori =  np.rad2deg(t_orient[istim-1]) # change orientation for each stim
                    basic_stim['grating'].phase = np.random.rand()
        
                    for frame in range(stim['stim_frames']): # drawing stim frames
                        if frame == 0:
                            if sst: win.callOnFlip(p_port.write, tg_stim.encode()) # send trigger
                        if frame == 1:     
                            if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
                        if frame < stim['stim_frames'] - 1:  # the last frame should be empty
                            basic_stim['grating'].draw()
                            st.fixation(win, basic_stim)
                            st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                            st.draw_contour(win,basic_stim)
                            if frame == 0:
                                basic_stim['white_test'].draw()  
                            t = win.flip()
                            if (frame ==0): 
                                basic_stim['white_test'].draw()  
                                trial_times = np.append(trial_times, t)
                        else: # flip empty frame
                            st.fixation(win, basic_stim)
                            st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                            st.draw_contour(win,basic_stim)
                            if frame == 0:
                                basic_stim['white_test'].draw()  
                            t = win.flip()               
            
            respClockStart = Clock.getTime()                       
            st.draw_mask(win, basic_stim)
            st.fixation(win, basic_stim)
            st.draw_contour(win,basic_stim)
            st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)     
            t = win.flip()
            trial_times = np.append(trial_times, t)
                           
            # Get responses
            while thisResp == None:
                    thisResp = exp.getResponse(win, ["z", "m"], Clock) # you have to pass a clock function
        
            rt_deci = thisResp[0] -  respClockStart
            rt_deci = np.around(rt_deci, decimals = 3)
            print(rt_deci) # you can make a function of this to assign the correctness
            
            if thisResp[1] == 'z':
                resp_ang = expInfo['resp_maps'][0] # assign response selected
                
            elif thisResp[1] == 'm':
                resp_ang = expInfo['resp_maps'][1]

            if (x > 0 and resp_ang == 45) or (x < 0 and  resp_ang == 0):
                correct = 1  # correct
                print("correct")
                if i_rep == 0: # update staircase only if is the first stim presentation
                    steps = staircase.stepSizes # change stepsize as suggested in Miguel Angel Perez paper
                    steps = np.array(steps)
                    staircase.stepSizes = 0.871*steps
                    
            else:
                print("incorrect")
                correct = -1  # incorrect
            trial_perform = np.append(trial_perform, correct) 
            
                
            for i_si in range(stim['wait_feedback_frames']): # wait time for feedback and send a couple of triggers
                if sst:
                    if i_si == 0: 
                        win.callOnFlip(p_port.write, tg_resp.encode()) # send trigger
                        basic_stim['white_test'].draw()  
                    if i_si == 1: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins           
                basic_stim['fixation_point'].draw()
                st.draw_contour(win,basic_stim)
                t = win.flip()
                if (i_si ==0): trial_times = np.append(trial_times, t)
                
              
            basic_stim['fixation_point'].draw()
            st.draw_contour(win,basic_stim)
            win.flip()   
            #print(correct)
            print(trial_perform)
            if i_rep == 0: staircase.addResponse(correct) # adding information to staircase
            if i_rep == 2:
                for i_si in range(stim['feedback_frames']): # time added between trials
                    st.draw_mask(win, basic_stim)
                    st.fixation(win, basic_stim)

                    if trial_perform[0] == 1:
                        basic_stim['feedback1'].color = [0.0, 1.0, 0.0]
                    else:
                        basic_stim['feedback1'].color = [1.0, 0.0, 0.0]
                    basic_stim['feedback1'].autoDraw = False
                    basic_stim['feedback1'].draw() 
                    
                    if trial_perform[1] == 1:
                        basic_stim['feedback2'].color = [0.0, 1.0, 0.0]
                    else:
                        basic_stim['feedback2'].color = [1.0, 0.0, 0.0]
                    basic_stim['feedback2'].autoDraw = False
                    basic_stim['feedback2'].draw()   
                    
                    if trial_perform[2] == 1:
                        basic_stim['feedback3'].color = [0.0, 1.0, 0.0]
                    else:
                        basic_stim['feedback3'].color = [1.0, 0.0, 0.0]
                    basic_stim['feedback3'].autoDraw = False
                    basic_stim['feedback3'].draw() 
                    
                    basic_stim['fixation_point'].draw()
                    #st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                    st.draw_contour(win,basic_stim)
                    t = win.flip()
                    #print('Overall, %i frames were dropped.' % win.nDroppedFrames)
            
            trial_times = np.array(trial_times) 
            trialClocktimes = np.vstack([trialClocktimes, trial_times]) if trialClocktimes.size else trial_times  
            # if there is an outlier trial, stop experiment
            this_times = np.diff(trial_times)
            mean_stim_time = np.mean(this_times[4:8]) * 1000
            if mean_stim_time > stim['stim_time'] + 75 or mean_stim_time < stim['stim_time'] - 75:
                print('Bad timing!! Cerrando programa')
                win.close()
                core.quit
            # Storing data in variables           
            thr_trials_var.append([subj_id, thisBlock , iTrial, i_rep, cond, x, thisResp[1], expInfo['resp_maps'][0], correct, rt_deci])# saving conditions here
            thr_trials_ori.append(t_orient.tolist())
            
    # Get datafiles in pandas format and attack to main Exp.variable
    headers =  thr_trials_var.pop(0)
    block['data'] = DataFrame(thr_trials_var, columns=headers)
    headers =  thr_trials_ori.pop(0)
    block['trial_orientations'] = DataFrame(thr_trials_ori , columns=headers)
    block['block_duration'] =  Clock.getTime() -  BlockClockStart
    
    main_frame_log = {}
    main_frame_log['droppedframes'] = win.nDroppedFrames
    main_frame_log['timings'] = np.diff(trialClocktimes,axis = 1)
    block['main_frame_log'] = main_frame_log 
    
    main_exp['Exp_blocks'][thisBlock] =  block  # saving block data
    toFile(resultspath, expInfo) #saving file to disk
            #trialClocktimes.append(trial_times)
        # Get datafiles in pandas format and attack to main Exp.variable

approxThreshold = np.average(staircase.reversalIntensities[-5:])

main_exp['monitor_features'] = monitor_features
main_exp['stim'] = stim

expInfo['subjInfo']['guess'] = approxThreshold  # save threshold for next experiment

headers =  thr_trials_ori.pop(0)

expInfo['main_exp'] =  main_exp

print('Overall, %i frames were dropped.' % win.nDroppedFrames)


toFile(resultspath, expInfo) #saving file to disk

instr.end_experiment(win)

# closing everything
win.close()

core.quit

