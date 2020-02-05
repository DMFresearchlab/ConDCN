#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 11:15:02 2019
Main experiment predictive cardinal diagonal orientations. 
In this second part of the experiment we should load the expInfo variable from first staircase part to update threshold and participant data.
Everything should be saved together
# Verify visual angles from here http://elvers.us/perception/visualAngle/

@author: Alexis Pérez-Bellido 2019
"""

version = 1.3

import numpy as np
import os
import matplotlib.pyplot as plt
import exp_func as exp
import stimuli as st
import instructions as instr
import serial

from psychopy import visual, logging, core, event,  gui, data, monitors
from psychopy.tools.filetools import fromFile, toFile # wrappers to save pickles
from scipy import signal, stats
from psychopy.preferences import prefs
from pandas import DataFrame

sst = False # Using parallel port to send triggers

if sst: 
    p_port = serial.Serial('COM3', 115200, timeout = 0) # this is for windows
    p_port.write(b'00')
    core.wait(0.2)
    p_port.write(b'RR')
    
    
# Some general presets
event.globalKeys.clear() # implementing a global event to quit the program at any time by pressing ctrl+q
event.globalKeys.add(key='q', modifiers=['ctrl'], func=core.quit)


prefs.hardware['audioLib']=['pyo'] # use Pyo audiolib for good temporal resolution
from psychopy.sound import Sound # This should be placed after changing the audio library
# monitors.getAllMonitors()
#from general_settings import * # reads the variables in one script

#subj_id    = input("Participant ID: ")
# Collect subject data
expInfo, resultspath = exp.mainexp_subject_info(version)
subj_id = expInfo['subjInfo']['observer']

# Loading monitor definitions
monitores = st.monitor_def() 

mon, expInfo['monitor'] = exp.define_monitor(monitores[3]) # select the correct monitor


# Creating a new experimental window
monitor_features = {}
monitor_features['monitor'] = mon
monitor_features['units'] = 'deg' # units to define your stimuli
monitor_features['screen_id'] = 0
monitor_features['full']  = True
monitor_features['Hz'] = 60  # 'auto'this can be set to "auto" to estimate the refreshing rate of the monitor, although it can fail often

tg_mblock =  '07'
tg_mtrial = '01'
tg_stim = '13'
tg_resp = '05'
tg_zero = '00'

win, monitor_features = exp.create_window(monitor_features)
ifi = monitor_features['ifi']

stim = st.stim_config(ifi) #Loading stim characteristics


basic_stim = st.draw_basic(win, stim)

# Experimental condition preparation
# response mapping (shuffled in each experiment ) LOAD FROM STAIRCASE PART
expInfo['resp_maps']  = np.array([0, 45]) # 0 -> cardinal; 45 -> diagonal
np.random.shuffle(expInfo['resp_maps'] ) # first response option (cardignal or diagonal) will be placed at left, and second at right 
 

# set path of stim file. Load the oriented trials dataset
stimfile = os.path.join(os.getcwd(),'stim_matrix') # os.sep

if os.path.isfile(stimfile + ".npy"): # check if file exist, otherwise, create a new one.
    orientations = np.load(stimfile + ".npy")
else: # This is the procedure to generate the trials matrix. I RECOMEND TO RUN THIS MANUALLY TO VISUALIZE THE RESULTS ONLINE
    from create_stimuli import stim_creation
    orientations = stim_creation(stim['nstim'], stimfile)


## Create the sounds that I am going to use
    
lowf_s = Sound(400, sampleRate=44100, secs=0.1, stereo=True ,loops=0) # ,hamming=True
medf_s = Sound(800, sampleRate=44100, secs=0.1, stereo=True ,loops=0)
highf_s = Sound(1200, sampleRate=44100, secs=0.1, stereo=True ,loops=0)

black_resps = np.array([[-1, -1, -1],[-1, -1, -1]]) # default color resp options


main_exp  = {}
main_exp['nblocks']     = 4 # 4 # totaltime = 90 * 6 * 5
main_exp['trial_reps']  = 5 # 1 # 1 bloque de 1 repeticion son 90 s aprox. 5
main_exp['Exp_blocks']  = [None] * main_exp['nblocks'] # assigning memory for storing block data

# shuffling beep and prediction 
predsounds  = [None] * 2
predsounds[0] = lowf_s 
predsounds[1] = highf_s 

main_exp['sound_map']  = np.array([0,1]) # Always 0 pos is C and 1 pos is D: [0 1] -> low(c) High(d) / 1 0 -> high(c) low(d) 
np.random.shuffle(main_exp['sound_map']) #
# Apply order 
predsounds = [predsounds[i] for i in main_exp['sound_map']]

# collapse the circular angles on one side of the circle to make the analyses easier (see that 45 is equal to 225
# degrees when drawing grating orientations)
orientations[orientations < 0] += np.deg2rad(180)
orientations = np.around(orientations, decimals = 3)
decision_var = signal.sawtooth(4 * ((orientations)), 0.5)  # DU decision variable
decision_var_T = np.mean(decision_var, 1)
decision_var_std = np.std(decision_var, 1)
orientation_var_std = stats.circstd((orientations), axis=1)

# instructions
instr.main_instructions(win)


# Create an experiment clock
Clock = core.Clock()

stepsize = [0.2/3]
# I create a staircase that I will update everytime after each cycle of conditions during the whole experiment
staircase = data.StairHandler(startVal = expInfo['subjInfo']['guess'],
                          stepType = 'lin', stepSizes=stepsize, # this determines the number of reversals and therefore the number of trials
                          nUp=1, nDown=2,  # will home in on the 80% threshold
                          minVal = 0.0, maxVal = 0.6,
                          nTrials= 1000)


# Experiment design
stimList = []
for cond in ['CP', 'NP', 'DP']:#
    for props in [1, 2, 3, 4]: # 75% vs 25% larger values than 1 in the CP & DP conditions correspond to one orientation
            stimList.append({'Pred':cond, 'prob': props}) #, 'n_reps': n_reps



# win.getMovieFrame()   # to print the screen in terminal and save it
# win.saveMovieFrames('


corr_lotery = []



for thisBlock in range(main_exp['nblocks']): # iterate over blocks
    instr.block_start(win)    
    # lets trigger the beggining of the experiment
    if sst: win.callOnFlip(p_port.write, tg_mblock.encode())
    win.flip()
    if sst: win.callOnFlip(p_port.write, tg_zero.encode())
    win.flip()
    
    BlockClockStart = Clock.getTime() # block experiment time
    block = {} # dummy variable to save block data
    trialClocktimes = np.array([]) # saving whole times here
    correct_seq = np.array([]) # saving seq. of correct responses per trial sequence
    trial_rep = 0 # update staricase after each iteraction
    thr_trials_var = [['subj','nblock', 'ntrial', 'pred','cond', 'prop', 'DV', 'resp', 'correct', 'RT']] # saving conditions here (variables names must match number of saved variables)
    thr_trials_ori = [['o1','o2','o3','o4','o5','o6','o7','o8']]

    trials = data.TrialHandler(stimList, main_exp['trial_reps'], method='random') # when calling next trial it continues to the next randomly ordered trial
    #staircase.stepSizes = stepsize  # restart stepsize in each block
    
    for thisTrial in trials: # iterate over trials
        #thisTrial = next(trials)
        trial_times = np.array([]) # logging trial event time stamps
        basic_stim['fixation_point'].color = [1, 1, 1] # changing back fixation to white (after feedback color)
        
        # Updating staircase if it is a new trial_rep
        if trials.thisRepN == trial_rep:  # update threshold ONLY after all the trials combinations have been displayed 
            if trial_rep > 0: # update after the first trial repetion
                stc_corr_upt = -1 # lets decide whether we update as correct or incorrect based on averaged responses during n trial repetion           
                
                if np.mean(correct_seq) > 0.5:
                    stc_corr_upt = 1
                    
                if np.mean(correct_seq) == 0.5: # randomize if it is 0.5
                    stc_corr_upt = np.random.choice([-1, 1])
                            
                if stc_corr_upt == 1: # correct response, change stepsize accordingly with m.a. garcia-perez Y-N staircases paper
                    steps = staircase.stepSizes
                    steps = np.array(steps)
                    staircase.stepSizes = 0.871*steps
         
                staircase.addResponse(stc_corr_upt) # adding information to staircase
                print('updating staircase with ', str(stc_corr_upt))
                correct_seq = np.array([])  # restart correct counter
                
            thisIncrement = np.around(next(staircase),decimals = 5)
            trial_rep = trial_rep + 1  # move flag to next trial rep
            
        # assigning orientation conditional to pred variable
        if thisTrial['Pred'] == 'NP':
            resp_cue = np.array([[0.3, 0.3, 0.3],[0.3, 0.3, 0.3]])
            if thisTrial['prob'] > 2:
                cond = 1
            else:
                cond = -1
                
        if thisTrial['Pred'] == 'CP':
            # assigning color cues conditional to probability- 1st row assign the color to left option and 2nd right option color(white is the most likely option. )
            resp_cue = np.array([[-1, -1, -1],[0.3, 0.3, 0.3]]) if expInfo['resp_maps'][0] == 45 else np.array([[0.3, 0.3, 0.3],[-1, -1, -1]])
            if thisTrial['prob'] > 1:
                cond = -1
            else:
                cond = 1
                
        if thisTrial['Pred'] == 'DP':
            resp_cue = np.array([[0.3, 0.3, 0.3],[-1, -1, -1]]) if expInfo['resp_maps'][0] == 45 else np.array([[-1, -1, -1],[0.3, 0.3, 0.3]])
            if thisTrial['prob'] > 1:
                cond = 1            
            else:
                cond = -1
     
        condlab = 'D' if cond == 1 else 'C' # assign label to cond
        
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
    
        stim['ISI1_frames'] = round(np.random.randint(650,850)/ifi)
        for i_si in range(stim['ISI1_frames']): # first period before the first beep
            st.fixation(win, basic_stim)
            t = win.flip()
            if (i_si ==0): trial_times = np.append(trial_times, t)
    
        medf_s.play() # reproduce 1st auditory cue

        
        for i_si in range(stim['ISI2_frames']): # second period before the second beep
            st.fixation(win, basic_stim)            
            st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
            st.draw_contour(win,basic_stim)
            t = win.flip()
            if (i_si ==0): trial_times = np.append(trial_times, t)
        
        # Auditory cue
        if thisTrial['Pred'] == 'NP':
            medf_s.play() # 
        if thisTrial['Pred'] == 'CP':
            predsounds[0].play()
        if thisTrial['Pred'] == 'DP':
            predsounds[1].play()
        
        
        for i_si in range(stim['ISI3_frames']): # blinking spatial cue
            if i_si == 0:     
                if sst: win.callOnFlip(p_port.write, tg_mtrial.encode()) # send trigger
            if i_si == 1:     
                if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins            
            st.fixation(win, basic_stim)
            #st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
            st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'],resp_cue ) 
            st.draw_contour(win,basic_stim)
            t = win.flip()
            if (i_si ==0): trial_times = np.append(trial_times, t)
            
    
        for i_si in range(stim['ISI4_frames']): # third period before the trial starts (response mapping cues)
            st.fixation(win, basic_stim)
            st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
            st.draw_contour(win,basic_stim)
            t = win.flip()
            if (i_si ==0): trial_times = np.append(trial_times, t)
    
        # Draw the stimuli
        for istim in range(stim['nstim']+2): # 2 stim for the masks sandwiching
            if i_si == 0:     
                if sst: win.callOnFlip(p_port.write, tg_stim.encode()) # send trigger
            if i_si == 1:     
                if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
            event.clearEvents()
            if (istim == 0) | (istim == stim['nstim']+1): # if first or last mask
               for frame in range(stim['stim_frames']):  # drawing stim frames
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
                        t = win.flip()
                    if (frame ==0): trial_times = np.append(trial_times, t)
            else:
                basic_stim['grating'].ori =  np.rad2deg(t_orient[istim-1]) # change orientation for each stim
                basic_stim['grating'].phase = np.random.rand()
    
                for frame in range(stim['stim_frames']): # drawing stim frames
                    if frame < stim['stim_frames'] - 1:  # the last frame should be empty
                        basic_stim['grating'].draw()
                        st.fixation(win, basic_stim)
                        st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                        st.draw_contour(win,basic_stim)
                        t = win.flip()
                        if (frame ==0): trial_times = np.append(trial_times, t)
                    else: # flip empty frame
                        st.fixation(win, basic_stim)
                        st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                        st.draw_contour(win,basic_stim)
                        t = win.flip()
                #st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)       
                #st.draw_contour(win,basic_stim)
                #st.fixation(win, basic_stim)
                #t = win.flip()
        respClockStart = Clock.getTime()                                       
        st.draw_mask(win, basic_stim)
        st.fixation(win, basic_stim)
        st.draw_contour(win,basic_stim)
        st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)      
        t = win.flip()
        trial_times = np.append(trial_times, t)
        
        
        # Get responses
        while thisResp == None:
                thisResp = exp.getResponse(win, ["c", "m"], Clock) # you have to pass a clock function
    
        rt_deci = thisResp[0] -  respClockStart
        rt_deci = np.around(rt_deci, decimals = 3)
        print(rt_deci) # you can make a function of this to assign the correctness
        
        if thisResp[1] == 'c':
            resp_ang = expInfo['resp_maps'][0] # assign response selected
            
        elif thisResp[1] == 'm':
            resp_ang = expInfo['resp_maps'][1]
            
        if (x > 0 and resp_ang == 45) or (x < 0 and  resp_ang == 0):
            correct = 1  # correct
            print("correct")                        
        else:
            print("incorrect")
            correct = 0  # incorrect
         
    
        for i_si in range(stim['feedback_frames']): # second period before the trial starts (response mapping cues)
            st.draw_mask(win, basic_stim)
            st.fixation(win, basic_stim)
            if correct == 1:
                 basic_stim['fixation_point'].color = [-1, 1,-1] # colors are expressed as deviation from grey red -> [1,-1,-1]
            else:
                 basic_stim['fixation_point'].color = [1, -1,-1]
                
            basic_stim['fixation_point'].draw()
            st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)  
            st.draw_contour(win,basic_stim)
            if i_si == 0:     
                if sst: win.callOnFlip(p_port.write, tg_resp.encode()) # send trigger
            if i_si == 1:     
                if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
            t = win.flip()
            if (i_si == 0): trial_times = np.append(trial_times, t)
            
        
        trial_times = np.array(trial_times) 
        trialClocktimes = np.vstack([trialClocktimes, trial_times]) if trialClocktimes.size else trial_times  
        
        # if there is an outlier trial, stop experiment
        this_times = np.diff(trial_times)
        mean_stim_time = np.mean(this_times[4:14]) * 1000
        if mean_stim_time > stim['stim_time'] + 75 or mean_stim_time < stim['stim_time'] - 75:
            print('Bad timing!! Cerrando programa')
            win.close()
            core.quit
            
        if thisTrial['Pred'] == 'NP': # only append correct values if the trial is Neutral
            correct_seq = np.append(correct_seq, correct)
        # Storing data in variables
        
        thr_trials_var.append([subj_id, thisBlock ,trials.thisN,thisTrial['Pred'], cond, thisTrial['prob'], x, thisResp[1], correct, rt_deci])# saving conditions here
        thr_trials_ori.append(t_orient.tolist())
        #trialClocktimes.append(trial_times)
    # Get datafiles in pandas format and attack to main Exp.variable
    headers =  thr_trials_var.pop(0)
    block['data'] = DataFrame(thr_trials_var, columns=headers)
    headers =  thr_trials_ori.pop(0)
    block['trial_orientations'] = DataFrame(thr_trials_ori , columns=headers)
    block['block_duration'] =  Clock.getTime() -  BlockClockStart
    block['BlockClockStart'] = BlockClockStart 
    
    main_frame_log = {}
    main_frame_log['droppedframes'] = win.nDroppedFrames
    main_frame_log['timings'] = np.diff(trialClocktimes,axis = 1)
    block['main_frame_log'] = main_frame_log 
    
    main_exp['Exp_blocks'][thisBlock] =  block  # saving block data
    
    
    cr_lot = instr.lotery(win, block, ifi)
    corr_lotery.append(cr_lot) # lotery win or lost
    toFile(resultspath, expInfo) #saving file to disk
    
main_exp['monitor_features'] = monitor_features
main_exp['stim'] = stim

expInfo['main_exp'] =  main_exp



print('Overall, %i frames were dropped.' % win.nDroppedFrames)

expInfo['lotery'] = corr_lotery

toFile(resultspath, expInfo) #saving file to disk

instr.end_experiment(win)

# closing everything
win.close()
if sst: p_port.close()

print('El participante ha ganado ' + str(sum(corr_lotery)) + ' puntos de ' + str(main_exp['nblocks']))

core.quit





    