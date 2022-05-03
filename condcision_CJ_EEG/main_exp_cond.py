#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Alexis Perez-Bellido
# Code for the confirmation bias experiment (using predcision code)


version = 1.5

#cd('/home/node2/Experiments/PreDCN-master/predcision')
import numpy as np

import matplotlib.pyplot as plt
import os, sys
import exp_func as exp
import stimuli as st
import instructions as instr
import serial
import pyxid2 #cedrux lib

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


sst = False # Using serial port with TTL to send triggers
cedrux = True # Using serial port with cedrux to send triggers
eyeT = True


if sst: 
    #p_port = serial.Serial('/dev/tty.usbserial-BBTKUSBTTL', 115200, timeout = 0) # 'COM3', 115200, timeout = 0 this is for windows
    p_port = serial.Serial('COM3', 115200, timeout = 0)
    p_port.write(b'00')
    core.wait(0.2)
    p_port.write(b'RR')
    
if cedrux:
    # get a list of all attached XID devices
    devices = pyxid2.get_xid_devices()
    dev = devices[0] # get the first device to use
    dev.reset_base_timer()
    dev.set_pulse_duration(10) # duration of trigger
        
# Trigger keys

tg_mblock =  '07'
tg_mtrial = '27'
tg_newtrial = '31'
tg_mask = '41'
tg_stim = '13'
tg_resp_start = '21'
tg_resp = '09'
tg_confi = '05'

tg_zero = '00'
#prefs.hardware['audioLib']=['pyo'] # use Pyo audiolib for good temporal resolution
#from psychopy.sound import Sound # This should be placed after changing the audio library
# monitors.getAllMonitors()
#from general_settings import * # reads the variables in one script

# Collect subject data
expInfo, resultspath = exp.mainexp_subject_info(version)
subj_id = expInfo['subjInfo']['observer']



if eyeT: # initializing eyetracker
    import pupilabs_func as pup
    import zmq #socket pupil labs
    import time
    # 1. Setup network connection
    ip_address = "127.0.0.1"
    port= 50020
    
    pup.check_capture_exists(ip_address, port)
    requester, pub_socket = pup.setup_pupil_remote_connection(ip_address, port)
         
    #ctx = zmq.Context()
    #requester = ctx.socket(zmq.REQ)
    #requester.connect('tcp://127.0.0.1:50020')
    # Request 'SUB_PORT' for reading data
    requester.send_string('SUB_PORT') 
    sub_port = requester.recv_string()
    # 2. Setup local clock function
    local_clock = time.perf_counter

    # 3. Measure clock offset accounting for network latency
    stable_offset_mean = pup.measure_clock_offset_stable(
    requester, clock_function=local_clock, n_samples=10 )
    
    pupil_time_actual = pup.request_pupil_time(requester)
    local_time_actual = local_clock()
    pupil_time_calculated_locally = local_time_actual + stable_offset_mean
    print(f"Pupil time actual: {pupil_time_actual}")
    print(f"Local time actual: {local_time_actual}")
    print(f"Stable offset: {stable_offset_mean}")
    print(f"Pupil time (calculated locally): {pupil_time_calculated_locally}")
    # Start the annotations plugin
    pup.notify(
        requester,
        {"subject": "start_plugin", "name": "Annotation_Capture", "args": {}},
    )
    
    # starting to record eye tracker data for this sesssion
    requester.send_string('R '+  subj_id )   
    requester.recv_string()
    # Calibrating eyetracker Check
    myDlg = gui.Dlg(title="Calibrate the eyetracker")
    myDlg.addText('Have you already calibrated the eyetracker?')
    ok_data = myDlg.show()   
    



# Loading monitor definitions
monitores = st.monitor_def() 

#mon, expInfo['monitor'] = exp.define_monitor(monitores[1]) # select the correct monitor
mon, expInfo['monitor'] = exp.define_monitor(monitores[1]) # select the correct monitor

# Creating a new experimental window
monitor_features = {}
monitor_features['monitor'] = mon
monitor_features['units'] = 'deg' # units to define your stimuli
monitor_features['screen_id'] = 0 # when using a extended display 
monitor_features['full']  = True
monitor_features['Hz'] =  60 #'auto' #144 #60 this can be set to "auto" to estimate the refreshing rate of the monitor, although it can fail often

## opening window   
win, monitor_features = exp.create_window(monitor_features)
ifi = monitor_features['ifi']

# Experiment timings and characteristics 
stim = st.stim_config(ifi) #Loading stim characteristics      
basic_stim = st.draw_basic(win,stim)


# Experimental condition preparation
# response mapping (shuffled in each experiment )
expInfo['resp_maps']  = np.array([0, 45]) # 0 -> cardinal; 45 -> diagonal
expInfo['oddblock_repeat']  = np.random.randint(2) # determine whether odd N blocks are the ones with more repeated trials (1 repeat in odd, 0 repear in even blocks)


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

# Experiment design
stimList = []
for cond in ['repeat', 'nonrepeat']: # 75% vs 25% larger values than 1 in the CP & DP conditions correspond to one orientation
            stimList.append({'Conds':cond}) #, 'n_reps': n_reps
            
            
            
main_exp  = {}
main_exp['nblocks']     = 4 # 4 # totaltime = 90 * 6 * 5
main_exp['Exp_blocks']  = [None] * main_exp['nblocks'] # assigning memory for storing block data
main_exp['trial_reps']  = 3 #12 is equal to 1 minute (2 trials x 3 repetitions)

# instr.block_ID(win, block_type)

corr_lotery = []
mouse = event.Mouse(visible = True)
mouse.setPos([0,0])

win.mouseVisible = False 

for thisBlock in range(main_exp['nblocks']): # iterate over blocks
      
    instr.block_start(win)
    # lets trigger the beggining of the experiment
    
    #instr.block_ID(win, block_type, True) # explicit information about the block type
    inst = visual.TextStim(win, pos = [0,8])
    inst.wrapWidth = 20
    inst.height = 1
    
    

    if sst: win.callOnFlip(p_port.write, tg_mblock.encode())
    if cedrux: dev.activate_line(bitmask=7)
    win.flip()
    if sst: win.callOnFlip(p_port.write, tg_zero.encode())
    win.flip()
    
    BlockClockStart = Clock.getTime() # block experiment time
    block = {} # dummy variable to save block data
    correct_seq = np.array([]) # saving seq. of correct responses per trial sequence
    trial_rep = 0 # update staricase after each iteraction
    thr_trials_var = [['subj','nblock', 'ntrial', 'nrep',  'trial_type', 'cond', 'DV', 'resp', 'r_map', 'correct', 'confi', 'RT']] # saving conditions here
    thr_trials_ori = [['o1','o2','o3','o4','o5','o6']]
    
    TS_eye = [['fp', 'onset', 'resp', 'confi', 'end']] # saving conditions here
    
    trialClocktimes = np.array([]) # saving whole times here
    correct_seq = np.array([]) # saving seq. of correct responses per trial sequence
    
    trials = data.TrialHandler(stimList, main_exp['trial_reps'], method='random') 
    
    for thisTrial in trials:  # will continue the staircase until it terminates!        
         
        thisIncrement = np.around(next(staircase),decimals = 5)
        cond = np.random.choice([-1, 1])  # will be cardinal or diagonal
        ExpClockTrial = Clock.getTime()
        
        # decision variable in this trial
        x =  cond*(thisIncrement)
        print(x)       
        fixation_color = [1, 1, 1]
        trial_perform = np.array([])
        instr.new_trial(win)
        
        for i_rep in range(stim['nreps']): 
            
            repeat = thisTrial['Conds'] 
             
            if  i_rep == 0:
                sel_trials = np.where((decision_var_T > x-0.025) & (decision_var_T < x+0.025))
                trial_sel_idx = np.random.choice(sel_trials[0]) # selecting orientation vector for this trial.
                t_orient = orientations[trial_sel_idx]
                decision_avg = decision_var_T[trial_sel_idx] 
               
            if  i_rep == 2:
                if repeat == 'nonrepeat': # whether in the third repetition we show the same or different information
                      t_orient = t_orient + pi/2 # getting orthogonal orientation
                      t_orient[t_orient>=pi] =  t_orient[t_orient>=pi] - pi # reorienting orientations in the range 0 to 360
 

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
                
                if i_si == 0:     
                    if sst: win.callOnFlip(p_port.write, tg_newtrial.encode()) # send trigger
                    if cedrux: dev.activate_line(bitmask=31)
                if i_si == 1:     
                    if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins   
                
                if i_rep == 1:
                     basic_stim['feedback2'].color = [0.5, 0.5, 0.5]
                     basic_stim['feedback2'].autoDraw = True
                     basic_stim['feedback2'].draw() 
                     
                if i_rep == 2:
                     basic_stim['feedback3'].color = [0.5, 0.5, 0.5]
                     basic_stim['feedback3'].autoDraw = True
                     basic_stim['feedback3'].draw()
                                    
                t = win.flip()
                
                if (i_si == 0): trial_times = np.append(trial_times, t)
        
            # medf_s.play() # reproduce 1st auditory cue
        
            for i_si in range(stim['ISI2_frames']): # second period before the second beep
                st.fixation(win, basic_stim)
               # st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                st.draw_contour(win,basic_stim)
                if i_si == 0:     
                    if sst: win.callOnFlip(p_port.write, tg_mtrial.encode()) # send trigger
                    if cedrux: dev.activate_line(bitmask=27)
                    if eyeT:
                       # time saved in results
                       requester.send_string('t')
                       fp_TS=requester.recv_string() 
                       # trigger saved pupil labs
                       local_time = local_clock()  
                       if i_rep == 0: label = "onset1"
                       if i_rep == 1: label = "onset2"
                       if i_rep == 2: label = "onset3"
                       duration = 0.0
                       eye_tgr = pup.new_trigger(label, duration, local_time + stable_offset_mean)
                       pup.send_trigger(pub_socket, eye_tgr)
                       
                if i_si == 1 and sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  

                t = win.flip()
                if (i_si ==0): trial_times = np.append(trial_times, t)
 
            # Draw the stimuli
            for istim in range(stim['nstim']+2): # 2 stim for the masks sandwiching
                event.clearEvents()
                if (istim == 0) | (istim == stim['nstim']+1): # if first or last mask
                   for frame in range(stim['stim_frames']):  # drawing stim frames
                        if frame == 0:     
                           if sst: win.callOnFlip(p_port.write, tg_mask.encode()) # send trigger
                           if cedrux: dev.activate_line(bitmask=41)
                           if eyeT:
                                requester.send_string('t')
                                onset_TS=requester.recv_string()
                                                                
                        if frame == 1:
                            if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
                        if (frame == stim['stim_frames']) and (istim == 0):  # the last frame should be emptyin the first mask
                           # -1 and istim == 0
                            st.fixation(win, basic_stim)
                            st.draw_contour(win,basic_stim)
                            t = win.flip()                    
                        else: # flip empty frame
                            st.draw_mask(win, basic_stim)
                            st.fixation(win, basic_stim)
                            st.draw_contour(win,basic_stim)
                            t = win.flip()
                        if (frame ==0): trial_times = np.append(trial_times, t)
                else:
                    basic_stim['grating'].ori =  np.rad2deg(t_orient[istim-1]) # change orientation for each stim
                    basic_stim['grating'].phase = np.random.rand()
        
                    for frame in range(stim['stim_frames']): # drawing stim frames
                        if frame == 0:
                            if sst: win.callOnFlip(p_port.write, tg_stim.encode()) # send trigger
                            if cedrux: dev.activate_line(bitmask=13)
                        if frame == 1: 
                            if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
                        if frame < stim['stim_frames'] - 1:  # the last frame should be empty
                            basic_stim['grating'].draw()
                            st.fixation(win, basic_stim)
                            st.draw_contour(win,basic_stim)
                            t = win.flip()
                            if (frame == 0): trial_times = np.append(trial_times, t)
                        else: # flip empty frame
                            st.fixation(win, basic_stim)
                            st.draw_contour(win,basic_stim)
                            t = win.flip()               
            
 
            
            for i_si in range(stim['ISI3_frames']): # second period before the second beep
                st.draw_mask(win, basic_stim)
                st.fixation(win, basic_stim)
               # st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                st.draw_contour(win,basic_stim)
                t = win.flip()
                trial_times = np.append(trial_times, t)
                
             
            if expInfo['resp_maps'][0] == 0:
                resp_pos = 1
            else:
                resp_pos = -1
                
                
            ticks = [-1,1] #['1','2','3','4','5']
            labels = ['Cardinal','','Diagonal'] # labels = labels,
            deci_slider = visual.Slider(win,ticks= ticks,  pos = [0,0],  
                                      size = (6,1), style = 'radio', granularity = 1)
            
            i_si = 0 # initialize counter
            mouse.setPos([0,0])
            while not deci_slider.rating:
                if i_si == 0:
                    if sst: win.callOnFlip(p_port.write, tg_resp_start.encode()) # send trigger
                    if cedrux: dev.activate_line(bitmask=21)
                    if eyeT:
                        requester.send_string('t')
                        deci_TS=requester.recv_string()
                        # trigger saved pupil labs
                        local_time = local_clock()  
                        if i_rep == 0: label = "deci1"
                        if i_rep == 1: label = "deci2"
                        if i_rep == 2: label = "deci3"
                        duration = 0.0
                        eye_tgr = pup.new_trigger(label, duration, local_time + stable_offset_mean)
                        pup.send_trigger(pub_socket, eye_tgr)
                if i_si == 1:
                    if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins
                st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps* -1)     
                deci_slider.draw() 
                st.fixation(win, basic_stim)
                #visual.Line(win=win,lineWidth = 20,units="deg",lineColor=[1, 1, 1],start = [0, -1],end = [0, 1]).draw()
                i_si += 1
                t = win.flip()
              
            resp_ang = deci_slider.rating * resp_pos
            rt_deci = deci_slider.getRT()
     
            # Confidence rating    
            ticks = [-1,1] #[-1,-0.5,0,0.5,1]
            labels = ['Dudosa','Segura'] # labels = labels,['Dudosa','','','','','Segura']
            confi_slider = visual.Slider(win,ticks= ticks, labels= labels, pos = [0,0],  
                                      size = (8,1), style = 'slider', granularity = 0.05, styleTweaks = ['triangleMarker'])
            
            i_si = 0 # initialize counter      
            while not confi_slider.rating:
               if i_si == 0:
                    if sst: win.callOnFlip(p_port.write, tg_resp.encode()) # send trigger
                    if cedrux: dev.activate_line(bitmask=9)
                    if eyeT:
                        requester.send_string('t')
                        confi_TS=requester.recv_string()
                        # trigger saved pupil labs
                        local_time = local_clock()  
                        if i_rep == 0: label = "confi1"
                        if i_rep == 1: label = "confi2"
                        if i_rep == 2: label = "confi3"
                        duration = 0.0
                        eye_tgr = pup.new_trigger(label, duration, local_time + stable_offset_mean)
                        pup.send_trigger(pub_socket, eye_tgr)
               if i_si == 1:
                    if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins
               i_si += 1
               confi_slider.draw() 
               t = win.flip()               
            
            confi = confi_slider.rating
            win.mouseVisible = False # hide mouse
            mouse.setPos([0,0])
            
            if (x > 0 and resp_ang > 0) or (x < 0 and  resp_ang < 0):
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
            correct_seq = np.append(correct_seq, correct)

                
            for i_si in range(stim['wait_feedback_frames']): # wait time for feedback and send a couple of triggers
                if i_si == 0:
                     if sst:  win.callOnFlip(p_port.write, tg_confi.encode()) # send trigger
                     if cedrux: dev.activate_line(bitmask=5)
                     if eyeT:
                        requester.send_string('t')
                        end_TS=requester.recv_string()
                        # trigger saved pupil labs
                        local_time = local_clock()  
                        if i_rep == 0: label = "end1"
                        if i_rep == 1: label = "end2"
                        if i_rep == 2: label = "end3"
                        duration = 0.0
                        eye_tgr = pup.new_trigger(label, duration, local_time + stable_offset_mean)
                        pup.send_trigger(pub_socket, eye_tgr)
                if i_si == 1 and sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins           
                basic_stim['fixation_point'].draw()
                t = win.flip()
                if (i_si ==0): trial_times = np.append(trial_times, t)
                
            
            basic_stim['fixation_point'].draw()
            win.flip()   
            #print(correct)
            print(trial_perform)
           # if i_rep == 0: staircase.addResponse(correct) # adding information to staircase
           # st.resp_option(win, basic_stim, expInfo['resp_maps'][0], black_resps[1], np.array([-2,-7]))
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
                    
 #st.resp_mapping(win, st.resp_option, basic_stim,  expInfo['resp_maps'], black_resps)
                    
            
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
            thr_trials_var.append([subj_id, thisBlock , trials.thisN, i_rep,  repeat, cond, x, resp_ang, expInfo['resp_maps'][0], correct, confi, rt_deci])# saving conditions here
            thr_trials_ori.append(t_orient.tolist())
            
            if eyeT: # appending pupil labs timestamps
                TS_eye.append( np.array( [float(fp_TS),  float(onset_TS),  float(deci_TS),  float(confi_TS),  float(end_TS)]))
    
    inst.autoDraw = False
    basic_stim['fixation_point'].draw()
    #st.draw_contour(win,basic_stim)
    t = win.flip()
    # Get datafiles in pandas format and attack to main Exp.variable   
    headers =  thr_trials_var.pop(0)
    block['data'] = DataFrame(thr_trials_var, columns=headers)
    headers =  thr_trials_ori.pop(0)
    block['trial_orientations'] = DataFrame(thr_trials_ori , columns=headers)
    block['block_duration'] =  Clock.getTime() -  BlockClockStart
    
    if eyeT:
        headers =  TS_eye.pop(0)            
        block['eye_TS'] = DataFrame(TS_eye, columns=headers) 
          
    #main_exp['Exp_blocks'][1]['trial_orientations']
    main_frame_log = {}
    main_frame_log['droppedframes'] = win.nDroppedFrames
    main_frame_log['timings'] = trialClocktimes
    block['main_frame_log'] = main_frame_log 
    main_exp['Exp_blocks'][thisBlock] =  block  # saving block data
    
    
  
    cr_lot = instr.lotery(win, block, ifi)
    corr_lotery.append(cr_lot) # lotery win or lost    
    
    toFile(resultspath, expInfo) #saving file to disk
            #trialClocktimes.append(trial_times)
        # Get datafiles in pandas format and attack to main Exp.variable
    if eyeT and (thisBlock < main_exp['nblocks']-1): # closing window to recalibrate eyetracker (not in the last block)
        win.close()
        myDlg = gui.Dlg(title="Calibrate the eyetracker")
        myDlg.addText('Have you already calibrated the eyetracker?')
        ok_data = myDlg.show()
        ## RE-opening window   
        win, monitor_features = exp.create_window(monitor_features)
        ifi = monitor_features['ifi']    
        # Experiment timings and characteristics 
        stim = st.stim_config(ifi) #Loading stim characteristics      
        basic_stim = st.draw_basic(win,stim)
        mouse = event.Mouse(visible = False)
        mouse.setPos([0,0])

approxThreshold = np.average(staircase.reversalIntensities[-5:])

main_exp['monitor_features'] = monitor_features
main_exp['stim'] = stim
main_exp['lotery'] = corr_lotery

expInfo['subjInfo']['guess'] = approxThreshold  # save threshold for next experiment

headers =  thr_trials_ori.pop(0)

expInfo['main_exp'] =  main_exp

print('Overall, %i frames were dropped.' % win.nDroppedFrames)


toFile(resultspath, expInfo) #saving file to disk
nblocks = main_exp['nblocks']

corr_lotery = np.array(corr_lotery)
corr_lotery[corr_lotery == - 1] = 0;



instr.end_experiment_lot(win, corr_lotery, nblocks)

# closing everything
win.close()

print('El participante ha ganado ' + str(sum(corr_lotery)) + ' puntos de ' + str(main_exp['nblocks']))

if eyeT:   # stopping the eyetracker recording    
    requester.send_string('r '+ subj_id )   
    requester.recv_string()

core.quit

