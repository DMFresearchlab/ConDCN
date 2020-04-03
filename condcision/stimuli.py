from psychopy import visual
import numpy as np
#import orient_decision.win
 
def stim_config(ifi):
    # Experiment timings design
    stim = {}
    stim['nreps'] = 3 # number of gratings per sequence
    stim['nstim'] = 6 # number of gratings per sequence
    stim['ntrials_per_cond'] = 1 # number of repetitions per trial type
    stim['stim_time'] = 250 # in miliseconds
    stim['stim_frames'] = round(stim['stim_time']/ifi)
    # stim['ISI1_frames'] IS determined randomly in each trial
    stim['ISI2_frames'] = round(750/ifi) # fixation  + circle + resp_options
    stim['wait_feedback_frames'] = round(50/ifi) # 
    stim['feedback_frames'] = round(750/ifi) # 
    stim['feedback_quest_frames'] = round(200/ifi) 
       
        # Create all the stimuli
    # Stimuli parameters
    stim['size_stim'] = 14 # degs # este es el lado del cuadrado del grating
    stim['grating_contrast'] = 0.5
    stim['SF'] = 2

    return stim

def draw_basic(win, stim): # creates a dictionary with basic stimulus features
    
    # Fication point
    basic_stim = {}
    basic_stim['fixation_point'] =  visual.PatchStim(win, color= [1, 1, 1], tex=None,mask='circle', size=0.2, units = "deg")
    basic_stim['fixation_point_c'] = visual.PatchStim(win, color= [-1, -1, -1], tex=None,mask='circle', size=0.3, units = "deg")
    basic_stim['fixation_circle'] =  visual.PatchStim(win, color= [0, 0, 0], tex=None, mask='circle', size=0.5, units = "deg")
      
    # Grating stim
    basic_stim['grating'] = visual.GratingStim(win=win, mask="raisedCos" , size=stim['size_stim'], pos=[0,0], sf=stim['SF'],
                                 units = "deg", contrast = stim['grating_contrast'], maskParams = {'sd': 3} ) # , color = [1,0,1]
    
    # Contour boundary
    basic_stim['stim_contour'] = visual.Circle(win=win,lineWidth = 8, units="deg", radius=stim['size_stim']/2, lineColor=[-1, -1, -1],edges=128)
    basic_stim['stim_contour_in'] = visual.Circle(win=win,lineWidth = 2, units="deg", radius=stim['size_stim']/2, lineColor=[1, 1, 1],edges=128)

    # Mask components
    basic_stim['g1'] = visual.GratingStim(win=win, mask="raisedCos" , opacity=0.5, size=stim['size_stim'], ori = 135, pos=[0,0], sf=stim['SF'], units = "deg", contrast = stim['grating_contrast'], maskParams = {'sd': 3} ) # , color = [1,0,1]
    basic_stim['g2'] = visual.GratingStim(win=win, mask="raisedCos" , opacity=0.5, size=stim['size_stim'], ori = 45, pos=[0,0], sf=stim['SF'], units = "deg", contrast = stim['grating_contrast'], maskParams = {'sd': 3} ) # , color = [1,0,1]
    basic_stim['g3'] = visual.GratingStim(win=win, mask="raisedCos" , opacity=0.5, size=stim['size_stim'], ori = 90, pos=[0,0], sf=stim['SF'], units = "deg", contrast = stim['grating_contrast'], maskParams = {'sd': 3} ) # , color = [1,0,1]
    basic_stim['g4'] = visual.GratingStim(win=win, mask="raisedCos" , opacity=0.5, size=stim['size_stim'], ori = 135, pos=[0,0], sf=stim['SF'], units = "deg", contrast = stim['grating_contrast'], maskParams = {'sd': 3} ) # , color = [1,0,1]
       
    # Circ. Mask
    basic_stim['circ_grating_mask'] = visual.RadialStim(win=win, mask="raisedCos" , opacity=0.5, size=stim['size_stim'], ori = 90, pos=[0,0], radialCycles=8, angularCycles=0, units = "deg", contrast = 0.8, maskParams = {'sd': 3} ) # , color = [1,0,1]
    
    # Response options Diagonal or Cardinal
    basic_stim['circle'] = visual.Circle(win=win,lineWidth = 10, units="deg", radius=2, lineColor=[1, 1, 1],edges=64)
    basic_stim['linev'] = visual.Line(win=win,lineWidth = 10,units="deg",lineColor=[1, 1, 1],start = [0, -2],end = [0, 2])
    basic_stim['lineh'] = visual.Line(win=win,lineWidth = 10,units="deg",lineColor=[1, 1, 1],start = [-2, 0],end = [2, 0])
    basic_stim['feedback1'] = visual.Circle(win=win,units="deg", radius=1, edges=64, color = [0.5, 0.5, 0.5], pos = [-4, -10])
    basic_stim['feedback2'] = visual.Circle(win=win,units="deg", radius=1, edges=64, color = [0.5, 0.5, 0.5], pos = [0, -10])
    basic_stim['feedback3'] = visual.Circle(win=win,units="deg", radius=1, edges=64, color = [0.5, 0.5, 0.5], pos = [4, -10])
    return basic_stim


def draw_mask(win, basic_stim):
    basic_stim['g1'].draw() 
    basic_stim['g2'].draw() 
    basic_stim['g3'].draw() 
    basic_stim['g4'].draw()
    return

def fixation(win, basic_stim):
    basic_stim['fixation_circle'].draw()
    basic_stim['fixation_point_c'].draw()
    basic_stim['fixation_point'].draw()
    return

def draw_contour(win, basic_stim):
    basic_stim['stim_contour'].draw()
    basic_stim['stim_contour_in'].draw()
    return

# Response options (cardinal/diagonal)
def resp_option(win, basic_stim, ori, rgb, pos): # response option object
    basic_stim['lineh'].ori = ori
    basic_stim['lineh'].pos = pos
    basic_stim['linev'].ori = ori
    basic_stim['linev'].pos = pos
    basic_stim['lineh'].lineColor = rgb
    basic_stim['linev'].lineColor = rgb
    basic_stim['circle'].lineColor = rgb
    basic_stim['circle'].pos = pos
    basic_stim['circle'].draw()
    basic_stim['linev'].draw()
    basic_stim['lineh'].draw()
    return

def resp_mapping(win, resp_option_f, basic_stim,  resp_maps, rgbs): # pass resp_option drawinf functions maps the response options spatially 
    resp_option_f(win, basic_stim, resp_maps[0], rgbs[0], np.array([-10,0])) 
    resp_option_f(win, basic_stim, resp_maps[1], rgbs[1], np.array([10,0]))
    return


def monitor_def():
    # lets define the monitor if necessary in a dictionary (you can define other monitors def2, def3)
    monitor_def1 = {}
    monitor_def1['monitor_name'] = 'mundet_screen' # monitor to use (make sure to define monitors in exp_monitors center.
    monitor_def1['monitor_pixels']  = (1920, 1080)
    monitor_def1['monitor_width'] = 54
    monitor_def1['distance2monitor'] = 50
    
    monitor_def2 = {}
    monitor_def2['monitor_name'] = 'multiple2' # monitor to use (make sure to define monitors in exp_monitors center.
    monitor_def2['monitor_pixels']  = (1024, 768)
    monitor_def2['monitor_width'] = 30
    monitor_def2['distance2monitor'] = 40
    
    monitor_def3 = {}
    monitor_def3['monitor_name'] = 'asus_home' # monitor to use (make sure to define monitors in exp_monitors center.
    monitor_def3['monitor_pixels']  = (1920,1080)
    monitor_def3['monitor_width'] = 60
    monitor_def3['distance2monitor'] = 40
    
    monitor_def4 = {}
    monitor_def4['monitor_name'] = 'dell_xps13' # monitor to use (make sure to define monitors in exp_monitors center.
    monitor_def4['monitor_pixels']  = (3200,1800)
    monitor_def4['monitor_width'] = 30
    monitor_def4['distance2monitor'] = 60
    
    monitors = [monitor_def1, monitor_def2, monitor_def3, monitor_def4]
    return monitors
        
        