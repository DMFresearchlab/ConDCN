from psychopy import visual, event, core
import random
import numpy as np
# Experiment instructions



def init_instructions(win, basic_stim, ifi):
       
    gratexample = visual.GratingStim(win=win, mask="raisedCos" , size=5, pos=[0,5], sf=3,
                                 units = "deg", contrast =0.5, maskParams = {'sd': 3} ) # , color = [1,0,1]
      
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.text = "En cada trial te presentaré una sequencia rápida de 6 estimulos con diferentes orientaciones... "
    inst.height = 0.7
    nextt = visual.TextStim(win, pos = [0,-6])
    nextt.wrapWidth = 20
    nextt.height = 0.7
    nextt.color = 'black'
    nextt.text = "Pulsa spacio para continuar"
    
    time_grating = round(250/ifi) # time in frames that each grating is displayed on the screen
    
    allKeys = []
    while allKeys != ['space']:
        #angle = np.random.uniform(0,1,1)
        gratexample.ori =  random.randint(0,359)
      #  gratexample.pos =  [0,5]
        
        for iframe in  range(time_grating):
            inst.draw()
            nextt.draw()
            gratexample.draw()
            win.flip()
        
        allKeys = event.getKeys()
        #print(allKeys)
    
    # Centering example wheel
    basic_stim['lineh'].pos = np.array([6, 5])
    basic_stim['linev'].pos = np.array([6, 5])
    basic_stim['circle'].pos = np.array([6, 5])
    #basic_stim['linev'].lineColor = np.array([-1,-1,-1])
    
    inst.text = "Tu tarea consiste en estimar si la orientación media de los estimulos presentados está más cerca de los ejes cardinales (vertical u horizonal)\
    pulsando la tecla ESPACIALMENTE congruente con el símbolo cardinal... ('z' izquierda, 'm' derecha)"
    basic_stim['lineh'].ori = 0
    basic_stim['linev'].ori = 0
    

    allKeys = []
    while allKeys != ['space']:
        #angle = np.random.uniform(0,1,1)
        gratexample.ori =  random.randint(0,359)
       # gratexample.pos =  [0,5]
        inst.draw()
        nextt.draw()
        
        basic_stim['linev'].draw()
        basic_stim['lineh'].draw()
        basic_stim['circle'].draw()
        gratexample.draw()
        win.flip()
        core.wait(0.2)
        allKeys = event.getKeys()
        #print(allKeys)       
    
    inst.text = "... o en los ejes diagonales pulsando la tecla ESPACIALMENTE congruente con el símbolo diagonal"
    nextt.text = "Pulsa spacio para continuar"
    basic_stim['lineh'].ori = 45
    basic_stim['linev'].ori = 45
    basic_stim['lineh'].pos = np.array([-6, 5])
    basic_stim['linev'].pos = np.array([-6, 5])
    basic_stim['circle'].pos = np.array([-6, 5])
    
    
    allKeys = []
    while allKeys != ['space']:
        #angle = np.random.uniform(0,1,1)
        gratexample.ori =  random.randint(0,359)
      #  gratexample.pos =  [0,5]
        inst.draw()
        nextt.draw()
        basic_stim['linev'].draw()
        basic_stim['lineh'].draw()
        basic_stim['circle'].draw()
        gratexample.draw()
        win.flip()
        core.wait(0.2)
        allKeys = event.getKeys()
        #print(allKeys)
           
    return



def localizer_instructions(win, basic_stim, ifi):
       
    gratexample = visual.GratingStim(win=win, mask="raisedCos" , size=5, pos=[0,5], sf=3,
                                 units = "deg", contrast =0.5, maskParams = {'sd': 3} ) # , color = [1,0,1]      
    #gratexample.draw()
    #win.flip()
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.text = "En esta tarea, en cada trial te presentaré una sequencia rápida de 8 estimulos con diferentes orientaciones... "
    inst.height = 0.7
    nextt = visual.TextStim(win, pos = [0,-6])
    nextt.wrapWidth = 20
    nextt.height = 0.7
    nextt.color = 'black'
    nextt.text = "Pulsa spacio para continuar"
    
    time_grating = round(250/ifi) # time in frames that each grating is displayed on the screen
    
    allKeys = []
    while allKeys != ['space']:
        #angle = np.random.uniform(0,1,1)
        gratexample.ori =  random.randint(0,359)
       #gratexample.pos =  [0,5]
        if random.randint(0, 14)>1:
            gratexample.sf =  3
        else:
            gratexample.sf = 2
        
        for iframe in  range(time_grating):
            inst.draw()
            nextt.draw()
            gratexample.draw()
            win.flip()
        
        allKeys = event.getKeys()
        #print(allKeys)
       
    inst.text = "Tu tarea consiste en pulsar la barra espaciadora al final de\
     cada trial si detectas que la anchura de las barras de alguno de los estimulos es diferente al resto"


    allKeys = []
    while allKeys != ['space']:
        #angle = np.random.uniform(0,1,1)
        gratexample.ori =  random.randint(0,359)
       # gratexample.pos =  [0,5]
        if random.randint(0, 14)>1:
            gratexample.sf =  3
        else:
            gratexample.sf = 2
        inst.draw()
        nextt.draw()
        gratexample.draw()
        win.flip()
        core.wait(0.2)
        allKeys = event.getKeys()
        #print(allKeys)       
    


def main_instructions(win):
       
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.text = "Empezamos el experimento. Seguimos con la misma tarea que antes\
    con la diferencia de que ahora te repetiré tres veces cada secuencia de estímulos y tendrás tres oportunidades de hacerlo bien"
    inst.height = 0.7
    nextt = visual.TextStim(win, pos = [0,-6])
    nextt.wrapWidth = 20
    nextt.height = 0.7
    nextt.color = 'black'
    nextt.text = "Pulsa spacio para continuar "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    
    inst.text = "Cuantos mas puntos verdes consigas durante el experimento,\
    MEJOR!"
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    
    return

def block_start(win):
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.height = 0.7
    nextt = visual.TextStim(win, pos = [0,-6])
    nextt.height = 0.7
    nextt.color = "black"
    inst.text = "Comenzando bloque! Coloca los dedos sobre botones de respuesta"
    nextt.text = "Cuando estes listo/a, pulsa espacio para comenzar "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys( keyList=['space'])
    core.wait(1)


def end_experiment(win):
    inst = visual.TextStim(win, pos=[0,0], height = 1.2)
    inst.text = "Final de esta parte del experimento!! Avisa al investigador y pulsa -espacio- para continuar."     
    inst.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    return
    
    

def lotery(win, block, ifi):
    inst = visual.TextStim(win, pos=[0,0], height = 1.2)
    inst.text = "Loteria de trials. Pulsa -espacio- para comenzar el sorteo."     
    inst.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    allKeys = []
    number_text = visual.TextStim(win, pos=[0,4], height = 2)
    inst.text = "Loteria de trials. \
    Pulsa -espacio- para parar el sorteo."     
    n_trials = len(block['data'])
    
    time_number = round(100/ifi) # time in frames that each lotery number is displayed on the screen
    
    
    allKeys = []
    while allKeys != ['space']:
        #angle = np.random.uniform(0,1,1)
        trial_ix =  random.randint(0,n_trials-1)
        corr = block['data'].loc[trial_ix,'correct']
        col = 'red' if corr == 0 else 'green'
        number_text.text = trial_ix
        number_text.color = col
        
        for iframe in  range(time_number):
             inst.draw()
             number_text.draw()
             win.flip()
        
        allKeys = event.getKeys()
        
    allKeys = []
    if corr:
        inst.text = "Genial, has ganado! Pulsa -espacio- para continuar"
    else:
        inst.text = "Oh! Qué mala suerte. Pulsa -espacio- para continuar"

    inst.draw()
    number_text.draw()
    win.flip()
        
    event.waitKeys(keyList = ["space"]) 
    return corr 
        
        
        