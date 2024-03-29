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
    
    inst.text = "Tu tarea consiste en estimar si la orientación media de los estimulos presentados está más cerca de los ejes cardinales (vertical u horizonal)"
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
    
    inst.text = "... o en los ejes diagonales, seleccionando con el rato el simbolo correcto"
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
       
    inst.text = "Tu tarea consiste en detectar si la anchura de las barras de alguno de los estimulos es diferente al resto. Justo al terminar cada secuencia tendrás unos pocos segundos para pulsar la barra espaciadora SOLO cuando hayas detectado un estímulo con barras más anchas"


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
    inst.text = "Empezamos el experimento de integración de orientaciones. \
    En esta parte te presentaré TRES secuencias de estímulos repetidas. Tendrás TRES oportunidades de hacerlo bien.\
    Además tras cada secuencia tendrás que indicar en una escala continua como de segur@ estás de tu decision (desde dudoso a seguro)"
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
    
    inst.text = "Tras terminar cada trial de dos repeticiones verás unos circulos verdes o rojos indicando si has acertado o fallado. \
    Cuantas más circulos verdes completes MEJOR!"
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    
    return

def block_start(win):
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.height = 0.7
    inst.color = 'white'
    nextt = visual.TextStim(win, pos = [0,-6])
    nextt.height = 0.7
    nextt.color = "black"
    
    inst.text = "Comenzamos nuevo bloque"
        
    nextt.text = "Cuando estes listo/a, pulsa espacio para comenzar "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys( keyList=['space'])
    core.wait(1)
    return

def block_loc_start(win):
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.height = 0.7
    inst.color = 'white'
    nextt = visual.TextStim(win, pos = [0,-6])
    nextt.height = 0.7
    nextt.color = "black"

    inst.text = "COMENZANDO BLOQUE"
    
        
    nextt.text = "Cuando estes listo/a, pulsa espacio para comenzar "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys( keyList=['space'])
    core.wait(1)
    return
    
    
def new_trial(win):
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.height = 0.7
    inst.color = 'white'
    inst.text = "Nuevo trial"
    inst.draw()
    win.flip()
    core.wait(0.75)
    return



    

def block_ID(win, rep, autodraw):
    inst = visual.TextStim(win, pos = [0,8])
    inst.wrapWidth = 20
    inst.height = 0.9
    inst.color = 'red'
    if rep == 'repeat':
        inst.text = "Secuencias repetidas"
    else:
        inst.text = "Secuencias parecidas"
    inst.autoDraw = autodraw   
    inst.draw()  
    
    
    


def end_experiment(win):
    inst = visual.TextStim(win, pos=[0,0], height = 1.2)
    inst.text = "Final de esta parte del experimento!! Avisa al investigador y pulsa -espacio- para continuar."     
    inst.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    return


def end_experiment_lot(win, corr_lotery, nblocks):
    inst = visual.TextStim(win, pos=[0,0], height = 1.2)
    inst.text = 'Final de esta parte del experimento!! El participante ha ganado ' + str(sum(corr_lotery)) + ' puntos de ' + str(nblocks) + '. Avisa al investigador'     
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
        col = 'red' if corr == -1 else 'green'
        number_text.text = trial_ix
        number_text.color = col
        
        for iframe in  range(time_number):
             inst.draw()
             number_text.draw()
             win.flip()
        
        allKeys = event.getKeys()
        
    allKeys = []
    
    inst.text = "Genial, has ganado! Pulsa -espacio- para continuar"
    if corr == -1:
        inst.text = "Oh! Qué mala suerte. Pulsa -espacio- para continuar"
        

    inst.draw()
    number_text.draw()
    win.flip()
        
    event.waitKeys(keyList = ["space"]) 
    return corr 
        
        
        