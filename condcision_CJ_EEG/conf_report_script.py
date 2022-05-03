#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
confidence reports

(you can add this to the display stimulus loop in order to collect confidence in each trial)
@author: alex
"""
confidence_rating = ["AL AZAR", "POCO", "MUCHO"]
confidence_color = np.array([[1, 0, 0] , [ 0.5, 0.5 , 0 ], [ 0, 1 , 0]]) #rojo, amarillo y verde


conf_resp = np.random.randint(3) # lets initialize the confidence randomly in each trial
    resp_option(resp_ang,  confidence_color[conf_resp,])
    fixation()  # text - Choose your level of confidence
    #exp.draw_this_text(win, "Cómo estás de segur@ en tu respuesta?", -3)
    exp.draw_this_text(win, confidence_rating[conf_resp], -2)
    win.flip()

    #press buttons to select confidence level
    thisKey = None

    while thisKey != 'space':
        allKeys = event.waitKeys()
        resp_cond = 0
        for thisKey in allKeys:
            if thisKey == 'left':
                resp_cond = -1
            elif thisKey == 'right':
                resp_cond = +1
            elif thisKey in ['q', 'escape']:                  
                win.quit()  # abort experiment
                core.quit()  # abort experiment
        event.clearEvents()  # clear other (eg mouse) events - they clog the buffer

        conf_resp = conf_resp + resp_cond  # lets initialize the confidence randomly in each trial
        if conf_resp < 0:
            conf_resp = 2
        elif conf_resp > 2:
            conf_resp = 0

        resp_option(resp_ang, confidence_color[conf_resp,])
        fixation()  # text - Choose your level of confidence
        # exp.draw_this_text(win, "Cómo estás de segur@ en tu respuesta?", -3)
        exp.draw_this_text(win, confidence_rating[conf_resp], -2)
        print(conf_resp)
        win.flip()

        trials.data.add('conf_resp', conf_resp)