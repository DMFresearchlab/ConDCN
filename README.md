# Condcision v1.0

This is the first version of the condcision project (confirmation bias project)
It includes a function to measure the participant threshold using a staircase, a functional localizer (it should be reviewed to check whether it is compatible with this new experiment) and the main experiment.

__Main experiment logic__

Participants will be presented with 3 similar grating sequences per trial. Their task is to judge whether the mean gratings orientation is closer to the diagonal or the cardinal axis. They have 3 opportunities to respond correctly. At the end of each triplet trial, participants will receive feedback about their performance.
Participants are supposed to earn a small ammount of money for each correct response (max 3 euros).
There is a staircase constantily updating the difficulty as a function of the response in the first trial sequence of each triplet.
The response mapping change randomly in each trial sequence presentation to minimize response biases that could explain the decision biases.
The fixation point duration previous to the surround circle presentation is jittered at the tr3ial level.

# Condcision_IEEG v1.2

__Condcision v1.0__ was used for the first behavioral experiment. In __Condcision_IEEG__ (v1.2) I have included the triggers for the iEEG/EEG versions of the experiment. I have also included two behavioral experiments (compatible with EEG). Control_exp, where three different stimuli sequences are presented 3 times in each trial (with the same average DV). The objective of this experiment is to assess whether the pattern of stimuli evaluation (disconfirmatory strategy), takes places only when exactly the same stimuli are encoded in memory.

and the Blocked_exp where there are intermixed blocks with 75% trials of random sign of DV and stimuli values (same average ABSOLUTE DV) and blocked trials blocks, where the same sequence (75%)is presented three trimes in a row (as in the original condcision). This experiment has the objective of seeing whether people can switch the information integration strategy flexibly, as a function of the "volatility" of the environment,

Coded by Pérez-Bellido 2020




