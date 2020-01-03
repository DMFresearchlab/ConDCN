from scipy import signal, stats
from numpy.random import vonmises
import numpy as np
import matplotlib.pyplot as plt

# What his function does is to create millions of X stimuli trials with randomly sampled orientations from a vonmises circular distribution.
# Then I use a triangle function to calculate the decision variable for each orientation. Then I eliminate those trials that have more or less than
# 1.5 standard deviation in orientation variability. Finally, to reduce the size of the matrix file (that is very space consuming), I use a for loop to
# subsample 1000 trials per level of decision variable from -1 to 1 in jumps of 0.1.
# There are multiple plotting functions that I used to visualize the distribution of the data when I create the stmuli for the
# first time. In principle you dont need to run this function unless the file that I generated the first time is gone.

def stim_creation(nstim, outfile):

    print("This process might take a couple of minutes. It can collapse if you dont have enough space in your computer")
    sample_size = 10000000 # I am going to generate thousands of randomly oriented trials
    orientations = vonmises(np.deg2rad(0), 0, [sample_size,nstim]) # generate a lot of random orientations in groups of 8
    orientations = np.around(orientations, decimals=3) # to save some space just keep 3 decimals in memory
    #plt.figure()
    #plt.hist(np.rad2deg(orientations))

    #max_ang = 179
    #min_ang = 0
    #orientations = (max_ang - min_ang ) * np.random.random_sample((sample_size, nstim)) + min_ang

    # Using vonmisses distributions
    # from scipy import signal
    #kappa_angle = 20
    #orientations = np.hstack( (vonmises(np.deg2rad(conditions[cond, 0]), kappa_angle, [sample_size,4]), vonmises(np.deg2rad(conditions[cond, 1]), kappa_angle, [sample_size,4])))
    #orientations = orientations[np.random.permutation(nstim)] # shuffle the order
    # Deviation of orientations
    #conditions = np.matrix('0 90; 45 135')
    #cond = 0
    #orientations = np.hstack( (vonmises(np.deg2rad(conditions[cond, 0]), kappa_angle, 4), vonmises(np.deg2rad(conditions[cond, 1]), kappa_angle, 4)))
    #orientations = orientations[np.random.permutation(nstim)] # shuffle the order
    #orientations = vonmises(np.deg2rad(conditions[cond, 0]), kappa_angle, 4), vonmises(np.deg2rad(conditions[cond, 1]), kappa_angle, 4)))

    # Visualize triangle function to estabkish the decision value of each stim
    #degs = np.deg2rad(np.linspace(1, 180, 180))
    #triangle = signal.sawtooth(4*degs, 0.5)
    #plt.figure()
    #plt.plot(np.rad2deg(degs ), triangle)
    #deg_x = np.arange(0, 181, step=45)
    #plt.xticks(deg_x)
    #[plt.axvline(_x, linewidth=1, color='black') for _x in deg_x]
    #plt.axhline(y=0.0, color='r')


    decision_var = signal.sawtooth(4 * ((orientations)), 0.5)  # DU decision variable
    decision_var_T = np.mean(decision_var, 1)
    decision_var_std = np.std(decision_var, 1)
    orientation_var_std = stats.circstd((orientations),axis = 1)

    # Visualize parameters distributions (decision, orientation, variability, etc)
    plt.figure()
    plt.plot(decision_var_T, decision_var_std, 'ro', markersize=0.1)
    plt.xlabel("Mean decision of trial")
    plt.ylabel("Stimuli trial decision variability")
    plt.figure()
    plt.plot(decision_var_T, orientation_var_std, 'go', markersize=0.1) # the shape is weird but because I have back transformed the negative orientations

    #LEts balance the standard deviation
    hbound = np.mean(orientation_var_std) + 1.5*np.std(orientation_var_std)
    lbound = np.mean(orientation_var_std) - 1.5*np.std(orientation_var_std)

    #Lets filter out low and very high variable trials in orientations (Based on 2*STD of the STD)
    indexes = np.where((orientation_var_std > lbound) & (orientation_var_std < hbound))
    orientation_var_std = orientation_var_std[indexes]

    decision_var = decision_var[indexes]

    decision_var_T = decision_var_T[indexes]
    decision_var_std = decision_var_std[indexes]
    orientations = orientations[indexes,]
    orientations = np.squeeze(orientations)


    # Visualize parameters distributions (decision, orientation, variability, etc)
    plt.figure()
    plt.plot(decision_var_T, decision_var_std, 'ro', markersize=0.5)
    plt.figure()
    plt.plot(decision_var_T, orientation_var_std, 'go', markersize=0.5)


    # Here I am going to create a smaller matrix in order to reduce the number of space taken by the computer.
    # The idea is to take randomly 1000 trials around each decision variable level.
    ntrials_per_dec = 1000
    range_of_dec = np.arange(-0.7, 0.8, 0.025)
    range_of_dec = np.round(range_of_dec, decimals = 3)


    orient_matrix = np.array([])

    for idec in range_of_dec:
        sel_trials = np.where((decision_var_T > idec-0.025) & (decision_var_T < idec+0.025))  # select a matrix with decisions in groups of difficult
        sel = orientations[sel_trials[0]][0:ntrials_per_dec,]
        orient_matrix = np.vstack([orient_matrix, sel]) if orient_matrix.size else sel # concatenate matrixes when not knowing the number of cells in empty array


    orient_matrix = np.around(orient_matrix, decimals = 3)
    decision_var = signal.sawtooth(4*(orient_matrix), 0.5) # DU decision variable
    decision_var_T = np.mean(decision_var , 1)
    orientation_var_std = stats.circstd((orient_matrix),axis = 1)
    decision_var_std = np.std(decision_var, 1)

    plt.figure()
    plt.hist(np.rad2deg(orient_matrix))
    plt.xlabel("Orientation (degrees)")

    plt.hist(decision_var_std)
    plt.figure()
    plt.plot(decision_var_T,decision_var_std,'ro', markersize=0.5)
    plt.xlabel("Mean decision of trial")
    plt.ylabel("Stimuli trial decision variability")
    plt.figure()
    plt.plot(decision_var_T,orientation_var_std,'go', markersize=0.5)
    plt.xlabel("Mean decision of trial")
    plt.ylabel("Stimuli trial orientation variability (in 360 degrees)")
    np.save(outfile, orient_matrix) # save matrix to data to be loaded in future experiments



    decision_var = signal.sawtooth(4*(orientations), 0.5) # DU decision variable
    decision_var_T = np.mean(decision_var , 1)
    orientation_var_std = stats.circstd(orientations,axis = 1)
    decision_var_std = np.std(decision_var, 1)

    # plt.figure()
    # plt.plot(decision_var_T,decision_var_std,'ro', markersize=0.5)
    # plt.xlabel("Mean decision of trial")
    # plt.ylabel("Stimuli trial decision variability")
    # plt.figure()
    # plt.plot(decision_var_T,orientation_var_std,'go', markersize=0.5)
    # plt.xlabel("Mean decision of trial")
    # plt.ylabel("Stimuli trial orientation variability (in 360 degrees)")
    return orient_matrix