'''
Created on 09 Jan 2019

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import numpy as np
import sys
import os
import matplotlib.pyplot as plt

DEBUG=False

def gillespieStep(state, N, gammas, alphas, rhos, sigmas, zetas, vectorsOfChange):
    # Computing the probabilities of change
    probabilitiesOfChange = []
    for i in range(0, len(state)-1 ): # PAY ATTENTION THAT i IS NOT THE CORRECT state[i] BUT IT MUST BE state[i+1], WHILE WORKS FOR gammas[i], alphas[i], rhos[i] 
        # Discovery
        probabilitiesOfChange.append( state[0]*gammas[i] )
        # Abandonment
        probabilitiesOfChange.append( state[i+1]*alphas[i] )
        # Recruitment
        probabilitiesOfChange.append( state[0]*state[i+1]*rhos[i]/N )
        #Cross-inhibition
        for j in range(0, len(state)-1):
            if (i==j):
                continue
            probabilitiesOfChange.append( state[i+1]*state[j+1]*sigmas[j]/N )
        # Self-inhibition
        probabilitiesOfChange.append( state[i+1]*state[i+1]*zetas[i]/N )
    
    probSum = sum(probabilitiesOfChange)
    if probSum == 0:
        return 10 #returning a 'large' timeInterval without state change
    timeInterval = np.random.exponential( 1/probSum ) 

    # Selecting the occurred reaction in a randomly, proportionally to their probabilities
    bottom = 0.0
    # Get a random between [0,1) (but we don't want 0!)
    reaction = 0.0
    while (reaction == 0.0):
        reaction = np.random.random_sample()
    # Normalising probOfChange in the range [0,1]
    probabilitiesOfChange = [pc/probSum for pc in probabilitiesOfChange]
    index = -1
    for i, prob in enumerate(probabilitiesOfChange):
        if ( reaction >= bottom and reaction < (bottom + prob)):
            index = i
            break
        bottom += prob
        #print(i, prob)
                
    if (index == -1):
        print("Transition not found. Error in the algorithm execution.")
        sys.exit()
    state += np.array(vectorsOfChange[index])
    return(timeInterval)

def runGillespie(state, T, N, gammas, alphas, rhos, sigmas, zetas, rnd_seed, finalStateFile, temporalEvolution, plot_evo, extraLog, quorum=-1):
    np.random.seed(rnd_seed)
    n = len(gammas)
    state = np.array(state)
    t = 0
    if DEBUG:
        print("t: ", t, "state: ", state)
    
    # Opening output file if needed
    if (temporalEvolution != "none"):
        os.makedirs(os.path.dirname(temporalEvolution), exist_ok=True)
        evoStream = open(temporalEvolution, "w+")
        out = '{:.20f}'.format(t) + "\t" + '\t'.join(str(x) for x in state) + "\n"
        evoStream.write(out)
    
    # Creating the list of vector of change
    vectorsOfChange = []
    for i in range(n):
        # Positive change
        plus = [-1] + [0]*n
        plus[i+1] = 1
        # Negative change
        negative = [1] + [0]*n
        negative[i+1] = -1
        
        # Discovery
        vectorsOfChange.append( plus )
        # Abandonment
        vectorsOfChange.append( negative )
        # Recruitment
        vectorsOfChange.append( plus )
        #Cross-inhibition
        for _ in range(n-1):
            vectorsOfChange.append( negative )
        # Self-inhibition
        vectorsOfChange.append( negative )      
        
    while t < T:
        t += gillespieStep(state, N, gammas, alphas, rhos, sigmas, zetas, vectorsOfChange)
        if DEBUG:
            print("t: ", t, "state: ", state)
        if (temporalEvolution != "none"):
            out = str(t) + "\t" + '\t'.join(str(x) for x in state) + "\n"
            evoStream.write(out)
        
        ## Checking each timestep if the quorum is reached
        if (quorum > 0):
            quorum_reached = False
            for i in np.arange(1,len(state)):
                if (state[i] > N*quorum):
                    quorum_reached = True
                    break
            if (quorum_reached):
                break
        
    if (finalStateFile != "none"):
        os.makedirs(os.path.dirname(finalStateFile), exist_ok=True)
        with open(finalStateFile, "a") as f:
            out = '\t'.join(str(x) for x in extraLog)
            if (len(out)>0): out += '\t'
            out += str(t) + "\t" + '\t'.join(str(x) for x in state) + "\n"
            f.write(out)
    
    if (plot_evo):
        if (temporalEvolution == "none"):
            print("WARNING! - to plot the temporal evolution, please specify a temporalEvolution file (e.g., a temp-file)")
        if DEBUG:
            print("Plotting temporal evolution...")
        # Moving at the beginning of the file
        evoStream.seek(0)
        matStr = evoStream.read()
        matStr = matStr.replace("\t"," ")
        matStr = matStr.replace("\n",";")
        matStr = "" + matStr[0:len(matStr)-1] + ""
        mat = np.matrix(matStr)
        
        colours = ['k', 'b', 'r', 'g', 'c', 'm', 'y', 'fuchsia', 'aqua', 'peru', 'lime']
        for c in range(len(state)):
            plt.plot(mat[0:,0],mat[0:,c+1],colours[c],linewidth=2)
        plt.xlim((0,T))
        plt.ylim((0,N))
        plt.xlabel('time')
        plt.ylabel('populations')
        plt.draw()
        return(plt)
        
    if DEBUG:
        print("Gillespie run ended")
    
