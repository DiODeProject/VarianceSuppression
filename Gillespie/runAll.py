'''
Created on 09 Jan 2019

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import selfInhiGill
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import axhline
import math 
import sys, os

DEBUG=False
batch=False

# The script accepts 3 optional arguments:
# * output directory (default is 'results')
# * positive feedback strength (i.e. rho)  (default is all values in range [1,200])
# * system size S  (default is S=200)
if __name__ == '__main__':
    if DEBUG: 
        print("Process Started")
    
    write_evo = True
    plot_evo = False
    
    outputdir = sys.argv[1] if sys.argv[1] else 'results' 
    optimalZetaFilename = os.path.join(os.path.dirname(__file__), '../published_data/optimalZeta_withV_wTimeS.txt')
    optimalZetas = np.loadtxt(optimalZetaFilename, skiprows=0, delimiter='\t', usecols=(0,1,2,3))
    #print(optimalZetas)
    
    repetitions = 1000
    
    # Input params
    g=1
    rhos= [ float(sys.argv[2]) ] if sys.argv[2] else np.arange(1, 200.1, 2)
        
    h=0.5
    quorum = -1
    
    discoMax = 1.0
    
    # Number of agents
    N = int( sys.argv[3] ) if sys.argv[3] else 200
    # Experiment time length
    T = 10
    
    #options = [2, 4, 6]
    options = [3]
    models = ['IFD', 'VASI','NOSOC']
    
    for r in rhos:
        for m,model in enumerate(models):
            for n in options:
                # Set initial state
                state = [N] + [0]*n
                
                # Set the quality values
                values=[0.75,0.5]
                if n==3: values = [0.75,0.5,0.25]
                if n==4: values = [0.8,0.6,0.4,0.2]
                
                print('model: ' + str(model))
                print('values: ' + str(values))
            
                # Define the transition rates 
                gammas = [g*values[i] for i in range(n)]
                #alphas = [k/v for v in values]
                alphas = [0.001]*n
                if model == 'CDCI':
                    rhos = [r*v for v in values]
                    sigmas = [r*v for v in values]
                    zetas  = [0]*n
                elif model == 'IFD':
                    #rhos = [0]*n
                    rhos = [r]*n
                    sigmas = [0]*n
                    zetas  = [0]*n
                elif model == 'VASI':
                    rhoVal=2*r
                    rhos = [rhoVal*v for v in values]
                    sigmas = [0]*n
                    #zetas  = [h]*n
                    #optZeta = [v[1] for v in optimalZetas if v[0]==r ][0]
                    valueB = 0.5
                    optZeta = [v[3] for v in optimalZetas if math.isclose(v[0],rhoVal) and math.isclose(v[1],values[0]) and math.isclose(v[2],valueB) ][0]
                    zetas  = [optZeta]*n
                elif model == 'NOSOC':
                    rhos = [0]*n
                    sigmas = [0]*n
                    zetas  = [0]*n
                
                # Create filename
                finalStateFile = outputdir + '/fs_N-' + str(N) + '_g-' + str(g) + '_r-' + str(r) + '_n-' + str(n) + '_v-' + str(values[0]) + '-' + str(values[1]) + '_' + str(model) + '.txt'
                    
                print('gammas: ' + str(gammas))
                print('rhos: ' + str(rhos))
                print('sigmas: ' + str(sigmas))
                print('zetas:  ' + str(zetas))
                
                if plot_evo: write_evo = True
                for i in range(0,repetitions):
                    temporalEvolution = outputdir + '/evo_N-' + str(N) + '_g-' + str(g) + '_r-' + str(r) + '_n-' + str(n) + '_v-' + str(values[0]) + '-' + str(values[1]) + '_' + str(model) + '_' + str(i) + '.txt' if write_evo else 'none'
                    rnd_seed = 601 + i*1000
                    plt = selfInhiGill.runGillespie(state, T, N, gammas, alphas, rhos, sigmas, zetas, rnd_seed, finalStateFile, temporalEvolution, plot_evo, extraLog=[n], quorum=quorum)
            
                if (plot_evo):
                    if model == 'IFD':
                        colours = ['k', 'b', 'r', 'g', 'c', 'm', 'y', 'fuchsia', 'aqua', 'peru', 'lime']
                        for i,v in enumerate(values):
                            axhline(y=v*N/sum([values[i] for i in range(n)]), color=colours[i+1], linestyle='--', linewidth=1.4)
                    elif model == 'VASI':
                        colours = ['k', 'b', 'r', 'g', 'c', 'm', 'y', 'fuchsia', 'aqua', 'peru', 'lime']
                        for i,v in enumerate(values):
                            axhline(y=v*N/sum(values), color=colours[i+1], linestyle='--', linewidth=1.5)
                    plt.show()
        
    if DEBUG:
        print("Process Ended")
    
