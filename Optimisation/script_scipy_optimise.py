'''
Created on 19 March 2021

@author: Andreagiovanni Reina.
IRIDIA, Universite' Libre de Bruxelles, Belgium.
'''

from scipy import integrate
import scipy.optimize
import numpy as np
#import matplotlib.pyplot as plt
import os, sys, getopt
import csv

def sys2opts(y, t, gamma_a, gamma_b, alpha_a, alpha_b, rho_a, rho_b, zeta_a, zeta_b):
    A, B = y
    dydt = [ (1-A-B)*gamma_a - A*alpha_a + A*(1-A-B)*rho_a - A*A*zeta_a,
             (1-A-B)*gamma_b - B*alpha_b + B*(1-A-B)*rho_b - B*B*zeta_b]
    return dydt

def sys3opts(y, t, gamma_a, gamma_b, gamma_c, alpha_a, alpha_b, alpha_c, rho_a, rho_b, rho_c, zeta_a, zeta_b, zeta_c):
    A, B, C = y
    dydt = [ (1-A-B-C)*gamma_a - A*alpha_a + A*(1-A-B-C)*rho_a - A*A*zeta_a,
             (1-A-B-C)*gamma_b - B*alpha_b + B*(1-A-B-C)*rho_b - B*B*zeta_b,
             (1-A-B-C)*gamma_c - C*alpha_c + C*(1-A-B-C)*rho_c - C*C*zeta_c]
    return dydt

def logisticFunc(value, slope, rate):
    return 2.0 * rate / (1.0 + np.exp(-slope*(value-0.5)))

def funcToOpitmise(x, *data):
    (tmax, zeta, rho) = data
        
    if rho == 'free':
        rhoRate = x[0]
        ai=1
    else:
        rhoRate = float(rho)
        ai=0
             
    if zeta == 'free':
        zetaRate = x[ai+1]
        zetaSlope = x[ai+5]
        i=ai+2
    else:
        zetaRate = 0
        zetaSlope = 0
        i=ai+1

    gammaRate=1
    alphaRate=x[ai]
    
    gammaSlope=x[i]
    alphaSlope=x[i+1]
    rhoSlope=x[i+2]
  
    all_errors=[]
    t = np.linspace(0, tmax, tmax*10+1)
    
    ## Two Options 
    vb=0.5
    
    ## Vary quality of A (va) in [0,1] 
    for va in np.linspace(0, 1, 11):
    
        gamma_a = logisticFunc(va, gammaSlope, gammaRate)
        gamma_b = logisticFunc(vb, gammaSlope, gammaRate)
        alpha_a = logisticFunc( 1-va, alphaSlope, alphaRate)
        alpha_b = logisticFunc( 1-vb, alphaSlope, alphaRate)
        rho_a   = logisticFunc( va, rhoSlope, rhoRate)
        rho_b   = logisticFunc( vb, rhoSlope, rhoRate)
        zeta_a  = logisticFunc( va, zetaSlope, zetaRate)
        zeta_b  = logisticFunc( vb, zetaSlope, zetaRate)
    
        ## Compute all initial states    
        all_initial_conditions = []
        for a in np.linspace(0, 1, 11):
            for b in np.linspace(0, 1, 11):
                if a+b <= 1.0001:
                    all_initial_conditions += [[a,b]]
        
        for y0 in all_initial_conditions:
            
            sol = integrate.odeint(sys2opts, y0, t, args=(gamma_a, gamma_b, alpha_a, alpha_b, rho_a, rho_b, zeta_a, zeta_b))
            
            error = 0 
            target_a = va/(va+vb)
            target_b = vb/(va+vb)
            for k in range(1,int(tmax)+1):
                idx = np.where(t==k)[0][0]
                error += abs(target_a-sol[idx, 0])
                error += abs(target_b-sol[idx, 1])
            error /= 2.0
            all_errors += [error]     
    
    
    ## Three Options 
    vb=0.6
    vc=0.3
    
    ## Vary quality of A (va) in [0,1] 
    for va in np.linspace(0, 1, 11):
    
        gamma_a = logisticFunc(va, gammaSlope, gammaRate)
        gamma_b = logisticFunc(vb, gammaSlope, gammaRate)
        gamma_c = logisticFunc(vc, gammaSlope, gammaRate)
        alpha_a = logisticFunc( 1-va, alphaSlope, alphaRate)
        alpha_b = logisticFunc( 1-vb, alphaSlope, alphaRate)
        alpha_c = logisticFunc( 1-vc, alphaSlope, alphaRate)
        rho_a   = logisticFunc( va, rhoSlope, rhoRate)
        rho_b   = logisticFunc( vb, rhoSlope, rhoRate)
        rho_c   = logisticFunc( vc, rhoSlope, rhoRate)
        zeta_a  = logisticFunc( va, zetaSlope, zetaRate)
        zeta_b  = logisticFunc( vb, zetaSlope, zetaRate)
        zeta_c  = logisticFunc( vc, zetaSlope, zetaRate)
        
        ## Compute all initial states
        all_initial_conditions = []
        for a in np.linspace(0, 1, 11):
            for b in np.linspace(0, 1, 11):
                for c in np.linspace(0, 1, 11):
                    if a+b+c <= 1.0001:
                        all_initial_conditions += [[a,b,c]]
        
        for y0 in all_initial_conditions:
            
            sol = integrate.odeint(sys3opts, y0, t, args=(gamma_a, gamma_b, gamma_c, alpha_a, alpha_b, alpha_c, rho_a, rho_b, rho_c, zeta_a, zeta_b, zeta_c))
            
            error = 0 
            target_a = va/(va+vb+vc)
            target_b = vb/(va+vb+vc)
            target_c = vc/(va+vb+vc)
            for k in range(1,int(tmax)+1):
                idx = np.where(t==k)[0][0]
                error += abs(target_a-sol[idx, 0])
                error += abs(target_b-sol[idx, 1])
                error += abs(target_c-sol[idx, 2])
            error /= 3.0
            all_errors += [error]     
    
    return np.mean(all_errors)
        
    # plt.plot(t, sol[:, 0], 'b', label='A(t)')
    # plt.plot(t, sol[:, 1], 'g', label='B(t)')
    # plt.plot(t, sol[:, 2], 'r', label='C(t)')
    # plt.axhline(target_a)
    # plt.axhline(target_b)
    # plt.ylim(0,1)
    # plt.legend(loc='best')
    # plt.xlabel('t')
    # plt.grid()
    # plt.show()
    
def myCallback(xk,convergence):
    print(xk)
    print(convergence)
    return(False)

if __name__ == "__main__":
    #print(funcToOpitm( [1,2,3,4,5,6,7] ))
    argv = sys.argv
    print(argv)
    
    zeta = "free"
    rho = "free"
    maxTime= None
    for i,arg in enumerate(argv):
        if arg == '-h':
            print(argv[0] + ' -t <maxTime> -z <("free"/0) zeta> -r <rho> -o <output directory>') 
            sys.exit()
        elif arg in ('-t', '-time'):
            maxTime = argv[i+1]
        elif arg in ('-r', '-rho'):
            rho = argv[i+1]
        elif arg in ('-z', '-zeta'):
            zeta = argv[i+1]
        elif arg in ('-o', '-output'):
            outdir = argv[i+1]

    if maxTime == None:
        print(argv[0] + ' -t <maxTime> -z <("free"/0) zeta> -r <rho> -o <output directory>')
        sys.exit(2)
    
    print( "Params. maxTime=", str(maxTime), " zeta=", str(zeta), " rho=", str(rho), " output=", str(outdir) )
    print("Optimisation process started...")
    args=(float(maxTime), zeta, rho)
    if rho == "free":
        if zeta == "free":
            ####### with inhibition -- free recruitment ######
            bounds = [(0,1000),(0.0000000001,10),(0.0000000001,1000),(0,5),(0,5),(0,5),(0,5)]
            results=scipy.optimize.differential_evolution(funcToOpitmise, bounds, maxiter=200, popsize=10, callback=myCallback, args=args) #, workers=-1)
            line = [maxTime, results.x[0], results.x[1], results.x[2], results.x[3], results.x[4], results.x[5], results.x[6], results.fun, int(results.success) ]
        else:
            ####### without inhibition -- free recruitment ######
            bounds = [(0,1000),(0.0000000001,10),(0,5),(0,5),(0,5)]
            results=scipy.optimize.differential_evolution(funcToOpitmise, bounds, maxiter=200, popsize=10, callback=myCallback, args=args) #, workers=-1)
            line = [maxTime, results.x[0], results.x[1], 0, results.x[2], results.x[3], results.x[4], 0, results.fun, int(results.success)]
    else:
        if zeta == "free":
            ####### with inhibition -- fixed recruitment ######
            bounds = [(0.0000000001,10),(0.0000000001,1000),(0,5),(0,5),(0,5),(0,5)]
            results=scipy.optimize.differential_evolution(funcToOpitmise, bounds, maxiter=200, popsize=5, callback=myCallback, args=args) #, workers=-1)
            line = [maxTime, float(rho), results.x[0], results.x[1], results.x[2], results.x[3], results.x[4], results.x[5], results.fun, int(results.success) ]
        else:
            ####### without inhibition -- fixed recruitment ######
            bounds = [(0.0000000001,10),(0,5),(0,5),(0,5)]
            results=scipy.optimize.differential_evolution(funcToOpitmise, bounds, maxiter=200, popsize=10, callback=myCallback, args=args) #, workers=-1)
            line = [maxTime, float(rho), results.x[0], 0, results.x[1], results.x[2], results.x[3], 0, results.fun, int(results.success)]

    print(results)
    print(line)
        
    # Save file
    rhoStr = "f" if rho=="free" else str(rho) 
    zetaStr = "f" if zeta=="free" else str(zeta) 
    with open(outdir + '/optimalRates_response-free_rates-1DfA' + rhoStr + 'R' + zetaStr + 'I' + '.txt', 'a+') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(line)
        
    print("Process ended.")


