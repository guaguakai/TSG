import util
from gurobipy import *
import time
import numpy as np
import random
import math
from StaffResourceAllocation import LPsolver
import csv
from Method1 import solve as solve_m1
from Method2 import solve as solve_m2

from DesignYNcombined import KStrategiesYNcomb
from DesignYNcombined import KStrategiesYNBnew

from Results import randomSetting

import pickle

if __name__ == "__main__":
    
    instance = 3
    if (instance == 1):
        W = 10 # number of time windows
        AI = 3 # interval in which passengers are arriving
        K = 20 # number of passenger types
        R = 5 # number of resources
        mR = 3 # max number of reosurces
        M = 3 # number of attack methods
        P = 15 # number of staff
        shift = 3 # d
        Q = 2 # number of strategies
        maxT = 10

    if (instance == 2):
        W = 5 # number of time windows
        AI = 3 # interval in which passengers are arriving
        K = 10 # number of passenger types
        R = 5 # number of resources
        mR = 3 # max number of reosurces
        M = 2 # number of attack methods
        P = 10 # number of staff
        shift = 1 # d
        Q = 2 # number of strategies
        maxT = 10
    
    if (instance == 3):
        W = 3 # number of time windows
        AI = 1 # interval in which passengers are arriving
        K = 10 # number of passenger types
        R = 3 # number of resources
        mR = 2 # max number of reosurces
        M = 3 # number of attack methods
        P = 10 # number of staff
        shift = 1 # d
        Q = 2 # number of strategies
        maxT = 10    
        
    teams = util.generateAllTeams(R, mR)
    
    Z = 20 # number of runs
    ZQ = 10 # max number of Q
    
    obj_relax = np.zeros((Z,ZQ))
    obj_yn = np.zeros((Z,ZQ))
    obj_finalm2= np.zeros((Z,ZQ))
    obj_finalm1 = np.zeros((Z,ZQ))

    time_relax = np.zeros((Z,ZQ))
    time_yn = np.zeros((Z,ZQ))
    time_final = np.zeros((Z,ZQ))
    
    for z in range(Z):
        seed = random.randint(1,1000)
        for zq in range(ZQ):
            Q = zq+1;  
            
            # Construct random instances
            resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, AI, K ,R, mR, M, P, teams, shift)
        
            minr = np.zeros((W,R))   
            
            
            obj, rt = solve_m1(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi)
            obj_relax, objyn1, obj, rt = solve_m2(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, TeamConstr)
            
            obj_final[z][zq], rt, t3,ni,oi_value,q_tem,o_temp  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)
            
            time_final[z][zq] = time.time() - start_time_final
            
            sol = [obj_relax[z][zq], obj_yn[z][zq], obj_final[z][zq], time_relax[z][zq], time_yn[z][zq], time_final[z][zq], seed]
            with open("resultsdiffQ.csv", 'a') as csvfile:
                w = csv.writer(csvfile)
                strlst = [str(s) for s in sol]
                w.writerow(strlst)  
            
       
    print obj_relax, obj_yn, obj_final
    print time_relax, time_yn, time_final
    
    pickle.dump((obj_relax, obj_yn, obj_final, time_relax, time_yn, time_final, seed), open("resultsDifferentQ_0125.data", "wb"))
    # obj_relax, obj_yn, obj_final, time_relax, time_yn, time_final = pickle.load(open("resultsDifferentQ_0125.data", "rb"))