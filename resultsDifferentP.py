import util
from gurobipy import *
import time
import numpy as np
import random
import math
from StaffResourceAllocation import LPsolver

from DesignYNcombined import KStrategiesYNcomb
from DesignYNcombined import KStrategiesYNBnew

from DesignYNcombinedInteger import KStrategiesYNcombInteger

from Results import randomSetting

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
    
    Z = 10# number of runs
    ZP = 5 # max number of P
    
    obj_relax = np.zeros((Z,ZP))
    obj_yn = np.zeros((Z,ZP))
    obj_final = np.zeros((Z,ZP))
    time_relax = np.zeros((Z,ZP))
    time_yn = np.zeros((Z,ZP))
    time_final = np.zeros((Z,ZP))
    
    for z in range(Z):
        seed = random.randint(1,1000)
        for zp in range(ZP):
            P = 5*zp 
            
            # Construct random instances
            resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, AI, K ,R, mR, M, P, teams, shift)
        
            minr = np.zeros((W,R))   
            
            
            start_time_relax = time.time()
            # Solve full LP relaxation
            obj_relax[z][zp], n_value0, overflow_value, y_value, s_value0, p,z_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=0, OverConstr = 0)
            
            time_relax[z][zp] = time.time() - start_time_relax
               
            for w in range(W):
                for r in range(R):
                    minr[w][r] = math.floor(y_value[w][r])
                    #print y_value[w][r]
                    
            
            q = np.zeros(Q)
            for i in range(Q):
                q[i] = float(1)/Q    
            
            # find Q strategies for y, using binaries and lower bound minr, n integer for each ys.    
            start_time_yn = time.time()
            
            obj_yn[z][zp], n, ns,ys,z_value,p,s,y = KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q, ar, phi, integer=0, OverConstr=False, OverConstr2=False)
            
            time_yn[z][zp] = time.time() - start_time_yn
            minn = np.zeros((Q,W,T,K))
            
            for i in range(Q):
                for w in range(W):
                    for k in range(K):
                        sum = 0 
                        for t in range(T):
                            minn[i][w][t][k] = math.floor(ns[i][w][t][k])
                            sum += math.floor(ns[i][w][t][k])
            
            # Find integer solution (using binaries and rounde) for each ns, given ys, p and s
            
            start_time_final = time.time()
            
            obj_final[z][zp], rt, t3,ni,oi_value  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)
            
            time_final[z][zp] = time.time() - start_time_final
       
    print obj_relax, obj_yn, obj_final
    print time_relax, time_yn, time_final