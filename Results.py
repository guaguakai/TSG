import util
from gurobipy import *
import time
import numpy as np
import random
import math
from StaffResourceAllocation import LPsolver

from DesignYNcombined import KStrategiesYNcomb
from DesignYNcombined import KStrategiesYNBnew
from DesignYNcombined import randomSetting

from DesignYNcombinedInteger import KStrategiesYNcombInteger

if __name__ == "__main__":
    
    print "======================== main ======================================"
    
    W = 10 # number of time windows
    K = 15 # number of passenger types
    R = 5 # number of resources
    mR = 3 # max number of reosurces
    M = 2 # number of attack methods
    P = 20 # number of staff
    shift = 3 # d
    Q = 2
    teams = util.generateAllTeams(R, mR)
    maxT = 25


    # ================= random generate game setting ===========================
    seed = 2345
    
    resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)

    minr = np.zeros((W,R))
    start_time = time.time()
    
    print "============================ FULL LP relaxation =============================="
    
    
    obj_relax, n_value0, overflow_value, y_value, s_value0, p,z_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=0, OverConstr = 0)
    

    for w in range(W):
        for r in range(R):
            minr[w][r] = math.floor(y_value[w][r])
            #print y_value[w][r]
    
    QT = 5
    q = [[0.1, 0.9],[0.2, 0.8],[0.25,0.75],[0.3333,0.6667],[0.5,0.5]]
    
    objyn = np.zeros(QT)
    obj = np.zeros(QT)
    for j in range(QT):
        objyn[j], n, ns,ys,z_value,p,s,y = KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q[j][:], ar, phi, integer=0, OverConstr=False, OverConstr2=False)
    
        minn = np.zeros((Q,W,T,K))
        
        for i in range(Q):
            for w in range(W):
                for k in range(K):
                    sum = 0 
                    for t in range(T):
                        minn[i][w][t][k] = math.floor(ns[i][w][t][k])
                        sum += math.floor(ns[i][w][t][k])
        
        obj[j], rt, t3  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)
    
    
    minr = np.zeros((W,R))
    
    objynint, n, ns,ys,z_value,p,s,y = KStrategiesYNcombInteger(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q[QT-1][:], ar, phi, integer=0, OverConstr=False, OverConstr2=False)
    
    minn = np.zeros((Q,W,T,K))
    
    for i in range(Q):
        for w in range(W):
            for k in range(K):
                sum = 0 
                for t in range(T):
                    minn[i][w][t][k] = math.floor(ns[i][w][t][k])
                    sum += math.floor(ns[i][w][t][k])
    
    objint, rt, t3  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)   
    #print obj_relax, objyn1, obj1

