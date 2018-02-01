import util
from gurobipy import *
import time
import numpy as np
import random
import math
from StaffResourceAllocation import LPsolver

from DesignYNcombined import KStrategiesYNcomb
from DesignYNcombined import KStrategiesYNBnew
from DesignYNcombined import randomSetting as rs

from DesignYNcombinedInteger import KStrategiesYNcombInteger

def randomSetting(seed, W, AI, K ,R, mR, M, P, teams, shift):
    # ========================== Random Seed ===================================
    #seed = 3234
    #seed = random.randint(1,10000)
    np.random.seed(seed)
    random.seed(seed)

    T = len(teams)

    resource2team = util.resourceTeamDict(R, T, teams)


    Er, C = util.genResources(R, M, 200)
    E = util.computeTeamsRate(R, M, T, teams, Er)
    print E     


    # suppose this is a zero-sum game
    U_plus = [] # covered (plus) utility of the defender
    U_minus = [] # uncovered (minus) utility of the defender
    for i in range(K):
        tmp_plus_utility = 0 #random.randint(100,500)
        U_plus.append(tmp_plus_utility)
        tmp_minus_utility = -random.randint(500,2000)
        U_minus.append(tmp_minus_utility)


    N_wk = [[ 0 for k in range(K)] for w in range(W)] # N_wk[w][k] represents the number of people getting in time window w with type k
    

    N_wk = np.zeros((W,K))        
    for k in range(K):
        startK = random.randint(AI,W)
        for w in range(startK-AI,startK):
            large_or_small = random.random()
            if large_or_small > 0.5:
                tmp_N = random.randint(150, 200)
            else:
                tmp_N = random.randint(10, 30)
            N_wk[w][k] = tmp_N


    mr = np.random.randint(5, 15, R)

    ar = np.random.randint(1, 3, R)


    phi = np.random.rand(R) # phi[r] overflow penalty


    return resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi

if __name__ == "__main__":
    
    print "======================== main ======================================"
     
    instance = 2
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
        P = 5 # number of staff
        shift = 2 # d
        Q = 1 # number of strategies
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
        Q = 1 # number of strategies
        maxT = 10    
   
    teams = util.generateAllTeams(R, mR)
    
    Z = 10 # number of runs
    ZQ = 2 # max number of Q
    
    obj_relax = np.zeros((Z,ZQ))
    obj_yn = np.zeros((Z,ZQ))
    obj_final = np.zeros((Z,ZQ))
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
            
            
            start_time_relax = time.time()
            # Solve full LP relaxation
            obj_relax[z][zq], n_value0, overflow_value, y_value, s_value0, p,z_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=0, OverConstr = 0)
            
            time_relax[z][zq] = time.time() - start_time_relax
               
            for w in range(W):
                for r in range(R):
                    minr[w][r] = math.floor(y_value[w][r])
                    #print y_value[w][r]
                    
            
            q = np.zeros(Q)
            for i in range(Q):
                q[i] = float(1)/Q    
            
            # find Q strategies for y, using binaries and lower bound minr, n integer for each ys.    
            start_time_yn = time.time()
            
            obj_yn[z][zq], n, ns,ys,z_value,p,s,y = KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q, ar, phi, integer=0, OverConstr=False, OverConstr2=False)
            
            time_yn[z][zq] = time.time() - start_time_yn
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
            
            obj_final[z][zq], rt, t3,ni,oi_value,q_temp,o_temp  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)
            
            time_final[z][zq] = time.time() - start_time_final
            
            
    """
    # Solve second step for different q-vectors (if Q=2 or Q=3)
    q = [[0.1667,0.8333],[0.2, 0.8],[0.4,0.6],[0.25,0.75],[0.3333,0.6667],[0.5,0.5]] # if Q=2
    q = [[0.1667,0.1667,0.6666],[0.1667,0.3333,0.5],[0.2,0.2,0.6],[0.2,0.4,0.4],[0.25,0.25,0.5],[0.3333,0.3333,0.3334]] # if Q=3
    
    QT = len(q)
    
    
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
        
        obj[j], rt, t3,ni,oi_value  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)
    
    
    # solve for MIP with ys integer (no rounding) and n marginal. q is given. 
    minr = np.zeros((W,R))
    
    q = np.zeros(Q)
    for i in range(Q):
        q[i] = float(1)/Q 
        
    objynint, n, ns,ys,z_value,p,s,y = KStrategiesYNcombInteger(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q[QT-1][:], ar, phi, integer=0, OverConstr=False, OverConstr2=False)
    
    minn = np.zeros((Q,W,T,K))
    
    for i in range(Q):
        for w in range(W):
            for k in range(K):
                sum = 0 
                for t in range(T):
                    minn[i][w][t][k] = math.floor(ns[i][w][t][k])
                    sum += math.floor(ns[i][w][t][k])
    
    objint, rt, t3,ni  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)   

    """
    print obj_relax, obj_yn, obj_final
    print time_relax, time_yn, time_final
