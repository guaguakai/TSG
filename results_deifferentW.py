import util
from gurobipy import *
import time
import numpy as np
import random
import math
from StaffResourceAllocation import LPsolver

from DesignYNcombined import KStrategiesYNcomb
from DesignYNcombined import KStrategiesYNBnew

from Method2 import randomSetting

from Method1 import solve

import pickle

if __name__ == "__main__":
    
    instance = 1
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
        maxT = 5

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
        maxT = 100
    
    if (instance == 3):
        W = 3 # number of time windows
        AI = 1 # interval in which passengers are arriving
        K = 10 # number of passenger types
        R = 3 # number of resources
        mR = 2 # max number of reosurces
        M = 3 # number of attack methods
        P = 10 # number of staff
        shift = 1 # d
        Q = 4 # number of strategies
        maxT = 10    
        
    teams = util.generateAllTeams(R, mR)
    
    Z = 10 # number of runs
    ZQ = 1 # max number of Q
    ZC = 10
    maxQT = 6
    
    obj_relax = np.zeros((Z,ZC))
    obj_uniform = np.zeros((Z,ZC))
    obj_yn = np.zeros((Z,ZC))
    obj_final = np.zeros((Z,ZC))
    time_uniform = np.zeros((Z,ZC))
    time_final = np.zeros((Z,ZC))
    
    K = 10
    P = 10
    
    for z in range(Z):
        seed = random.randint(1,1000)
        for zc in range(ZC):
            W = zc+6
            # Construct random instances
            resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
        
            minr = np.zeros((W,R))   
            
            
            start_time_relax = time.time()
            # Solve full LP relaxation
            obj_relax[z][zc], n_value0, overflow_value, y_value, s_value0, p,z_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=0, OverConstr = 0)
            
            #time_relax[z][zq] = time.time() - start_time_relax
               
            for w in range(W):
                for r in range(R):
                    minr[w][r] = math.floor(y_value[w][r])
                    #print y_value[w][r]
            
            
            # find Q strategies for y, using binaries and lower bound minr, n integer for each ys.    
            start_time_yn = time.time()
            
            
            
            q = np.zeros(Q)
            for i in range(Q):
                q[i] = float(1)/Q 
                    
            obj_yn[z][zc], n, ns,ys,z_value,p,s,y = KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q, ar, phi, integer=0, OverConstr=False, OverConstr2=False)
                 
                   
            
            #time_yn[z][zq][j] = time.time() - start_time_yn
            minn = np.zeros((Q,W,T,K))
            
            for i in range(Q):
                for w in range(W):
                    for k in range(K):
                        sum = 0 
                        for t in range(T):
                            minn[i][w][t][k] = math.floor(ns[i][w][t][k])
                            sum += math.floor(ns[i][w][t][k])
            
            # Find integer solution (using binaries and rounde) for each ns, given ys, p and s
            
            #start_time_final = time.time()
            
            obj_final[z][zc], rt, t3,ni,oi_value,q_tem,o_temp  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)
            
            time_final[z][zc] = time.time() - start_time_yn
            
            start_time_uniform = time.time()
            
            obj_uniform[z][zc], rt = solve(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi)
                    
            time_uniform[z][zc] = time.time() - start_time_uniform
            file = open('DifferentQ_10.txt','w')
            file.write('run%s\n' %str(z+1))
            file.write("obj:\n" + str(obj_relax)+'\n\n'+str(obj_final)+'\n\n'+str(obj_uniform))
            file.write("\n\ntime:\n" +str(time_final) +'\n\n'+str(time_uniform))
            file.close()
            
    
    obj_relax = np.zeros((Z,ZC))
    obj_uniform = np.zeros((Z,ZC))
    obj_yn = np.zeros((Z,ZC))
    obj_final = np.zeros((Z,ZC))
    time_uniform = np.zeros((Z,ZC))
    time_final = np.zeros((Z,ZC))
    
    K=20
    P=20
    for z in range(Z):
        seed = random.randint(1,1000)
        for zc in range(ZC):
            W = zc+6
            # Construct random instances
            resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
        
            minr = np.zeros((W,R))   
            
            
            start_time_relax = time.time()
            # Solve full LP relaxation
            obj_relax[z][zc], n_value0, overflow_value, y_value, s_value0, p,z_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=0, OverConstr = 0)
            
            #time_relax[z][zq] = time.time() - start_time_relax
               
            for w in range(W):
                for r in range(R):
                    minr[w][r] = math.floor(y_value[w][r])
                    #print y_value[w][r]
            
            
            # find Q strategies for y, using binaries and lower bound minr, n integer for each ys.    
            start_time_yn = time.time()
            
            
            
            q = np.zeros(Q)
            for i in range(Q):
                q[i] = float(1)/Q 
                    
            obj_yn[z][zc], n, ns,ys,z_value,p,s,y = KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q, ar, phi, integer=0, OverConstr=False, OverConstr2=False)
                 
                   
            
            #time_yn[z][zq][j] = time.time() - start_time_yn
            minn = np.zeros((Q,W,T,K))
            
            for i in range(Q):
                for w in range(W):
                    for k in range(K):
                        sum = 0 
                        for t in range(T):
                            minn[i][w][t][k] = math.floor(ns[i][w][t][k])
                            sum += math.floor(ns[i][w][t][k])
            
            # Find integer solution (using binaries and rounde) for each ns, given ys, p and s
            
            #start_time_final = time.time()
            
            obj_final[z][zc], rt, t3,ni,oi_value,q_tem,o_temp  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)
            
            time_final[z][zc] = time.time() - start_time_yn
            
            start_time_uniform = time.time()
            
            obj_uniform[z][zc], rt = solve(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi)
                    
            time_uniform[z][zc] = time.time() - start_time_uniform
            file = open('DifferentW_20_0131.txt','w')
            file.write('run%s\n' %str(z+1))
            file.write("obj:\n" + str(obj_relax)+'\n\n'+str(obj_final)+'\n\n'+str(obj_uniform))
            file.write("\n\ntime:\n" +str(time_final) +'\n\n'+str(time_uniform))
            file.close()
            
    obj_relax = np.zeros((Z,ZC))
    obj_uniform = np.zeros((Z,ZC))
    obj_yn = np.zeros((Z,ZC))
    obj_final = np.zeros((Z,ZC))
    time_uniform = np.zeros((Z,ZC))
    time_final = np.zeros((Z,ZC))
    
    K=30
    P=30
    for z in range(Z):
        seed = random.randint(1,1000)
        for zc in range(ZC):
            W = zc+6
            # Construct random instances
            resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
        
            minr = np.zeros((W,R))   
            
            
            start_time_relax = time.time()
            # Solve full LP relaxation
            obj_relax[z][zc], n_value0, overflow_value, y_value, s_value0, p,z_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=0, OverConstr = 0)
            
            #time_relax[z][zq] = time.time() - start_time_relax
               
            for w in range(W):
                for r in range(R):
                    minr[w][r] = math.floor(y_value[w][r])
                    #print y_value[w][r]
            
            
            # find Q strategies for y, using binaries and lower bound minr, n integer for each ys.    
            start_time_yn = time.time()
            
            
            
            q = np.zeros(Q)
            for i in range(Q):
                q[i] = float(1)/Q 
                    
            obj_yn[z][zc], n, ns,ys,z_value,p,s,y = KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q, ar, phi, integer=0, OverConstr=False, OverConstr2=False)
                 
                   
            
            #time_yn[z][zq][j] = time.time() - start_time_yn
            minn = np.zeros((Q,W,T,K))
            
            for i in range(Q):
                for w in range(W):
                    for k in range(K):
                        sum = 0 
                        for t in range(T):
                            minn[i][w][t][k] = math.floor(ns[i][w][t][k])
                            sum += math.floor(ns[i][w][t][k])
            
            # Find integer solution (using binaries and rounde) for each ns, given ys, p and s
            
            #start_time_final = time.time()
            
            obj_final[z][zc], rt, t3,ni,oi_value,q_tem,o_temp  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)
            
            time_final[z][zc] = time.time() - start_time_yn
            
            start_time_uniform = time.time()
            
            obj_uniform[z][zc], rt = solve(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi)
                    
            time_uniform[z][zc] = time.time() - start_time_uniform
            file = open('DifferentW_30.txt','w')
            file.write('run%s\n' %str(z+1))
            file.write("obj:\n" + str(obj_relax)+'\n\n'+str(obj_final)+'\n\n'+str(obj_uniform))
            file.write("\n\ntime:\n" +str(time_final) +'\n\n'+str(time_uniform))
            file.close()
    
    
    #print obj_relax, obj_yn, obj_final
    #print time_relax, time_yn, time_final
    
    
    #pickle.dump((obj_relax, obj_yn, obj_final, time_relax, time_yn, time_final), open("MediumresultsDifferentDistr_q_0128.data", "wb"))
    #obj_relax, obj_yn, obj_final, time_relax, time_yn, time_final = pickle.load(open("MediumresultsDifferentDistr_q_0128.data", "rb))
"""
    instance = 1
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
    
    Z = 0 # number of runs
    ZQ = 3 # max number of Q
    
    maxQT = 6
    obj_relax = np.zeros((Z,ZQ))
    obj_yn = np.zeros((Z,ZQ,maxQT))
    obj_final = np.zeros((Z,ZQ,maxQT))
    time_relax = np.zeros((Z,ZQ))
    time_yn = np.zeros((Z,ZQ,maxQT))
    time_final = np.zeros((Z,ZQ,maxQT))
    
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
                    
            # Solve second step for different q-vectors (if Q=2 or Q=3)
            if (Q==1):
                q = 1*[Q]
                QT = 1
            if (Q==2):
                q = [[0.1667,0.8333],[0.2, 0.8],[0.4,0.6],[0.25,0.75],[0.3333,0.6667],[0.5,0.5]] # if Q=2
                QT = 6
            if (Q==3):
                q = [[0.1667,0.1667,0.6666],[0.1667,0.3333,0.5],[0.2,0.2,0.6],[0.2,0.4,0.4],[0.25,0.25,0.5],[0.3333,0.3333,0.3334]] # if Q=3
                QT = 6 
            if (Q==4):
                q = [[0.1667,0.1667,0.1666,0.5],[0.1667,0.1667,0.3333,0.3333],[0.2,0.2,0.2,0.4],[0.25,0.25,0.25,0.25]]
                QT = 4
            
            # find Q strategies for y, using binaries and lower bound minr, n integer for each ys.    
            for j in range(QT):
                start_time_yn = time.time()
                
                if (Q == 1):
                    obj_yn[z][zq][j], n, ns,ys,z_value,p,s,y = KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q, ar, phi, integer=0, OverConstr=False, OverConstr2=False)
                else:
                    obj_yn[z][zq][j], n, ns,ys,z_value,p,s,y = KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q[j][:], ar, phi, integer=0, OverConstr=False, OverConstr2=False)    
                
                time_yn[z][zq][j] = time.time() - start_time_yn
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
                
                obj_final[z][zq][j], rt, t3,ni,oi_value,q_tem,o_temp  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)
                
                time_final[z][zq][j] = time.time() - start_time_final
                
        file = open('LargeresultsDifferentDistr_q_0128.txt','w')
        file.write('run%s\n' %str(z+1))
        file.write("obj:\n" + str(obj_relax)+'\n\n'+str(obj_yn)+'\n\n'+str(obj_final))
        file.write("\n\ntime:\n" + str(time_relax)+'\n\n'+str(time_yn)+'\n\n'+str(time_final))
        file.close()
    
    
    #print obj_relax, obj_yn, obj_final
    #print time_relax, time_yn, time_final
    
    
    pickle.dump((obj_relax, obj_yn, obj_final, time_relax, time_yn, time_final), open("LargeresultsDifferentDistr_q_0128.data", "wb"))
    #obj_relax, obj_yn, obj_final, time_relax, time_yn, time_final = pickle.load(open("LargeresultsDifferentDistr_q_0128.data", "rb))
    """