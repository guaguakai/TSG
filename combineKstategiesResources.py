import util
from gurobipy import *
import numpy as np
import random
import time
import math
from StaffResourceAllocation import LPsolver
from KStrategies import Ksolver
from KStrategies import randomSetting

if __name__ == "__main__":
    # ============================= main =======================================
    print "======================== main ======================================"
    # ========================= Game Setting ===================================
    W = 15 # number of time windows
    K = 10 # number of passenger types
    R = 5 # number of resources
    mR = 10 # max number of reosurces
    M = 2 # number of attack methods
    P = 10 # number of staff
    Q = 5 # number of pure strategies
    shift = 3 # d

    nT = 20
    teams = util.generateAllTeams(R, mR)
    
    seed = 2345
    start_time = time.time()
    resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
    #resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = testSetting()
    #resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = bigSetting( W, K , mR,  P, teams, shift)
    
    
    minr=np.zeros((W,R))
    
    
    obj_relaxed_all15, n_value0, overflow_value0, y_value0, s_value0, p_value0 = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=0, OverConstr=True)
    
    for w in range(W):
        for r in range(R):
            minr[w][r] = math.floor(y_value0[w][r])
            
    obj_relaxed_n15, n_value, overflow_value, y_value, s_value, p_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=1, OverConstr=True)

    
    print "============================ MIP SOLVER =============================="

    obj15, rt, q15, n2, o, att_set = Ksolver(W, K, R, mR, M, P, Q, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, n_value, overflow_value, y_value, s_value, p_value, integer=0, OverConstr=False)
    walltime = time.time() - start_time
    
    
    print "relaxation all objective: ", obj_relaxed_all15
    print "relaxation n objective : ", obj_relaxed_n15
    print "MIP objective: ", obj15
    print "Runtime/walltime ", rt, " ", walltime