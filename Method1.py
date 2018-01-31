import util
from gurobipy import *
import numpy as np
import random
from KStrategies import randomSetting
from relaxed_feed import LPsolver
import math
import time

def FixedQ(Q, W, K, R, M, P, y_start, n_start,q, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, OverConstr=False, OverConstr2=False, verbose=True): # integer indicates different relaxation method
    # ======================= Gurobi Setting ===================================
    model = Model("MIP")
    if not verbose:
        model.params.OutputFlag=0
        model.params.TuneOutput=0
    #model.params.DualReductions = 0
    #model.params.TimeLimit=600 # at most 10 min
    model.params.MIPGap=0.01;

    team = [[ model.addVar(lb=0.0, ub = 1.0, vtype=GRB.BINARY, name="team_t{0}_s{1}".format(t,i)) for t in range(T)] for i in range(Q)]

    theta = model.addVar(vtype=GRB.CONTINUOUS, lb=-10000, name="theta")
    z = [] # z[w][k][m]
    for w in range(W):
        z.append([])
        for k in range(K):
            z[w].append([])
            for m in range(M):
                tmp_z_var = model.addVar(vtype=GRB.CONTINUOUS, name="z_w{0}_k{1}_m{2}".format(w, k, m))
                z[w][k].append(tmp_z_var)

    pi = [] # pi[w][t][k]
    for w in range(W):
        pi.append([])
        for t in range(T):
            pi[w].append([])
            for k in range(K):
                tmp_pi_var = model.addVar(vtype=GRB.CONTINUOUS, name="pi_w{0}_t{1}_k{2}".format(w, t, k))
                pi[w][t].append(tmp_pi_var)

    n_wtk = [] # n_wtk[w][t][k] # integer value of N_wk[w][k] * pi[w][t][k]
    for w in range(W):
        n_wtk.append([])
        for t in range(T):
            n_wtk[w].append([])
            for k in range(K):
                tmp_pi_var = model.addVar(vtype=GRB.CONTINUOUS, name="n_w{0}_t{1}_k{2}".format(w, t, k))
                n_wtk[w][t].append(tmp_pi_var)
                
    ni = [[[[model.addVar(vtype=GRB.INTEGER, name="ni_s{0}_w{1}_t{2}_k{3}".format(i, w, t, k)) for k in range(K)] for t in range(T)]for w in range(W)]for i in range(Q)]

    overflow = [] # overflow[w][r]
    for w in range(W):
        overflow.append([])
        for r in range(R):
            tmp_overflow_var = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="o_w{0}_r{1}".format(w, r)) # this is original integer but can be approximated by continuous value
            overflow[w].append(tmp_overflow_var)
    
    O = [[[model.addVar(vtype=GRB.CONTINUOUS, name="O_{0}_w{1}_r{2}".format(i, w, r)) for r in range(R)] for w in range(W)] for i in range(Q)]
    
    #y = [[model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="y_w{0}_r{1}".format(w, r)) for r in range(R)] for w in range(W)]
    yi = [[[model.addVar(vtype=GRB.INTEGER, lb=0, name="y_w{0}_r{1}_s{2}".format(w, r,i)) for r in range(R)] for w in range(W)] for i in range(Q)]
     # y[i][w][r]: number of operating resources r at time w
    

    p = [] # available staff
    s = [] # working staff
    for w in range(W):
        tmp_staff = model.addVar(vtype=GRB.CONTINUOUS, name="p_w{0}".format(w))
        tmp_working_staff = model.addVar(vtype=GRB.INTEGER, name="s_w{0}".format(w))
        p.append(tmp_staff)
        s.append(tmp_working_staff)
        
    model.update()
    
    for i in range(Q):
        for w in range(W):
            for r in range(R):
                yi[i][w][r].start = y_start[w][r]
            for k in range(K):
                for t in range(T):
                    ni[i][w][t][k].start = n_start[w][t][k]
    
    model.update()

    # ========================= Gurobi Objective ===============================
    objective_variables = [theta] + [overflow[w][r] for w in range(W) for r in range(R)]
    objective_coefficients = [1] + [-phi[r] for r in range(R)]*W
    objective_value = LinExpr(objective_coefficients, objective_variables)

    #objective_value = theta
    #for w in range(W):
    #    for r in range(R):
    #        objective_value += phi[r]*overflow[w][r]
    model.setObjective(objective_value, GRB.MAXIMIZE)

    # ======================= Gurobi Constraints ===============================
    for w in range(W):
        for k in range(K):
            if N_wk[w][k] > 0 :          
                for m in range(M):
                    model.addConstr(theta - z[w][k][m]*(U_plus[k] - U_minus[k]) - U_minus[k] <= 0, "(1)_w{0}_k{1}_m{2}".format(w,k,m))

    for w in range(W):
        for k in range(K):
            for m in range(M):
                tmp_sum = LinExpr([E[t][m] for t in range(T)], [pi[w][t][k] for t in range(T)])
                model.addConstr(z[w][k][m] - tmp_sum == 0, name="(2)_w{0}_k{1}_m{2}".format(w, k, m))

    for w in range(W):
        for k in range(K):
            tmp_sum = LinExpr([1]*T, [pi[w][t][k] for t in range(T)])
            model.addConstr(tmp_sum == 1, name="(3)_w{0}_k{1}".format(w, k))

    for w in range(W):
        for t in range(T):
            for k in range(K):
                #if N_wk[w][k] > 0 :
                model.addConstr(pi[w][t][k] * N_wk[w][k] - n_wtk[w][t][k] == 0, name="(3.5)_w{0}_t{1}_k{2}".format(w,t,k))
                #else:
                #    model.addConstr(pi[w][t][k] == 0, name="(3.5)_w{0}_t{1}_k{2}".format(w,t,k))
                #    model.addConstr(n_wtk[w][t][k] == 0, name="(3.6)_w{0}_t{1}_k{2}".format(w,t,k))

    #pre_overflow = np.random.randint(0, 100, R)
    pre_overflow = [0] * R
   
    for i in range(Q):
        for w in range(W):
            tmp_sum = LinExpr(ar, [yi[i][w][r] for r in range(R)])
            model.addConstr(tmp_sum - p[w] <= 0, name="(6)_w{0}_s{1}".format(w,i))

        for w in range(W):
            for r in range(R):
                model.addConstr(yi[i][w][r] - mr[r] <= 0, name="(7)_w{0}_r{1}_s{2}".format(w, r,i)) 

    for w in range(W):
        start_index = max(0, w - shift + 1)
        tmp_sum = LinExpr([1]*(w - start_index + 1), [s[j] for j in range(start_index, w+1)])
        model.addConstr(tmp_sum - p[w] == 0, name="(8)_w{0}".format(w))

    tmp_sum = LinExpr([1]*W, [s[w] for w in range(W)])
    model.addConstr(tmp_sum - P <= 0, name="(9)")

   
    #pure strategy constraints
    pre_overflow = [0] * R
    for i in range(Q):
        for w in range(W):
            for k in range(K):
                model.addConstr( quicksum(ni[i][w][t][k] for t in range(T)) == N_wk[w][k])  
                
    for i in range(Q):
        for w in range(W):
            for r in range(R):  
                tmp_sum = LinExpr([1 for k in range(K)]*len(resource2team[r]), [ni[i][w][t][k] for t in resource2team[r] for k in range(K)])
                if w == 0:
                    model.addConstr(tmp_sum + pre_overflow[r] - yi[i][w][r] * C[r] - O[i][w][r] <= 0, name="(10)_w{0}_r{1}_s{2}".format(w, r, i))
                else:
                    model.addConstr(tmp_sum + O[i][w-1][r] - yi[i][w][r] * C[r] - O[i][w][r] <= 0, name="(10)_w{0}_r{1}_s{2}".format(w, r, i ))
    
                model.addConstr( tmp_sum <= yi[i][w][r]*10000, name="(10.5)_w{0}_r{1}_s{2}".format(w, r,i))
    
    if OverConstr2:        
        for r in range(R): # OPTIONAL
            model.addConstr(O[i][W-1][r] == 0, name="(.95)_r{0}".format(r))
    
        for w in range(W):
            for r in range(R):
                if w > 0:
                    model.addConstr(yi[i][w][r] * C[r] - O[i][w-1][r] >= 0, name="(5.9)_w{0}_r{1}".format(w, r))

   
    for w in range(W):
        for t in range(T):
            for k in range(K):
                tmp_sum = LinExpr([q[i] for i in range(Q)], [ni[i][w][t][k] for i in range(Q)])
                model.addConstr(tmp_sum == n_wtk[w][t][k])
        for r in range(R):
            tmp_sum = LinExpr([q[i] for i in range(Q)], [O[i][w][r] for i in range(Q)])
            model.addConstr(tmp_sum == overflow[w][r])

    

    #TEAM CONSTRAINTS 
    if True:
        for i in range(Q):
            for w in range(W):
                for k in range(K):
                    for t in range(T):
                        model.addConstr( ni[i][w][t][k] <= team[i][t]*N_wk[w][k], name="team{0}{1}{2}{3}".format(t,i,k,w))         
            model.addConstr(quicksum(team[i][t] for t in range(T)) <= maxT)
   
   
   
   
    model.update()

    model.write("tsgkynteam.lp")

    model.optimize()
    model.write("tsgkynteam.sol")

    

    n_value = np.zeros((W,T,K))
    for w in range(W):
        for t in range(T):
            for k in range(K):
                n_value[w][t][k] = n_wtk[w][t][k].x
                #if n_value[w][t][k] != int(n_value[w][t][k]):
                    #print n_value[w][t][k]
    ns_value = np.zeros((Q,W,T,K))
    team_val = []
    for i in range(Q):
        team_val_i = []
        for w in range(W):
            for t in range(T):
                for k in range(K):
                    ns_value[i][w][t][k] = ni[i][w][t][k].x
                    if ns_value[i][w][t][k] > 0 and t not in team_val_i:
                        team_val_i.append(t)
        team_val.append(team_val_i)
    
    overflow_value = np.zeros((W,R))
    for w in range(W):
        for r in range(R):
            overflow_value[w][r] = overflow[w][r].x
    
    #team_val = np.zeros((Q,T))
    #for i in range(Q):
    #    for t in range(T):
    #        team_val[i][t] = team[i][t].x
    
    obj = model.getAttr('ObjVal')

    model.terminate()

    return obj, n_value, ns_value, team_val

def solve(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, verbose=True):
   
    start_time = time.time()
    
    minr = np.zeros((W,R))
    
    obj_relax, n_value0, overflow_value0, y_value0, s_value0, p, attset, f, tl = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=1, OverConstr=False, TeamConstr=False, MaxT=maxT, Q=Q, verbose=verbose)
    
    minr = np.zeros((W,R))
    minn = np.zeros((W,T,K))
    for w in range(W):
        for r in range(R):
            minr[w][r] = math.floor(y_value0[w][r])
        for k in range(K):
            for t in range(T):
                minn[w][t][k] = math.floor(n_value0[w][t][k])
    q_val = np.zeros(Q)
    for i in range(Q):
        q_val[i] = float(1)/Q
        
    obj, n_value, ns_value, team_val = FixedQ(Q, W, K, R, M, P, minr, minn, q_val, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, False, False, verbose=verbose)
    rt = time.time() - start_time
        
    return obj, rt
          
if __name__ == "__main__":
    # ============================= main =======================================
    print "======================== main ======================================"
    # ========================= Game Setting ===================================
    W = 2 # number of time windows
    K = 2 # number of passenger types
    R = 2 # number of resources
    mR = 1 # max number of reosurces
    M = 1 # number of attack methods
    P = 10 # number of staff
    Q= 2
    shift = 1 # d
    maxT=2

    nT = 22
    teams = util.generateAllTeams(R, mR)
    #teams = util.randomGenerateTeams(R, mR, nT)
    Nmax = 150

    #print teams

    # ================= random generate game setting ===========================
    seed = 2
    resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
    
    obj, rt = solve(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi)
    
    print obj, rt

