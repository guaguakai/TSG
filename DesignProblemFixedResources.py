import util
from gurobipy import *
import time
import numpy as np
import random
import math
from relaxed_feed import LPsolver, LPsolverR
from KStrategies import randomSetting as rs
from KStrategiesFixedYRoundN import Ksolver

def KStrategiesYNB(Q, W, K, R, M, resource2team, T, E, C, U_plus, U_minus, N_wk, y, yi, n, p, s, phi, integer=0, OverConstr=False, OverConstr2=False): # integer indicates different relaxation method
    # ======================= Gurobi Setting ===================================
    model = Model("MIP")
    model.params.DualReductions = 0
    model.params.MIPGap=0.0005;

    q = [ model.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="q_s{0}".format(i)) for i in range(Q)]

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
                if (integer == 2) or (integer == 3):
                    tmp_pi_var = model.addVar(vtype=GRB.INTEGER, name="n_w{0}_t{1}_k{2}".format(w, t, k))
                else:
                    tmp_pi_var = model.addVar(vtype=GRB.CONTINUOUS, name="n_w{0}_t{1}_k{2}".format(w, t, k))
                n_wtk[w][t].append(tmp_pi_var)
                
    nb = [[[[model.addVar(vtype=GRB.BINARY, name="ni_s{0}_w{1}_t{2}_k{3}".format(i, w, t, k)) for k in range(K)] for t in range(T)]for w in range(W)]for i in range(Q)]
    X = [[[[model.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="X(s%d,w%d,k%d,t%d)" %(i, w,k,t))  for k in range(K)] for t in range(T)] for w in range(W)] for i in range(Q)]

    # n_wtk[w][t][k] # integer value of N_wk[w][k] * pi[w][t][k]

    overflow = [] # overflow[w][r]
    for w in range(W):
        overflow.append([])
        for r in range(R):
            tmp_overflow_var = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="o_w{0}_r{1}".format(w, r)) # this is original integer but can be approximated by continuous value
            overflow[w].append(tmp_overflow_var)
    
    O = [[[model.addVar(vtype=GRB.CONTINUOUS, name="O_{0}_w{1}_r{2}".format(i, w, r)) for r in range(R)] for w in range(W)] for i in range(Q)]
    y = y
    yi = yi 
     # y[i][w][r]: number of operating resources r at time w
    

    p = p # available staff
    s = s # working staff
    
        
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
                model.addConstr(pi[w][t][k] * N_wk[w][k] - n_wtk[w][t][k] == 0, name="(3.5)_w{0}_t{1}_k{2}".format(w,t,k))

    #pre_overflow = np.random.randint(0, 100, R)
    pre_overflow = [0] * R
    for w in range(W):
        for r in range(R):

            tmp_sum = LinExpr([1 for k in range(K)]*len(resource2team[r]), [n_wtk[w][t][k] for t in resource2team[r] for k in range(K)])
            if w == 0:
                model.addConstr(tmp_sum + pre_overflow[r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(4)_w{0}_r{1}".format(w, r))
            else:
                model.addConstr(tmp_sum + overflow[w-1][r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(4)_w{0}_r{1}".format(w, r))

            model.addConstr( tmp_sum <= y[w][r]*10000 ,  name="(5.6)_w{0}_r{1}".format(w, r))
    
    if OverConstr:        
        for r in range(R): # OPTIONAL
            model.addConstr(overflow[W-1][r] == 0, name="(5)_r{0}".format(r))
    
        for w in range(W):
            for r in range(R):
                if w > 0:
                    model.addConstr(y[w][r] * C[r] - overflow[w-1][r] >= 0, name="(5.5)_w{0}_r{1}".format(w, r))

            

   
    #pure strategy constraints
    pre_overflow = [0] * R
    for i in range(Q):
        for w in range(W):
            for r in range(R):  
                marginal_sum = quicksum(n[i][w][t][k] for t in resource2team[r] for k in range(K))
                tmp_sum = LinExpr([1 for k in range(K)]*len(resource2team[r]), [nb[i][w][t][k] for t in resource2team[r] for k in range(K)])
                if w == 0:
                    model.addConstr(marginal_sum + tmp_sum + pre_overflow[r] - yi[i][w][r] * C[r] - O[i][w][r] <= 0, name="(10)_w{0}_r{1}_{2}".format(w, r, i))
                else:
                    model.addConstr(marginal_sum + tmp_sum + overflow[w-1][r] - yi[i][w][r] * C[r] - O[i][w][r] <= 0, name="(10)_w{0}_r{1}_{2}".format(w, r, i))
    
                model.addConstr( marginal_sum+ tmp_sum <= y[w][r]*10000, name="(10.5)_w{0}_r{1}".format(w, r))
        
    for w in range(W):
        for t in range(T):
            for k in range(K):
                tmp_sum = LinExpr([1 for i in range(Q)], [X[i][w][t][k] for i in range(Q)])
                marginal_sum = LinExpr([n[i][w][t][k] for i in range(Q)], [q[i] for i in range(Q)])
                model.addConstr(marginal_sum + tmp_sum == n_wtk[w][t][k])
    
    tmp_sum = LinExpr([1]*Q, [q[i] for i in range(Q)])
    model.addConstr(tmp_sum == 1, name="sumQ")   

    # Linearization Constraints
    for i in range(Q):
        for w in range(W):
            for k in range(K):
                model.addConstr( quicksum(n[i][w][t][k] + nb[i][w][t][k] for t in range(T)) == N_wk[w][k])  
                for t in range(T):
                    model.addConstr(X[i][w][t][k] <= q[i]) 
                    model.addConstr(X[i][w][t][k] <= nb[i][w][t][k] ) 
                    model.addConstr(X[i][w][t][k] >= q[i] -(1-nb[i][w][t][k])) 
    
    if OverConstr2:  
        for i in range(Q):      
            for r in range(R): # OPTIONAL
                model.addConstr(O[i][W-1][r] == 0, name="(5)_r{0}".format(r))
        
            for w in range(W):
                for r in range(R):
                    if w > 0:
                        model.addConstr(yi[i][w][r] * C[r] - O[i][w-1][r] >= 0, name="(5.5)_w{0}_r{1}".format(w, r))

            


   
    model.update()

    model.write("tsgkpMIP.lp")
    start_time = time.time()
    model.optimize()
    runtime = time.time() - start_time

    model.write("tsgkpMIP.sol")

    

    n_value = np.zeros((W,T,K))
    for w in range(W):
        for t in range(T):
            for k in range(K):
                n_value[w][t][k] = n_wtk[w][t][k].x
                #if n_value[w][t][k] != int(n_value[w][t][k]):
                    #print n_value[w][t][k]
    
    
    
    obj = model.getAttr('ObjVal')


    return obj, runtime
    


def KStrategiesYN(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, y, yi, q, p, s, ar, phi, integer=0, OverConstr=False, OverConstr2=False): # integer indicates different relaxation method
    # ======================= Gurobi Setting ===================================
    model = Model("MIP")
    model.params.DualReductions = 0
    model.params.MIPGap=0.0005;

    team = [[ model.addVar(lb=0.0, ub = 1.0, vtype=GRB.INTEGER, name="team{0}".format(t)) for t in range(T)] for i in range(Q)]

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
                if (integer == 2) or (integer == 3):
                    tmp_pi_var = model.addVar(vtype=GRB.INTEGER, name="n_w{0}_t{1}_k{2}".format(w, t, k))
                else:
                    tmp_pi_var = model.addVar(vtype=GRB.CONTINUOUS, name="n_w{0}_t{1}_k{2}".format(w, t, k))
                n_wtk[w][t].append(tmp_pi_var)
                
    ni = [[[[model.addVar(vtype=GRB.CONTINUOUS, name="ni_s{0}_w{1}_t{2}_k{3}".format(i, w, t, k)) for k in range(K)] for t in range(T)]for w in range(W)]for i in range(Q)]
    # n_wtk[w][t][k] # integer value of N_wk[w][k] * pi[w][t][k]

    overflow = [] # overflow[w][r]
    for w in range(W):
        overflow.append([])
        for r in range(R):
            tmp_overflow_var = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="o_w{0}_r{1}".format(w, r)) # this is original integer but can be approximated by continuous value
            overflow[w].append(tmp_overflow_var)
    
    O = [[[model.addVar(vtype=GRB.CONTINUOUS, name="O_{0}_w{1}_r{2}".format(i, w, r)) for r in range(R)] for w in range(W)] for i in range(Q)]
    y = y
    yi = yi 
     # y[i][w][r]: number of operating resources r at time w
    

    p = p # available staff
    s = s # working staff
    
        
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
                model.addConstr(pi[w][t][k] * N_wk[w][k] - n_wtk[w][t][k] == 0, name="(3.5)_w{0}_t{1}_k{2}".format(w,t,k))

    #pre_overflow = np.random.randint(0, 100, R)
    pre_overflow = [0] * R
    for w in range(W):
        for r in range(R):

            tmp_sum = LinExpr([1 for k in range(K)]*len(resource2team[r]), [n_wtk[w][t][k] for t in resource2team[r] for k in range(K)])
            if w == 0:
                model.addConstr(tmp_sum + pre_overflow[r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(4)_w{0}_r{1}".format(w, r))
            else:
                model.addConstr(tmp_sum + overflow[w-1][r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(4)_w{0}_r{1}".format(w, r))

            model.addConstr( tmp_sum <= y[w][r]*10000 ,  name="(5.6)_w{0}_r{1}".format(w, r))
    
    if OverConstr:        
        for r in range(R): # OPTIONAL
            model.addConstr(overflow[W-1][r] == 0, name="(5)_r{0}".format(r))
    
        for w in range(W):
            for r in range(R):
                if w > 0:
                    model.addConstr(y[w][r] * C[r] - overflow[w-1][r] >= 0, name="(5.5)_w{0}_r{1}".format(w, r))

            

   
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
                    model.addConstr(tmp_sum + pre_overflow[r] - yi[i][w][r] * C[r] - O[i][w][r] <= 0, name="(10)_w{0}_r{1}".format(w, r))
                else:
                    model.addConstr(tmp_sum + O[i][w-1][r] - yi[i][w][r] * C[r] - O[i][w][r] <= 0, name="(10)_w{0}_r{1}".format(w, r))
    
                model.addConstr( tmp_sum <= y[w][r]*10000, name="(10.5)_w{0}_r{1}".format(w, r))
    
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
    for i in range(Q):
        for t in range(T):
            tmp = LinExpr([1/(np.maximum(N_wk[w][k],1)) for w in range(W) for k in range(K)], [ni[i][w][t][k] for w in range(W) for k in range(K)])
            model.addConstr(tmp  <= team[i][t], name="team{0}".format(t))         
        model.addConstr(quicksum(team[i][t] for t in range(T)) <= maxT)
   
    model.update()

    model.write("tsgkyn.lp")

    model.optimize()
    model.write("tsgkyn.sol")

    

    n_value = np.zeros((W,T,K))
    for w in range(W):
        for t in range(T):
            for k in range(K):
                n_value[w][t][k] = n_wtk[w][t][k].x
                #if n_value[w][t][k] != int(n_value[w][t][k]):
                    #print n_value[w][t][k]
    ns_value = np.zeros((Q,W,T,K))
    for i in range(Q):
        for w in range(W):
            for t in range(T):
                for k in range(K):
                    ns_value[i][w][t][k] = ni[i][w][t][k].x
    
    
    overflow_value = np.zeros((W,R))
    for w in range(W):
        for r in range(R):
            overflow_value[w][r] = overflow[w][r].x

    
    obj = model.getAttr('ObjVal')


    return obj, n_value, ns_value
    


def KStrategiesY(Q, W, K, R, mR, M, P, teams, resource2team, T, MaxT, E, C, U_plus, U_minus, N_wk, shift, mr, miny, ar, phi, integer=0, OverConstr=False): # integer indicates different relaxation method
    # ======================= Gurobi Setting ===================================
    model = Model("MIP")
    model.params.DualReductions = 0
    model.params.MIPGap=0.0005;

    team = [ model.addVar(lb=0.0, ub = 1.0, vtype=GRB.INTEGER, name="team{0}".format(t)) for t in range(T)]
    q = [ model.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="q_s{0}".format(i)) for i in range(Q)]

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
                if (integer == 2) or (integer == 3):
                    tmp_pi_var = model.addVar(vtype=GRB.INTEGER, name="n_w{0}_t{1}_k{2}".format(w, t, k))
                else:
                    tmp_pi_var = model.addVar(vtype=GRB.CONTINUOUS, name="n_w{0}_t{1}_k{2}".format(w, t, k))
                n_wtk[w][t].append(tmp_pi_var)

    overflow = [] # overflow[w][r]
    for w in range(W):
        overflow.append([])
        for r in range(R):
            tmp_overflow_var = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="o_w{0}_r{1}".format(w, r)) # this is original integer but can be approximated by continuous value
            overflow[w].append(tmp_overflow_var)

    var = GRB.CONTINUOUS
    y = [[model.addVar(vtype=var, lb=0, name="y_w{0}_r{1}".format(w, r)) for r in range(R)] for w in range(W)]
     # y[i][w][r]: number of operating resources r at time w
    

    p = [] # available staff
    s = [] # working staff
    for w in range(W):
        tmp_staff = model.addVar(vtype=GRB.CONTINUOUS, name="p_w{0}".format(w))
        tmp_working_staff = model.addVar(vtype=GRB.CONTINUOUS, name="s_w{0}".format(w))
        p.append(tmp_staff)
        s.append(tmp_working_staff)
        
    yb = [[[model.addVar(vtype=GRB.BINARY, lb=0, name="yb_w{0}_r{1}_s{2}".format(w, r,i)) for r in range(R)] for w in range(W)] for i in range(Q)]
    X = [[[model.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="X(s%d,w%d,r%d)" %(i,w,r))  for r in range(R)] for w in range(W)] for i in range(Q)]
 
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
                model.addConstr(pi[w][t][k] * N_wk[w][k] - n_wtk[w][t][k] == 0, name="(3.5)_w{0}_t{1}_k{2}".format(w,t,k))

    #TEAM CONSTRAINTS  
    for t in range(T):
        model.addConstr( quicksum(pi[w][t][k] for w in range(W) for k in range(K))/(W*K) <= team[t], name="team{0}".format(t))         
    model.addConstr(quicksum(team[t] for t in range(T)) <= MaxT)
    
    
    #pre_overflow = np.random.randint(0, 100, R)
    pre_overflow = [0] * R
    for w in range(W):
        for r in range(R):
            #tmp_sum = LinExpr([N_wk[w][k] for k in range(K)]*len(resource2team[r]), [pi[w][t][k] for t in resource2team[r] for k in range(K)])
            tmp_sum = LinExpr([1 for k in range(K)]*len(resource2team[r]), [n_wtk[w][t][k] for t in resource2team[r] for k in range(K)])
            if w == 0:
                model.addConstr(tmp_sum + pre_overflow[r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(4)_w{0}_r{1}".format(w, r))
            else:
                model.addConstr(tmp_sum + overflow[w-1][r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(4)_w{0}_r{1}".format(w, r))

            model.addConstr(y[w][r]*10000 >= tmp_sum, name="(5.6)_w{0}_r{1}".format(w, r))
    if OverConstr:        
        for r in range(R): # OPTIONAL
            model.addConstr(overflow[W-1][r] == 0, name="(5)_r{0}".format(r))
    
        for w in range(W):
            for r in range(R):
                if w > 0:
                    model.addConstr(y[w][r] * C[r] - overflow[w-1][r] >= 0, name="(5.5)_w{0}_r{1}".format(w, r))

    for w in range(W):
        tmp_sum = LinExpr(ar, [y[w][r] for r in range(R)])
        model.addConstr(tmp_sum - p[w] <= 0, name="(6)_w{0}".format(w))

    for w in range(W):
        for r in range(R):
            model.addConstr(y[w][r] - mr[r] <= 0, name="(7)_w{0}_r{1}".format(w, r))
            
    
    #for w in range(W):
    #    for r in range(R):
    #        model.addConstr(-y[w][r] - minr[w][r] <= 0, name="(10)_w{0}_r{1}".format(w, r))
            
    for w in range(W):
        for r in range(R):
            model.addConstr(y[w][r] - quicksum(X[i][w][r] for i in range(Q)) - miny[w][r] == 0, name="(11)_w{0}_r{1}".format(w, r))
    
    tmp_sum = LinExpr([1]*Q, [q[i] for i in range(Q)])
    model.addConstr(tmp_sum == 1, name="sumQ")         
#   # Linearization Constraints
    for i in range(Q):
        for w in range(W):
            for r in range(R):
                model.addConstr(X[i][w][r] <= q[i]) 
                model.addConstr(X[i][w][r] <= yb[i][w][r] ) 
                model.addConstr(X[i][w][r] >= q[i] -(1-yb[i][w][r])) 
    
    
                
    for w in range(W):
        start_index = max(0, w - shift + 1)
        tmp_sum = LinExpr([1]*(w - start_index + 1), [s[i] for i in range(start_index, w+1)])
        model.addConstr(tmp_sum - p[w] == 0, name="(8)_w{0}".format(w))

    tmp_sum = LinExpr([1]*W, [s[w] for w in range(W)])
    model.addConstr(tmp_sum - P <= 0, name="(9)")
    
    tmp_sum = LinExpr([1]*Q, [q[i] for i in range(Q)])
    model.addConstr(tmp_sum == 1, name="sumQ")
   
    model.update()

    model.write("tsgkp.lp")

    model.optimize()
    model.write("tsgkp.sol")

    

    n_value = np.zeros((W,T,K))
    for w in range(W):
        for t in range(T):
            for k in range(K):
                n_value[w][t][k] = n_wtk[w][t][k].x
                #if n_value[w][t][k] != int(n_value[w][t][k]):
                    #print n_value[w][t][k]

    overflow_value = np.zeros((W,R))
    for w in range(W):
        for r in range(R):
            overflow_value[w][r] = overflow[w][r].x

    y_value = np.zeros((W,R))
    ys_value = np.zeros((Q,W,R))

    for w in range(W):
        for r in range(R):
            y_value[w][r] = y[w][r].x
            for i in range(Q):
                ys_value[i][w][r] = miny[w][r] + yb[i][w][r].x
                

    s_value = np.zeros(W)
    for w in range(W):
        s_value[w] = s[w].x
    
    p_value = np.zeros(W)
    for w in range(W):
        p_value[w] = p[w].x    
    q_val = np.zeros(Q)
    for i in range(Q):
        q_val[i]=q[i].x
    obj = model.getAttr('ObjVal')


    return obj, n_value, overflow_value, y_value, ys_value, s_value, p_value, q_val
    

def randomSetting(seed, W, K ,R, mR, M, P, teams, shift):
    # ========================== Random Seed ===================================
    #seed = 3234
    #seed = random.randint(1,10000)
    np.random.seed(seed)
    random.seed(seed)

    T = len(teams)
    resource2team = util.resourceTeamDict(R, T, teams)
    #print resource2team

    Er = np.random.rand(R, M)/2 + 0.5 # Er[m][r]
    Er = Er / 2
    #print "Er"
    #print Er
    Er = [[0.3, 0.5, 0.2], 
          [0.6, 0.3, 0.4], 
          [0.4, 0.6, 0.5],
          [0.6, 0.3, 0.8], 
          [0.7, 0.4, 0.7], 
          [0.7, 0.6, 0.9]]

    E = util.computeTeamsRate(R, M, T, teams, Er)
    #print E     

    # RatesEr = [0.357, 0.916, 0.511, 0.916, 1.204, 1.204,
    #     0.693, 0.357, 0.916, 0.357, 0.511, 0.916,
    #     0.223, 0.511, 0.693, 1.609, 1.204, 2.303]

    # suppose this is a zero-sum game
    U_plus = [] # covered (plus) utility of the defender
    U_minus = [] # uncovered (minus) utility of the defender
    for i in range(K):
        tmp_plus_utility = 0 #random.randint(100,500)
        U_plus.append(tmp_plus_utility)
        tmp_minus_utility = -random.randint(500,2000)
        U_minus.append(tmp_minus_utility)
    #print "\nutilities"
    #print "plus"
    #print U_plus
    #print "minus"
    #print U_minus

#    N_wk = [] # N_wk[w][k] represents the number of people getting in time window w with type k
#    for w in range(W):
#        N_wk.append([])
#        for k in range(K):
#            large_or_small = random.random()
#            if large_or_small > 0.5:
#                tmp_N = random.randint(100, 300)
#            else:
#                tmp_N = random.randint(10, 30)
#            N_wk[w].append(tmp_N)
    N_wk = np.zeros((W,K))        
    for k in range(K):
        startK = random.randint(3,W)
        for w in range(startK-3,startK):
            large_or_small = random.random()
            if large_or_small > 0.5:
                tmp_N = random.randint(100, 300)
            else:
                tmp_N = random.randint(10, 30)
            N_wk[w][k] = tmp_N

    #print "number of passengers: {0}".format([sum(N_wk[w]) for w in range(W)])


    C = np.random.randint(200, 500, R) # C[r] is the capacity of resource r
    C = [100,80,75,50,60,60]
    #print "\nC"
    #print C
    #C = [100,80,75,50,30,15]
    mr = np.random.randint(10, 30, R)
    #print "\nmr"
    #print mr
    #mr = [5, 5, 5, 3, 2, 4] # maximum number of people to operate resource r
    mr = [10,10,10,15,10,5]
    
    minr = np.zeros((W,R))

    ar = np.random.randint(1, 5, R)
    #print "\nar"
    #print ar
    #ar = [2, 1, 2, 1, 1, 2] # ar[r] number of people required to operate resource r

    phi = np.random.rand(R)/10 # phi[r] overflow penalty
    print "\nphi"
    print phi
    #phi = np.random.rand(W, R) # phi[w][r] overflow penalty

    return resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi

if __name__ == "__main__":
    # ============================= main =======================================
    print "======================== main ======================================"
    # ========================= Game Setting ===================================
    W = 15 # number of time windows
    K = 1 # number of passenger types
    R = 5 # number of resources
    mR = 10 # max number of reosurces
    M = 2 # number of attack methods
    P = 10 # number of staff
    shift = 5 # d
    Q = 20
    nT = 20
    teams = util.generateAllTeams(R, mR)
    #teams = util.randomGenerateTeams(R, mR, nT)


    # ================= random generate game setting ===========================
    seed = 2345
    #resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
    resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = rs(seed, W, K ,R, mR, M, P, teams, shift)

    minr = np.zeros((W,R))
    start_time = time.time()

    print "============================ FULL LP relaxation =============================="
    obj_relax, n_value0, overflow_value0, y_value0, s_value0, p, attset, f = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=1, OverConstr=True)
    for w in range(W):
        for r in range(R):
            minr[w][r] = math.floor(y_value0[w][r])
            print y_value0[w][r]
            
    print "============================ Binary Y / single relaxed n_wtk (allocated arrivals) MIP ======================"
    objy, n_value, overflow_value, y_value, ys,  s_value, p_value, q = KStrategiesY( Q, W, K, R, mR, M, P, teams, resource2team, T, 5, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, OverConstr=True)
    
    print "============================ multiple relaxed n_wtk (allocated arrivals) MIP ==============================="

    objyn, n_val, n = KStrategiesYN(Q, W, K, R, M, P, resource2team, T,5,  E, C, U_plus, U_minus, N_wk, shift, mr, y_value, ys, q, p_value, s_value, ar, phi, integer=0, OverConstr=True, OverConstr2=False)
    
    minn = np.zeros((Q,W,T,K))
    
    for i in range(Q):
        for w in range(W):
            for k in range(K):
                sum = 0 
                for t in range(T):
                    minn[i][w][t][k] = math.floor(n[i][w][t][k])
                    sum += math.floor(n[i][w][t][k])
                print "n",i,w,k," ", sum, " ", N_wk[w][k],  q[i]
    #n =[]
    #objs =[]
    #for i in range(Q):
    #    obj_relaxi, n_value0, overflow_value0, attset, f = LPsolverR(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, ys[i], s_value, p_value, mr, ar, phi, integer=0, OverConstr=True)
    #    n.append(n_value0)
    #    objs.append(obj_relaxi)
    #obj, rt, q, n2, o, att_set = Ksolver(W, K, R, mR, M, P, Q,1, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, minn, overflow_value, ys, s_value, p_value, integer=0, OverConstr=False)
    
    print "============================ Integer n_wtk  MIP ==============================="

    obj, rt  = KStrategiesYNB(Q, W, K, R, M, resource2team, T, E, C, U_plus, U_minus, N_wk, y_value, ys, minn, p, s_value, phi, integer=0, OverConstr=True, OverConstr2=False)
    walltime = time.time() - start_time

    print "Runtime/walltime ", rt, " ", walltime

    #print y_value0

    #print y_value
    
    print obj_relax, objy, objyn, obj, 
    #print y_value
    #for i in range(Q):
    #    if q[i]>0:
    #        print q[i]
    #        print ys[i]


