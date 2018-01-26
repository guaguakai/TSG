import util
from gurobipy import *
import time
import numpy as np
import random
import math
#from relaxed_feed import LPsolver, LPsolverR
from StaffResourceAllocation import LPsolver
from KStrategiesFixedYRoundN import Ksolver
#from DesignProblemFixedResources import KStrategiesYNB

def KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, n, p, s, phi, integer=0, OverConstr=False, OverConstr2=False): # integer indicates different relaxation method
    # ======================= Gurobi Setting ===================================
    model = Model("MIP")
    model.params.DualReductions = 0
    model.params.MIPGap=0.01;

    team = [[ model.addVar(lb=0.0, ub = 1.0, vtype=GRB.BINARY, name="team_t{0}_s{1}".format(t,i)) for t in range(T)] for i in range(Q)]

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
    

    ni = [[[[model.addVar(vtype=GRB.INTEGER, name="ni_s{0}_w{1}_t{2}_k{3}".format(i, w, t, k)) for k in range(K)] for t in range(T)]for w in range(W)]for i in range(Q)]           
    nb = [[[[model.addVar(vtype=GRB.BINARY, name="nb_s{0}_w{1}_t{2}_k{3}".format(i, w, t, k)) for k in range(K)] for t in range(T)]for w in range(W)]for i in range(Q)]
    X = [[[[model.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="X(s%d,w%d,k%d,t%d)" %(i, w,k,t))  for k in range(K)] for t in range(T)] for w in range(W)] for i in range(Q)]

    # n_wtk[w][t][k] # integer value of N_wk[w][k] * pi[w][t][k]

    overflow = [] # overflow[w][r]
    for w in range(W):
        overflow.append([])
        for r in range(R):
            tmp_overflow_var = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="o_w{0}_r{1}".format(w, r)) # this is original integer but can be approximated by continuous value
            overflow[w].append(tmp_overflow_var)
    
    O = [[[model.addVar(vtype=GRB.CONTINUOUS, name="O_{0}_w{1}_r{2}".format(i, w, r)) for r in range(R)] for w in range(W)] for i in range(Q)]
    
    y = [] # y[w][r]: number of operating resources r at time w
    for w in range(W):
        y.append([])
        for r in range(R):
            tmp_resource_r = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="y_w{0}_r{1}".format(w, r))
            y[w].append(tmp_resource_r)
            
    yi = ys
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

    #TEAM CONSTRAINTS 
    for i in range(Q):
        for w in range(W):
            for k in range(K):
                for t in range(T):
                    #tmp = LinExpr([1/(np.maximum(N_wk[w][k],1)) for w in range(W) for k in range(K)], [ni[i][w][t][k] for w in range(W) for k in range(K)])
                    model.addConstr( n[i][w][t][k]+nb[i][w][t][k] <= team[i][t]*N_wk[w][k], name="team{0}{1}{2}{3}".format(t,i,k,w))         
                    model.addConstr( n[i][w][t][k]+nb[i][w][t][k]-ni[i][w][t][k]==0)
        model.addConstr(quicksum(team[i][t] for t in range(T)) <= maxT)

    
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
                    model.addConstr(marginal_sum + tmp_sum + overflow[w-1][r] - yi[i][w][r] * C[r] - O[i][w][r] <= 0, name="(11)_w{0}_r{1}_{2}".format(w, r, i))
    
                model.addConstr( marginal_sum+ tmp_sum <= yi[i][w][r]*10000, name="(10.5)_w{0}_r{1}".format(w, r))
        
    for w in range(W):
        for t in range(T):
            for k in range(K):
                tmp_sum = LinExpr([1 for i in range(Q)], [X[i][w][t][k] for i in range(Q)])
                marginal_sum = LinExpr([n[i][w][t][k] for i in range(Q)], [q[i] for i in range(Q)])
                model.addConstr(marginal_sum + tmp_sum == n_wtk[w][t][k])
                
    for w in range(W):
        for r in range(R):
            tmp_sum = LinExpr([yi[i][w][r] for i in range(Q)],[q[i] for i in range(Q)])
            model.addConstr(tmp_sum == y[w][r])
    
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

    model.write("tsgkpMIPteam.lp")
    start_time = time.time()
    model.optimize()
    runtime = time.time() - start_time

    model.write("tsgkpMIPteam.sol")

    ni_value = np.zeros((Q,W,T,K))
    for i in range(Q):
        for w in range(W):
            for t in range(T):
                for k in range(K):
                    ni_value[i][w][t][k] = ni[i][w][t][k].x

    n_value = np.zeros((W,T,K))
    for w in range(W):
        for t in range(T):
            for k in range(K):
                n_value[w][t][k] = n_wtk[w][t][k].x
                #if n_value[w][t][k] != int(n_value[w][t][k]):
                    #print n_value[w][t][k]

    O_value = np.zeros((Q, W, R))
    for i in range(Q):
        for w in range(W):
            for r in range(R):
                O_value[i][w][r] = O[i][w][r].x

    q_value = np.zeros(Q)
    for i in range(Q):
        q_value[i] = q[i].x

    team_val = np.zeros((Q,T))
    for i in range(Q):
        for t in range(T):
            team_val[i][t] = team[i][t].x

    overflow_value = np.zeros((W, R))
    for w in range(W):
        for r in range(R):
            overflow_value[w][r] = overflow[w][r].x
    
    obj = model.getAttr('ObjVal')

    model.terminate()

    return obj, runtime, team_val, ni_value, O_value, q_value, overflow_value


def KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q, ar, phi, integer=0, OverConstr=False, OverConstr2=False): # integer indicates different relaxation method
    # ======================= Gurobi Setting ===================================
    model = Model("MIP")
    model.params.DualReductions = 0
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
                if (integer == 2) or (integer == 3):
                    tmp_pi_var = model.addVar(vtype=GRB.INTEGER, name="n_w{0}_t{1}_k{2}".format(w, t, k))
                else:
                    tmp_pi_var = model.addVar(vtype=GRB.CONTINUOUS, name="n_w{0}_t{1}_k{2}".format(w, t, k))
                n_wtk[w][t].append(tmp_pi_var)
                
    ni_wtk = [[[[model.addVar(vtype=GRB.CONTINUOUS, name="ni_s{0}_w{1}_t{2}_k{3}".format(i, w, t, k)) for k in range(K)] for t in range(T)]for w in range(W)]for i in range(Q)]
    # n_wtk[w][t][k] # integer value of N_wk[w][k] * pi[w][t][k]

    overflow = [] # overflow[w][r]
    for w in range(W):
        overflow.append([])
        for r in range(R):
            tmp_overflow_var = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="o_w{0}_r{1}".format(w, r)) # this is original integer but can be approximated by continuous value
            overflow[w].append(tmp_overflow_var)
    
    #O = [[[model.addVar(vtype=GRB.CONTINUOUS, name="O_{0}_w{1}_r{2}".format(i, w, r)) for r in range(R)] for w in range(W)] for i in range(Q)]
    oi = [[[model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="oi_w{0}_r{1}_s{2}".format(i,w, r)) for r in range(R)] for w in range(W)] for i in range(Q)]
    
    y = [[model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="y_w{0}_r{1}".format(w, r)) for r in range(R)] for w in range(W)]
    
    yi = [[[model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="yi_w{0}_r{1}_s{2}".format(i,w, r)) for r in range(R)] for w in range(W)] for i in range(Q)]
    
    yb = [[[model.addVar(vtype=GRB.BINARY, lb=0, name="yb_w{0}_r{1}_s{2}".format(i,w, r)) for r in range(R)] for w in range(W)] for i in range(Q)]
    
    yb2 = [[[model.addVar(vtype=GRB.BINARY, lb=0, name="yb_w2{0}_r{1}_s{2}".format(i,w, r)) for r in range(R)] for w in range(W)] for i in range(Q)]
     # y[i][w][r]: number of operating resources r at time w
    

    p = [] # available staff
    s = [] # working staff
    for w in range(W):
        tmp_staff = model.addVar(vtype=GRB.CONTINUOUS, name="p_w{0}".format(w))
        tmp_working_staff = model.addVar(vtype=GRB.CONTINUOUS, name="s_w{0}".format(w))
        p.append(tmp_staff)
        s.append(tmp_working_staff)
        
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
    for w in range(W):
        for r in range(R):
            for i in range(Q):
                tmp_sum = LinExpr([1 for k in range(K)]*len(resource2team[r]), [ni_wtk[i][w][t][k] for t in resource2team[r] for k in range(K)])
                if w == 0:
                    model.addConstr(tmp_sum + pre_overflow[r] - yi[i][w][r] * C[r] - oi[i][w][r]<= 0, name="(4)_w{0}_r{1}".format(w, r,i))
                else:
                    model.addConstr(tmp_sum + oi[i][w-1][r] - yi[i][w][r] * C[r] - oi[i][w][r] <= 0, name="(4)_w{0}_r{1}".format(w, r,i))
    
                model.addConstr( tmp_sum <= yi[i][w][r]*10000 ,  name="(5.6)_w{0}_r{1}".format(w, r))
    
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
                model.addConstr( quicksum(ni_wtk[i][w][t][k] for t in range(T)) == N_wk[w][k])  
                
    #for i in range(Q):
    #    for w in range(W):
    #        for r in range(R):  
    #            tmp_sum = LinExpr([1 for k in range(K)]*len(resource2team[r]), [ni[i][w][t][k] for t in resource2team[r] for k in range(K)])
    #            if w == 0:
    #                model.addConstr(tmp_sum + pre_overflow[r] - yi[i][w][r] * C[r] - O[i][w][r] <= 0, name="(10)_w{0}_r{1}_s{2}".format(w, r, i))
    #            else:
    #                model.addConstr(tmp_sum + O[i][w-1][r] - yi[i][w][r] * C[r] - O[i][w][r] <= 0, name="(10)_w{0}_r{1}_s{2}".format(w, r, i ))
    #
    #            model.addConstr( tmp_sum <= y[w][r]*10000, name="(10.5)_w{0}_r{1}".format(w, r))
    
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
                tmp_sum = LinExpr([q[i] for i in range(Q)], [ni_wtk[i][w][t][k] for i in range(Q)])
                model.addConstr(tmp_sum == n_wtk[w][t][k])
        for r in range(R):
            tmp_sum = LinExpr([q[i] for i in range(Q)], [oi[i][w][r] for i in range(Q)])
            model.addConstr(tmp_sum == overflow[w][r])
            tmp_sum = LinExpr([q[i] for i in range(Q)], [yi[i][w][r] for i in range(Q)])
            model.addConstr(tmp_sum == y[w][r])


    for w in range(W):
        for r in range(R):
            for i in range(Q):
                model.addConstr(yi[i][w][r] - mr[r] <= 0, name="(7)_w{0}_r{1}_{2}".format(w, r,i))
            
            
    for w in range(W):
        for r in range(R):
            for i in range(Q):
                model.addConstr(yi[i][w][r] - yb[i][w][r]- minr[w][r] == 0, name="(11)_w{0}_r{1}_{2}".format(w, r,i))
                          
    #for w in range(W):
    #    for r in range(R):
    #        for i in range(Q):
    #            model.addConstr(yi[i][w][r] + yb2[i][w][r]- yb[i][w][r]- minr[w][r] == 0, name="(11)_w{0}_r{1}_{2}".format(w, r,i))
    #            model.addConstr(yb2[i][w][r]+yb[i][w][r]<= 1, name="(11)_w{0}_r{1}_{2}".format(w, r,i))
                
    #for i in range(Q):         
    #    model.addConstr(quicksum(yb2[i][w][r] for r in range(R) for w in range(W)) <= 10)
              
    for w in range(W):
        for i in range(Q):
            tmp_sum = LinExpr(ar, [yi[i][w][r] for r in range(R)])
            model.addConstr(tmp_sum - p[w] <= 0, name="(6)_w{0}".format(w))

                
    for w in range(W):
        start_index = max(0, w - shift + 1)
        tmp_sum = LinExpr([1]*(w - start_index + 1), [s[i] for i in range(start_index, w+1)])
        model.addConstr(tmp_sum - p[w] == 0, name="(8)_w{0}".format(w))

    tmp_sum = LinExpr([1]*W, [s[w] for w in range(W)])
    model.addConstr(tmp_sum - P <= 0, name="(9)")
    
    

    #TEAM CONSTRAINTS 
    if True:
        for i in range(Q):
            for w in range(W):
                for k in range(K):
                    for t in range(T):
                        #tmp = LinExpr([1/(np.maximum(N_wk[w][k],1)) for w in range(W) for k in range(K)], [ni[i][w][t][k] for w in range(W) for k in range(K)])
                        model.addConstr( ni_wtk[i][w][t][k] <= team[i][t]*N_wk[w][k], name="team{0}{1}{2}{3}".format(t,i,k,w))         
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
    for i in range(Q):
        for w in range(W):
            for t in range(T):
                for k in range(K):
                    ns_value[i][w][t][k] = ni_wtk[i][w][t][k].x
                    
    ys_value = np.zeros((Q,W,R))

    for w in range(W):
        for r in range(R):
            for i in range(Q):
                ys_value[i][w][r] = yi[i][w][r].x
                
    z_value = np.zeros((W,K,M))
    for w in range(W):
        for k in range(K):
            for m in range(M):
                z_value[w][k][m] = z[w][k][m].x
    
    
    overflow_value = np.zeros((W,R))
    for w in range(W):
        for r in range(R):
            overflow_value[w][r] = overflow[w][r].x
            
    s_value = np.zeros(W)
    for w in range(W):
        s_value[w] = s[w].x
        
    p_value = np.zeros(W)
    for w in range(W):
        p_value[w] = p[w].x
        
    y_value = np.zeros((W,R))
    for w in range(W):
        for r in range(R):
            y_value[w][r] = y[w][r].x
    
    team_val = np.zeros((Q,T))
    for i in range(Q):
        for t in range(T):
            team_val[i][t] = team[i][t].x
    
    obj = model.getAttr('ObjVal')

    model.terminate()

    return obj, n_value, ns_value, ys_value,z_value,p_value,s_value,y_value
    

def solve(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, TeamConstr=False):
   
    start_time = time.time()

    print "============================ FULL LP relaxation =============================="
    obj_relax, n_value0, overflow_value0, y_value, s_value0, p, attset, f = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=0, OverConstr=False, TeamConstr=TeamConstr, MaxT=maxT, Q=Q)
    
    
    minr = np.zeros((W,R))

    for w in range(W):
        for r in range(R):
            minr[w][r] = math.floor(y_value[w][r])
            #print y_value[w][r]
    
    q = np.zeros(Q)
    for i in range(Q):
        q[i] = float(1)/Q
        
    
    objyn, n, ns,ys = KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q, p, s_value0, ar, phi, integer=0, OverConstr=False, OverConstr2=False)
            
    #print "============================ Binary Y / single relaxed n_wtk (allocated arrivals) MIP ======================"
    #[x] = KStrategiesY( Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, OverConstr=True,  TeamConstr=TeamConstr), 
    #objy, n_value, overflow_value, y_value, ys,  s_value, p_value, q, t1
    
    #[objy, y_value, ys,  s_value, p_value, q, t1] = [x[0],x[3],x[4],x[5],x[6],x[7],x[8]]
    #print "============================ multiple relaxed n_wtk (allocated arrivals) MIP ==============================="

    #objyn, n_val, n, t2 = KStrategiesYN(Q, W, K, R, M, P, resource2team, T,maxT,  E, C, U_plus, U_minus, N_wk, shift, mr, y_value, ys, q, p_value, s_value, ar, phi, integer=0, OverConstr=True, OverConstr2=False)
    
    #minn = np.zeros((Q,W,T,K))
    
    #for i in range(Q):
    #    for w in range(W):
    #        for k in range(K):
    #            sum = 0 
    #            for t in range(T):
    #                minn[i][w][t][k] = math.floor(n[i][w][t][k])
    #                sum += math.floor(n[i][w][t][k])
    
    print "============================ Integer n_wtk  MIP ==============================="

    #obj, rt, t3  = KStrategiesYNB(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, y_value, ys, minn, p, s_value, phi, integer=0, OverConstr=True, OverConstr2=False)
    #walltime = time.time() - start_time
    
    #print "Runtime/walltime ", rt, " ", walltime
    
    print obj_relax, objyn#, objyn, obj
    
    return [obj_relax, objyn]#, walltime, [t1, t2, t3]

def randomSetting(seed, W, K ,R, mR, M, P, teams, shift):
    # ========================== Random Seed ===================================
    #seed = 3234
    #seed = random.randint(1,10000)
    np.random.seed(seed)
    random.seed(seed)

    T = len(teams)
    print T
    resource2team = util.resourceTeamDict(R, T, teams)
    print resource2team

    Er = np.random.rand(R, M)/2 + 0.5 # Er[m][r]
    Er = Er / 2
    print "Er"
    print Er

    Er, C = util.genResources(R, M, 600)
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
    print "\nutilities"
    print "plus"
    print U_plus
    print "minus"
    print U_minus

    N_wk = [[ 0 for k in range(K)] for w in range(W)] # N_wk[w][k] represents the number of people getting in time window w with type k
    
    #arrival_window= np.minimum(3,W)

    N_wk = np.zeros((W,K))        
    for k in range(K):
        startK = random.randint(3,W)
        for w in range(startK-3,startK):
            large_or_small = random.random()
            if large_or_small > 0.5:
                tmp_N = random.randint(50, 200)
            else:
                tmp_N = random.randint(10, 30)
            N_wk[w][k] = tmp_N
    """        
    for w in range(W):
        N_wk.append([])
        for k in range(K):
            large_or_small = random.random()
            if large_or_small > 0.5:
                tmp_N = random.randint(100, 300)
            else:
                tmp_N = random.randint(10, 30)
            N_wk[w].append(tmp_N)
    """
#    print "number of passengers: {0}".format([sum(N_wk[w]) for w in range(W)])

    #C = np.random.randint(200, 500, R) # C[r] is the capacity of resource r
 #   print "\nC"
  #  print C
    #C = [100,80,75,50,30,15]
    mr = np.random.randint(5, 15, R)
 #   print "\nmr"
  #  print mr
    #mr = [5, 5, 5, 3, 2, 4] # maximum number of people to operate resource r
    ar = np.random.randint(1, 5, R)
 #   print "\nar"
  #  print ar
    #ar = [2, 1, 2, 1, 1, 2] # ar[r] number of people required to operate resource r

    phi = np.random.rand(R) # phi[r] overflow penalty
   # print "\nphi"
    #print phi
    #phi = np.random.rand(W, R) # phi[w][r] overflow penalty

    return resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi

def fullYNcombined(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, Q, maxT):
    minr = np.zeros((W,R))
    obj_relax, n_value0, overflow_value, y_value, s_value0, p,z_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=0, OverConstr = 0)
    for w in range(W):
        for r in range(R):
            minr[w][r] = math.floor(y_value[w][r])

    q = np.zeros(Q)
    for i in range(Q):
        q[i] = float(1)/Q
    #q = [0.25, 0.75]
    
    print "============================ K strategies Y, N combined =============================="
    objyn1, n, ns, ys, z_value, p, s, y = KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q, ar, phi, integer=0, OverConstr=False, OverConstr2=False)
    minn = np.zeros((Q,W,T,K))
    
    for i in range(Q):
        for w in range(W):
            for k in range(K):
                #sum = 0 
                for t in range(T):
                    minn[i][w][t][k] = math.floor(ns[i][w][t][k])
                    #sum += math.floor(ns[i][w][t][k])
    
    print "============================ K strategies Y, N, B new =============================="
    obj1, rt, t3, ni_value, O_value, q_value, overflow_value  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)
    print obj_relax, objyn1, obj1

    new_O_value = np.zeros((Q, W, R))
    for i in range(Q):
        for w in range(W):
            for r in range(R):
                if w > 0:
                    new_O_value[i][w][r] = max(new_O_value[i][w-1][r] + np.sum([ni_value[i][w][t][k] for k in range(K) for t in resource2team[r]]) - ys[i][w][r] * C[r], 0)
                else:
                    new_O_value[i][w][r] = max(np.sum([ni_value[i][w][t][k] for k in range(K) for t in resource2team[r]]) - ys[i][w][r] * C[r], 0)

    new_overflow = np.zeros((W, R))
    for i in range(Q):
        new_overflow += new_O_value[i] * q_value[i]
    print "old overflow: {0}".format(overflow_value)
    print "new overflow: {0}".format(new_overflow)
                

    strategySet = []
    for i in range(Q):
        tmp_strategy = {}
        tmp_strategy["n"] = ni_value[i]
        tmp_strategy["overflow"] = new_O_value[i]
        strategySet.append(tmp_strategy)

    return strategySet, obj_relax, objyn1, obj1, q_value

if __name__ == "__main__":

    W = 5 # number of time windows
    K = 3 # number of passenger types
    R = 6 # number of resources
    mR = 3 # max number of reosurces
    M = 2 # number of attack methods
    P = 15 # number of staff
    shift = 2 # d
    Q = 5
    nT = 25
    teams = util.generateAllTeams(R, mR)
    maxT = 5


    # ================= random generate game setting ===========================
    seed = 3345
    #resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
    resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)

    print "============================ FULL LP relaxation =============================="
    start_time = time.time()
    minr = np.zeros((W,R))
    
    obj_relax, n_value0, overflow_value, y_value, s_value0, p,z_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=0, OverConstr = 0)

    for w in range(W):
        for r in range(R):
            minr[w][r] = math.floor(y_value[w][r])
            #print y_value[w][r]
    
    #objyn1, n_value0, overflow_value, y_value, s_value0, p,z_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=1, OverConstr = 0)
    
    q = np.zeros(Q)
    for i in range(Q):
        q[i] = float(1)/Q
    #q = [0.25, 0.75]
    
    print "============================ K strategies Y, N combined =============================="
    objyn1, n, ns,ys,z_value,p,s,y = KStrategiesYNcomb(Q, W, K, R, M, P, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, minr, q, ar, phi, integer=0, OverConstr=False, OverConstr2=False)
    minn = np.zeros((Q,W,T,K))
    
    for i in range(Q):
        for w in range(W):
            for k in range(K):
                #sum = 0 
                for t in range(T):
                    minn[i][w][t][k] = math.floor(ns[i][w][t][k])
                    #sum += math.floor(ns[i][w][t][k])
    
    print "============================ K strategies Y, N, B new =============================="
    obj1, rt, t3, ni_value, O_value, q_value  = KStrategiesYNBnew(Q, W, K, R, M, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, ys, minn, p, s, phi, integer=0, OverConstr=False, OverConstr2=False)
    
    print obj_relax, objyn1, obj1
    
