import util
from gurobipy import *
import numpy as np
import random
import math
from KStrategies import randomSetting
from StaffResourceAllocation import LPsolver

def integerSolution(Nmax, Q, W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=0, binary_y=0, OverConstr=False): # integer indicates different relaxation method
    # ======================= Gurobi Setting ===================================
    model = Model("MIP")
    model.params.DualReductions = 0
    model.params.MIPGap=0.0001;

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

    y = [] # y[w][r]: number of operating resources r at time w
    for w in range(W):
        y.append([])
        for r in range(R):
            if ((integer == 1) or (integer == 3)) and (binary_y == 0):
                tmp_resource_r = model.addVar(vtype=GRB.INTEGER, lb=0, name="y_w{0}_r{1}".format(w, r))
            else:
                tmp_resource_r = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="y_w{0}_r{1}".format(w, r))
            y[w].append(tmp_resource_r)

    p = [] # available staff
    s = [] # working staff
    for w in range(W):
        tmp_staff = model.addVar(vtype=GRB.CONTINUOUS, name="p_w{0}".format(w))
        tmp_working_staff = model.addVar(vtype=GRB.CONTINUOUS, name="s_w{0}".format(w))
        p.append(tmp_staff)
        s.append(tmp_working_staff)
        
    yi = [] # binary_y
    for w in range(W):
        yi.append([])
        for r in range(R):
            yi[w].append([])
            for i in range(Q):
                tmp_resource_binary_r = model.addVar(vtype=GRB.INTEGER, lb=0, name="yb_w{0}_r{1}_s{2}".format(w, r, i))
                yi[w][r].append(tmp_resource_binary_r)
            
    ni_wtk = [] # n_wtk[w][t][k] # integer value of N_wk[w][k] * pi[w][t][k]
    for w in range(W):
        ni_wtk.append([])
        for t in range(T):
            ni_wtk[w].append([])
            for k in range(K):
                ni_wtk[w][t].append([])
                for i in range(Q):
                    tmp_pi_var = model.addVar(vtype=GRB.INTEGER, name="n_w{0}_t{1}_k{2}_s{3}".format(w, t, k, i))
                    ni_wtk[w][t][k].append(tmp_pi_var)
                    
    oi = [] # overflow[w][r]
    for w in range(W):
        oi.append([])
        for r in range(R):
            oi[w].append([])
            for i in range(Q):
                tmp_overflow_var = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="o_w{0}_r{1}_s{2}".format(w, r,i)) # this is original integer but can be approximated by continuous value
                oi[w][r].append(tmp_overflow_var)

    
    q = [ model.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="q_s{0}".format(i)) for i in range(Q)]
    
    
    X = [[[[[model.addVar(vtype=GRB.BINARY, name="X_{0}_w{1}_t{2}_{3}_{4}".format(w,t,k,i,j)) for j in range(Nmax)] for i in range(Q)]for k in range(K)] for t in range(T)] for w in range(W)]                    
    O = [[[[model.addVar(vtype=GRB.BINARY, name="O_w{0}_r{1}_{2}_{3}".format(w, r,i,j)) for j in range(Nmax)] for i in range(Q)] for r in range(R)] for w in range(W)]
    XC = [[[[[model.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="XC_{0}_w{1}_t{2}_{3}_{4}".format(w,t,k,i,j)) for j in range(Nmax)] for i in range(Q)]  for k in range(K)] for t in range(T)] for w in range(W)]
    OC = [[[[model.addVar(vtype=GRB.CONTINUOUS, name="OC_{0}_w{1}_r{2}_{3}".format(w, r,i,j)) for j in range(Nmax)] for i in range(Q)] for r in range(R)] for w in range(W)]
    # ========================= Gurobi Objective ===============================
    model.update()

    objective_variables = [theta] + [overflow[w][r] for w in range(W) for r in range(R)]
    objective_coefficients = [1] + [-phi[r] for r in range(R)]*W
    objective_value = LinExpr(objective_coefficients, objective_variables)

    #objective_value = theta
    #for w in range(W):
    #    for r in range(R):
    #        objective_value += phi[r]*overflow[w][r]
    model.setObjective(objective_value, GRB.MAXIMIZE)
    model.update()

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
            for i in range(Q):
                #tmp_sum = LinExpr([N_wk[w][k] for k in range(K)]*len(resource2team[r]), [pi[w][t][k] for t in resource2team[r] for k in range(K)])
                tmp_sum = LinExpr([1 for k in range(K)]*len(resource2team[r]), [ni_wtk[w][t][k][i] for t in resource2team[r] for k in range(K)])
                model.addConstr(yi[w][r][i]*10000 >= tmp_sum, name="(5.6)_w{0}_r{1}".format(w, r,i))
                if w == 0:
                    model.addConstr(tmp_sum + pre_overflow[r] - yi[w][r][i] * C[r] - oi[w][r][i] <= 0, name="(4)_w{0}_r{1}".format(w, r,i))
                else:
                    model.addConstr(tmp_sum + oi[w-1][r][i] - yi[w][r][i] * C[r] - oi[w][r][i] <= 0, name="(4)_w{0}_r{1}".format(w, r,i))

    if OverConstr:
        for r in range(R): # OPTIONAL
            model.addConstr(overflow[W-1][r] == 0, name="(5)_r{0}".format(r))
    
        for w in range(W):
            for r in range(R):
                if w > 0:
                    model.addConstr(y[w][r] * C[r] - overflow[w-1][r] >= 0, name="(5.5)_w{0}_r{1}".format(w, r))

    for w in range(W):
        for t in range(T):
            for k in range(K):
                for i in range(Q):
                    tmp_sum = quicksum(X[w][t][k][i][j] for j in range(Nmax))
                    model.addConstr(tmp_sum-ni_wtk[w][t][k][i] == 0, name="(5.5)_w{0}_t{1}_k{2}_s{3}".format(w, t,k,i))
                    for j in range(1,Nmax):
                        model.addConstr(X[w][t][k][i][j-1]-X[w][t][k][i][j] >= 0, name="(5.5)_w{0}_t{1}_k{2}_s{3}_j{4}".format(w, t,k,i,j))    
                       
    for w in range(W):
        for r in range(R):
            for i in range(Q):
                tmp_sum = quicksum(O[w][r][i][j] for j in range(Nmax))
                model.addConstr(tmp_sum-oi[w][r][i] == 0, name="(5.5)_w{0}_r{1}_s{2}".format(w, r,i))
                for j in range(1,Nmax):
                    model.addConstr(O[w][r][i][j-1]-O[w][r][i][j] >= 0, name="(5.5)_w{0}_r{1}_s{2}_j{3}".format(w, r,i,j))    
            
    for w in range(W):
        for t in range(T):
            for k in range(K):
                tmp_sum = quicksum(XC[w][t][k][i][j] for j in range(Nmax) for i in range(Q))
                model.addConstr(tmp_sum-n_wtk[w][t][k] == 0, name="(5.5)_w{0}_t{1}_k{2}".format(w, t,k))
                for i in range(Q):
                    for j in range(Nmax):
                        model.addConstr(XC[w][t][k][i][j] <= q[i])
                        model.addConstr(XC[w][t][k][i][j] <= X[w][t][k][i][j])
                        model.addConstr(XC[w][t][k][i][j] >= q[i]-(1-X[w][t][k][i][j]))
                        
    for w in range(W):
        for r in range(R):
            tmp_sum = quicksum(OC[w][r][i][j] for j in range(Nmax) for i in range(Q))
            model.addConstr(tmp_sum-overflow[w][r] == 0, name="(5.5)_w{0}_r{1}".format(w, r))
            for i in range(Q):
                for j in range(Nmax):
                    model.addConstr(OC[w][r][i][j] <= q[i])
                    model.addConstr(OC[w][r][i][j] <= O[w][r][i][j])
                    model.addConstr(OC[w][r][i][j] >= q[i]-(1-O[w][r][i][j]))
                    
    model.addConstr(quicksum(q[i] for i in range(Q))==1)
            

    for w in range(W):
        for i in range(Q):
            tmp_sum = LinExpr(ar, [yi[w][r][i] for r in range(R)])
            model.addConstr(tmp_sum - p[w] <= 0, name="(6)_w{0}".format(w))

    for w in range(W):
        for r in range(R):
            for i in range(Q):
                model.addConstr(yi[w][r][i] - mr[r] <= 0, name="(7)_w{0}_r{1}".format(w, r,i))
                            
#                model.addConstr(yb[w][r] == 0, name="(11)_w{0}_r{1}".format(w,r))
#            else:
                
    for w in range(W):
        start_index = max(0, w - shift + 1)
        tmp_sum = LinExpr([1]*(w - start_index + 1), [s[i] for i in range(start_index, w+1)])
        model.addConstr(tmp_sum - p[w] == 0, name="(8)_w{0}".format(w))

    tmp_sum = LinExpr([1]*W, [s[w] for w in range(W)])
    model.addConstr(tmp_sum - P <= 0, name="(9)")

    """
    for w in range(W):
        for t in range(T-1):
            for k in range(K):
                model.addConstr(pi[w][t][k] == 0, name="(10-1)_w{0}_t{1}_k{2}".format(w, t, k))

    for w in range(W):
        for r in range(R):
            model.addConstr(overflow[w][r] == 0, name="(10-2)_w{0}_r{1}".format(w, r))

    for w in range(W):
        model.addConstr(s[w] == 0, name="(11-1)_w{0}".format(w))

    for w in range(W):
        for r in range(R):
            model.addConstr(y[w][r] == 0, name="(11-2)_w{0}_r{1}".format(w, r))
    """
    model.update()

    model.write("tsg.lp")

    model.optimize()

    #if True:
    #    for v in model.getVars():
    #        if v.x > 0:
    #            print "{0} {1}".format(v.varName, v.x)

    #defender_utility = np.inf
    #for w in range(W):
    #    for k in range(K):
    #        for m in range(M):
    #            tmp_z_wkm = model.getVarByName("z_w{0}_k{1}_m{2}".format(w, k, m)).x
    #            tmp_utility = U_plus[k] * tmp_z_wkm + (1 - tmp_z_wkm) * U_minus[k]
    #            if tmp_utility < defender_utility:
    #               defender_utility = tmp_utility
    #            print "utility of w={0}, k={1}, m={2} is: {3}".format(w,k,m, tmp_utility)

    #print "defender utility: {0}".format(defender_utility)

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
    for w in range(W):
        for r in range(R):
            y_value[w][r] = y[w][r].x

    s_value = np.zeros(W)
    for w in range(W):
        s_value[w] = s[w].x
        
    p_value = np.zeros(W)
    for w in range(W):
        p_value[w] = p[w].x
    
    obj = model.getAttr('ObjVal')
    
    oi_value = np.zeros((W,R,Q))
    for w in range(W):
        for r in range(R):
            for i in range(Q):
                oi_value[w][r] = oi[w][r][i].x

    yi_value = np.zeros((W,R,Q))
    for w in range(W):
        for r in range(R):
            for i in range(Q):
                yi_value[w][r][i] = yi[w][r][i].x
                
    ni_value = np.zeros((W,T,K,Q))
    for w in range(W):
        for t in range(T):
            for k in range(K):
                for i in range(Q):
                    ni_value[w][t][k][i] = ni_wtk[w][t][k][i].x
                
    q_value = np.zeros(Q)
    for i in range(Q):
        q_value[i] = q[i].x
        
    z_value = np.zeros((W,K,M))
    for w in range(W):
        for k in range(K):
            for m in range(M):
                z_value[w][k][m] = z[w][k][m].x
    

    return obj, n_value, overflow_value, y_value, s_value, p_value,ni_value,oi_value,yi_value,q_value,z_value

if __name__ == "__main__":
    # ============================= main =======================================
    print "======================== main ======================================"
    # ========================= Game Setting ===================================
    W = 2 # number of time windows
    K = 3 # number of passenger types
    R = 2 # number of resources
    mR = 2 # max number of reosurces
    M = 1 # number of attack methods
    P = 3 # number of staff
    Q= 1
    shift = 1 # d

    nT = 22
    teams = util.generateAllTeams(R, mR)
    #teams = util.randomGenerateTeams(R, mR, nT)
    Nmax = 150

    #print teams

    # ================= random generate game setting ===========================
    seed = 2345
    resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
    
    minr = np.zeros((W,R))
    
    #obj, n_value, overflow_value, y_value, s_value, p_value,z_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=3, binary_y=0, OverConstr = 0)
    #obj1, n_value, overflow_value, y_value, s_value, p_value,z_value = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=3, binary_y=0, OverConstr = 0)
    from DesignProblemFixedResources import solve
    
    obj0, n_value0, overflow_value0, y_value0, s_value0, p_value0,ni_value,oi_value,yi_value,q_value,z_value0 = integerSolution(Nmax,Q,W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=0, binary_y=0, OverConstr = 0)
    util1, rt1, q1, teams1 = solve(Q, W, K, R, mR, M, P, teams, resource2team, T, 5, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi,  TeamConstr=False, YConstr=False)
    
    
    print obj0, util1
    
    
    
    
