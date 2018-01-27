import util
from gurobipy import *
import numpy as np
import random
import math
from relaxed_feed import LPsolver, LPsolverR
from KStrategiesFixedYRoundN import Ksolver

def KStrategiesY(Q, W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, miny, ar, phi, integer=0, OverConstr=False): # integer indicates different relaxation method
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
    YX = [[[model.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="X(s%d,w%d,r%d)" %(i,w,r))  for r in range(R)] for w in range(W)] for i in range(Q)]
 
    nb = [[[[model.addVar(lb=0.0, ub = 1.0, vtype=GRB.BINARY, name="Np(s%d,w%d,k%d,t%d)" %(i,w,k,t))  for k in range(K)]  for t in range(T)] for w in range(W)] for i in range(Q)]
    NX = [[[[model.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="X(s%d,w%d,k%d,t%d)" %(i,w,k,t))  for k in range(K)] for t in range(T)] for w in range(W)] for i in range(Q)]
 
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
            model.addConstr(y[w][r] - quicksum(YX[i][w][r] for i in range(Q)) - miny[w][r] == 0, name="(11)_w{0}_r{1}".format(w, r))
    
    # Mixed Strategy Constraints
    for w in range(W):
        for k in range(K):
            for t in range(T):
                model.addConstr(n_wtk[w][t][k] == nmin[w][t][k] + quicksum(NX[i][w][t][k] for i in range(Q)))
        
        for r in range(R):
            nmarginal = quicksum(nmin[w][t][k]*q[i] for t in resource2team[r] for k in range(K) for i in range(Q))
            nX = [ NX[i][w][t][k] for t in resource2team[r] for k in range(K) for i in range(Q)]
                
            tmp_sum = LinExpr([1 for k in range(K) for i in range(Q)]*len(resource2team[r]), nX)
            
            model.addConstr( tmp_sum + nmarginal - y[w][r]*10000 <=0 , name="(5.6)_w{0}_r{1}".format(w, r))    
            if w == 0:
                model.addConstr(tmp_sum + nmarginal - y[w][r] * C[r] - overflow[w][r] <= 0)
            else:
                model.addConstr(tmp_sum + nmarginal + overflow[w-1][r] - y[w][r] * C[r] - overflow[w][r] <= 0)
            
    
    tmp_sum = LinExpr([1]*Q, [q[i] for i in range(Q)])
    model.addConstr(tmp_sum == 1, name="sumQ")         
 
    # Linearization Constraints
    for i in range(Q):
        for w in range(W):
            for r in range(R):
                model.addConstr(YX[i][w][r] <= q[i]) 
                model.addConstr(YX[i][w][r] <= yb[i][w][r] ) 
                model.addConstr(YX[i][w][r] >= q[i] -(1-yb[i][w][r])) 
    
    # Linearization Constraints
    for i in range(Q):
        for w in range(W):
            for k in range(K):
                model.addConstr( quicksum( nmin[w][t][k] + nb[i][w][t][k] for t in range(T)) == N_wk[w][k])  
                for t in range(T):
                    model.addConstr(NX[i][w][t][k] <= q[i]) 
                    model.addConstr(NX[i][w][t][k] <= nb[i][w][t][k] ) 
                    model.addConstr(NX[i][w][t][k] >= q[i] -(1-nb[i][w][t][k])) 
            
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
    K = 10 # number of passenger types
    R = 3 # number of resources
    mR = 10 # max number of reosurces
    M = 1 # number of attack methods
    P = 10 # number of staff
    shift = 5 # d
    Q = 15
    nT = 20
    teams = util.generateAllTeams(R, mR)
    #teams = util.randomGenerateTeams(R, mR, nT)


    # ================= random generate game setting ===========================
    seed = 2345
    resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)

    print "============================ LP relaxation =============================="
    obj_relax, n_value0, overflow_value0, y_value0, s_value0, p, attset, f = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=1, OverConstr=True)
    for w in range(W):
        for r in range(R):
            minr[w][r] = math.floor(y_value0[w][r])
            print y_value0[w][r]
            
    print "============================ relaxed n_wtk (allocated arrivals) MIP ==============================="
    objy, n_value, overflow_value, y_value, ys,  s_value, p_value, q = KStrategiesY( Q, W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0)
    
    n =[]
    objs =[]
    for i in range(Q):
        obj_relaxi, n_value0, overflow_value0, attset, f = LPsolverR(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, ys[i], s_value, p_value, mr, ar, phi, integer=0, OverConstr=False)
        n.append(n_value0)
        objs.append(obj_relaxi)
    obj, rt, q, n2, o, att_set = Ksolver(W, K, R, mR, M, P, Q,5, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, n,overflow_value, ys, s_value, p_value, integer=0, OverConstr=False)
    print y_value0

    print y_value
    
    print obj_relax, objy, obj, 
    print objs
    print y_value
    for i in range(Q):
        if q[i]>0:
            print q[i]
            print ys[i]



