import util
from gurobipy import *
import numpy as np
import random
import time
from relaxed_feed import LPsolver
import pickle

def Ksolver(W, K, R, mR, M, P, Q, QN, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, N_marginal, O_marginal, y_marginal, s_marginal, p_marginal, integer=0, OverConstr=False): # integer indicates different relaxation method
    # ======================= Gurobi Variables =================================
    #
    #
    # ==========================================================================
    model = Model("MIP")
    model.params.DualReductions = 0
    model.params.MIPGap=0.0001

    theta = model.addVar(vtype=GRB.CONTINUOUS, lb=-10000, name="theta")
    
    # Pure strategy probability
    q = [[ model.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="q_s{0}_sn{1}".format(i,i2)) for i in range(Q)] for i2 in range(QN)]
  
    # binary variables for pure strateies
    OS = [[[[model.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="Op(s%d,sn%d,w%d,r%d)" %(i,i2,w,r)) for r in range(R)] for w in range(W)] for i in range(Q)] for i2 in range(QN)]
    test = [[[model.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="ONP(s%d,w%d,r%d)" %(i,w,r)) for r in range(R)] for w in range(W)] for i in range(Q)]

   
    Nplus = [[[[[model.addVar(lb=0.0, ub = 1.0, vtype=GRB.BINARY, name="Np(s%d,sn%d,w%d,k%d,t%d)" %(i,i2,w,k,t))  for k in range(K)]  for t in range(T)] for w in range(W)] for i in range(Q)] for i2 in range(QN)]

    # linearization variables

    X = [[[[[model.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="X(s%d,sn%d,w%d,k%d,t%d)" %(i,i2, w,k,t))  for k in range(K)] for t in range(T)] for w in range(W)] for i in range(Q)] for i2 in range(QN)]
 
    # z[w][k][m]   
    z = [[[ model.addVar(vtype=GRB.CONTINUOUS, name="z_w{0}_k{1}_m{2}".format(w, k, m)) for m in range(M)] for k in range(K)] for w in range(W)]
    
    # pi[w][t][k]
    pi = [[[model.addVar(vtype=GRB.CONTINUOUS, name="pi_w{0}_t{1}_k{2}".format(w, t, k)) for k in range(K)] for t in range(T)] for w in range(W)]
    
    vtype=GRB.CONTINUOUS
    if (integer == 2) or (integer == 3):
        vtype=GRB.INTEGER
        
    n_wtk = [[[model.addVar(vtype=vtype, name="n_w{0}_t{1}_k{2}".format(w, t, k)) for k in range(K)] for t in range(T)] for w in range(W)] 
    # n_wtk[w][t][k] # integer value of N_wk[w][k] * pi[w][t][k]
    

    overflow = [[model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="o_w{0}_r{1}".format(w, r)) for r in range(R)] for w in range(W)] 
    # overflow[w][r]
    

    y = y_marginal#[[ model.addVar(vtype=GRB.INTEGER, lb=0, name="y_w{0}_r{1}".format(w, r)) for r in range(R)] for w in range(W) ]
    # y[w][r]: number of operating resources r at time w
   

    p = p_marginal #[ model.addVar(vtype=GRB.INTEGER, name="p_w{0}".format(w)) for w in range(W)] # available staff
    s = s_marginal #[model.addVar(vtype=GRB.INTEGER, lb=np.floor(s_marginal[w]), name="s_w{0}".format(w)) for w in range(W)] # working staff
        
    model.update()
    print "Variables Added"
    # ========================= Gurobi Objective ===============================
    #
    #
    # ==========================================================================
    
    objective_variables = [theta] + [overflow[w][r] for w in range(W) for r in range(R)]
    objective_coefficients = [1] + [-phi[r] for r in range(R)]*W
    objective_value = LinExpr(objective_coefficients, objective_variables)

    
    model.setObjective(objective_value, GRB.MAXIMIZE)
    print "Objective Set"
    # ======================= Gurobi Constraints ===============================
    #
    #
    # ==========================================================================
    
  
    def Overflow(i, i2,  w, r):
        #return quicksum(Oplus[i][w][r][i] for t in range(len(resource2team[r]))) + quicksum(-Ominus[i][w][r][i] for t in range(len(resource2team[r]))) +np.floor(O_marginal[w][r])
        return OS[i2][i][w][r]

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

    # Throughput constraints
    
    #pre_overflow = np.random.randint(0, 100, R)
    pre_overflow = [0] * R
    for i2 in range(QN):
       for i in range(Q):
            for w in range(W):
                for r in range(R):
                    nmarginal = quicksum(np.floor(N_marginal[i][w][t][k]) for t in resource2team[r] for k in range(K))
                    nbinary = [ Nplus[i2][i][w][t][k] for t in resource2team[r] for k in range(K)]
                    
                    tmp_sum = LinExpr([1 for k in range(K)]*len(resource2team[r]), nbinary)
                    
                    overflow_w = Overflow(i, i2, w, r)    
                    if w == 0:
                        model.addConstr(tmp_sum + nmarginal + pre_overflow[r] - y[i][w][r] * C[r] - overflow_w <= 0, name="(4)_w{0}_r{1}_s{2}{3}".format(w, r,i,i2))
                    else:
                        overflow_w2 = Overflow(i,i2, w-1, r)
                        model.addConstr(tmp_sum + nmarginal + overflow_w2 - y[i][w][r] * C[r] - overflow_w <= 0, name="(4)_w{0}_r{1}_s{2}{3}".format(w, r,i,i2))      
    if OverConstr :
        for i in range(Q):
            for r in range(R): # OPTIONAL
                overflow_w = Overflow(i, i2, W-1, r)
               # model.addConstr(overflow_w == 0, name="(5)_r{0}_s{1}".format(r,i))
        
        for i in range(Q):
            for w in range(W):
                for r in range(R):
                    overflow_w = Overflow(i, i2, w, r)
                    if w > 0:
                        model.addConstr(y[i][w][r] * C[r] - overflow_w >= 0, name="(5.5)_w{0}_r{1}_s{2}".format(w, r,i))

    if False: 
        for w in range(W):
            #tmp_sum = LinExpr(ar, [y[w][r] for r in range(R)])
            tmp_sum = quicksum( ar[r]*y[i][w][r] for r in range(R) )
            model.addConstr(tmp_sum - p[w] <= 0, name="(6)_w{0}".format(w))
    
        for w in range(W):
            for r in range(R):
                model.addConstr(y[i][w][r] - mr[r] <= 0, name="(7)_w{0}_r{1}".format(w, r))
    
        for w in range(W):
            start_index = max(0, w - shift + 1)
            #tmp_sum = LinExpr([1]*(w - start_index + 1), [s[i] for i in range(start_index, w+1)])
            tmp_sum = quicksum(s[i] for i in range(start_index, w+1))
            model.addConstr(tmp_sum - p[w] == 0, name="(8)_w{0}".format(w))
    
        #tmp_sum = LinExpr([1]*W, [s[w] for w in range(W)])
        tmp_sum = quicksum(s[w] for w in range(W))
        model.addConstr(tmp_sum - P <= 0, name="(9)")


    tmp_sum = LinExpr([1]*Q, [q[i2][i] for i in range(Q) for i2 in range(QN)])
    model.addConstr(tmp_sum == 1, name="sumQ")
    
    
    # Linearization Constraints
    for i in range(Q):
        for i2 in range(QN):
            for w in range(W):
                for k in range(K):
                    model.addConstr( quicksum(np.floor(N_marginal[i][w][t][k]) + Nplus[i2][i][w][t][k] for t in range(T)) == N_wk[w][k])  
                    for t in range(T):
                        model.addConstr(X[i2][i][w][t][k] <= q[i2][i]) 
                        model.addConstr(X[i2][i][w][t][k] <= Nplus[i2][i][w][t][k] ) 
                        model.addConstr(X[i2][i][w][t][k] >= q[i2][i] -(1-Nplus[i2][i][w][t][k])) 

            """
            for r in range(R):
                for t in range(len(resource2team[r])):
                    model.addConstr(Up[i][w][r][t] <= q[i]) 
                    model.addConstr(Up[i][w][r][t] <= Oplus[i][w][r][t]) 
                    model.addConstr(Up[i][w][r][t] >= q[i] -(1-Oplus[i][w][r][t])) 
                    
                    model.addConstr(Um[i][w][r][t] <= -q[i]) 
                    model.addConstr(Um[i][w][r][t] <= -Ominus[i][w][r][t]) 
                    model.addConstr(Um[i][w][r][t] >= -q[i] -(1+Ominus[i][w][r][t])) 
                    """
    
    # Mixed Strategy Constraints
    for w in range(W):
        for k in range(K):
            for t in range(T):
                model.addConstr(n_wtk[w][t][k] == quicksum(quicksum(q[i2][i]*np.floor(N_marginal[i][w][t][k]) for i in range(Q)) + quicksum(X[i2][i][w][t][k] for i in range(Q)) for i2 in range(QN)))
        for r in range(R):
            nmarginal = quicksum(np.floor(N_marginal[i][w][t][k])*q[i2][i] for t in resource2team[r] for k in range(K) for i in range(Q) for i2 in range(QN))
            nX = [ X[i2][i][w][t][k] for t in resource2team[r] for k in range(K) for i in range(Q) for i2 in range(QN)]
                
            tmp_sum = LinExpr([1 for k in range(K) for i in range(Q)]*len(resource2team[r]), nX)
 
            model.addConstr( tmp_sum + nmarginal - quicksum(q[i2][i]*y[i][w][r] for i in range(Q) for i2 in range(QN))*10000 <=0 , name="(5.6)_w{0}_r{1}".format(w, r))    
            if w == 0:
                model.addConstr(tmp_sum + nmarginal - quicksum(q[i2][i]*y[i][w][r] for i in range(Q) for i2 in range(QN)) * C[r] - overflow[w][r] <= 0)
            else:
                model.addConstr(tmp_sum + nmarginal + overflow[w-1][r] - quicksum(q[i2][i]*y[i][w][r] for i in range(Q) for i2 in range(QN)) * C[r] - overflow[w][r] <= 0)
            #model.addConstr(overflow[w][r] == quicksum(Overflow(i,w,r) for i in range(Q))) 
            
    
    model.update()
    
    print "Constraints Added"  
    model.write("tsg_MIP.lp")
    start_time = time.time()
    model.optimize()
    runtime = time.time() - start_time
    model.write("tsg_MIP.sol")


    #======================GET SOLUTION===========================#
    
    n_value = np.zeros((W,T,K))
    for w in range(W):
        for t in range(T):
            for k in range(K):
                n_value[w][t][k] = n_wtk[w][t][k].x
                

    overflow_value = np.zeros((W,R))
    for w in range(W):
        for r in range(R):
            overflow_value[w][r] = overflow[w][r].x

    """
    y_value = np.zeros((W,R))
    for w in range(W):
        for r in range(R):
            y_value[w][r] = y[w][r].x

    s_value = np.zeros(W)
    for w in range(W):
        s_value[w] = s[w].x
    """
    n_binary = np.zeros((Q,QN,W,T,K))
    X_sol = np.zeros((Q,QN,W,T,K))

    mixed_strategy = np.zeros((Q,QN))
    for i in range(Q):
        for i2 in range(QN):
            mixed_strategy[i][i2] = q[i2][i].x
    
    for w in range(W):
        for t in range(T):
            for k in range(K):
                sum_q = 0
                for i in range(Q):
                    for i2 in range(QN):

                        X_sol [i][i2][w][t][k] = X[i2][i][w][t][k].x
    
    attack_set = 0
    AV = theta.x
    for w in range(W):
        for m in range(M):
            for k in range(K):
                att_val = z[w][k][m].x*(U_plus[k] - U_minus[k]) - U_minus[k]
                diff = np.abs(AV - att_val)
                if(diff>0.0001): attack_set+=1
    
    
    obj = model.getAttr('ObjVal')
    
    return obj, runtime, mixed_strategy, n_value, overflow_value, attack_set #, y_value, s_value

def testSetting(): 
    K=1
    R=2
    M=2
    teams = [[0],[1]] 
    T = len(teams)
    resource2team = util.resourceTeamDict(R, T, teams)
    print resource2team

    Er = [[1,0],[0,1]]
    print "Er"
    print Er

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

    N_wk = [[1],[1]]
   
    C = [1,1]
    print "\nC"
    print C
    #C = [100,80,75,50,30,15]
    mr = np.random.randint(10, 30, R)
    print "\nmr"
    print mr
    #mr = [5, 5, 5, 3, 2, 4] # maximum number of people to operate resource r
    ar = [1,1]
    print "\nar"
    print ar
    #ar = [2, 1, 2, 1, 1, 2] # ar[r] number of people required to operate resource r

    phi = np.random.rand(R)/10 # phi[r] overflow penalty
    print "\nphi"
    print phi
    #phi = np.random.rand(W, R) # phi[w][r] overflow penalty

    return resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi  

def randomSetting(seed, W, K ,R, mR, M, P, teams, shift):
    # ========================== Random Seed ===================================
    #seed = 3234
    #seed = random.randint(1,10000)
    np.random.seed(seed)
    random.seed(seed)

    T = len(teams)
    resource2team = util.resourceTeamDict(R, T, teams)
    print resource2team

    Er = np.random.rand(R, M)/2 + 0.5 # Er[m][r]
    Er = Er / 2
    print "Er"
    print Er
    #Er = [[0.3, 0.5, 0.2], [0.6, 0.3, 0.4], [0.4, 0.6, 0.5]
    #     ,[0.6, 0.3, 0.8], [0.7, 0.4, 0.7], [0.7, 0.6, 0.9]]

    Er, C = util.genResources(R, M, 100)
    E = util.computeTeamsRate(R, M, T, teams, Er)
    print E     

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
    print "\nutilities"
    print "plus"
    print U_plus
    print "minus"
    print U_minus

    N_wk = [[ 0 for k in range(K)] for w in range(W)] # N_wk[w][k] represents the number of people getting in time window w with type k
    
    arrival_window=3

    for k in range(K):
        flight_time = random.randint(0, W-arrival_window)+arrival_window
        for w in range(flight_time-arrival_window, flight_time):
            large_or_small = random.random()
            if large_or_small > 0.5:
                tmp_N = random.randint(100, 300)
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
    print "number of passengers: {0}".format([sum(N_wk[w]) for w in range(W)])

    #C = np.random.randint(200, 500, R) # C[r] is the capacity of resource r
    print "\nC"
    print C
    #C = [100,80,75,50,30,15]
    mr = np.random.randint(10, 30, R)
    print "\nmr"
    print mr
    #mr = [5, 5, 5, 3, 2, 4] # maximum number of people to operate resource r
    ar = np.random.randint(1, 5, R)
    print "\nar"
    print ar
    #ar = [2, 1, 2, 1, 1, 2] # ar[r] number of people required to operate resource r

    phi = np.random.rand(R)/10 # phi[r] overflow penalty
    print "\nphi"
    print phi
    #phi = np.random.rand(W, R) # phi[w][r] overflow penalty

    return resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi

def bigSetting(W, K, mR, P, teams, shift):
    # ========================== Random Seed ===================================
    #seed = 3234
    #seed = random.randint(1,10000)

    T = len(teams)
    resource2team = util.resourceTeamDict(R, T, teams)
    print resource2team

    print "Er"
    
    Er = [[0.5, 0.5], [0.6, 0.4], [0.4, 0.6]
         ,[0.3, 0.7], [0.7, 0.3]]
    print Er
    C = [100,50,50,10,10]
    #Er, C = util.genResources(R, M, 100)
    E = util.computeTeamsRate(R, M, T, teams, Er)
    print E     

    # RatesEr = [0.357, 0.916, 0.511, 0.916, 1.204, 1.204,
    #     0.693, 0.357, 0.916, 0.357, 0.511, 0.916,
    #     0.223, 0.511, 0.693, 1.609, 1.204, 2.303]

    # suppose this is a zero-sum game
    U_plus = [0 for k in range(K)] # covered (plus) utility of the defender
    U_minus = [-1000 for k in range(K)] # uncovered (minus) utility of the defender
    
    print "\nutilities"
    print "plus"
    print U_plus
    print "minus"
    print U_minus

    N_wk = [[ 0 for k in range(K)] for w in range(W)] # N_wk[w][k] represents the number of people getting in time window w with type k
    
    arrival_window=3

    for k in range(K):
        flight_time = random.randint(0, W-arrival_window)+arrival_window
        for w in range(flight_time-arrival_window, flight_time):
            N_wk[w][k] = 60
   
    print "number of passengers: {0}".format([sum(N_wk[w]) for w in range(W)])

    #C = np.random.randint(200, 500, R) # C[r] is the capacity of resource r
    print "\nC"
    #print C
    
    #mr = np.random.randint(10, 30, R)
    print "\nmr"
    #print mr
    mr = [5, 5, 5, 3, 2, 4] # maximum number of people to operate resource r
    #ar = np.random.randint(1, 5, R)
    print "\nar"
    #print ar
    ar = [2, 1, 1, 1, 1] # ar[r] number of people required to operate resource r

    phi = [0.6 for r in range(5)] # phi[r] overflow penalty
    print "\nphi"
    print phi
    #phi = np.random.rand(W, R) # phi[w][r] overflow penalty

    return resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi    
if __name__ == "__main__":
    # ============================= main =======================================
    print "======================== main ======================================"
    # ========================= Game Setting ===================================
    W = 10 # number of time windows
    K = 5 # number of passenger types
    R = 5 # number of resources
    mR = 5 # max number of reosurces
    M = 2 # number of attack methods
    P = 3 # number of staff
    Q=5
    shift = 3 # d

    nT = 20
    teams = util.generateAllTeams(R, mR)
    #teams = util.randomGenerateTeams(R, mR, nT)
    #teams = [[0],[1]]
    print teams

    # ================= random generate game setting ===========================
    seed = 2345
    #resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
    #resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = testSetting()
    #resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = bigSetting( W, K , mR,  P, teams, shift)
    
    #lst = [W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi]
    #pickle.dump( lst, open( "data.p", "wb" ) )
    [W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi] = pickle.load( open( "data.p", "rb" ) )
    print "============================ LP relaxation =============================="
    obj_relaxed, n, o, y, s, p, att_set1, f = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=0, OverConstr=False)
    
    start_time = time.time()
    print "============================ MIP SOLVER =============================="

    obj, rt, q, n2, o, att_set = Ksolver(W, K, R, mR, M, P, Q, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, n, o, y, s,p,integer=0,OverConstr=False)
    walltime = time.time() - start_time
    print "relaxation objective: ", obj_relaxed
    print "MIP objective: ", obj
    print "Runtime/walltime ", rt, " ", walltime
    print "attack set LP", att_set1," ",f, " MIP ", att_set


