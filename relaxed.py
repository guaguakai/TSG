import util
from gurobipy import *
import numpy as np
import random

def LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=True): # integer indicates different relaxation method
    # ======================= Gurobi Setting ===================================
    model = Model("MIP")
    model.params.DualReductions = 0

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
            if (integer == 1) or (integer == 3):
                tmp_resource_r = model.addVar(vtype=GRB.INTEGER, lb=0, name="y_w{0}_r{1}".format(w, r))
            else:
                tmp_resource_r = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="y_w{0}_r{1}".format(w, r))
            y[w].append(tmp_resource_r)

    p = [] # available staff
    s = [] # working staff
    for w in range(W):
        tmp_staff = model.addVar(vtype=GRB.CONTINUOUS, name="p_w{0}".format(w))
        tmp_working_staff = model.addVar(vtype=GRB.INTEGER, name="s_w{0}".format(w))
        p.append(tmp_staff)
        s.append(tmp_working_staff)

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
            #tmp_sum = LinExpr([N_wk[w][k] for k in range(K)]*len(resource2team[r]), [pi[w][t][k] for t in resource2team[r] for k in range(K)])
            tmp_sum = LinExpr([1 for k in range(K)]*len(resource2team[r]), [n_wtk[w][t][k] for t in resource2team[r] for k in range(K)])
            if w == 0:
                model.addConstr(tmp_sum + pre_overflow[r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(4)_w{0}_r{1}".format(w, r))
            else:
                model.addConstr(tmp_sum + overflow[w-1][r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(4)_w{0}_r{1}".format(w, r))

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

    model.write("tsg.mps")
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
    return model
    

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

    N_wk = [] # N_wk[w][k] represents the number of people getting in time window w with type k
    for w in range(W):
        N_wk.append([])
        for k in range(K):
            large_or_small = random.random()
            if large_or_small > 0.5:
                tmp_N = random.randint(1000, 3000)
            else:
                tmp_N = random.randint(100, 300)
            N_wk[w].append(tmp_N)

    print "number of passengers: {0}".format([sum(N_wk[w]) for w in range(W)])

    C = np.random.randint(200, 500, R) # C[r] is the capacity of resource r
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

if __name__ == "__main__":
    # ============================= main =======================================
    print "======================== main ======================================"
    # ========================= Game Setting ===================================
    W = 15 # number of time windows
    K = 20 # number of passenger types
    R = 10 # number of resources
    mR = 3 # max number of reosurces
    M = 10 # number of attack methods
    P = 100 # number of staff
    shift = 8 # d

    nT = 20
    teams = util.generateAllTeams(R, mR)
    #teams = util.randomGenerateTeams(R, mR, nT)

    print teams

    # ================= random generate game setting ===========================
    seed = 2345
    resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)

    print "============================ LP relaxation =============================="
    model = LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=0)
    #print "============================ relaxed n_wtk (allocated arrivals) MIP ==============================="
    #LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=1)
    #print "============================ relaxed y (number of resources) MIP ==============================="
    #LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=2)
    #print "============================ full MIP ==============================="
    #LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=3)


