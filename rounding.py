from gurobipy import *
import util
import relaxed
import numpy as np
import random

def randomizedRounding(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=0):
    n_value, overflow_value, y_value, s_value = relaxed.LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, 1)

    print "========================================== randomized rounding ================================================"

    pure_strategy_limit = 10
    large_number = 100

    model = Model("MIP")
    model.params.DualReductions = 0

    theta = model.addVar(vtype=GRB.CONTINUOUS, name="theta")

    q = []
    for i in range(pure_strategy_limit):
        tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="q_i{0}".format(i))
        q.append(tmp_var)

    n_value_floor = np.array(n_value)
    n_value_integer_indicator = np.zeros((W,T,K))
    n_binary = [] # n[w][t][k][i]
    n_effective = [] # n[w][t][k][i]
    n_mixed = [] # n[w][t][k]
    for w in range(W):
        n_binary.append([])
        n_effective.append([])
        n_mixed.append([])
        for t in range(T):
            n_binary[w].append([])
            n_effective[w].append([])
            n_mixed[w].append([])
            for k in range(K):
                n_binary[w][t].append([])
                n_effective[w][t].append([])
                for i in range(pure_strategy_limit):
                    tmp_var = model.addVar(vtype=GRB.BINARY, name="n_binary_w{0}_t{1}_k{2}_i{3}".format(w,t,k,i))
                    n_binary[w][t][k].append(tmp_var)
                    tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="n_effective_w{0}_t{1}_k{2}_i{3}".format(w,t,k,i))
                    n_effective[w][t][k].append(tmp_var)
                tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="n_mixed_w{0}_t{1}_k{2}".format(w,t,k))
                n_mixed[w][t].append(tmp_var)

                n_value_floor[w][t][k] = int(n_value[w][t][k]) # since all of them are positive, so this is fine
                if n_value_floor[w][t][k] == n_value[w][t][k]:
                    n_value_integer_indicator[w][t][k] = 1
                    for i in range(pure_strategy_limit):
                        model.addConstr(n_binary[w][t][k][i] == 0, name="(0.1)_w{0}_t{1}_k{2}_i{3} (n is already integer)".format(w,t,k,i))
                #else:
                #    print n_value[w][t][k]
                    

    overflow_value_floor = np.array(overflow_value)
    overflow_value_integer_indicator = np.zeros((W,R))
    overflow_binary = [] # o[w][r][i]
    overflow_effective = [] # o[w][r][i]
    overflow_mixed = [] # o[w][r]
    for w in range(W):
        overflow_binary.append([])
        overflow_effective.append([])
        overflow_mixed.append([])
        for r in range(R):
            overflow_binary[w].append([])
            overflow_effective[w].append([])
            for i in range(pure_strategy_limit):
                tmp_var = model.addVar(vtype=GRB.BINARY, name="overflow_binary_w{0}_r{1}_i{2}".format(w,r,i))
                overflow_binary[w][r].append(tmp_var)
                tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="overflow_effective_w{0}_r{1}_i{2}".format(w,r,i))
                overflow_effective[w][r].append(tmp_var)
            tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="overflow_mixed_w{0}_r{1}".format(w,r))
            overflow_mixed[w].append(tmp_var)

            overflow_value_floor[w][r] = int(overflow_value[w][r])
            if overflow_value_floor[w][r] == overflow_value[w][r]:
                overflow_value_integer_indicator[w][r] = 1
                for i in range(pure_strategy_limit):
                    model.addConstr(overflow_binary[w][r][i] == 0, name="(0.2)_w{0}_r{1}_i{2} (overflow is already integer)".format(w,r,i))

    y_value_floor = np.array(y_value)
    y_value_integer_indicator = np.zeros((W,R))
    y_binary = [] # y[w][r][i]
    y_effective = [] # y[w][r][i]
    y_mixed = [] # y[w][r]
    for w in range(W):
        y_binary.append([])
        y_effective.append([])
        y_mixed.append([])
        for r in range(R):
            y_binary[w].append([])
            y_effective[w].append([])
            for i in range(pure_strategy_limit):
                tmp_var = model.addVar(vtype=GRB.BINARY, name="y_binary_w{0}_r{1}_i{2}".format(w,r,i))
                y_binary[w][r].append(tmp_var)
                tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="y_effective_w{0}_r{1}_i{2}".format(w,r,i))
                y_effective[w][r].append(tmp_var)
            tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="y_mixed_w{0}_r{1}".format(w,r))
            y_mixed[w].append(tmp_var)

            y_value_floor[w][r] = int(y_value[w][r])
            if y_value_floor[w][r] == y_value[w][r]:
                y_value_integer_indicator[w][r] = 1
                for i in range(pure_strategy_limit):
                    model.addConstr(y_binary[w][r][i] == 0, name="(0.3)_w{0}_r{1}_i{2} (y is already integer)".format(w,r,i))

    s_value_floor = np.array(s_value)
    s_value_integer_indicator = np.zeros(W)
    s_binary = [] # s[r][i]
    s_effective = [] # s[r][i]
    s_mixed = [] # s[r]
    for w in range(W):
        s_binary.append([])
        s_effective.append([])
        for i in range(pure_strategy_limit):
            tmp_var = model.addVar(vtype=GRB.BINARY, name="s_binary_w{0}_i{1}".format(w,i))
            s_binary[w].append(tmp_var)
            tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="s_effective_w{0}_i{1}".format(w,i))
            s_effective[w].append(tmp_var)
        tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="s_mixed_w{0}".format(w))
        s_mixed.append(tmp_var)

        s_value_floor[w] = int(s_value[w])
        if s_value_floor[w] == s_value[w]:
            s_value_integer_indicator[w] = 1
            for i in range(pure_strategy_limit):
                model.addConstr(s_binary[w][i] == 0, name="(0.4)_w{0}_i{1} (s is already integer)".format(w,i))

    p = []
    for w in range(W):
        p.append([])
        for i in range(pure_strategy_limit):
            tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="p_w{0}_i{1}".format(w,i))
            p[w].append(tmp_var)

    # =================================== objective value ===========================================
    model.setObjective(n_mixed[0][0][0], GRB.MAXIMIZE)

    # =============================== gurobi constraints ============================================
    for w in range(W):
        for t in range(T):
            for k in range(K):
                for i in range(pure_strategy_limit):
                    model.addConstr(n_effective[w][t][k][i] <= q[i], name="(n 1)_w{0}_t{1}_k{2}_i{3}".format(w,t,k,i))
                    model.addConstr(n_effective[w][t][k][i] <= n_binary[w][t][k][i] * large_number, name="(n 2)_w{0}_t{1}_k{2}_i{3}".format(w,t,k,i)) # if n_binary == 1, then n_effective can be q[i]; if n_binary = 0, then n_effective = 0.
                tmp_sum = LinExpr([1]*pure_strategy_limit, n_effective[w][t][k])
                model.addConstr(n_mixed[w][t][k] == n_value_floor[w][t][k] + tmp_sum, name="(n 3)_w{0}_t{1}_k{2}".format(w,t,k))

    for w in range(W):
        for r in range(R):
            for i in range(pure_strategy_limit):
                model.addConstr(overflow_effective[w][r][i] <= q[i], name="(o 1)_w{0}_r{1}_i{2}".format(w,r,i))
                model.addConstr(overflow_effective[w][r][i] <= overflow_binary[w][r][i] * large_number, name="(o 2)_w{0}_r{1}_i{2}".format(w,r,i))
            tmp_sum = LinExpr([1]*pure_strategy_limit, overflow_effective[w][r])
            model.addConstr(overflow_mixed[w][r] == overflow_value_floor[w][r] + tmp_sum, name="(o 3)_w{0}_r{1}".format(w,r))

    for w in range(W):
        for r in range(R):
            for i in range(pure_strategy_limit):
                model.addConstr(y_effective[w][r][i] <= q[i], name="(y 1)_w{0}_r{1}_i{2}".format(w,r,i))
                model.addConstr(y_effective[w][r][i] <= y_binary[w][r][i] * large_number, name="(y 2)_w{0}_r{1}_i{2}".format(w,r,i))
            tmp_sum = LinExpr([1]*pure_strategy_limit, y_effective[w][r])
            model.addConstr(y_mixed[w][r] == y_value_floor[w][r] + tmp_sum, name="(y 3)_w{0}_r{1}".format(w,r))

    for w in range(W):
        for i in range(pure_strategy_limit):
            model.addConstr(s_effective[w][i] <= q[i], name="(s 1)_w{0}_i{1}".format(w,i))
            model.addConstr(s_effective[w][i] <= s_binary[w][r] * large_number, name="(s 2)_w{0}_i{1}".format(w,i))
        tmp_sum = LinExpr([1]*pure_strategy_limit, s_effective[w])
        model.addConstr(s_mixed[w] == s_value_floor[w] + tmp_sum, name="(s 3)_w{0}".format(w))

    # ============================== pure strategy constraints ======================================
    for w in range(W):
        for k in range(K):
            for i in range(pure_strategy_limit):
                tmp_sum = LinExpr([1]*T, [n_binary[w][t][k][i] for t in range(T)])
                tmp_sum2 = sum([n_value_floor[w][t][k] for t in range(T)])
                model.addConstr(tmp_sum + tmp_sum2 == N_wk[w][k], name="(N_wk)_w{0}_k{1}_i{2}".format(w,k,i))
    for w in range(W):
        for r in range(R):
            for i in range(pure_strategy_limit):
                tmp_sum = LinExpr([1 for t in resource2team[r] for k in range(K)], [n_binary[w][t][k][i] for t in resource2team[r] for k in range(K)])
                tmp_sum2 = sum([n_value_floor[w][t][k] for t in resource2team[r] for k in range(K)])
                if w == 0:
                    model.addConstr(tmp_sum + tmp_sum2 <= (y_value_floor[w][r] + y_binary[w][r][i]) * C[r] - 0 + (overflow_value_floor[w][r] + overflow_binary[w][r][i]))
                else:
                    model.addConstr(tmp_sum + tmp_sum2 <= (y_value_floor[w][r] + y_binary[w][r][i]) * C[r] - (overflow_value_floor[w-1][r] + overflow_binary[w-1][r][i]) + (overflow_value_floor[w][r] + overflow_binary[w][r][i]))

    for w in range(W):
        for r in range(R):
            for i in range(pure_strategy_limit):
                if w >= 1:
                    model.addConstr((y_value_floor[w][r] + y_binary[w][r][i]) * C[r] >= (overflow_value_floor[w-1][r] + overflow_binary[w-1][r][i]))
    for w in range(W):
        for r in range(R):
            for i in range(pure_strategy_limit):
                model.addConstr((overflow_value_floor[w][r] + overflow_binary[w][r][i]) == 0)

    for w in range(W):
        for i in range(pure_strategy_limit):
            tmp_sum = LinExpr([ar[r] for r in range(R)], [y_binary[w][r][i] for r in range(R)])
            tmp_sum2 = sum([y_value_floor[w][r] * ar[r] for r in range(R)])
            model.addConstr(tmp_sum + tmp_sum2 <= p[w][i])

    for w in range(W):
        for r in range(R):
            for i in range(pure_strategy_limit):
                model.addConstr(y_value_floor[w][r] + y_binary[w][r][i] <= mr[r])

    for w in range(W):
        for i in range(pure_strategy_limit):
            start_index = max(0, w - shift + 1)
            tmp_sum = LinExpr([1]*(w - start_index + 1), [s_binary[w_0][i] for w_0 in range(start_index, w+1)])
            tmp_sum2 = sum([s_value_floor[w_0] for w_0 in range(start_index, w+1)])
            model.addConstr(tmp_sum + tmp_sum2 - p[w][i] == 0, name="(8)_w{0}_i{1}".format(w,i))

    for i in range(pure_strategy_limit):
        tmp_sum = LinExpr([1]*W, [s_binary[w][i] for w in range(W)])
        tmp_sum2 = sum([s_value_floor[w] for w in range(W)])
        model.addConstr(tmp_sum + tmp_sum2 <= P)


    model.write("lp/k_tsg.lp")

    model.optimize()

    print "TODO"

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
    seed = random.randint(1,10000)
    resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = relaxed.randomSetting(seed, W, K ,R, mR, M, P, teams, shift)

    #print "============================ LP relaxation =============================="
    #relaxed.LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=0)

    print "============================ rounding method ============================"
    randomizedRounding(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=0)

