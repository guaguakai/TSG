import util
from gurobipy import *
import numpy as np
import random
import relaxed

# =========================================== column generation =================================================

def slaveProblem(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, gamma): # slave problem solver
    model = Model("MIP")
    model.params.DualReductions = 0

    n_wtk = [] # n_wtk[w][t][k] # integer value of N_wk[w][k] * pi[w][t][k]
    for w in range(W):
        n_wtk.append([])
        for t in range(T):
            n_wtk[w].append([])
            for k in range(K):
                tmp_pi_var = model.addVar(vtype=GRB.INTEGER, name="n_w{0}_t{1}_k{2}".format(w, t, k))
                n_wtk[w][t].append(tmp_pi_var)

    overflow = [] # overflow[w][r]
    for w in range(W):
        overflow.append([])
        for r in range(R):
            tmp_overflow_var = model.addVar(vtype=GRB.INTEGER, lb=0, name="o_w{0}_r{1}".format(w, r))
            overflow[w].append(tmp_overflow_var)

    y = [] # y[w][r]: number of operating resources r at time w
    for w in range(W):
        y.append([])
        for r in range(R):
            tmp_resource_r = model.addVar(vtype=GRB.INTEGER, lb=0, name="y_w{0}_r{1}".format(w, r))
            y[w].append(tmp_resource_r)

    p = [] # available staff
    s = [] # working staff
    for w in range(W):
        tmp_staff = model.addVar(vtype=GRB.CONTINUOUS, name="p_w{0}".format(w)) # continuous is fine here since p = sum of s
        tmp_working_staff = model.addVar(vtype=GRB.INTEGER, name="s_w{0}".format(w))
        p.append(tmp_staff)
        s.append(tmp_working_staff)

    # ========================= Gurobi Objective ===============================
    print "objective setting"
    objective_variables = [n_wtk[w][t][k] for w in range(W) for t in range(T) for k in range(K)] + [overflow[w][r] for w in range(W) for r in range(R)]
    objective_coefficients = [gamma[w][t][k] for w in range(W) for t in range(T) for k in range(K)] + [-phi[r] for r in range(R)] * W
    print len(objective_variables), len(objective_coefficients)
    objective_value = LinExpr(objective_coefficients, objective_variables)

    #objective_value = theta
    #for w in range(W):
    #    for r in range(R):
    #        objective_value += phi[r]*overflow[w][r]
    model.setObjective(objective_value, GRB.MAXIMIZE)

    # ======================= Gurobi Constraints ===============================
    for w in range(W):
        for k in range(K):
            tmp_sum = LinExpr([1]*T, [n_wtk[w][t][k] for t in range(T)])
            model.addConstr(tmp_sum == N_wk[w][k], name="(1)_w{0}_k{1}".format(w, k))

    pre_overflow = [0] * R
    for w in range(W):
        for r in range(R):
            #tmp_sum = LinExpr([N_wk[w][k] for k in range(K)]*len(resource2team[r]), [pi[w][t][k] for t in resource2team[r] for k in range(K)])
            tmp_sum = LinExpr([1 for k in range(K)]*len(resource2team[r]), [n_wtk[w][t][k] for t in resource2team[r] for k in range(K)])
            if w == 0:
                model.addConstr(tmp_sum + pre_overflow[r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(2)_w{0}_r{1}".format(w, r))
            else:
                model.addConstr(tmp_sum + overflow[w-1][r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(2)_w{0}_r{1}".format(w, r))

    for w in range(W):
        for r in range(R):
            if w > 0:
                model.addConstr(y[w][r] * C[r] - overflow[w-1][r] >= 0, name="(3)_w{0}_r{1}".format(w, r))

    for r in range(R): # OPTIONAL
        model.addConstr(overflow[W-1][r] == 0, name="(4)_r{0}".format(r))

    for w in range(W):
        tmp_sum = LinExpr(ar, [y[w][r] for r in range(R)])
        model.addConstr(tmp_sum - p[w] <= 0, name="(5)_w{0}".format(w))

    for w in range(W):
        for r in range(R):
            model.addConstr(y[w][r] - mr[r] <= 0, name="(6)_w{0}_r{1}".format(w, r))

    for w in range(W):
        start_index = max(0, w - shift + 1)
        tmp_sum = LinExpr([1]*(w - start_index + 1), [s[i] for i in range(start_index, w+1)])
        model.addConstr(tmp_sum - p[w] == 0, name="(7)_w{0}".format(w))

    tmp_sum = LinExpr([1]*W, [s[w] for w in range(W)])
    model.addConstr(tmp_sum - P <= 0, name="(8)")

    model.write("lp/slave.lp")

    model.optimize()

    # ====================== retrieve the values ===============================
    n_value = np.zeros((W, T, K))
    for w in range(W):
        for t in range(T):
            for k in range(K):
                n_value[w][t][k] = model.getVarByName("n_w{0}_t{1}_k{2}".format(w, t, k)).x
                #print "n_w{0}_t{1}_k{2}: {3}".format(w, t, k, n_value[w][t][k])

    overflow_value = np.zeros((W, R))
    for w in range(W):
        for r in range(R):
            overflow_value[w][r] = model.getVarByName("o_w{0}_r{1}".format(w, r)).x
            #print "o_w{0}_r{1}: {2}".format(w, r, overflow_value[w][r])

    return n_value, overflow_value


def columnGenerationSolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, Q, iteration=10): # integer indicates different relaxation method
    model = Model("LP")
    model.params.DualReductions = 0

    theta = model.addVar(vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, name="theta")

    q = []
    for i in range(len(Q)):
        tmp_q = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="q_{0}".format(i))
        q.append(tmp_q)

    z = [] # z[w][k][m]
    for w in range(W):
        z.append([])
        for k in range(K):
            z[w].append([])
            for m in range(M):
                tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="z_w{0}_k{1}_m{2}".format(w, k, m))
                z[w][k].append(tmp_var)

    tilde_n = [] # tilde_n[w][t][k]
    for w in range(W):
        tilde_n.append([])
        for t in range(T):
            tilde_n[w].append([])
            for k in range(K):
                tmp_var = model.addVar(vtype=GRB.CONTINUOUS, name="tilde_n_w{0}_t{1}_k{2}".format(w, t, k))
                tilde_n[w][t].append(tmp_var)

    # ========================= Gurobi Objective ===============================
    objective_variables = [theta] + [q[i] for i in range(len(Q))]
    objective_coefficients = [1] + [-sum([Q[i]["overflow"][w][r] * phi[r] for w in range(W) for r in range(R)]) for i in range(len(Q))]
    objective_value = LinExpr(objective_coefficients, objective_variables)

    #objective_value = theta
    #for w in range(W):
    #    for r in range(R):
    #        objective_value += phi[r]*overflow[w][r]
    model.setObjective(objective_value, GRB.MAXIMIZE)

    # ======================= Gurobi Constraints ===============================
    c_alpha = []
    for w in range(W):
        c_alpha.append([])
        for k in range(K):
            c_alpha[w].append([])
            for m in range(M):
                tmp_c = model.addConstr(theta - z[w][k][m]*(U_plus[k] - U_minus[k]) - U_minus[k] <= 0, name="(1)_w{0}_k{1}_m{2}".format(w,k,m))
                c_alpha[w][k].append(tmp_c)

    for w in range(W):
        for k in range(K):
            for m in range(M):
                tmp_sum = LinExpr([E[t][m] for t in range(T)], [tilde_n[w][t][k] for t in range(T)])
                model.addConstr(z[w][k][m]*N_wk[w][k] == tmp_sum, name="(2)_w{0}_k{1}_m{2}".format(w, k, m))

    c_gamma = []
    for w in range(W):
        c_gamma.append([])
        for t in range(T):
            c_gamma[w].append([])
            for k in range(K):
                tmp_sum = LinExpr([Q[i]['n'][w][t][k] for i in range(len(Q))], [q[i] for i in range(len(Q))] )
                tmp_c = model.addConstr(tilde_n[w][t][k] == tmp_sum, name="(3)_w{0}_t{1}_k{2}".format(w, t, k))
                c_gamma[w][t].append(tmp_c)

    tmp_sum = LinExpr([1 for i in range(len(Q))], [q[i] for i in range(len(Q))])
    c_q_sum = model.addConstr(tmp_sum == 1, name="(4)")

    model.write("lp/CG.lp")

    model.optimize()

    objective_value = 0
    alpha_value = np.zeros((W, K, M))
    for w in range(W):
        for k in range(K):
            for m in range(M):
                alpha_value[w][k][m] = c_alpha[w][k][m].Pi
                objective_value += alpha_value[w][k][m] * U_minus[k]

    gamma_value = np.zeros((W, T, K))
    for w in range(W):
        for t in range(T):
            for k in range(K):
                gamma_value[w][t][k] = c_gamma[w][t][k].Pi
                #if gamma_value[w][t][k] != 0:
                #    print gamma_value[w][t][k]

    delta_value = c_q_sum.Pi
    #print "delta value: {0}".format(delta_value)
    #print "dual objective value: {0}".format(objective_value)

    q_value = np.zeros(len(Q))
    for i in range(len(Q)):
        #print q[i].varName, q[i].x
        q_value[i] = q[i].x

    return model, gamma_value, delta_value, q_value


if __name__ == "__main__":
    # ============================= main =======================================
    print "======================== main ======================================"
    # ========================= Game Setting ===================================
    W = 5 # number of time windows
    K = 10 # number of passenger types
    R = 6 # number of resources
    mR = 3 # max number of reosurces
    M = 2 # number of attack methods
    P = 30 # number of staff
    shift = 3 # d
    Q = 4 # here no use
    nT = 25
    teams = util.generateAllTeams(R, mR)
    maxT = 5 # here no use either
    #teams = util.randomGenerateTeams(R, mR, nT)

    print teams

    # ================= random generate game setting ===========================
    seed = 2345

    resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = relaxed.randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
    total_arrivals = 0
    for w in range(W):
        for k in range(K):
            total_arrivals += N_wk[w][k]
    print "total arrivals: {0}".format(total_arrivals)

    Q = []
    gamma = np.ones((W, T, K))

    for j in range(1000):
        print "==================================== solving slave problem ===================================="
        slave_n_value, slave_overflow_value = slaveProblem(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, gamma)
        tmpQ = {}
        tmpQ["n"] = slave_n_value
        tmpQ["overflow"] = slave_overflow_value

        Q.append(tmpQ)
        #print "==================================== q sanity check ==========================================="
        #for i in range(len(Q)):
        #    tmp_sum = sum([gamma[w][t][k] * Q[i]['n'][w][t][k] for w in range(W) for t in range(T) for k in range(K)]) - sum([phi[r] * Q[i]["overflow"][w][r] for w in range(W) for r in range(R)])
        #    print "{0} sanity check: {1}".format(i, tmp_sum)

        print "================================== column generation testing =================================="
        cg_model, gamma, delta_value, q_value = columnGenerationSolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, Q, iteration=10)
        #new_Q = [Q[i] for i in range(len(Q)) if q_value[i] > 0]
        #Q = new_Q

        #print "delta value: {0}".format(delta_value)


    print "============================ LP relaxation =============================="
    model = relaxed.LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, integer=0)


