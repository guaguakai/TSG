import util
from gurobipy import *
import numpy as np
import random
import relaxed
import DesignYNcombined
import StaffResourceAllocation
import time
import itertools

# =========================================== column generation =================================================

def slaveProblem(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, gamma, fix_s=None): # slave problem solver
    model = Model("MIP")
    model.params.OutputFlag = 0
    model.params.TuneOutput = 0

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
    #print "objective setting"
    objective_variables = [n_wtk[w][t][k] for w in range(W) for t in range(T) for k in range(K)] + [overflow[w][r] for w in range(W) for r in range(R)]
    objective_coefficients = [gamma[w][t][k] for w in range(W) for t in range(T) for k in range(K)] + [-phi[r] for r in range(R)] * W
    #print len(objective_variables), len(objective_coefficients)
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
            model.addConstr(y[w][r]*10000 >= tmp_sum, name="(5.6)_w{0}_r{1}".format(w, r))
            if w == 0:
                model.addConstr(tmp_sum + pre_overflow[r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(2)_w{0}_r{1}".format(w, r))
            else:
                model.addConstr(tmp_sum + overflow[w-1][r] - y[w][r] * C[r] - overflow[w][r] <= 0, name="(2)_w{0}_r{1}".format(w, r))

    #for w in range(W): # DON'T NEED IT ANYMORE
    #    for r in range(R):
    #        if w > 0:
    #            model.addConstr(y[w][r] * C[r] - overflow[w-1][r] >= 0, name="(3)_w{0}_r{1}".format(w, r))

    #for r in range(R): # OPTIONAL # DON'T NEED IT ANYMORE
    #    model.addConstr(overflow[W-1][r] == 0, name="(4)_r{0}".format(r))

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

    if fix_s != None:
        for w in range(W):
            model.addConstr(s[i] == fix_s[i], name="fixed_s_constraint_w{0}".format(w))

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

    slave_optimal_value = model.getAttr("ObjVal")

    return n_value, overflow_value, slave_optimal_value

def checkFeasibility(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, Q, strategySet):
    print " =========================== checking feasibility ================================ "
    for i in range(Q):
        n = strategySet[i]['n']
        overflow = strategySet[i]['overflow']
        y = strategySet[i]['y']
        s = strategySet[i]['s'] # b

        for (w,k) in itertools.product(range(W), range(K)):
            n_wk_sum = np.sum([n[w][t][k] for t in range(T)])
            if n_wk_sum != N_wk[w][k]:
                print i, w, k, "passenger equation"
        for (w,r) in itertools.product(range(W), range(R)):
            if w == 0:
                income = np.sum([n[w][t][k] for k in range(K) for t in resource2team[r] ])
            else:
                income = np.sum([n[w][t][k] for k in range(K) for t in resource2team[r] ]) + overflow[w-1][r]
            outcome = y[w][r] * C[r] + overflow[w][r]
            if outcome < income:
                print i, w, r, "outcome"
        p = np.zeros(W)
        for w in range(W):
           p[w] = np.sum([s[ww] for ww in range(max(w+1-shift, 0),w+1)])
        if (P > np.sum(s) + 0.01) or (P < np.sum(s) - 0.01):
            print "staff constraint P", P, np.sum(s)
        for w in range(W):
            if np.sum([y[w][r] * ar[r] for r in range(R)]) > p[w] + 0.01:
                print i, w, r, "resource staff constraint", np.sum([y[w][r] * ar[r] for r in range(R)]), p[w]
            if y[w][r] > mr[r]:
                print i, w, r, "max resource constraint"

def columnGenerationSolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, Q, slave_optimal_value): # integer indicates different relaxation method
    model = Model("LP")
    model.params.OutputFlag = 0
    model.params.TuneOutput = 0

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

    objective_value_gurobi = model.objVal
    print "# of strategies: {0}, objective value (gurobi): {1}, dual objective value: {2}, delta value {3}, slave optimal value {4}".format(len(Q), objective_value_gurobi, objective_value, delta_value, slave_optimal_value)

    return model, gamma_value, delta_value, q_value, objective_value_gurobi

def columnGeneration(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, Q, maxT, column_generation_iterations=1000, warm_start=True):
    # ========================= column generation ==============================
    start_time = time.time()
    # ================= warm start by precomputing =============================
    if warm_start:
        strategySet, obj_relax, objyn, obj_our, q_value = DesignYNcombined.fullYNcombined(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, Q, maxT)
        mixed_n = np.zeros((W, T, K))
        mixed_z = np.zeros((W, K, M))
        mixed_overflow = np.zeros((W, R))
        for i in range(Q):
            mixed_n += strategySet[i]["n"] * q_value[i]
            mixed_overflow += strategySet[i]["overflow"] * q_value[i]
        for w in range(W):
            for k in range(K):
                for m in range(M):
                    if N_wk[w][k] > 0:
                        mixed_z[w][k][m] = np.sum([E[t][m] * mixed_n[w][t][k] / float(N_wk[w][k]) for t in range(T)])
                    else:
                        mixed_z[w][k][m] = 100 # a large cover

        utility_matrix = np.zeros((W, K, M))
        for w in range(W):
            for k in range(K):
                for m in range(M):
                    utility_matrix[w][k][m] = U_plus[k] * mixed_z[w][k][m] + U_minus[k] * (1 - mixed_z[w][k][m])
        overflow_penalty = np.sum([mixed_overflow[w][r] * phi[r] for r in range(R) for w in range(W)])
        min_utility = np.min(utility_matrix) - overflow_penalty
        print "our solution with usability constraints: {0}".format(min_utility)
        fix_s = strategySet[0]['s']

        checkFeasibility(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, Q, strategySet)

    else:
        strategySet = [] # NO WARM START
        fix_s = None

    gamma = np.ones((W, T, K))
    slave_optimal_value = 0

    # ========================= column generation ==============================

    print "\n\n ======================== column generation ============================="

    pre_obj = None
    obj_cg = None
    cumulative_obj = 0
    decay_factor = 0.6
    cutoff_value = 0.00000001
    for j in range(column_generation_iterations):
        #print "================================== column generation testing =================================="
        if len(strategySet) > 0:
            pre_obj = obj_cg
            cg_model, gamma, delta_value, q_value, obj_cg = columnGenerationSolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, strategySet, slave_optimal_value)
            if not pre_obj:
                pre_obj = obj_cg # no previous objective value
            cumulative_obj = cumulative_obj * decay_factor + (obj_cg - pre_obj) # cumulative improvement
            if (cumulative_obj < cutoff_value) and (j > 100):
                break
            #print "cumulative objective value improvement: {0}".format(cumulative_obj)

        #print "==================================== solving slave problem ===================================="
        slave_n_value, slave_overflow_value, slave_optimal_value = slaveProblem(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, gamma, fix_s)
        tmpQ = {}
        tmpQ["n"] = slave_n_value
        tmpQ["overflow"] = slave_overflow_value

        strategySet.append(tmpQ)
        #print "==================================== q sanity check ==========================================="
        #for i in range(len(Q)):
        #    tmp_sum = sum([gamma[w][t][k] * Q[i]['n'][w][t][k] for w in range(W) for t in range(T) for k in range(K)]) - sum([phi[r] * Q[i]["overflow"][w][r] for w in range(W) for r in range(R)])
        #    print "{0} sanity check: {1}".format(i, tmp_sum)

    elapsed_time = time.time() - start_time
    #print "elapsed time: {0}".format(elapsed_time)
    print "true optimal: {0}, our method: {1}, relaxed solution: {2}".format(obj_cg, obj_our, obj_relax)
    num_iterations = j

    return obj_cg, elapsed_time, num_iterations

    #print "============================ LP relaxation =============================="
    #obj_relax, n_value0, overflow_value, y_value, s_value0, p, z_value = StaffResourceAllocation.LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=0, OverConstr = 0)




if __name__ == "__main__":
    # ============================= main =======================================
    print "======================== main ======================================"
    # ========================= Game Setting ===================================
    W = 5 # number of time windows
    K = 10 # number of passenger types
    R = 5 # number of resources
    mR = 2 # max number of reosurces
    M = 3 # number of attack methods
    P = 15 # number of staff
    shift = 3 # d
    Q = 5
    nT = 25
    teams = util.generateAllTeams(R, mR)
    maxT = 5
    #teams = util.randomGenerateTeams(R, mR, nT)

    print teams

    # ================= random generate game setting ===========================
    #seed = 2345
    seed = random.randint(1, 10000)
    print "random seed: {0}".format(seed)

    resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = DesignYNcombined.randomSetting(seed, W, K ,R, mR, M, P, teams, shift)


    total_arrivals = 0
    for w in range(W):
        for k in range(K):
            total_arrivals += N_wk[w][k]
    print "total arrivals: {0}".format(total_arrivals)

    obj_cg, time_cg, iterations_cg = columnGeneration(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, Q, maxT, column_generation_iterations=1000, warm_start=True)


    """
    # ========================= column generation ==============================
    start_time = time.time()
    # ================= warm start by precomputing =============================
    if warm_start:
        strategySet, obj_relax, objyn, obj_our, q_value = DesignYNcombined.fullYNcombined(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, Q, maxT)
        mixed_n = np.zeros((W, T, K))
        mixed_z = np.zeros((W, K, M))
        mixed_overflow = np.zeros((W, R))
        for i in range(Q):
            mixed_n += strategySet[i]["n"] * q_value[i]
            mixed_overflow += strategySet[i]["overflow"] * q_value[i]
        for w in range(W):
            for k in range(K):
                for m in range(M):
                    if N_wk[w][k] > 0:
                        mixed_z[w][k][m] = np.sum([E[t][m] * mixed_n[w][t][k] / float(N_wk[w][k]) for t in range(T)])
                    else:
                        mixed_z[w][k][m] = 100 # a large cover

        utility_matrix = np.zeros((W, K, M))
        for w in range(W):
            for k in range(K):
                for m in range(M):
                    utility_matrix[w][k][m] = U_plus[k] * mixed_z[w][k][m] + U_minus[k] * (1 - mixed_z[w][k][m])
        overflow_penalty = np.sum([mixed_overflow[w][r] * phi[r] for r in range(R) for w in range(W)])
        min_utility = np.min(utility_matrix) - overflow_penalty
        print "our solution with usability constraints: {0}".format(min_utility)
        fix_s = strategySet[0]['s']

    else:
        strategySet = [] # NO WARM START
        fix_s = None

    checkFeasibility(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, Q, strategySet)

    gamma = np.ones((W, T, K))
    slave_optimal_value = 0

    # ========================= column generation ==============================

    print "\n\n ======================== column generation ============================="

    column_generation_iterations = 2000
    for j in range(column_generation_iterations):
        #print "================================== column generation testing =================================="
        if len(strategySet) > 0:
            cg_model, gamma, delta_value, q_value, obj_cg = columnGenerationSolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, strategySet, slave_optimal_value)
        #new_Q = [Q[i] for i in range(len(Q)) if q_value[i] > 0]
        #Q = new_Q

        #print "delta value: {0}".format(delta_value)

        #print "==================================== solving slave problem ===================================="
        slave_n_value, slave_overflow_value, slave_optimal_value = slaveProblem(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, gamma, fix_s)
        tmpQ = {}
        tmpQ["n"] = slave_n_value
        tmpQ["overflow"] = slave_overflow_value

        strategySet.append(tmpQ)
        #print "==================================== q sanity check ==========================================="
        #for i in range(len(Q)):
        #    tmp_sum = sum([gamma[w][t][k] * Q[i]['n'][w][t][k] for w in range(W) for t in range(T) for k in range(K)]) - sum([phi[r] * Q[i]["overflow"][w][r] for w in range(W) for r in range(R)])
        #    print "{0} sanity check: {1}".format(i, tmp_sum)

    elapsed_time = time.time() - start_time
    print "elapsed time: {0}".format(elapsed_time)
    print "true optimal: {0}, our method: {1}, relaxed solution: {2}".format(obj_cg, obj_our, obj_relax)

    #print "============================ LP relaxation =============================="
    #obj_relax, n_value0, overflow_value, y_value, s_value0, p, z_value = StaffResourceAllocation.LPsolver(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi, integer=0, binary_y=0, OverConstr = 0)

    """

