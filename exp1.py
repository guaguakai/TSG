import Method1
import Method2
import column
import random
import util
import numpy as np
import pickle

if __name__ == "__main__":
    W = 5 # number of time windows
    K = 3 # number of passenger types
    R = 6 # number of resources
    mR = 3 # max number of reosurces
    M = 2 # number of attack methods
    P = 15 # number of staff
    shift = 2 # d
    nT = 25
    teams = util.generateAllTeams(R, mR)
    maxT = 5

    iterations = 5
    Q_start = 1
    Q_end = 3

    objective_values = np.zeros((Q_end - Q_start + 1, 3, iterations))
    running_time = np.zeros((Q_end - Q_start + 1, 3, iterations))
    cg_objective = np.zeros(iterations)
    cg_time = np.zeros(iterations)

    obj_method1 = obj_method2 = obj_relax = rt_method1 = rt_method2 = rt_relax = 0 # initialization

    for i in range(iterations):
        # ================= random generate game setting ===========================
        seed = random.randint(1, 10000)
        #resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
        resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = Method2.randomSetting(seed, W, K ,R, mR, M, P, teams, shift)

        tmp_obj_cg, tmp_time_cg, tmp_iterations_cg = column.columnGeneration(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, 1, maxT, column_generation_iterations=1000, warm_start=True)
        print "Column Generation, i = {0}, obj = {1}, running time = {2}".format(i, tmp_obj_cg, tmp_time_cg)
        cg_objective[i] = tmp_obj_cg
        cg_time[i] = tmp_time_cg

        for Q in range(Q_start, Q_end + 1):
            print " ============================================ Q: {0}, i: {1} ==============================================".format(Q, i)

            obj_method1, rt_method1 = Method1.solve(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, verbose=False)
            print "Method 1, Q = {0}, i = {1}, obj = {2}, running time = {3}".format(Q, i, obj_method1, rt_method1)
            obj_relax, objyn_method2, obj_method2, rt_method2, rt_relax = Method2.solve(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, TeamConstr=True, verbose=False)
            print "Method 2, Q = {0}, i = {1}, obj = {2}, running time = {3}".format(Q, i, obj_method2, rt_method2)
            print "Relaxed obj = {0}, YN combined obj = {1}, running time = {2}".format(obj_relax, objyn_method2, rt_relax)

            objective_values[Q-Q_start][0][i] = obj_method1
            objective_values[Q-Q_start][1][i] = obj_method2
            objective_values[Q-Q_start][2][i] = obj_relax

            running_time[Q-Q_start][0][i] = rt_method1
            running_time[Q-Q_start][1][i] = rt_method2
            running_time[Q-Q_start][2][i] = rt_relax

    print "objective values: {0}".format(objective_values)
    print "running time: {0}".format(running_time)
    print "CG objective: {0}".format(cg_objective)
    print "CG time: {0}".format(cg_time)

    pickle.dump((objective_values, running_time, cg_objective, cg_time), open("data/0128_kai_1_5.pl", "wb"))
