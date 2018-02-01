import Method1
import Method2
import column
import random
import util
import numpy as np
import pickle
import time
import Results

if __name__ == "__main__":
    instance = 1
    if (instance == 1):
        W = 10 # number of time windows
        AI = 3 # interval in which passengers are arriving
        K = 30 # number of passenger types
        R = 7 # number of resources
        mR = 3 # max number of reosurces
        M = 3 # number of attack methods
        P = 50 # number of staff
        shift = 3 # d
        Q = 3 # number of strategies
        maxT = 10

    if (instance == 2):
        W = 7 # number of time windows
        AI = 3 # interval in which passengers are arriving
        K = 15 # number of passenger types
        R = 5 # number of resources
        mR = 3 # max number of reosurces
        M = 2 # number of attack methods
        P = 30 # number of staff
        shift = 2 # d
        Q = 2 # number of strategies
        maxT = 10

    if (instance == 3):
        W = 3 # number of time windows
        AI = 1 # interval in which passengers are arriving
        K = 5 # number of passenger types
        R = 3 # number of resources
        mR = 2 # max number of reosurces
        M = 3 # number of attack methods
        P = 10 # number of staff
        shift = 1 # d
        Q = 1 # number of strategies
        maxT = 10

    teams = util.generateAllTeams(R, mR)

    iterations = 1
    maxT_start = 1
    maxT_list = [10]
    maxT_end = max(maxT_list)

    # ========================= file storage ==============================
    f_q = open("exp/exp4/exp4_0130_0200.csv", "a")
    f_cg = open("exp/exp4/exp4_cg_0130_0200.csv", "a")

    objective_values = np.zeros((maxT_end - maxT_start + 1, 3, iterations))
    running_time = np.zeros((maxT_end - maxT_start + 1, 3, iterations))
    cg_objective = np.zeros(iterations)
    cg_time = np.zeros(iterations)

    obj_method1 = obj_method2 = obj_relax = rt_method1 = rt_method2 = rt_relax = 0 # initialization

    for i in range(iterations):
        # ================= random generate game setting ===========================
        seed = random.randint(1, 10000)
        #resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
        resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = Method2.randomSetting(seed, W, K ,R, mR, M, P, teams, shift)

        #tmp_obj_cg, tmp_time_cg, tmp_iterations_cg, tmp_all_objectives = column.columnGeneration(W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, 1, maxT, column_generation_iterations=1000, warm_start=True)
        #print "Column Generation, i = {0}, obj = {1}, running time = {2}".format(i, tmp_obj_cg, tmp_time_cg)
        #cg_objective[i] = tmp_obj_cg
        #cg_time[i] = tmp_time_cg

        #f_cg.write("i, {0}, maxT, {1}, method, {2}, obj, {3}, running time, {4}, \n".format(i, "NA", "cg", tmp_obj_cg, tmp_time_cg))
        #f_cg.write(", ".join([str(j) for j in tmp_all_objectives]) + "\n")

        for maxT in maxT_list:
            print " ============================================ maxT: {0}, i: {1} ==============================================".format(maxT, i)

            print " ============================================= Method 1 =================================================="
            start_time = time.time()
            obj_method1, rt_method1 = Method1.solve(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, verbose=True)
            print "Method 1, maxT = {0}, i = {1}, obj = {2}, running time = {3}".format(maxT, i, obj_method1, rt_method1)
            print "time: {0}".format(time.time() - start_time)
            f_q.write("i, {0}, maxT, {1}, method, {2}, obj, {3}, running time, {4}, \n".format(i, maxT, 1, obj_method1, rt_method1))

            print " ============================================= Method 2 =================================================="
            start_time = time.time()
            obj_relax, objyn_method2, obj_method2, rt_method2, rt_relax = Method2.solve(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, TeamConstr=True, verbose=True)
            print "Method 2, maxT = {0}, i = {1}, obj = {2}, running time = {3}".format(maxT, i, obj_method2, rt_method2)
            print "time: {0}".format(time.time() - start_time)
            print "Relaxed obj = {0}, YN combined obj = {1}, running time = {2}".format(obj_relax, objyn_method2, rt_relax)
            f_q.write("i, {0}, maxT, {1}, method, {2}, obj, {3}, running time, {4}, \n".format(i, maxT, 2, obj_method2, rt_method2))
            f_q.write("i, {0}, maxT, {1}, method, {2}, obj, {3}, running time, {4}, \n".format(i, maxT, "relaxed", obj_relax, rt_relax))

            objective_values[maxT-maxT_start][0][i] = obj_method1
            objective_values[maxT-maxT_start][1][i] = obj_method2
            objective_values[maxT-maxT_start][2][i] = obj_relax

            running_time[maxT-maxT_start][0][i] = rt_method1
            running_time[maxT-maxT_start][1][i] = rt_method2
            running_time[maxT-maxT_start][2][i] = rt_relax

    print objective_values
    print running_time

    pickle.dump((objective_values, running_time, cg_objective, cg_time), open("data/0128_kai_testing.pl", "wb"))
