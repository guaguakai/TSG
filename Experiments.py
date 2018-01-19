'''
Created on Jan 17, 2018

@author: Sara
'''
from DesignProblem import solve as solve1
from DesignProblemFixedResources import solve as solve2
from KStrategies import randomSetting as rs
import util

if __name__ == "__main__":
    # ============================= main =======================================
    print "======================== main ======================================"
    # ========================= Game Setting ===================================
    W = 15 # number of time windows
    K = 10 # number of passenger types
    R = 5 # number of resources
    mR = 10 # max number of reosurces
    M = 2 # number of attack methods
    P = 5 # number of staff
    shift = 5 # d
    Q = 6
    nT = 20
    maxT=5
    teams = util.generateAllTeams(R, mR)
    #teams = util.randomGenerateTeams(R, mR, nT)


    # ================= random generate game setting ===========================
    seed = 2345
    #resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
    resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = rs(seed, W, K ,R, mR, M, P, teams, shift)

    
    util1, rt1, tt = solve1(Q, W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi)
    util2, rt2, teams = solve2(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi,  TeamConstr=True)
    util3, rt3, teams3 = solve2(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi,  TeamConstr=False)

    print ""
    print util1, rt1
    print tt[0]
    print tt[1]
    print util2, rt2
    print util3, rt3

    def to_team(v):
        team = []
        j=0
        for i in v:
            if i > 0:
                team.append(j)
            j+=1
        return team
    print to_team(teams[0])
    print
    for i in range(Q):
        print to_team(teams[1][i])
    print
    for i in range(Q):
        print to_team(teams[2][i])


