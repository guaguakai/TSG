'''
Created on Jan 17, 2018

@author: Sara
'''
from DesignProblem import solve as solve1
from DesignProblemFixedResources import solve as solve2
from KStrategies import randomSetting as rs
import util
def to_team(v):
        team = []
        j=0
        for i in v:
            if i > 0:
                team.append(j)
            j+=1
        return team
def printTeams(teams, q):
    print "full relaxation"
    print teams[0]
    print "IntegerY Single N"
    print teams[1]
    print "IntegerY multiple N"
    for i in range(len(q)):
        if q[i]>0:
            print teams[2][i]
    print "Full MIP"
    for i in range(len(q)):
        if q[i]>0:
            print teams[3][i]
def usedStrategies(q):
    s=0
    for i in q:
        if i>0: s+=1
    return s
if __name__ == "__main__":
    # ============================= main =======================================
    print "======================== main ======================================"
    # ========================= Game Setting ===================================
    W = 5 # number of time windows
    K = 5 # number of passenger types
    R = 5 # number of resources
    mR = 10 # max number of resources
    M = 2 # number of attack methods
    P = 25 # number of staff
    shift = 5 # d
    Q = 5
    nT = 20
    maxT=5
    teams = util.generateAllTeams(R, mR)
    #teams = util.randomGenerateTeams(R, mR, nT)


    # ================= random generate game setting ===========================
    seed = 2345
    #resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, minr, ar, phi = randomSetting(seed, W, K ,R, mR, M, P, teams, shift)
    resource2team, T, Er, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi = rs(seed, W, K ,R, mR, M, P, teams, shift)

    
    #util1, rt1, teams1 = solve1(3, W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, YConstr=True)
    #util2, rt2, teams2 = solve1(3, W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, YConstr=False)
    #util3, rt3, teams3 = solve1(5, W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, YConstr=True)
    #util4, rt4, teams4 = solve1(5, W, K, R, mR, M, P, teams, resource2team, T, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, YConstr=False)


    from DesignProblemFixedResources import testsolve
    #testsolve(Q, W, K, R, mR, M, P, teams, resource2team, T, maxT, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi, TeamConstr=False)
    
    util1, rt1, q1, teams1 = solve2(3, W, K, R, mR, M, P, teams, resource2team, T, 5, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi,  TeamConstr=False, YConstr=False)
    util2, rt2, q2, teams2 = solve2(5, W, K, R, mR, M, P, teams, resource2team, T, 5, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi,  TeamConstr=False, YConstr=False)
    #util3, rt3, q3, teams3 = solve2(5, W, K, R, mR, M, P, teams, resource2team, T, 15, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi,  TeamConstr=True, YConstr=False)
    #util4, rt4, q4, teams4 = solve2(5, W, K, R, mR, M, P, teams, resource2team, T, 15, E, C, U_plus, U_minus, N_wk, shift, mr, ar, phi,  TeamConstr=True, YConstr=True)


    print ""
    print util1, rt1, usedStrategies(q1)

    print util2, rt2, usedStrategies(q2)
    print util3, rt3, usedStrategies(q3)
    print util4, rt4, usedStrategies(q4)

    
    printTeams(teams1, q1)
    printTeams(teams2, q2)
    printTeams(teams3, q3)


