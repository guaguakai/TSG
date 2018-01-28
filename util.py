import numpy as np
from gurobipy import *
import itertools
import random
import math


def genResources(R, M, NP):
    ER = []
    capacity = []
    for r in range(R):
        E = np.random.rand(M)
        cap = math.floor((1-np.max(E))*NP)
        
        ER.append(E)
        capacity.append(cap)
    return ER, capacity
        
def generateAllTeams(R, mR):
    teams = []
    i =0
    for nR in range(mR):
        iterator = itertools.combinations(range(R), nR+1)
        for x in iterator:
            i+=1
            teams.append(x)
    teams.append([])
    #teams.append([])
    return teams

def computeTeamsRate(R, M, T, teams, Er):
    E = np.zeros((T, M))
    for t in range(T):
        for m in range(M):
            tmp_rate = 1
            for r in teams[t]:
                tmp_rate *= (1 - Er[r][m])
            E[t][m] = 1 - tmp_rate
    return E

def randomGenerateTeams(R, mR, nT):
    allTeams = []
    for nR in range(mR):
        iterator = itertools.combinations(range(R), nR+1)
        for x in iterator:
            allTeams.append(x)

    assert(len(allTeams) >= nT)
    teams = random.sample(allTeams, nT)
    #teams.append([])
    
    return teams

def resourceTeamDict(R, T, teams):
    #print T
    resource2team = []
    for r in range(R):
        resource2team.append([])
    for t in range(T):
        for tmp_r in teams[t]:
            resource2team[tmp_r].append(t)

    return resource2team
