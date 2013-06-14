#!/usr/bin/env python2
# -*- coding: utf8 -*-



# Skrypt wypisujący wartości sił działających na 2 cząstki.
# Potrzebuje pewnych funkcji 'print' w kodzie DPD_main!!!



from DPD_main import *

# Parametry imitujące prawdziwą symulację.
SimParams.SpaceSize       = 3
SimParams.StepsN          = 10
SimParams.Temperature     = 100
SimParams.StepTime        = 0.01

SimParams.ConsForceFactor = 10
SimParams.DissForceFactor = 1
SimParams.RandForceFactor = 1

SimParams.GravityFactor   = 0.0001
SimParams.CutoffRadius    = 5

# Funkcja wypisująca cząsteczkę, w ludzkim formacie.
def printPart(P):
    print('T: {}, M: {}, R: {}\tRx: {}, Ry: {}, Rz: {}\t\
            Vx: {}, Vy: {}, Vz: {}\tAx: {}, Ay: {}, Az: {}'.format(
            P[0], P[1], P[2], P[3], P[4], P[5], P[6],
            P[7], P[8], P[9], P[10], P[11]))

# Wylosowanie dwóch cząsteczek (tylko pozycji, prędkość i przyspieszenie są zerami)
P, Q = genInitCond2(1, 2, 1, 1)

# Początkowe wypisanie.
printPart(P)
printPart(Q)

# Pętla imitująca działanie symulacji, oblicza interakcje między dwiema
# powyższymi cząstkami i wypisuje ich stan.
for i in range(100):

    # Ten if jest przydatny: chcemy zobaczyć jak siły zachowują się po wielu
    # krokach symulacji, szczególnie cieakwa jest siła dyssypatywna, która zależy od prędkości
    if i>96:
        print('Iteration----------------------------------------')
        twoParticles(P,Q)

        printPart(P)
        printPart(Q)

        calculateNewPositions(P)
        calculateNewPositions(Q)

        printPart(P)
        printPart(Q)

