#!/usr/bin/python2

from DPD_main import *

def printPart(P):
    print('T: {}, M: {}, R: {}\tRx: {}, Ry: {}, Rz: {}\t\
            Vx: {}, Vy: {}, Vz: {}\tAx: {}, Ay: {}, Az: {}'.format(
            P[0], P[1], P[2], P[3], P[4], P[5], P[6],
            P[7], P[8], P[9], P[10], P[11]))

#def genInitCond2(ptype, numOfParticles, spaceSize, mass, radius):
uni = genInitCond2(1, 2, 10, 1, 1)

printPart(uni[0])
printPart(uni[1])

E = relativeVec(uni[0], uni[1])

print(E)

