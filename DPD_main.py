#!/usr/bin/env python2
# -*- coding: utf8 -*-
'''
Created on Mar 15, 2013

@author: easyrider

Non object oriented version
'''

import csv, argparse, os, glob
import ConfigParser
from random import randrange

#------------------------------------------------------------------------------.
# List contatining particles. Main data structure.
Universe = []

#------------------------------------------------------------------------------.
# Particle is a tuple consisting of (numbers are indexes):
Type = 0
Mass = 1
Rads = 2

Rx = 3
Ry = 4
Rz = 5

Vx = 6
Vy = 7
Vz = 8

Ax = 9
Ay = 10
Az = 11

#------------------------------------------------------------------------------.
class SimParams:
    SpaceSize       = 0
    StepsN          = 0
    Temperature     = 0
    StepTime        = 0

    ConsForceFactor = 0
    DissForceFactor = 0
    RandForceFactor = 0

    GravityFactor   = 0
    CutoffRadius    = 0

#------------------------------------------------------------------------------.
def initSimParams(cParser):
    ParamsSectionName = 'SimulationParameters'
    SimParams.SpaceSize   = cParser.getint(ParamsSectionName, 'spacesize')
    SimParams.StepsN      = cParser.getint(ParamsSectionName, 'steps')
    SimParams.Temperature = cParser.getfloat(ParamsSectionName, 'temperature')
    SimParams.StepTime    = cParser.getfloat(ParamsSectionName, 'time_per_step')

    SimParams.ConsForceFactor = cParser.getfloat(ParamsSectionName, 's_c')
    SimParams.DissForceFactor = cParser.getfloat(ParamsSectionName, 's_c')
    SimParams.RandForceFactor = cParser.getfloat(ParamsSectionName, 'xi')

    SimParams.GravityFactor = cParser.getfloat(ParamsSectionName, 'gravity_factor')
    SimParams.CutoffRadius  = cParser.getfloat(ParamsSectionName, 'radius_const')

    cParser.remove_section(ParamsSectionName)

#------------------------------------------------------------------------------.
def initPartTypes(cParser):
    # Section should be removed by initSimParams()
    assert 'SimulationParameters' not in cParser.sections()

    pTypes = []
    for partTypeSect in cParser.sections():
        partN  = cParser.getint(partTypeSect, 'numberofparticles')
        mass   = cParser.getfloat(partTypeSect, 'mass')
        radius = cParser.getfloat(partTypeSect, 'radius')
        pTypes.append((partN, mass, radius))
    return pTypes

#------------------------------------------------------------------------------.
def initUniverse(partTypes):
    for (pI, pT) in enumerate(partTypes, start=1):
        Universe.extend(genInitCond2(pI, pT[0], pT[1], pT[2]))

#------------------------------------------------------------------------------.
def writeHeader(OutFilePrefix, pTypes):
    with open(OutFilePrefix, 'wb') as f:
        wr = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
        
        header = [1, len(Universe), SimParams.StepsN, SimParams.StepTime, SimParams.SpaceSize]
        wr.writerow(header)
        wr.writerow("")
        
        for pT in pTypes:
            wr.writerow(pT)

#------------------------------------------------------------------------------.
def writeUniverseState(OutFilePrefix, stepnr):
    filename = OutFilePrefix + "_step_" + str(stepnr).zfill(6)
    with open(filename, 'wb') as f:
        wr = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
        
        for line in Universe:
            wr.writerow(line[:1] + line[3:9])

#------------------------------------------------------------------------------.

'''
Generates initial condition with based on the given particle. 
The function places randomly x Particles into the given space and returns a list of Particles.
'''
def genInitCond2(ptype, numOfParticles, mass, radius):
    particles=[]
    for _ in range(numOfParticles):
        particles.append([ptype, mass, radius,
            randrange(0, SimParams.SpaceSize), #x
            randrange(0, SimParams.SpaceSize), #y
            randrange(0, SimParams.SpaceSize), #z
            0, 0, 0, #V
            0, 0, 0]) #A
    return particles

#------------------------------------------------------------------------------.
# see below for definition of arguments
def conservativeForce(Sc, R, Ex, Ey, Ez):
    T = Sc * (1.0 - R)
    return (T*Ex, T*Ey, T*Ez)

def dissipativeForce(Sd, R, Vx, Vy, Vz, Ex, Ey, Ez):
    T = -1 * Sd*(1.0-R)**2 * (Vx*Ex + Vy*Ey + Vz*Ez)
    return (T*Ex, T*Ey, T*Ez)

def randomForce(Sr, Sd, kB, T, stepT, R, Ex, Ey, Ez):
    # TODO Nie ma zmiennej losowej!!!
    T = (2*Sd*kB*T/stepT)**0.5 * (1.0 - R) * Sr
    return (T*Ex, T*Ey, T*Ez)

#------------------------------------------------------------------------------.
def applyForces(P, pToQ, pqRelDist):
# args: P - particle
#       pToQ - normalized vector from P to Q (other particle affecting P)
#       pqRelDist - distance between P and Q relative to CutoffRadius
    k_b = 1.3806488 * 10.0 **(-23)
    
    fCons = conservativeForce(SimParams.ConsForceFactor, pqRelDist,
            pToQ[0], pToQ[1], pToQ[2])
    print('fCons: {}'.format(fCons))
    
    fDiss = dissipativeForce(SimParams.DissForceFactor, pqRelDist,
            P[Vx], P[Vy], P[Vz], pToQ[0], pToQ[1], pToQ[2])
    print('fDiss: {}'.format(fDiss))
    
    fRand = randomForce(SimParams.RandForceFactor, SimParams.DissForceFactor,
            k_b, SimParams.Temperature, SimParams.StepTime, pqRelDist,
            pToQ[0], pToQ[1], pToQ[2])
    print('fRand: {}'.format(fRand))
    
    massRev = 1.0 / P[Mass]
    accelChange = (
            massRev * (fCons[0] + fDiss[0] + fRand[0]),
            massRev * (fCons[1] + fDiss[1] + fRand[1]),
            massRev * (fCons[2] + fDiss[2] + fRand[2]))
    
    P[Ax] += accelChange[0]
    P[Ay] += accelChange[1]
    P[Az] += accelChange[2]
        
#------------------------------------------------------------------------------.
def relativeVec(P, Q, pqDist):
    #vector pointing from P to Q
    pToQ = [P[3]-Q[3], P[4]-Q[4], P[5]-Q[5]]

    #normalization
    if pqDist!=0:
        pToQ[0] /= pqDist
        pToQ[1] /= pqDist
        pToQ[2] /= pqDist
    return pToQ

#------------------------------------------------------------------------------.
def twoParticles(P, Q):
    
    pqDist=((P[Rx]-Q[Rx]) ** 2.0 +
            (P[Ry]-Q[Ry]) ** 2.0 +
            (P[Rz]-Q[Rz]) ** 2.0) ** 0.5
    print('P to Q dist: {}'.format(pqDist))
    
    if pqDist < SimParams.CutoffRadius:
        pqRelDist = pqDist / SimParams.CutoffRadius

        # TODO Ten wektor de facto wskazuje od Q do P. To daje dobre wyniki dla
        # sił, jeśli to dyssypatywna ma minus. Zamienić nazwę zmiennej???
        pToQ = relativeVec(P, Q, pqDist)
        
        qToP = (pToQ[0]*(-1), pToQ[1]*(-1), pToQ[2]*(-1))
        
        applyForces(P, pToQ, pqRelDist)
        applyForces(Q, qToP, pqRelDist)
    
#------------------------------------------------------------------------------.
def applyGravity(P):
    P[Az] -= 9.82*SimParams.GravityFactor
    
#------------------------------------------------------------------------------.
def calculateNewPositions(p):
    stepT = SimParams.StepTime
    spcSz = SimParams.SpaceSize

    currentVelocity=(p[Vx], p[Vy], p[Vz])
    
    #calculate current velocity + acceleration
    p[Vx] += p[Ax] * stepT
    p[Vy] += p[Ay] * stepT
    p[Vz] += p[Az] * stepT
    
    #calculate new position
    p[Rx] += ((currentVelocity[0]+p[Vx]) / 2) * stepT
    p[Ry] += ((currentVelocity[1]+p[Vy]) / 2) * stepT
    p[Rz] += ((currentVelocity[2]+p[Vz]) / 2) * stepT

    # if any coordinate jumps out of space, wrap it around
    if p[Rx] < 0 or p[Rx] > spcSz:
        p[Rx] %= spcSz

    if p[Ry] < 0 or p[Ry] > spcSz:
        p[Ry] %= spcSz

    if p[Rz] < 0 or p[Rz] > spcSz:
        p[Rz] %= spcSz
        
    p[Ax]=0
    p[Ay]=0
    p[Az]=0

#------------------------------------------------------------------------------.
def nextStep():
    
    '''apply forces...'''
    for i in range(0, len(Universe)):
        for j in range(i+1, len(Universe)):
            # TODO może tutaj liczyć odległość i wołać jeśli < CutOffRadius???
            twoParticles(Universe[i], Universe[j])
        # TODO może tutaj dodawać grawitację?
           
    '''apply gravity'''
    for p in Universe:
        applyGravity(p)  
    
    '''Calculate new positions'''
    for p in Universe:
        calculateNewPositions(p)

#------------------------------------------------------------------------------.
if __name__ == "__main__":
    print('Program starting.')
    
    #Program start
    aparser = argparse.ArgumentParser()
    aparser.add_argument('ConfigFile', help='input parameter')
    aparser.add_argument('outfile', help='output file')
    
    args = aparser.parse_args()
    
    # TODO Odkomentować!!!
    #OutDirContent = os.listdir(args.OutDir)
    #if OutDirContent:
    #    print('Output directory is not empty! Do you want to erase it? [y/n]')
    #    if raw_input() == 'y':
    #        for f in OutDirContent: os.remove(f)
    #clean output folder
    for fl in glob.glob(args.outfile + "*"):
        os.remove(fl)

    # Read configuration file.
    cParser = ConfigParser.RawConfigParser()
    cParser.read(args.ConfigFile)

    # Initialize simulation parameters...
    initSimParams(cParser)
    # ... and list of particle types.
    partTypes = initPartTypes(cParser)

    # Initialize particles in 'Universe' list (randomize positions, zero other params)
    initUniverse(partTypes)

    # Save header file - mostly simulation parameters.
    writeHeader(args.outfile, partTypes)
    
    # Save initial state
    writeUniverseState(args.outfile, 0)
    
    # Actual simulation loop.
    for i in range(1, SimParams.StepsN):
        nextStep()
        writeUniverseState(args.outfile, i)
        print ("%(i)i/%(total)i" % {"i":i+1, "total":SimParams.StepsN})
    
    print("Finished.")

