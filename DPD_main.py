'''
Created on Mar 15, 2013

@author: easyrider

Non object oriented version
'''

import csv, argparse, os, glob
import ConfigParser
from random import randrange

#------------------------------------------------------------------------------.
def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print ("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

#------------------------------------------------------------------------------.
def write_Header(simulatedSteps, totalNumberOfParticles, timePerStep, spacesize):
    f = open(args.outfile, 'wb')
    
    header = [1,totalNumberOfParticles, simulatedSteps, timePerStep, spacesize]
    
    wr = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
    
    wr.writerow(header)
    wr.writerow("")
    
    i=1
    for (n, m) in zip(numberOfParticles,mass):
        wr.writerow([i, n, m])
        i+=1

#------------------------------------------------------------------------------.
def save_state(state, stepnr):
    #filename = args.outfile + "_step_" + str(stepnr).zfill(6)
    filename = args.outfile + "_step_" + str(stepnr).zfill(6)
    f = open(filename, 'wb')
    wr = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
    
    for line in state:
        wr.writerow(line[:1] + line[3:9])
    
    f.close()
    pass

#------------------------------------------------------------------------------.
def conservativeForce(Sc, R, Ex, Ey, Ez):
    # temporary
    T = Sc * (1.0 - R)
    return (T*Ex, T*Ey, T*Ez)

def dissipativeForce(Sd, R, Vx, Vy, Vz, Ex, Ey, Ez):
    # temporary
    T = Sd*(1.0-R)**2 * (Vx*Ex + Vy*Ey + Vz*Ez)
    return (T*Ex, T*Ey, T*Ez)

def randomForce(Sr,Sd,kB,temp,delta_t, R, Ex, Ey, Ez):
    #temporary
    T = (2*Sd*kB*(temp/delta_t))**0.5 * (1.0 - R) * Sr
    return (T*Ex, T*Ey, T*Ez)

#------------------------------------------------------------------------------.

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

'''
Generates initial condition with based on the given particle. 
The function places randomly x Particles into the given space and returns a list of Particles.
'''
def genInitCond2(ptype, numOfParticles, spaceSize, mass, radius):
    particles=[]
    for _ in range(numOfParticles):
        '''currentParticle = [type, mass, radius, randrange(0, spaceSize), randrange(0, spaceSize), randrange(0, spaceSize), 0, 0, 0, 0, 0, 0]'''
        particles.append([ptype, mass, radius,
            randrange(0, spaceSize), #x
            randrange(0, spaceSize), #y
            randrange(0, spaceSize), #z
            0, 0, 0, #V
            0, 0, 0]) #A
    return particles


    
#------------------------------------------------------------------------------.
def gravity2(P, step):
    P[Az] -= 9.82*gravityFactor
    return P
    
#------------------------------------------------------------------------------.
def applyForces(P, R, step, e, Sc, Sd, xi, temperature):
    k_b=1.3806488 * 10.0**(-23)
    
    #conservativeForce(Sc, R, e_x, e_y, e_z)
    fCons=conservativeForce(Sc, R, e[0], e[1], e[2])
    
    #dissipativeForce(Sd, R, v_x, v_y, v_z, e_x, e_y, e_z)
    fDiss=dissipativeForce(Sd, R, P[6], P[7], P[8], e[0], e[1], e[2])
    
    #randomForce(xi,Sd,k_b,temperature,delta_t, R, e_x, e_y, e_z)
    fRand=randomForce(xi,Sd,k_b,temperature,step, R, e[0], e[1], e[2])
    
    massRev = 1.0 / P[Mass]
    accelChange=(
            massRev * (fCons[0] + fDiss[0] + fRand[0]),
            massRev * (fCons[1] + fDiss[1] + fRand[1]),
            massRev * (fCons[2] + fDiss[2] + fRand[2]))
    
    P[Ax] += accelChange[0]
    P[Ay] += accelChange[1]
    P[Az] += accelChange[2]
    return P
        
#------------------------------------------------------------------------------.
def calculateNewPositions(particles, step, spacesize):
    for p in particles:
        currentVelocity=(p[Vx], p[Vy], p[Vz])
        
        #calculate current velocity + acceleration
        p[Vx] += p[Ax] * step
        p[Vy] += p[Ay] * step
        p[Vz] += p[Az] * step
        
        #calculate new position
        p[Rx] += ((currentVelocity[0]+p[Vx]) / 2) * step
        p[Ry] += ((currentVelocity[1]+p[Vy]) / 2) * step
        p[Rz] += ((currentVelocity[2]+p[Vz]) / 2) * step

        # if any coordinate jumps out of space, wrap it around
        if p[Rx] < 0 or p[Rx] > spacesize:
            p[Rx] %= spacesize

        if p[Ry] < 0 or p[Ry] > spacesize:
            p[Ry] %= spacesize

        if p[Rz] < 0 or p[Rz] > spacesize:
            p[Rz] %= spacesize
            
        p[Ax]=0
        p[Ay]=0
        p[Az]=0
    return particles

#------------------------------------------------------------------------------.
def relativeVec(P, Q):
    #e vector
    Ei=[P[3]-Q[3], P[4]-Q[4], P[5]-Q[5]]

    #normalization
    L = (Ei[0]**2.0 + Ei[1]**2.0 + Ei[2]**2.0) ** 0.5
    print('L: {}'.format(L))
    if L!=0:
        Ei[0]/=L
        Ei[1]/=L
        Ei[2]/=L
    return Ei

#------------------------------------------------------------------------------.
def twoParticles(P, Q, step, Rcut, s_c, s_d, xi, temperature):
    
    distance=((P[Rx]-Q[Rx]) ** 2.0 +
            (P[Ry]-Q[Ry]) ** 2.0 +
            (P[Rz]-Q[Rz]) ** 2.0) ** 0.5
    
    if distance < Rcut:
        r = distance / Rcut

        Ei = relativeVec(P, Q)
        #Ej=(state[j][3]-state[i][3], state[j][4]-state[i][4], state[j][5]-state[i][5])
        
        Ej=(Ei[0]*(-1), Ei[1]*(-1), Ei[2]*(-1))
       
        
        P = applyForces(P, r, step, Ei, s_c, s_d, xi, temperature)
        Q = applyForces(Q, r, step, Ej, s_c, s_d, xi, temperature)
    return [P, Q]
    
#------------------------------------------------------------------------------.
def nextStep2(state, step, spacesize):
    
    '''apply forces...'''
    cnt=0
    for i in range(0, len(state)):
        for j in range(i+1, len(state)):
            result = twoParticles(state[i], state[j], step, radius_const, s_c, s_d, xi, temperature)
            state[i]=result[0]
            state[j]=result[1]
            cnt+=1
                
    #print cnt
           
            #print state[i][3] #3,4,5 is position. 2 is radius      
    '''apply gravity'''
    for p in state:
        p=gravity2(p, step)  
    
    '''Calculate new positions'''
    particles = calculateNewPositions(state, step, spacesize)
    
    return particles

#------------------------------------------------------------------------------.
if __name__ == "__main__":
    print('start')
    log = open("log", 'wb')
    
    wr_log = csv.writer(log, quoting=csv.QUOTE_NONNUMERIC)
        
    
    
    
    #Program start
    aparser = argparse.ArgumentParser()
    aparser.add_argument('infile', help='input parameter')
    aparser.add_argument('outfile', help='output file')
    
    args = aparser.parse_args()
    
    #generate from inpufile
    Config = ConfigParser.ConfigParser()
    Config.read(args.infile)
    
    #clean output folder
    for fl in glob.glob(args.outfile + "*"):
        os.remove(fl)
    
    #f = open(args.infile)
    #lines = f.readlines()
    #f.close()
    
    spacesize = int(ConfigSectionMap("SimulationParameters") ['spacesize'])
    simulatedSteps = int(ConfigSectionMap("SimulationParameters") ['steps'])
    s_c = float(ConfigSectionMap("SimulationParameters") ['s_c'])
    s_d = float(ConfigSectionMap("SimulationParameters") ['s_d'])
    xi = float(ConfigSectionMap("SimulationParameters") ['xi'])
    temperature = float(ConfigSectionMap("SimulationParameters") ['temperature'])
    timePerStep = float(ConfigSectionMap("SimulationParameters") ['time_per_step'])
    gravityFactor = float(ConfigSectionMap("SimulationParameters") ['gravity_factor'])
    radius_const = float(ConfigSectionMap("SimulationParameters") ['radius_const'])
    #gravityFactor = 0.001
    
    
    #timePerStep=0.1
    
    
    numberOfTypes = len(Config.sections()) - 1
    
    numberOfParticles=[]
    mass=[]
    radius=[]
    
    for p in range (1, len(Config.sections())):
        numberOfParticles.append(int(ConfigSectionMap(Config.sections()[p]) ['numberofparticles']))
        mass.append(int(ConfigSectionMap(Config.sections()[p]) ['mass']))
        radius.append(int(ConfigSectionMap(Config.sections()[p]) ['radius']))
    
    
    totalNumberOfParticles=0
    for p in numberOfParticles:
        totalNumberOfParticles+=p
    
    #timePerStep=0.005
    numberOfStates=simulatedSteps*totalNumberOfParticles
    
    state=[]
    for i in range(0, numberOfTypes):
        state.extend(genInitCond2(ptype=i+1, numOfParticles=numberOfParticles[i], spaceSize=spacesize, mass=mass[i], radius=radius[i]))
    
    
    #Deprecated: Used for classic one file output
    #world_history = []
    #world_history.extend(copy.deepcopy(state))
    
    write_Header(simulatedSteps, totalNumberOfParticles, timePerStep, spacesize)
    
    #save initial state
    save_state(state, 0)
    
    for i in range(1, simulatedSteps):
        state=nextStep2(state, timePerStep, spacesize)
        save_state(state, i)
        print ("%(i)i/%(total)i" % {"i":i+1, "total":simulatedSteps})
        #Deprecated: Used for classic one file output
        #world_history.extend(copy.deepcopy(state))
    
    #Deprecated: Used for classic one file output    
    #save also in one plike for compatibility purposes
    #f = open("output_classic", 'wb')
        
    #header = [1,totalNumberOfParticles, simulatedSteps, timePerStep, spacesize]
    
    #wr = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
    
    #wr.writerow(header)
    #wr.writerow("")
    
    #i=1
    #for (n, m) in zip(numberOfParticles,mass):
    #    wr.writerow([i, n, m])
    #    i+=1
    
    
    #wr.writerow("")
    #cnt = 0
    #for line in world_history:
    #    if cnt == totalNumberOfParticles:
    #        wr.writerow("")
    #        cnt=0
            
    #    cnt=cnt+1
    #    wr.writerow(line[:1] + line[3:9])
        
    #f.close()
    
        
    
    print("end")
