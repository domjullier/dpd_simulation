'''
Created on Mar 15, 2013

@author: easyrider

Non object oriented version
'''

import copy
import csv
from Vector import Vector
from Particle import Particle
from random import randrange

if __name__ == '__main__':
    pass



'''
Generates initial condition with based on the given particle. 
The function places randomly x Particles into the given space and returns a list of Particles.
'''
def genInitCond(numOfParticles, spaceSize, mass, radius):
    
    particles=[]
       
    for _ in range(numOfParticles):
        particles.append(Particle(mass, radius, Vector(randrange(0, spaceSize), randrange(0, spaceSize), randrange(0, spaceSize)), Vector(0, 0, 0), Vector(0, 0, 0)))
        
    return particles

def genInitCond2(ptype, numOfParticles, spaceSize, mass, radius):
    
    particles=[]
       
    for _ in range(numOfParticles):
        '''currentParticle = [type, mass, radius, randrange(0, spaceSize), randrange(0, spaceSize), randrange(0, spaceSize), 0, 0, 0, 0, 0, 0]'''
        particles.append([ptype, mass, radius, randrange(0, spaceSize), randrange(0, spaceSize), randrange(0, spaceSize), 0, 0, 0, 0, 0, 0])
        
    return particles

def gravity(particle, step):
    if(particle.position.z<=0):
        particle.position.z=0
        particle.acceleration.z=0
        particle.velocity.z=0

        return particle
    
    else:
        particle.acceleration+=Vector(0,0,(-9.82*step))
        return particle
    
def gravity2(particle, step):
    if(particle[5]<=0):
        particle[5]=0 #position z
        particle[11]=0 #acceleration z
        particle[8]=0 #velocity z

        return particle
    
    else:
        particle[11]+=(-9.82*step)
        return particle
        
def calculateNewPositions(particles, step):
    for p in particles:

        currentVelocity=(p[6], p[7], p[8])
        
        #calculate current velocity + acceleration
        p[6]+=p[9]
        p[7]+=p[10]
        p[8]+=p[11]
        
        #calculate new position
        p[3]+=((currentVelocity[0]+p[6])/2)*step
        p[4]+=((currentVelocity[1]+p[7])/2)*step
        p[5]+=((currentVelocity[2]+p[8])/2)*step

        
        if p[5]<=0: #pos z
            p[5]=0
            p[8]=0 #velo z
            

        p[9]=0
        p[10]=0
        p[11]=0

        
    return particles

def nextStep(particles, step):
    
    '''apply gravity'''
    for p in particles:
        p=gravity(p, step)
    
    '''apply more forces...'''
    
    '''Calculate new positions'''
    particles = calculateNewPositions(particles, step)
    
    return particles

def nextStep2(state, step):
    
    '''apply gravity'''
    for p in state:
        p=gravity2(p, step)
    
    '''apply more forces...'''
    
    '''Calculate new positions'''
    particles = calculateNewPositions(state, step)
    
    return particles

def getListsForVisualization(particleWorld):
    x = []
    y = []
    z = []
    
    newList = [[]]
    
    newList.append(x, y, z)
    
'''#Full structure
|Full structure        |number of particles    |number of states    |time step between states    |size of world
|type of particle    |mass    |radius        |position x|y|z    |velocity x|y|z    |acceleration x|y|z
|type of particle    |mass    |radius        |position x|y|z    |velocity x|y|z    |acceleration x|y|z
...
'''
    

state = genInitCond2(ptype=1, numOfParticles=20, spaceSize=1000, mass=1, radius=1)
state.extend(genInitCond2(ptype=2, numOfParticles=20, spaceSize=1000, mass=1, radius=1))

'''add initial state to history'''
world_history = []
world_history.extend(copy.deepcopy(state))




world_history.extend(copy.deepcopy(state))

#print world_history 
     

for i in range(1, 1500):
    state=nextStep2(state, 0.01)
    world_history.extend(copy.deepcopy(state))
    #print "Step %s" % (i)
    #for p in state:
    #    print p

#print world_history 

#|Full structure        |number of particles    |number of states    |time step between states    |size of world

f = open('workfile', 'wb')

#write fileheader
numberOfStates = 0
numberOfParticles = 40

for _ in world_history:
    numberOfStates+=1
    
numberOfStates/=numberOfParticles

header = [1,numberOfParticles, numberOfStates, 0.01, 1000]

wr = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)

wr.writerow(header)
wr.writerow("")
wr.writerow([1, 1, 1])
wr.writerow([2, 1, 1])

wr.writerow("")
cnt = 0
for line in world_history:
    if cnt == numberOfParticles:
        wr.writerow("")
        cnt=0
        
    cnt=cnt+1
    wr.writerow(line[:1] + line[3:9])
    
#wr.writerow(line[:1] + line[3:6])
print "end"
