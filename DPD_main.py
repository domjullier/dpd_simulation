'''
Created on Mar 15, 2013

@author: easyrider
'''
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

def gravity(particle, step):
    if(particle.position.z<=0):
        particle.position.z=0
        particle.acceleration.z=0
        particle.velocity.z=0

        return particle
    
    else:
        particle.acceleration+=Vector(0,0,(-9.82*step))
        return particle
        
def calculateNewPositions(particles, step):
    for p in particles:

        currentVelocity=Vector(p.velocity.x, p.velocity.y, p.velocity.z)
        p.velocity+=p.acceleration

        p.position+=((currentVelocity+p.velocity)/2)*step
        
        if p.position.z<=0:
            p.position.z=0
            p.velocity.z=0
            

        p.acceleration=Vector(0,0,0)

        
    return particles

def nextStep(particles, step):
    
    '''apply gravity'''
    for p in particles:
        p=gravity(p, step)
    
    '''apply more forces...'''
    
    '''Calculate new positions'''
    particles = calculateNewPositions(particles, step)
    
    return particles

def getListsForVisualization(particleWorld):
    x = []
    y = []
    z = []
    
    newList = [[]]
    
    newList.append(x, y, z)

particle_world = []
start = genInitCond(numOfParticles=10, spaceSize=1000, mass=1, radius=1)

particle_world.append(start)

print "Beginning"
for p in start:
    print p



for i in range(10):
    start=nextStep(start, 0.1)
    particle_world.append(start)
    print "Step %s" % (i)
    for p in start:
        print p

print "end"
