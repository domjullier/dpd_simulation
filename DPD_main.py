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

def gravity(particle):
    if(particle.position.z<=0):
        particle.position.z=0
        particle.acceleration.z=0
        particle.velocity.z=0

        return particle
    
    else:
        particle.acceleration+=Vector(0,0,-9.82)
        return particle
        
def calculateNewPositions(particles):
    for p in particles:

        currentVelocity=Vector(p.velocity.x, p.velocity.y, p.velocity.z)
        p.velocity+=p.acceleration

        p.position+=((currentVelocity+p.velocity)/2)
        
        if p.position.z<=0:
            p.position.z=0
            p.velocity.z=0
            

        p.acceleration=Vector(0,0,0)

        
    return particles

def nextStep(particles):
    
    '''apply gravity'''
    for p in particles:
        p=gravity(p)
    
    '''apply more forces...'''
    
    '''Calculate new positions'''
    particles = calculateNewPositions(particles)
    
    return particles

start = genInitCond(numOfParticles=1, spaceSize=1000, mass=1, radius=1)


print "Beginning"
for p in start:
    print p


for i in range(10):
    start=nextStep(start)
    print "Step %s" % (i)
    for p in start:
        print p
