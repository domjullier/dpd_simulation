'''
Created on Mar 15, 2013

@author: easyrider

Non object oriented version
'''

import copy
import csv, argparse
from random import randrange

if __name__ == '__main__':
    pass



'''
Generates initial condition with based on the given particle. 
The function places randomly x Particles into the given space and returns a list of Particles.
'''

def genInitCond2(ptype, numOfParticles, spaceSize, mass, radius):
    
    particles=[]
       
    for _ in range(numOfParticles):
        '''currentParticle = [type, mass, radius, randrange(0, spaceSize), randrange(0, spaceSize), randrange(0, spaceSize), 0, 0, 0, 0, 0, 0]'''
        particles.append([ptype, mass, radius, randrange(0, spaceSize), randrange(0, spaceSize), randrange(0, spaceSize), 0, 0, 0, 0, 0, 0])
        
    return particles


    
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


def nextStep2(state, step):
    
    '''apply gravity'''
    for p in state:
        p=gravity2(p, step)
    
    '''apply more forces...'''
    
    '''Calculate new positions'''
    particles = calculateNewPositions(state, step)
    
    return particles



aparser = argparse.ArgumentParser()
aparser.add_argument('infile', help='input parameter')
aparser.add_argument('outfile', help='output file')

args = aparser.parse_args()

#generate from inpufile

f = open(args.infile)
lines = f.readlines()
f.close()


spacesize=int(lines[0])
simulatedSteps=int(lines[1])
numberOfTypes=len(lines)-2
numberOfParticles=[]
mass=[]

for i in range(2, len(lines)):
    tmp=lines[i].split(',')
    numberOfParticles.append(int(tmp[0]))
    mass.append(int(tmp[1]))

totalNumberOfParticles=0
for p in numberOfParticles:
    totalNumberOfParticles+=p

simulatedSteps=int(lines[1])
timerPerStep=0.1
numberOfStates=simulatedSteps*totalNumberOfParticles

state=[]
for i in range(0, numberOfTypes):
    state.extend(genInitCond2(ptype=i+1, numOfParticles=numberOfParticles[i], spaceSize=spacesize, mass=mass[i], radius=1))


'''add initial state to history'''
world_history = []
world_history.extend(copy.deepcopy(state))


for i in range(1, simulatedSteps):
    state=nextStep2(state, timerPerStep)
    world_history.extend(copy.deepcopy(state))
 
f = open(args.outfile, 'wb')

    
numberOfStates/=totalNumberOfParticles

header = [1,totalNumberOfParticles, numberOfStates, timerPerStep, spacesize]

wr = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)

wr.writerow(header)
wr.writerow("")

i=1
for n, m  in numberOfParticles, mass:
    wr.writerow([i, n, m])
    i+=1


wr.writerow("")
cnt = 0
for line in world_history:
    if cnt == totalNumberOfParticles:
        wr.writerow("")
        cnt=0
        
    cnt=cnt+1
    wr.writerow(line[:1] + line[3:9])
    

print "end"
