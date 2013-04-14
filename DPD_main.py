'''
Created on Mar 15, 2013

@author: easyrider

Non object oriented version
'''

import copy
import csv, argparse, math
from random import randrange

if __name__ == '__main__':
    pass



def conservativeForce(s_c, r, e_x, e_y, e_z):
    #e = 1,1,1
    return (s_c*(1-r)*e_x, s_c*(1-r)*e_y, s_c*(1-r)*e_z)

def dissipativeForce(s_d, r, v_x, v_y, v_z, e_x, e_y, e_z): 
    return (s_d*(1-r)**2*(v_x*e_x)*e_x, s_d*(1-r)**2*(v_y*e_y)*e_y, s_d*(1-r)**2*(v_z*e_z)*e_z)

def randomForce(s,s_d,k_b,temp,delta_t, r, e_x, e_y, e_z):
    return (s*(2*s_d*k_b*(temp/delta_t))**(1/2)*(1-r)*e_x, s*(2*s_d*k_b*(temp/delta_t))**(1/2)*(1-r)*e_y, s*(2*s_d*k_b*(temp/delta_t))**(1/2)*(1-r)*e_z)


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
    
def applyForces(particle, r, step):
    k_b=1.3806488 * 10.0**(-23)
    e_x=1
    e_y=1
    e_z=1
    
    #conservativeForce(s_c, r, e_x, e_y, e_z)
    f_c=conservativeForce(10, r, e_x, e_y, e_z)
    
    #dissipativeForce(s_d, r, v_x, v_y, v_z, e_x, e_y, e_z)
    f_d=dissipativeForce(1, r, particle[6], particle[7], particle[8], e_x, e_y, e_z)
    
    #randomForce(s,s_d,k_b,temp,delta_t, r, e_x, e_y, e_z)
    f_r=randomForce(1,10,k_b,100,step, r, e_x, e_y, e_z)
    
    #9,10,11: acceleration vector 
    a=((1.0/particle[1])*(f_c[0]+f_d[0]+f_r[0]), (1.0/particle[1])*(f_c[1]+f_d[1]+f_r[1]), (1.0/particle[1])*(f_c[2]+f_d[2]+f_r[2]))
    
    particle[9]=a[0]
    particle[10]=a[1]
    particle[11]=a[2]
    
    
    return particle
        
def calculateNewPositions(particles, step, spacesize):
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

        #z position top and bottom is the end
        if p[5]<0: #pos z
            p[5]=0
            p[8]=0 #velo z
            
        if p[5]>spacesize: #pos z
            p[5]=spacesize
            p[8]=0 #velo z
            
        #x position is going around between right and left    
        if(p[3]<0):
            p[3]=p[3]+spacesize #put the particle on to the opposite side of room
            
        if(p[3]>spacesize):
            p[3]=p[3]-spacesize
            
        #y position
        if(p[4]<0):
            p[4]=p[4]+spacesize #put the particle on to the opposite side of room
            
        if(p[4]>spacesize):
            p[4]=p[4]-spacesize
            

        p[9]=0
        p[10]=0
        p[11]=0

        
    return particles


def nextStep2(state, step, spacesize):
    
    '''apply gravity'''
    #for p in state:
    #    p=gravity2(p, step)
    
    '''apply more forces...'''
    cnt=0
    for i in range(0, len(state)):
        for j in range(i+1, len(state)):
            radius_sum=state[i][2]+state[j][2]
            #radius_sum=100
            distance=((state[i][3]-state[j][3])**2.0+(state[i][4]-state[j][4])**2.0+(state[i][5]-state[j][5])**2.0)**(1.0/2.0)
            if(distance<=radius_sum):
                r=1.0-(distance/radius_sum)
                state[i]=applyForces(state[i], r, step)
                state[j]=applyForces(state[j], r, step)
                cnt+=1
                print "yes"
    #print cnt
           
            #print state[i][3] #3,4,5 is position. 2 is radius      
        
    
    '''Calculate new positions'''
    particles = calculateNewPositions(state, step, spacesize)
    
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
    state.extend(genInitCond2(ptype=i+1, numOfParticles=numberOfParticles[i], spaceSize=spacesize, mass=mass[i], radius=100))


'''add initial state to history'''
world_history = []
world_history.extend(copy.deepcopy(state))


for i in range(1, simulatedSteps):
    state=nextStep2(state, timerPerStep, spacesize)
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
