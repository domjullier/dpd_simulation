'''
Created on Mar 15, 2013

@author: easyrider

Non object oriented version
'''

import csv, argparse
import ConfigParser
from random import randrange

if __name__ == '__main__':
    pass

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

def save_state(state, stepnr):
    
    filename = args.outfile + "_step_" + str(stepnr).zfill(6)
    f = open(filename, 'wb')
    wr = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
    
    for line in state:
        wr.writerow(line[:1] + line[3:9])
    
    f.close()
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
    
def applyForces(particle, r, step, e):
    k_b=1.3806488 * 10.0**(-23)
 
    
    #conservativeForce(s_c, r, e_x, e_y, e_z)
    f_c=conservativeForce(s_c, r, e[0], e[1], e[2])
    
    #dissipativeForce(s_d, r, v_x, v_y, v_z, e_x, e_y, e_z)
    f_d=dissipativeForce(s_d, r, particle[6], particle[7], particle[8], e[0], e[1], e[2])
    
    #randomForce(xi,s_d,k_b,temperature,delta_t, r, e_x, e_y, e_z)
    f_r=randomForce(xi,s_d,k_b,temperature,step, r, e[0], e[1], e[2])
    
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
            p[3]=(p[3]%spacesize) #put the particle on to the opposite side of room
            
        if(p[3]>spacesize):
            p[3]=(p[3]%spacesize)
            
        #y position
        if(p[4]<0):
            p[4]=(p[4]%spacesize) #put the particle on to the opposite side of room
            
        if(p[4]>spacesize):
            p[4]=(p[4]%spacesize)
            

        p[9]=0
        p[10]=0
        p[11]=0

        
    return particles


def nextStep2(state, step, spacesize):
    
    
    
    '''apply more forces...'''
    cnt=0
    for i in range(0, len(state)):
        for j in range(i+1, len(state)):
            radius_sum=state[i][2]+state[j][2]
            #radius_sum=100
            distance=((state[i][3]-state[j][3])**2.0+(state[i][4]-state[j][4])**2.0+(state[i][5]-state[j][5])**2.0)**(1.0/2.0)
            if(distance<=radius_sum):
                r=1.0-(distance/radius_sum)
                #e vector
                e_i=[state[i][3]-state[j][3], state[i][4]-state[j][4], state[i][5]-state[j][5]]
                #normalization
                length=(e_i[0]**2.0+e_i[1]**2.0+e_i[2]**2.0)**(1.0/2)
                e_i[0]/=length
                e_i[1]/=length
                e_i[2]/=length
                
               
                #e_j=(state[j][3]-state[i][3], state[j][4]-state[i][4], state[j][5]-state[i][5])
                
                e_j=(e_i[0]*(-1), e_i[1]*(-1), e_i[2]*(-1))
               
                
                state[i]=applyForces(state[i], r, step, e_i)
                state[j]=applyForces(state[j], r, step, e_j)
                cnt+=1
                
    #print cnt
           
            #print state[i][3] #3,4,5 is position. 2 is radius      
    '''apply gravity'''
    for p in state:
        p=gravity2(p, step)  
    
    '''Calculate new positions'''
    particles = calculateNewPositions(state, step, spacesize)
    
    return particles


#Program start
aparser = argparse.ArgumentParser()
aparser.add_argument('infile', help='input parameter')
aparser.add_argument('outfile', help='output file')

args = aparser.parse_args()

#generate from inpufile
Config = ConfigParser.ConfigParser()
Config.read(args.infile)



#f = open(args.infile)
#lines = f.readlines()
#f.close()

spacesize = int(ConfigSectionMap("SimulationParameters") ['spacesize'])
simulatedSteps = int(ConfigSectionMap("SimulationParameters") ['steps'])
s_c = int(ConfigSectionMap("SimulationParameters") ['s_c'])
s_d = int(ConfigSectionMap("SimulationParameters") ['s_d'])
xi = int(ConfigSectionMap("SimulationParameters") ['xi'])
temperature = int(ConfigSectionMap("SimulationParameters") ['temperature'])

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

timePerStep=0.01
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

    

print "end"
