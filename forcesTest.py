'''
Created on May 17, 2013

@author: easyrider
'''

pm = __import__('DPD_main')

#parameters
s_c=10.
s_d=1.
xi=1
temperature=100
radius_const=200.
distance=50
k_b=1.3806488 * 10.0**(-23)
step=0.1
radius_const=200

#Particles
#particles.append([ptype, mass, radius, randrange(0, spaceSize), randrange(0, spaceSize), randrange(0, spaceSize), 0, 0, 0, 0, 0, 0])
particle_1 = ([1, 10, 100, 10, 10, 10, 0, 0, 0, 0, 0, 0])  
particle_2 = ([1, 10, 100, 20, 10, 10, 0, 0, 0, 0, 0, 0])     

print [particle_1, particle_2]
print pm.twoParticles(particle_1, particle_2, step, radius_const, s_c, s_d, xi, temperature)

print "Calculate new position"
print pm.calculateNewPositions([particle_1, particle_2], step, spacesize=100)