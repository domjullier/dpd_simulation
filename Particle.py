'''
Created on Mar 15, 2013

@author: easyrider
'''

class Particle(object):
    '''
    classdocs
    '''


    def __init__(self, mass, radius, position, velocity, acceleration):
        self.mass=mass
        self.radius=radius
        self.position=position
        self.velocity=velocity
        self.acceleration=acceleration

    def __repr__(self):
        return ('PARTICLE DETAILS\nPosition: %s\nVelocity: %s\nAcceleration: %s\n' % (self.position, self.velocity, self.acceleration))
       

    def get_mass(self):
        return self.__mass


    def get_radius(self):
        return self.__radius


    def get_position(self):
        return self.__position


    def get_velocity(self):
        return self.__velocity


    def get_acceleration(self):
        return self.__acceleration


    def set_mass(self, value):
        self.__mass = value


    def set_radius(self, value):
        self.__radius = value


    def set_position(self, value):
        self.__position = value


    def set_velocity(self, value):
        self.__velocity = value


    def set_acceleration(self, value):
        self.__acceleration = value

    mass = property(get_mass, set_mass, None, None)
    radius = property(get_radius, set_radius, None, None)
    position = property(get_position, set_position, None, None)
    velocity = property(get_velocity, set_velocity, None, None)
    acceleration = property(get_acceleration, set_acceleration, None, None)
        
        