'''
Created on Mar 15, 2013

@author: easyrider
'''

class Vector(object):
    '''
    classdocs
    '''


    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        
    def __repr__(self):
        return repr((self.x, self.y, self.z))

    
    # overload []
    def __getitem__(self, index):
        '''data = [self.x,self.y,self.z]'''
        return (self.x, self.y, self.z)[index]
    
    # overload set []
    def __setitem__(self, key, item):
        if (key == 0):
            self.x = item
        elif (key == 1):
            self.y = item
        elif (key == 2):
            self.z = item
        else:
            raise Exception("Index doesnt exist!")
        
    def __add__(self, other):
        self.x+=other.x
        self.y+=other.y
        self.z+=other.z
        
        return Vector(self.x, self.y, self.z)
    
    def __div__(self, other):
        if (isinstance(other, int)):
            self.x/=other
            self.y/=other
            self.z/=other
        else:
            self.x/=other.x
            self.y/=other.y
            self.z/=other.z
        
        return Vector(self.x, self.y, self.z)