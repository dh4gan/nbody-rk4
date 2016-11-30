import math
import numpy as np

class Vector3D(object):
    """3D cartesian vector object for dynamics"""
# Initialising Function
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        s= '3D Vector ( %e, %e, %e)' % (self.x, self.y, self.z)
        return s

    def clone(self):
        return Vector3D(self.x, self.y, self.z)
# Vector addition
            
    def add(self,other):
        """ Adds another vector"""
        return Vector3D (self.x + other.x, self.y + other.y, self.z + other.z)
            
# Vector subtraction
            
    def subtract(self,other):
        """ Subtracts another vector"""
        return Vector3D (self.x - other.x, self.y - other.y, self.z - other.z)
            
    def scalarmult(self, num):
        """ Multiplies vector by scalar"""
        return Vector3D(num*self.x, num*self.y,num*self.z)

    
# Magnitude of the Vector
    
    def mag(self):
        """ Takes magnitude of the vector"""
        mag = math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
        return mag

    def unitVector(self):
        """Returns the unit vector"""
        return self.scalarmult(1.0/self.mag())
    
# Scalar Product
            
    def dot(self,other):
        """Returns dot product of two vectors"""
        dotproduct = 0.0
        dotproduct += self.x*other.x
        dotproduct += self.y*other.y
        dotproduct += self.z*other.z
        return dotproduct
    
# Vector Product
    def cross(self,other):
        """Calculates cross product (self x other)"""

        cross = Vector3D(0.0,0.0,0.0)
        cross.x = self.y * other.z - self.z * other.y;
        cross.y = self.z * other.x - self.x * other.z;
        cross.z = self.x * other.y - self.y * other.x;
        return cross

    
# Rotate around the Z axis

    def rotateX(self,angle):
        """Rotates the vector around the x axis"""

        oldvec = self.clone()
        self.x - oldvec.x
        self.y = oldvec.y*np.cos(angle) - oldvec.z*np.sin(angle)
        self.z = oldvec.y*np.sin(angle) + oldvec.z*np.cos(angle)

    def rotateY(self,angle):
        """Rotates the vector around the y axis"""

        oldvec = self.clone()

        self.x = oldvec.x*np.cos(angle) + oldvec.z*np.sin(angle);
        self.y = oldvec.y;
        self.z = -oldvec.x*np.sin(angle) + oldvec.z*np.cos(angle);

    def rotateZ(self, angle):
        """Rotates the vector around the z axis"""
 
        oldvec = self.clone()

        self.x = oldvec.x*np.cos(angle) - oldvec.y*np.sin(angle);
        self.y = oldvec.x*np.sin(angle) + oldvec.y*np.cos(angle);
        self.z = oldvec.z;

    
    
