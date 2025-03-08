#!/usr/bin/env python
"""
Module Name: geom

Description:
This module provides geometry operation related classes and functions.
"""

import numpy as np

class Vector:
    """
    Vector class
    """
    def __init__(self, values=None):
        self.x, self.y, self.z = values if values else (0.0, 0.0, 0.0)

    def __add__(self, other):
        return Vector([self.x + other.x, self.y + other.y, self.z + other.z])

    def __sub__(self, other):
        return Vector([self.x - other.x, self.y - other.y, self.z - other.z])

    def __mul__(self, other):
        return Vector([self.x * other, self.y * other, self.z * other])

    def __truediv__(self, other):
        return Vector([self.x / other, self.y / other, self.z / other])

    def __str__(self):
        return f"({self.x}, {self.y}, {self.z})"

    def __repr__(self):
        return self.__str__()

    def to_np(self):
        return np.array([self.x, self.y, self.z])
    
    def norm(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)

    def normalize(self):
        n = self.norm()
        self.x /= n
        self.y /= n
        self.z /= n

    def dot(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other):
        return Vector([
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x
        ])

    def angle(self, other):
        return np.arccos(self.dot(other) / (self.norm() * other.norm()))
    
    def distance(self, other):
        return (self - other).norm()
    
    def copy(self):
        return Vector([self.x, self.y, self.z])



class Matrix:
    """
    Matrix class for 3D transformations
    """
    def __init__(self, values=None):
        if values is None:
            self.values = np.identity(4)
        else:
            self.values = np.array(values).reshape(4, 4)

    def __mul__(self, other):
        if isinstance(other, Matrix):
            return Matrix(np.dot(self.values, other.values))
        elif isinstance(other, Vector):
            vec = np.array([other.x, other.y, other.z, 1.0])
            result = np.dot(self.values, vec)
            return Vector(result[:3])
        else:
            raise TypeError("Unsupported multiplication")

    def __str__(self):
        return str(self.values)

    def __repr__(self):
        return self.__str__()

    def reset(self):
        self.values = np.identity(4)

    def translate(self, translation_vector):
        tx, ty, tz = translation_vector.x, translation_vector.y, translation_vector.z
        translation_matrix = np.array([
            [1, 0, 0, tx],
            [0, 1, 0, ty],
            [0, 0, 1, tz],
            [0, 0, 0, 1]
        ])
        self.values = np.dot(self.values, translation_matrix)

    def rotation_axis(self, point1_vector, point2_vector, angle):
        axis = point2_vector - point1_vector
        axis.normalize()
        x, y, z = axis.x, axis.y, axis.z
        c, s = np.cos(angle), np.sin(angle)
        t = 1 - c

        rotation_matrix = np.array([
            [t*x*x + c, t*x*y - s*z, t*x*z + s*y, 0],
            [t*x*y + s*z, t*y*y + c, t*y*z - s*x, 0],
            [t*x*z - s*y, t*y*z + s*x, t*z*z + c, 0],
            [0, 0, 0, 1]
        ])
        self.values = np.dot(self.values, rotation_matrix)

    def apply_to_vector(self, vector):
        vec = np.array([vector.x, vector.y, vector.z, 1.0])
        result = np.dot(self.values, vec)
        return Vector(result[:3])


if __name__ == "__main__":
    # test vector class
    print("Test Vector class:")
    v0 = Vector()
    v1 = Vector([1, 0, 0])
    v2 = Vector([0, 1, 0])
    print("v0:", v0)
    print("v1:", v1)
    print("v2:", v2)

    # test basic vector operations
    print("v1 + v2:", v1 + v2)
    print("v0 + v1 + v2:", v0 + v1 + v2)
    print("v1 - v2:", v1 - v2)
    print("v1 * 2:", v1 * 2)
    print("v1 / 2:", v1 / 2)

    # test advanced vector operations
    print("v1.norm():", v1.norm())
    v1.normalize()
    print("v1 after v1.normalize():", v1)
    print("v1.dot(v2):", v1.dot(v2))
    print("v1.cross(v2):", v1.cross(v2))
    angle = v1.angle(v2)
    print("v1.angle(v2) in radian:", angle)
    print("v1.angle(v2) in degree:", angle / np.pi * 180)
    print("v1.to_np():", v1.to_np())


