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
        self.x, self.y, self.z = values if values is not None else (0.0, 0.0, 0.0)

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
        #self.values = np.dot(self.values, translation_matrix)
        self.values = np.dot(translation_matrix, self.values)  # translation matrix is on the left, same as my original code

    def rotate_axis(self, axis_vector, angle):
        axis_vector = axis_vector.copy()
        axis_vector.normalize()
        x, y, z = axis_vector.x, axis_vector.y, axis_vector.z
        c, s = np.cos(angle), np.sin(angle)
        t = 1 - c

        rotation_matrix = np.array([
            [t*x*x + c, t*x*y - s*z, t*x*z + s*y, 0],
            [t*x*y + s*z, t*y*y + c, t*y*z - s*x, 0],
            [t*x*z - s*y, t*y*z + s*x, t*z*z + c, 0],
            [0, 0, 0, 1]
        ])
        #self.values = np.dot(self.values, rotation_matrix)
        self.values = np.dot(rotation_matrix, self.values) # rotation matrix is on the left, same as my original code

    def roll_line(self, p1, p2, angle):
        """
        Rotate the line defined by p1 and p2 around the line by angle
        """
        axis_vector = p2 - p1
        axis_vector.normalize()
        self.translate(p1 * -1)
        self.rotate_axis(axis_vector, angle)
        self.translate(p1)


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

    # test matrix class
    print("\nTest Matrix class:")
    m = Matrix()
    
    o = Vector([0, 0, 0])
    vt = Vector([1, 2, 3])
    m.translate(vt)
    print("Translate by vt %s:" % vt)
    print(m)
    print("Apply to origin:", m.apply_to_vector(o))
    print("Apply to vt:", m.apply_to_vector(vt))

    print("\nTranslet it forward and back:")
    m.reset()
    m.translate(vt)
    print("Translate by vt: %s" % vt)
    print(m)
    m.translate(vt * -1)
    print("Translate by -vt: %s" % (vt * -1))
    print(m)

    # test rotation
    print("\nTest rotation:")
    m.reset()
    p1 = Vector([0, 0, 0])
    p2 = Vector([1, 0, 0])
    angle = np.pi / 2
    m.roll_line(p1, p2, angle)
    print("Rotate 90 degree around x-axis:")
    vr = Vector([0, 1, 0])
    print("Apply to vr %s, expect to be (0, 0, 1):" % vr)
    print(m.apply_to_vector(vr))

    m.reset()
    p1 = Vector([0, 0, 0])
    p2 = Vector([0, 1, 0])
    angle = np.pi / 2
    m.roll_line(p1, p2, angle)
    print("\nRotate 90 degree around y-axis:")
    vr = Vector([2, 0, 0])
    print("Apply to vr %s, expect to be (0, 0, -2):" % vr)
    print(m.apply_to_vector(vr))

    m.reset()
    p1 = Vector([1, 0, 0])
    p2 = Vector([0, 1, 0])
    angle = np.pi / 2
    m.roll_line(p1, p2, angle)
    print("\nRotate 90 degree around %s -> %s:" % (p1, p2))
    vr = Vector([0, 0, 0])
    print("Apply to vr %s, expect to be (0.5, 0.5, 0.707):" % vr)
    print(m.apply_to_vector(vr))


