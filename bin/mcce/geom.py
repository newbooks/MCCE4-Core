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


