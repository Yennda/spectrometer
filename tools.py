from algebra import la
import numpy as np
import numpy.linalg as l
import math as m


class Tools:
    @staticmethod
    def mm_from_ang(x):
        return x * 1e-12

    @staticmethod
    def deg_from_rad(x):
        return x / (2 * m.pi) * 360

    @staticmethod
    def rad_from_deg(x):
        return x / 180 * m.pi

    @staticmethod
    def mol(x, y, e):
        return ((y - e) < x) and (x < (y + e))

    @staticmethod
    def mols(x, y, e: float):
        norm = la.norm(la.minus(x, y))
        return norm < e

    @staticmethod
    def normalize(vec):
        norm = l.norm(vec)
        if norm == 0:
            raise ValueError
        return vec / norm

    @staticmethod
    def rotate(vec: list, angles, u='rad'):
        if u == 'r':
            angles = la.x(angles, m.pi)
        elif u == 'd':
            angles = [Tools.rad_from_deg(a) for a in angles]

        Rx = [
            [1, 0, 0],
            [0, m.cos(angles[0]), -m.sin(angles[0])],
            [0, m.sin(angles[0]), m.cos(angles[0])]
        ]
        Ry = [
            [m.cos(angles[1]), 0, m.sin(angles[1])],
            [0, 1, 0],
            [-m.sin(angles[1]), 0, m.cos(angles[1])]
        ]
        Rz = [
            [m.cos(angles[2]), -m.sin(angles[2]), 0],
            [m.sin(angles[2]), m.cos(angles[2]), 0],
            [0, 0, 1]
        ]
        # return (Rx * Ry * Rz * vec.T).T
        return la.dotmv(Rx, la.dotmv(Ry, la.dotmv(Rz, vec)))

    @staticmethod
    def gauss(x, s, mi):
        output=m.exp(-(x - mi) ** 2 / (2 * s ** 2))
        if output<1e-9:
            return 0
        return output
