from algebra import la
import numpy.linalg as l
import math as m
from scipy import constants as c


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
    def ang_from_kev(e):
        return c.c / (e * 1e3 * c.eV / c.h) * 1e10

    @staticmethod
    def kev_from_ang(l):
        return c.c * c.h / c.e / (l * 1e-10) / 1e3

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
    def rotate(vec, angles):
        # if u == 'r':
        #     angles = la.x(angles, m.pi)
        # elif u == 'd':
        #     angles = [Tools.rad_from_deg(a) for a in angles]

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
        return la.dotmv(Rx, la.dotmv(Ry, la.dotmv(Rz, vec)))

    @staticmethod
    def gauss(x, s, mi):
        output = m.exp(-(x - mi) ** 2 / (2 * s ** 2))
        if output < 1e-9:
            return 0
        return output

    @staticmethod
    def qroot(a, b, c):
        # print(a)
        # print(b)
        # print(c)
        D = b ** 2 - 4 * a * c
        if D < 0:
            raise ValueError
        elif D == 0:
            return -b / (2 * a)
        else:
            x1 = (-b + D ** 0.5) / (2 * a)
            x2 = (-b - D ** 0.5) / (2 * a)
            return x1

    @staticmethod
    def imaging_equation(f, a):
        return 1 / (1 / f - 1 / a)

    @staticmethod
    def shoelace(x: list):
        sum_first = 0
        sum_second = 0
        for i in range(len(x) - 1):
            sum_first += x[i][0] * x[i + 1][1]
            sum_second += x[i + 1][0] * x[i][1]
        return 0.5 * (sum_first + x[-1][0] * x[0][1] - sum_second - x[0][0] * x[-1][1])

    @staticmethod
    def solid_angle(D, l, tb):
        def z(g):
            return (l + m.cos(tb) * D / 2) * m.tan(tb) / (m.tan(tb) - m.tan(g))

        g1 = m.atan(m.sin(tb) * D / (2 * l))
        g2 = m.atan(m.sin(tb) * D / (2 * (l + m.cos(tb) * D)))

        integral = 0
        n = 10000

        for i in range(n):
            d = (g1 + g2) / n
            g = g1 - d * i
            integral += m.tan(d) * 2 * m.atan(D / 2 / z(g))
        return integral

    @staticmethod
    def distances(D, d, l, tb):
        g1 = m.atan(m.sin(tb) * D / (2 * l))
        g2 = m.atan(m.sin(tb) * D / (2 * (l + m.cos(tb) * D)))

        tb1 = tb + g1
        tb2 = tb - g2

        AA = (m.cos(2 * tb) * d + m.sin(tb) * D / 2) / m.tan(tb1 + tb) + m.sin(tb) * D / 2 / m.tan(tb2 + tb) - m.cos(
            tb) * D + m.sin(2 * tb) * d
        BB = 1 / m.tan(tb2 + tb) - 1 / m.tan(tb1 + tb)

        h = AA / BB
        d = (h - m.sin(tb) * D / 2) / m.tan(tb2 + tb) + m.cos(tb) * D
        return [0, h, l + d]