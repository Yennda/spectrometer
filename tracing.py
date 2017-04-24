import os
import math as m
import numpy as np
import random
import scipy as sc
from copy import copy

"""
All the positions and dimensions are in centimetres. Except for the wavelengths and lattice constants, that are in Angstroems.
"""


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
    def normalize(vec):
        norm = np.linalg.norm(vec)
        if norm == 0:
            raise ValueError
        return vec / norm

    @staticmethod
    def rotate(vec, angles, u='rad'):
        if u == 'r':
            angles = angles * m.pi
        elif u == 'd':
            angles = [Tools.rad_from_deg(a) for a in angles]

        Rx = np.matrix(
            [
                [1, 0, 0],
                [0, m.cos(angles[0]), -m.sin(angles[0])],
                [0, m.sin(angles[0]), m.cos(angles[0])]
            ])
        Ry = np.matrix(
            [
                [m.cos(angles[1]), 0, m.sin(angles[1])],
                [1, 0, 0],
                [-m.sin(angles[1]), 0, m.cos(angles[1])]
            ])
        Rz = np.matrix(
            [
                [m.cos(angles[2]), -m.sin(angles[2]), 0],
                [m.sin(angles[2]), m.cos(angles[2]), 0],
                [0, 0, 1]
            ])
        print(vec * (Rz * Ry * Rx).T)
        return np.array(vec * (Rz * Ry * Rx).T)


class Detector:
    def __init__(self, dim: list, loc: list, resolution, ):
        """
        :param dimensions [cm]: 
        :param loc: 
        :param resolution: [micrometres]
        """
        self.dim = dim
        self.loc = loc
        self.res = resolution / 1e4  # to centimeters
        self.points
        self.n = None
        self.mesh = None

    def generate_mesh(self):
        if self.n is None:
            raise ValueError

        nx = int(self.dim[0] / self.res)
        ny = int(self.dim[1] / self.res)

        # generating mesh
        self.mesh = np.ndarray([nx, ny])
        for i in range(-nx // 2, nx // 2 + 1):
            for j in range(-ny // 2, ny // 2 + 1):
                self.mesh[i, j] = DetectorPoint(
                    self.res * np.array([i, j, 0]),
                )

    def rotate(self):
        pass


class DetectorPoint:
    def __init__(self, loc):
        self.loc = loc


class Source:
    def __init__(self, loc: list, wavelength):
        self.loc = np.array(loc)
        self.wl = wavelength


class CrystalPoint:
    def __init__(self, loc: np.array, n: np.array):
        """
        :param loc: coordinates [cm]
        :param n: normal vector
        """
        self.loc = loc
        self.n = n
        self.out = list()


class Crystal:
    def __init__(self, d, D, R, n):
        """
        :param d: lattice constant [Ang]
        :param D: diameter of the crystal [cm]
        :param R: radius of curvature [cm]
        :param n: number of points 
        :return: 
        """
        self.d = d
        self.R = R
        self.D = D
        self.theta_max = m.atan(d / R)
        self.n = n
        self.points = self.generate_points_random()

    def bragg_angle(self, wavelength, order: int):
        return m.asin((wavelength * order) / (2 * self.d))

    def abs_coordinates(self, theta, phi):
        x = self.R * m.sin(phi) * m.sin(theta)
        y = self.R * m.cos(phi) * m.sin(theta)
        z = self.R * m.cos(theta)
        return np.array([x, y, z])

    def generate_points(self):
        points = list()
        for i in range(-self.n, self.n + 1):
            theta = self.theta_max / self.n * i
            n_phi = int(self.n * theta / self.theta_max)
            for j in range(abs(n_phi)):
                phi = 2 * m.pi / n_phi * j
                points.append(CrystalPoint(
                    self.abs_coordinates(theta, phi),
                    self.abs_coordinates(theta, phi) / self.R)
                )
        return points

    def generate_points_random(self):
        points = list()

        for i in range(self.n):
            x = 2
            y = 2
            while ((x - 0.5) ** 2 + (y - 0.5) ** 2) > 0.5 ** 2:
                x = random.random()
                y = random.random()

            loc = np.array([
                x * self.D - self.D / 2,
                y * self.D - self.D / 2,
                (self.R ** 2 -
                 (x * self.D - self.D / 2) ** 2 -
                 (y * self.D - self.D / 2) ** 2) ** 0.5
            ]
            )
            points.append(CrystalPoint(loc, -loc / self.R))

        return points


class SetUp:
    def __init__(self, source: Source, crystal: Crystal, detector: Detector):
        self.source = source
        self.crystal = crystal
        self.detector = detector

    def compute_reflected(self):
        bragg = self.crystal.bragg_angle(self.source.wl, 2)
        print(Tools.deg_from_rad(bragg))
        tf = list()
        for p in self.crystal.points:
            s = Tools.normalize(p.loc - self.source.loc)
            # s = s / np.linalg.norm(s)

            # print(Tools.deg_from_rad(m.pi / 2 - m.acos(s.dot(-p.n))))

            tf.append(Tools.mol((m.pi / 2 - m.acos(s.dot(-p.n))), bragg, 0.02))

            if Tools.mol((m.pi / 2 - m.acos(s.dot(-p.n))), bragg, 0.04):
                p.out.append(2 * (p.n + s) - s)

        print('{}/{}'.format(tf.count(True), self.crystal.n))

    def compute_image(self):
        self.detector.n = Tools.normalize(
            np.array([0, 0, self.crystal.R]) - self.loc
        )
