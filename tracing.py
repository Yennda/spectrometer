import os
import math as m
import numpy as np
import numpy.linalg as l
import random

from skimage import img_as_ubyte
from scipy import misc

from scipy import ndimage

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
    def mols(x: np.matrix, y: np.matrix, e: float):
        norm = l.norm(x - y)

        return norm < e

    @staticmethod
    def normalize(vec):
        norm = l.norm(vec)
        if norm == 0:
            raise ValueError
        return vec / norm

    @staticmethod
    def rotate(vec: np.matrix, angles, u='rad'):
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
                [0, 1, 0],
                [-m.sin(angles[1]), 0, m.cos(angles[1])]
            ])
        Rz = np.matrix(
            [
                [m.cos(angles[2]), -m.sin(angles[2]), 0],
                [m.sin(angles[2]), m.cos(angles[2]), 0],
                [0, 0, 1]
            ])
        # print(Rx * Ry * Rz * vec.T)
        return (Rx * Ry * Rz * vec.T).T

        # @staticmethod
        # def iter_detector(detector: Detector, function):
        #     for i in range(-detector.nx // 2, detector.nx // 2 + 1):
        #         for j in range(-detector.ny // 2, detector.ny // 2 + 1):
        #             function


class Detector:
    def __init__(self, dim: list, loc: list, res):
        """
        :param dim [cm]: 
        :param loc: 
        :param res: [micrometres]
        """
        self.dim = dim
        self.loc = np.matrix(loc)
        self.res = res / 1e4  # to centimeters
        self.nx = int(self.dim[0] / self.res)
        self.ny = int(self.dim[1] / self.res)

        self.n, self.mesh = self.generate_mesh()
        self.translate(self.loc)

    def generate_mesh(self):
        n = np.matrix([0, 0, 1])

        # generating mesh
        mesh = [[i for i in range(self.ny)] for i in range(self.nx)]

        for i in range(self.nx):
            for j in range(self.ny):
                mesh[i][j] = DetectorPoint(
                    self.res * np.matrix([(i - self.nx) / 2, (j - self.ny) / 2, 0])
                )
        return n, mesh

    def translate(self, vector: np.matrix):
        for i in range(self.nx):
            for j in range(self.ny):
                self.mesh[i][j].loc += vector

    def rotate(self, angles, u='rad'):
        self.n = Tools.rotate(self.n, angles, u)

        for i in range(-self.nx // 2, self.nx // 2 + 1):
            for j in range(-self.ny // 2, self.ny // 2 + 1):
                self.mesh[i][j].loc = Tools.rotate(self.mesh[i][j].loc, angles, u)


class DetectorPoint:
    def __init__(self, loc):
        self.loc = loc
        self.intensity = 0


class Source:
    def __init__(self, loc: list, wavelength, intensity=1):
        self.loc = np.matrix(loc)
        self.wl = wavelength
        self.intensity = intensity


class CrystalPoint:
    def __init__(self, loc: np.matrix, n: np.matrix):
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

            loc = np.matrix([
                x * self.D - self.D / 2,
                y * self.D - self.D / 2,
                (self.R ** 2 -
                 (x * self.D - self.D / 2) ** 2 -
                 (y * self.D - self.D / 2) ** 2) ** 0.5
            ]
            )
            points.append(CrystalPoint(loc=loc, n=-loc / self.R))

        return points


class SetUp:
    def __init__(self, source: Source, crystal: Crystal, detector: Detector):
        self.source = source
        self.crystal = crystal
        self.detector = detector

    def compute_reflected(self):
        bragg = self.crystal.bragg_angle(self.source.wl, 2)
        # print(Tools.deg_from_rad(bragg))
        # print(bragg)

        tf = list()
        for p in self.crystal.points:
            # unit 'wavevector' s, distance r, output intensity io
            # necessary to implement reflectivitz of the crystal
            s = Tools.normalize(p.loc - self.source.loc)
            r = l.norm(p.loc - self.source.loc)
            io = self.source.intensity / r ** 2

            tf.append(Tools.mol((m.pi / 2 - m.acos((s * (-p.n).T)[0, 0])), bragg, 0.02))

            if Tools.mol((m.pi / 2 - m.acos((s * (-p.n).T)[0, 0])), bragg, 0.04):
                p.out.append((2 * (p.n + s) - s) * io)

                # print('{}/{}'.format(tf.count(True), self.crystal.n))

    def intensity_for_detector(self):
        for i in range(self.detector.nx):
            for j in range(self.detector.ny):
                self.detector.mesh[i][j].intensity = self.intensity_for_point(self.detector.mesh[i][j])

    def intensity_for_point(self, det_point: DetectorPoint):
        intensity = 0
        for p in self.crystal.points:
            for o in p.out:
                if Tools.mols(Tools.normalize(o), Tools.normalize(det_point.loc - p.loc), e=0.1):
                    intensity += l.norm(o)
                    # print('zasah!!!')
        return intensity

    def graph(self):
        f = open('name.txt', 'r+')
        number = int(f.read())
        f.close()

        image = np.ndarray([self.detector.nx, self.detector.ny])
        max = list()
        for i in self.detector.mesh:
            for j in i:
                max.append(j.intensity)
        max.sort()
        print(max[-1])

        for i in range(self.detector.nx):
            for j in range(self.detector.ny):
                image[i, j] = self.detector.mesh[i][j].intensity / max[-1] * 255
        print(image)
        misc.imsave('images/detector{:03d}.png'.format(number), image)
        f=open('name.txt','w')
        f.write(str(number+1))
        f.close()