import math as m
import numpy as np
import random
import time

from scipy import misc
from algebra import la
from tools import Tools as tl

"""
All the positions and dimensions are in centimetres. Except for the wavelengths and lattice constants, that are in Angstroems.
"""


class Detector:
    def __init__(self, dim: list, loc: list, res):
        self.dim = dim
        self.loc = loc
        self.res = res / 1e4  # micrometers to centimeters
        self.nx = int(self.dim[0] / self.res)
        self.ny = int(self.dim[1] / self.res)

        self.n, self.mesh = self.generate_mesh()
        self.translate(self.loc)

    def generate_mesh(self):
        n = np.matrix([0, 0, 1])

        mesh = [[i for i in range(self.ny)] for i in range(self.nx)]

        for i in range(self.nx):
            for j in range(self.ny):
                mesh[i][j] = DetectorPoint(
                    la.x([i - self.nx / 2, j - self.ny / 2, 0], self.res)
                )
        return n, mesh

    def translate(self, vector: list):
        for i in range(self.nx):
            for j in range(self.ny):
                self.mesh[i][j].loc = la.plus(self.mesh[i][j].loc, vector)

    def rotate(self, angles, u='rad'):
        self.n = tl.rotate(self.n, angles, u)

        self.translate(la.i(self.loc))

        for i in range(-self.nx // 2, self.nx // 2 + 1):
            for j in range(-self.ny // 2, self.ny // 2 + 1):
                self.mesh[i][j].loc = tl.rotate(self.mesh[i][j].loc, angles, u)

        self.translate(self.loc)


class DetectorPoint:
    def __init__(self, loc):
        self.loc = loc
        self.intensity = 0


class Source:
    def __init__(self, loc: list, wavelength, intensity=1):
        self.loc = np.array(loc)
        self.wl = wavelength
        self.intensity = intensity


class CrystalPoint:
    def __init__(self, loc: np.array, n: np.array):
        self.loc = loc
        self.n = n
        self.out = list()


class Crystal:
    def __init__(self, d, D, R, n):
        self.d = d
        self.R = R
        self.D = D
        self.theta_max = m.atan(d / R)
        self.n = n
        self.points = self.generate_points_random()
        self.points_em = []

    def bragg_angle(self, wavelength, order: int):
        return m.asin((wavelength * order) / (2 * self.d))

    def abs_coordinates(self, theta, phi):
        x = self.R * m.sin(phi) * m.sin(theta)
        y = self.R * m.cos(phi) * m.sin(theta)
        z = self.R * m.cos(theta)
        return [x, y, z]

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

            loc = [
                x * self.D - self.D / 2,
                y * self.D - self.D / 2,
                (self.R ** 2 -
                 (x * self.D - self.D / 2) ** 2 -
                 (y * self.D - self.D / 2) ** 2) ** 0.5
            ]

            points.append(CrystalPoint(loc=loc, n=la.x(loc, (-1 / self.R))))

        return points


class SetUp:
    def __init__(self, source: Source, crystal: Crystal, detector: Detector):
        self.source = source
        self.crystal = crystal
        self.detector = detector

    def compute_reflected(self):
        bragg = self.crystal.bragg_angle(self.source.wl, 2)

        def rock_curve(x):
            return tl.gauss(x, mi=bragg, s=0.0014544410433286077 / 3)

        for p in self.crystal.points:
            # unit 'wavevector' s, distance r, output intensity io
            # necessary to implement reflectivity of the crystal
            s = la.normalize(la.minus(p.loc, self.source.loc))
            r = la.norm(la.minus(p.loc, self.source.loc))
            io = self.source.intensity / (self.crystal.n * r ** 2) * rock_curve(m.pi / 2 - m.acos(la.dot(s, la.i(p.n))))

            if io != 0:
                p.out.append(la.x((la.minus(la.plus(p.n, s), s)), 2 * io))

        self.crystal.points_em = [p for p in self.crystal.points if p.out != []]

    def intensity_for_detector(self):
        suma = 0
        for i in range(self.detector.nx):
            for j in range(self.detector.ny):
                self.detector.mesh[i][j].intensity = self.intensity_for_point(self.detector.mesh[i][j])
                suma += self.detector.mesh[i][j].intensity
        print('Total photon count: {}'.format(suma))
        print('Photon fraction: {}'.format(suma/self.source.intensity))

    def intensity_for_point(self, det_point: DetectorPoint):
        intensity = 0
        for p in self.crystal.points_em:
            for o in p.out:
                if tl.mols(la.normalize(o), la.normalize(la.minus(det_point.loc, p.loc)), e=0.05):
                    intensity += la.norm(o)
                    # print('zasah!!!')
        return intensity

    def graph(self):
        # f = open('name.txt', 'r+')
        # number = int(f.read())
        # f.close()

        image = np.ndarray([self.detector.nx, self.detector.ny])
        maximum = list()
        for i in self.detector.mesh:
            for j in i:
                maximum.append(j.intensity)
        maximum.sort()
        print(maximum[-1])

        for i in range(self.detector.nx):
            for j in range(self.detector.ny):
                if maximum[-1] != 0:
                    image[i, j] = self.detector.mesh[i][j].intensity / maximum[-1] * 255
                else:
                    image[i, j] = self.detector.mesh[i][j].intensity
        number = time.gmtime()
        misc.imsave(
            'images/detector{:02d}{:02d}{:02d}{:02d}{:02d}.png'.format(number.tm_mon, number.tm_mday, number.tm_hour,
                                                                       number.tm_min, number.tm_sec), image)

        # misc.imsave('images/detector{:03d}.png'.format(number), image)
        # f = open('name.txt', 'w')
        # f.write(str(number + 1))
        # f.close()

    def do(self):
        t = time.time()
        print('reflection...')
        self.compute_reflected()
        print('\t done: {}'.format(time.time() - t))
        print('detector...')
        self.intensity_for_detector()
        print('\t done: {}'.format(time.time() - t))
        print('building image...')
        self.graph()
        print('\t done')

        print('elapsed time: {}'.format(time.time() - t))
