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
        self.rot = 3 * [0]
        self.ux = [1, 0, 0]
        self.uy = [0, 1, 0]

        self.res = res / 1e4  # micrometers to centimeters
        self.nx = int(self.dim[0] / self.res)
        self.ny = int(self.dim[1] / self.res)

        self.n, self.mesh = self.generate_mesh()
        self.translate(self.loc)

        self.suma_intensity = 0

    def generate_mesh(self):
        n = [0, 0, 1]

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
        if u == 'r':
            angles = la.x(angles, m.pi)
        elif u == 'd':
            angles = [tl.rad_from_deg(a) for a in angles]
        self.rot += angles

        self.n = tl.rotate(self.n, angles)
        self.ux = tl.rotate(self.ux, angles)
        self.uy = tl.rotate(self.uy, angles)

        self.translate(la.i(self.loc))

        for i in range(-self.nx // 2, self.nx // 2 + 1):
            for j in range(-self.ny // 2, self.ny // 2 + 1):
                self.mesh[i][j].loc = tl.rotate(self.mesh[i][j].loc, angles)

        self.translate(self.loc)


class DetectorPoint:
    def __init__(self, loc):
        self.loc = loc
        self.intensity = 0


class Source:
    def __init__(self, wavelength: float, intensity: int, number: int):
        self.wl = wavelength
        self.intensity = intensity
        self.intensity_per_photon = None
        self.number = number
        self.rays = list()
        self.count_reached_crystal = 0


class CrystalPoint:
    def __init__(self, loc: list, n: list):
        self.loc = loc
        self.n = n
        self.ray_in = list()
        self.ray_out = list()


class Crystal:
    def __init__(self, d, D, r, loc):
        self.d = d
        self.r = r
        self.D = D
        self.loc_centre = la.minus(loc, [0, 0, r])
        self.loc = loc
        self.points = list()
        self.count_reflected = 0

    def bragg_angle(self, wavelength, order: int):
        return m.asin((wavelength * order) / (2 * self.d))


class SetUp:
    def __init__(self, source: Source, crystal: Crystal, detector: Detector):
        self.source = source
        self.crystal = crystal
        self.detector = detector
        self.bragg = crystal.bragg_angle(source.wl, 2)
        self.direction = la.normalize(crystal.loc)
        self.max_angle = m.atan(crystal.D / 2 / la.norm(crystal.loc))
        self.source.intensity_per_photon = (1 - m.cos(self.max_angle)) / 2 * source.intensity / source.number

        print('Bragg angle: {:.4f}°'.format(tl.deg_from_rad(self.bragg)))

    def reflect_point_eff(self, s):
        cp_loc = la.x(s, tl.qroot(
            a=1,
            b=-2 * la.dot(s, self.crystal.loc_centre),
            c=la.norm(self.crystal.loc_centre) ** 2 - self.crystal.r ** 2
        ))
        if la.norm(la.minus(cp_loc, self.crystal.loc)[:2]) < self.crystal.D / 2:
            normal = la.i(la.normalize(la.minus(cp_loc, self.crystal.loc_centre)))
            cpoint = CrystalPoint(loc=cp_loc, n=normal)
            if not self.ray_on_point_eff(cpoint, la.normalize(cp_loc)):
                return False

            self.crystal.points.append(cpoint)
            self.source.count_reached_crystal += 1

            return True
        return False

    def ray_on_point_eff(self, point: CrystalPoint, ray: list):

        def rock_curve(x):
            return tl.gauss(x, mi=self.bragg, s=0.0014544410433286077 / 3)

        out_intensity = self.source.intensity_per_photon * rock_curve(m.pi / 2 - m.acos(la.cos(ray, la.i(point.n))))

        if out_intensity != 0:
            self.crystal.count_reflected += 1
            rayout = [out_intensity * (-ray[i] + 2 * (point.n[i] + ray[i])) for i in range(3)]

            point.ray_in.append(ray)
            point.ray_out.append(rayout)
            return True
        # if True, shines to the whole crystal surface, if False, shines only to reflecting area
        return False

    def shine_eff(self):
        for i in range(self.source.number):
            done = False
            while not done:
                # s = [0.5 - random.random() for j in range(3)]
                s = la.plus(
                    self.crystal.loc,
                    [self.crystal.D * (0.5 - random.random()) for j in range(2)] + [0]
                )

                if la.cos(self.direction, s) < m.cos(self.max_angle):
                    continue
                done = self.reflect_point_eff(la.normalize(s))
            self.source.rays.append(la.normalize(s))

    def intensity_for_detector(self):
        self.detector.suma_intensity = 0
        t = time.time()
        for i in range(self.detector.nx):
            if i != 0:
                print('{}/{}, {}s'.format(i, self.detector.nx, int((time.time() - t) * (self.detector.nx / i - 1))))
            for j in range(self.detector.ny):
                self.detector.mesh[i][j].intensity = self.intensity_for_point(self.detector.mesh[i][j])
                self.detector.suma_intensity += self.detector.mesh[i][j].intensity

    def intensity_for_point(self, det_point: DetectorPoint):
        intensity = 0
        for c in self.crystal.points:
            for o in c.ray_out:
                s = la.normalize(o)
                r = la.minus(det_point.loc, c.loc)
                a = m.fabs(la.dot(r, r)) / la.dot(s, r)

                proj_x = la.x(self.detector.ux, la.dot(la.x(s, a), self.detector.ux) - la.dot(self.detector.ux, r))
                proj_y = la.x(self.detector.uy, la.dot(la.x(s, a), self.detector.uy) - la.dot(self.detector.uy, r))

                if la.norm(proj_x) < self.detector.res and la.norm(proj_y) < self.detector.res:
                    intensity += la.norm(o)

        return intensity

    def relativize_image(self):
        image = np.ndarray([self.detector.nx, self.detector.ny])

        maximum = list()
        for i in self.detector.mesh:
            for j in i:
                maximum.append(j.intensity)
        maximum.sort()
        print('Maximal intensity {}'.format(maximum[-1]))

        for i in range(self.detector.nx):
            for j in range(self.detector.ny):
                if maximum[-1] != 0:
                    image[i, j] = self.detector.mesh[i][j].intensity / maximum[-1] * 255
                else:
                    image[i, j] = self.detector.mesh[i][j].intensity

        return image

    def graph(self):
        # f = open('name.txt', 'r+')
        # number = int(f.read())
        # f.close()
        image = self.relativize_image()
        # image = np.ndarray([self.detector.nx, self.detector.ny])
        #
        # for i in range(self.detector.nx):
        #     for j in range(self.detector.ny):
        #         image[i, j] = self.detector.mesh[i][j].intensity
        # if image[i, j] != 0:
        #     print(image[i, j])

        number = time.gmtime()
        misc.imsave(
            'images/detector{:02d}{:02d}{:02d}{:02d}{:02d}.tiff'.format(number.tm_mon, number.tm_mday, number.tm_hour,
                                                                        number.tm_min, number.tm_sec), image, 'tiff')

    def solid_angle(self):
        for i in range(1000)
    def statistics(self):
        print('--------------------\nStatisitcs')
        print('Total photon intensity: {}'.format(self.detector.suma_intensity))
        print('Photon fraction: {}'.format(self.detector.suma_intensity / self.source.intensity))
        print('Photons on crystal {}'.format(self.source.count_reached_crystal))
        print('Photons reflected {}'.format(self.crystal.count_reflected))

