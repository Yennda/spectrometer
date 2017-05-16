import math as m
import numpy as np
import random
import time

from algebra import la
from tools import Tools as tl
from curves import Curves

from PIL import Image

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

        self.image = [[0 for j in range(self.ny)] for i in range(self.nx)]
        self.n = [0, 0, 1]
        self.detected_intensity = list()
        self.detected_position = list()

    @property
    def integral(self):
        suma = 0
        for i in range(self.nx):
            for j in range(self.ny):
                suma += self.image[i][j]
        return suma

    def rotate(self, angles, u='rad'):
        if u == 'r':
            angles = la.x(angles, m.pi)
        elif u == 'd':
            angles = [tl.rad_from_deg(a) for a in angles]
        self.rot += angles

        self.n = tl.rotate(self.n, angles)
        self.ux = tl.rotate(self.ux, angles)
        self.uy = tl.rotate(self.uy, angles)


class DetectorPoint:
    def __init__(self, loc):
        self.loc = loc
        self.intensity = 0


class Source:
    def __init__(self, loc: list, wavelength: float, intensity: int, number: int):
        self.loc = loc
        self.wl = wavelength
        self.intensity = intensity

        if type(number) is list:
            self.number = number[0] * number[1]
            self.number_list = number
        else:
            self.number = number
        self.rays = list()
        self.intensity_per_photon = None
        self.photons_total = 0
        self.photons_reached = 0

        self.max_angle = float()
        self.solid_angle = float()


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
    def __init__(self, source, crystal: Crystal, detector: Detector):
        self.source = source
        self.crystal = crystal
        self.detector = detector
        self.curveSi = Curves('reflectivity.csv')

        self.bragg = [crystal.bragg_angle(s.wl, 2) for s in source]
        self.direction = la.normalize(crystal.loc)
        # self.max_angle = m.atan(crystal.D / 2 / la.norm(crystal.loc))

        for s in self.source:
            s.max_angle = m.atan(crystal.D / 2 / la.norm(la.minus(crystal.loc, s.loc)))
            s.max_angle = m.asin(crystal.D / 2 / la.norm(la.minus(crystal.loc, s.loc)))

            calpha = la.cos([0, self.crystal.loc[1] - s.loc[1], self.crystal.loc[2] - s.loc[2]], [0, 0, 1])
            cbeta = la.cos([self.crystal.loc[0] - s.loc[0], 0, self.crystal.loc[2] - s.loc[2]], [0, 0, 1])
            s.solid_angle = (m.pi * self.crystal.D ** 2 / 4 * calpha * cbeta) / (
                4 * m.pi * la.norm(self.crystal.loc) ** 2)
            # print('Bragg angle: {:.4f}°'.format(tl.deg_from_rad(self.bragg)))
            # print('Bragg angle: {:.4f}°'.format(self.bragg[0]))

    def reflection_crystal(self, s, source):
        cp_loc = la.x(s, tl.qroot(
            a=1,
            b=-2 * la.dot(s, la.minus(self.crystal.loc_centre, source.loc)),
            c=la.norm(self.crystal.loc_centre) ** 2 + la.norm(source.loc) ** 2 - 2 * la.dot(
                self.crystal.loc_centre, source.loc) - self.crystal.r ** 2
        ))
        cp_loc = la.plus(cp_loc, source.loc)
        if la.norm(la.minus(cp_loc, self.crystal.loc)[:2]) < self.crystal.D / 2:
            normal = la.i(la.normalize(la.minus(cp_loc, self.crystal.loc_centre)))

            if not self.reflection_point(cp_loc, normal, s):
                source.photons_total += 1
                return False

            source.photons_total += 1
            return True
        return False

    def reflection_point(self, loc, n, ray: list):
        # self.crystal.points.append(point)
        out_intensity = self.curveSi.curve(m.pi / 2 - m.acos(la.cos(ray, la.i(n))))

        if out_intensity != 0:
            self.crystal.count_reflected += 1
            o = [out_intensity * (-ray[i] + 2 * (n[i] + ray[i])) for i in range(3)]

            r = np.matrix(la.minus(self.detector.loc, loc))
            A = np.matrix([
                la.normalize(o),
                self.detector.ux,
                self.detector.uy
            ]).T
            coeff = A.I * r.T

            i = int(coeff[1, 0] // self.detector.res + self.detector.nx / 2)
            j = int(coeff[2, 0] // self.detector.res + self.detector.ny / 2)

            if 0 <= i < self.detector.nx and 0 <= j < self.detector.ny:
                self.detector.detected_intensity.append(la.norm(o))
                self.detector.detected_position.append([i, j])
            return True
        return False

    def pre_shine(self, source: Source):
        angles = [
            la.cos(la.plus([self.crystal.D / 2, 0, self.crystal.r], la.minus(self.crystal.loc_centre, source.loc)),
                   la.normalize([self.crystal.D / 2, 0, self.crystal.r])),
            la.cos(la.plus([-self.crystal.D / 2, 0, self.crystal.r], la.minus(self.crystal.loc_centre, source.loc)),
                   la.normalize([-self.crystal.D / 2, 0, self.crystal.r])),
            la.cos(la.plus([0, self.crystal.D / 2, self.crystal.r], la.minus(self.crystal.loc_centre, source.loc)),
                   la.normalize([0, self.crystal.D / 2, self.crystal.r])),
            la.cos(la.plus([0, -self.crystal.D / 2, self.crystal.r], la.minus(self.crystal.loc_centre, source.loc)),
                   la.normalize([0, -self.crystal.D / 2, self.crystal.r]))
        ]
        angles = [m.pi / 2 - m.acos(a) for a in angles]

        if angles[0] < self.curveSi.bragg < angles[1] or angles[0] > self.curveSi.bragg > angles[1]:
            return True
        if angles[2] < self.curveSi.bragg < angles[3] or angles[2] > self.curveSi.bragg > angles[3]:
            return True
        return False

    def shine_random(self, source: Source):

        for i in range(source.number):
            done = False
            while not done:
                s = la.plus(
                    la.minus(self.crystal.loc, source.loc),
                    [self.crystal.D * (0.5 - random.random()) for j in range(2)] + [0]
                )

                if la.cos(self.direction, s) < m.cos(source.max_angle):
                    continue

                done = self.reflection_crystal(la.normalize(s), source)

        source.intensity_per_photon = source.intensity * source.solid_angle * source.number / (
            4 * m.pi * source.photons_total * source.number)

    def shine_matrix(self, source: Source):
        num = source.number_list
        for i in range(num[0]):
            for j in range(num[1]):
                vec = [self.crystal.D * (0.5 - i / num[0]), self.crystal.D * (0.5 - j / num[1]), 0]
                s = la.plus(self.crystal.loc, vec)

                if la.cos(self.direction, s) < m.cos(source.max_angle):
                    continue
                self.reflection_crystal(la.normalize(s), source)

        source.intensity_per_photon = source.intensity * source.solid_angle * source.number / (
            4 * m.pi * source.photons_total * source.number)

    def detect(self, point: CrystalPoint):
        r = np.matrix(la.minus(self.detector.loc, point.loc))
        for o in point.ray_out:
            A = np.matrix([
                la.normalize(o),
                self.detector.ux,
                self.detector.uy
            ]).T
            coeff = A.I * r.T

            i = int(coeff[1, 0] // self.detector.res + self.detector.nx / 2)
            j = int(coeff[2, 0] // self.detector.res + self.detector.ny / 2)

            if 0 <= i < self.detector.nx and 0 <= j < self.detector.ny:
                self.detector.detected_intensity.append(la.norm(o))
                self.detector.detected_position.append([i, j])

    def mesh_to_image(self, source):
        pos = self.detector.detected_position
        inte = self.detector.detected_intensity
        for i in range(len(inte)):
            self.detector.image[pos[i][0]][pos[i][1]] += inte[i] * source.intensity_per_photon
        pos.clear()
        inte.clear()

    def relativize_image(self):
        image = np.ndarray([self.detector.nx, self.detector.ny, 3])

        maximum = list()
        for i in self.detector.image:
            for j in i:
                maximum.append(j)
        max_int = max(maximum)

        if max_int != 0:
            for i in range(self.detector.nx):
                for j in range(self.detector.ny):
                    image[i, j] = (
                        int(self.detector.image[i][j] / max_int * 2 ** 24) // 2 ** 16,
                        int(self.detector.image[i][j] / max_int * 2 ** 24) % 2 ** 16 // 2 ** 8,
                        int(self.detector.image[i][j] / max_int * 2 ** 24) % 2 ** 16 % 2 ** 8
                    )
        else:
            image[i, j] = self.detector.image[i][j]

        return image

    def graph(self):
        f = open('name.txt', 'r+')
        number = int(f.read())
        f.close()
        image = np.ndarray([self.detector.ny, self.detector.nx])

        for i in range(self.detector.nx):
            for j in range(self.detector.ny):
                image[j, i] = self.detector.image[i][j]

        pilimage = Image.fromarray(self.detector.image)
        pilimage.save('images/{}.tiff'.format(number))

        text_info = '''
###################################
{}.txt
---------------------------------
c = Crystal(d={}, D={}, r={}, loc={})
d = Detector(dim={}, loc={}, res={})


        '''.format(
            number,
            self.crystal.d,
            self.crystal.D,
            self.crystal.r,
            self.crystal.loc,
            self.detector.dim,
            self.detector.loc,
            self.detector.res * 1e4,

        )
        for s in self.source:
            text_info += '''
Source(loc={}, wavelength={}, intensity={}, number={})
            
            '''.format(
                s.loc,
                s.wl,
                s.intensity,
                s.number
            )
        text_info += '''
---------------------------------
Detector intensity: {}
Photon fraction on detector: {}
Photons on crystal {}
Photons reflected {}
        '''.format(
            self.detector.integral,
            self.detector.integral / sum([s.intensity for s in self.source]),
            sum([s.photons_reached for s in self.source]),
            self.crystal.count_reflected
        )

        info = open('info.txt', 'a')
        info.write(text_info)
        info.close()

        f = open('name.txt', 'w')
        f.write(str(number + 1))
        f.close()

    def statistics(self):
        print('--------------------\nStatisitcs')
        print('Detector intensity: {}'.format(self.detector.integral))
        print('Photon fraction on detector: {}'.format(
            self.detector.integral / sum([s.intensity for s in self.source])
        ))

    def work(self):
        t = time.time()

        for s in self.source:
            if self.source.index(s) != 0:
                print('{}/{}, {}s'.format(self.source.index(s), len(self.source),
                                          int((time.time() - t) * (len(self.source) / self.source.index(s) - 1))))
            if not self.pre_shine(s):
                continue
            self.shine_random(s)
            self.mesh_to_image(s)

        print('Elapsed time: {}s'.format(time.time() - t))
        # print('All photons to crystal: {}'.format(s.total))

        self.graph()
        self.statistics()
