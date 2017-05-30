import math as m
import numpy as np
import random
import time

from algebra import la
from tools import Tools as tl
from curves import Curves

from PIL import Image
from scipy import misc

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

        self.image = np.ndarray([self.nx, self.ny])
        self.image.fill(0)
        self.n = [0, 0, 1]
        self.detected_intensity = list()
        self.detected_position = list()

    @property
    def integral(self):
        suma = 0
        for i in range(self.nx):
            for j in range(self.ny):
                suma += self.image[i, j]
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


class Source:
    def __init__(self, loc: list, wavelength: float, intensity: int, number):
        self.loc = loc
        self.wl = wavelength
        self.intensity = intensity
        if type(number) is list:
            self.number = number[0] * number[1]
            self.number_list = number
        else:
            self.number = number
        # self.rays = list()
        self.intensity_per_photon = None
        self.photons_total = 0
        self.photons_reached = 0

        self.max_angle = float()
        self.direction = float()
        self.solid_angle = float()

    def reset(self):
        self.__init__(self.loc, self.wl, self.intensity, self.number)


class Crystal:
    def __init__(self, d, D, r, loc):
        self.d = d
        self.r = r
        self.D = D
        self.centre = la.minus(loc, [0, 0, (r ** 2 - (D / 2) ** 2) ** 0.5])
        self.loc = loc
        self.points = list()
        self.count_reflected = 0

    def bragg_angle(self, wavelength, order: int):
        return m.asin((wavelength * order) / (2 * self.d))


class SetUp:
    def __init__(self, source, crystal: Crystal, detector: Detector):
        for s in source:
            s.reset()
        self.source = source
        self.crystal = crystal
        self.detector = detector
        self.curveSi = Curves('reflectivity.csv')

        for s in self.source:
            s.direction = la.normalize(la.minus(crystal.loc, s.loc))

            alpha = la.cos([0, self.crystal.loc[1] - s.loc[1], self.crystal.loc[2] - s.loc[2]], [0, 0, 1])
            beta = la.cos([self.crystal.loc[0] - s.loc[0], 0, self.crystal.loc[2] - s.loc[2]], [0, 0, 1])

            s.solid_angle = (m.pi * self.crystal.D ** 2 / 4 * alpha * beta) / (
                la.norm(la.minus(crystal.loc, s.loc)) ** 2)

    def shine_random(self, source: Source):
        for i in range(source.number):
            done = False

            while not done:
                s = la.plus(
                    la.minus(self.crystal.loc, source.loc),
                    [self.crystal.D * (0.5 - random.random()) for j in range(2)] + [0]
                )
                done = self.reflection_crystal(la.normalize(s), source)

        source.intensity_per_photon = source.intensity * (source.solid_angle / (4 * m.pi)) / source.photons_total
        # print('ipp {}'.format(source.intensity_per_photon))

    def reflection_crystal(self, s, source):
        cp_loc = la.x(s, tl.qroot(
            a=1,
            b=-2 * la.dot(s, la.minus(self.crystal.centre, source.loc)),
            c=la.norm(self.crystal.centre) ** 2 + la.norm(source.loc) ** 2 - 2 * la.dot(
                self.crystal.centre, source.loc) - self.crystal.r ** 2
        ))
        cp_loc = la.plus(cp_loc, source.loc)

        if la.norm(la.minus(cp_loc, self.crystal.loc)[:2]) < self.crystal.D / 2:
            normal = la.normalize(la.minus(self.crystal.centre, cp_loc))
            source.photons_total += 1

            if not self.reflection_point(cp_loc, normal, s):
                return False

            return True
        return False

    def reflection_point(self, loc, n, ray: list):
        out_intensity = self.curveSi.curve(m.pi / 2 - m.acos(la.cos(ray, la.i(n))))

        if out_intensity != 0:
            self.crystal.points.append(loc)
            ray = la.x(ray, 1 / la.dot(ray, n))
            o = [+ray[i] + 2 * (n[i] - ray[i]) for i in range(3)]
            r = np.array(la.minus(self.detector.loc, loc))

            coeff = np.linalg.solve(np.array([la.normalize(o),
                                              self.detector.ux,
                                              self.detector.uy]).T, r)

            i = int(coeff[1] // self.detector.res + self.detector.nx / 2)
            j = int(coeff[2] // self.detector.res + self.detector.ny / 2)

            if 0 <= i < self.detector.nx and 0 <= j < self.detector.ny:
                self.detector.detected_intensity.append(out_intensity)
                self.detector.detected_position.append([i, j])
            return True
        return False

    def pre_shine(self, source: Source):
        angles = [
            la.cos(la.plus([self.crystal.D / 2, 0, self.crystal.r], la.minus(self.crystal.centre, source.loc)),
                   [self.crystal.D / 2, 0, self.crystal.r]),
            la.cos(la.plus([-self.crystal.D / 2, 0, self.crystal.r], la.minus(self.crystal.centre, source.loc)),
                   [-self.crystal.D / 2, 0, self.crystal.r]),
            la.cos(la.plus([0, self.crystal.D / 2, self.crystal.r], la.minus(self.crystal.centre, source.loc)),
                   [0, self.crystal.D / 2, self.crystal.r]),
            la.cos(la.plus([0, -self.crystal.D / 2, self.crystal.r], la.minus(self.crystal.centre, source.loc)),
                   [0, -self.crystal.D / 2, self.crystal.r])
        ]

        angles = [m.pi / 2 - m.acos(a) for a in angles]

        if angles[0] < self.curveSi.bragg < angles[1] or angles[0] > self.curveSi.bragg > angles[1]:
            return True
        if angles[2] < self.curveSi.bragg < angles[3] or angles[2] > self.curveSi.bragg > angles[3]:
            return True
        return False

    def shine_matrix(self, source: Source):
        # angles = list()
        # for i in range(100):
        #     radius = tl.rotate([self.crystal.D / 2, 0, 0], [0, 0, 2 * m.pi * i / 100])
        #
        #     ray = la.plus(la.minus(self.crystal.loc, source.loc), radius)
        #     normal = la.plus([0, 0, self.crystal.r], radius)
        #     # print(la.cos(ray, normal))
        #     angles.append([i, m.pi / 2 - m.acos(la.cos(ray, normal))])
        # angles_bragg = [a[0] for a in angles if self.curveSi.bragg_lim[0] < a[1] < self.curveSi.bragg_lim[1]]
        # limits = []
        # for i in range(len(angles_bragg) - 1):
        #     if (angles_bragg[i + 1] - angles_bragg[i]) > 1:
        #         limits.append(angles_bragg[i])
        # if len(limits) == 1:
        #     lim_x = [m.cos(angles_bragg[0] * m.pi / 100), m.cos(limits[0] * m.pi / 100)]
        #     lim_y= [m.cos(angles_bragg[0] * m.pi / 100), m.cos(limits[0] * m.pi / 100)]

        # print(self.curveSi.bragg)
        # print(angles_bragg)

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

    def mesh_to_image(self, source: Source):
        pos = self.detector.detected_position
        inte = self.detector.detected_intensity
        for i in range(len(inte)):
            self.detector.image[pos[i][0], pos[i][1]] += inte[i] * source.intensity_per_photon
        pos.clear()
        inte.clear()

    def graph(self):
        f = open('name.txt', 'r+')
        number = int(f.read())
        f.close()

        pilimage = Image.fromarray(self.detector.image.T)
        pilimage.save('images/{}.tiff'.format(number))
        # misc.imsave('images/{}.png'.format(number), self.detector.image.T)

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
        text_info += '''
Source(loc={}, wavelength={}, intensity={}, number={})
        
        '''.format(
            self.source[0].loc,
            self.source[0].wl,
            self.source[0].intensity,
            self.source[0].number
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
        integral = self.detector.integral
        print('Intensity on detector: {}'.format(integral))
        print('Photon fraction on detector: {}'.format(integral / sum([s.intensity for s in self.source])))

    def work(self):
        t = time.time()

        for s in self.source:
            if self.source.index(s) != 0:
                print('{}/{}, {}s'.format(self.source.index(s), len(self.source),
                                          int((time.time() - t) * (len(self.source) / self.source.index(s) - 1))))
            if not self.pre_shine(s):
                print('no')
                continue
            self.shine_random(s)
            self.mesh_to_image(s)

        print('Elapsed time: {}s'.format(time.time() - t))

        self.graph()
        self.statistics()