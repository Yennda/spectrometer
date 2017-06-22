import math as m
import numpy as np
import random
import time

from algebra import la
from tools import Tools as tl
from curves import Polycurve
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
    def __init__(self, loc: list, wavelength: list, intensity: int, number):
        self.loc = loc
        self.wl = wavelength
        self.intensity = intensity
        self.number = number
        self.rays = list()
        self.intensity_per_photon = None
        self.photons_total = 0
        self.photons_reached = 0

        self.max_angle = float()
        self.direction = float()
        self.solid_angle = float()

    def reset(self):
        self.__init__(self.loc, self.wl, self.intensity, self.number)


class Crystal:
    def __init__(self, d, D: list, l, tb, rc: int):
        self.d = d
        self.D = D
        self.l = l
        self.tb = tb
        self.rc = rc
        self.ux = tl.rotate([1, 0, 0], [-tb, 0, 0])
        self.uy = tl.rotate([0, 0, 1], [-tb, 0, 0])
        self.normal = tl.rotate([0, 1, 0], [-tb, 0, 0])

        self.points = list()
        self.count_reflected = 0

    def bragg_angle(self, wavelength, order: int):
        return m.asin((wavelength * order) / (2 * self.d))


class SetUp:
    def __init__(self, source, crystal: Crystal, detector: Detector):
        self.source = source
        self.crystal = crystal
        self.detector = detector
        self.curveSi = Curves()

    @property
    def leng(self):
        a, h, d = tl.distances(self.crystal.D[1], self.detector.dim[1], self.crystal.l, self.crystal.tb)
        d -= self.crystal.l
        leng = ((d - m.cos(self.crystal.tb) * self.crystal.D[0] / 2 - m.sin(2 * self.crystal.tb) *
                 self.detector.dim[0] / 2) ** 2 + (
                    h + m.cos(2 * self.crystal.tb) * self.detector.dim[0] / 2) ** 2) ** 0.5
        return leng

    def shine_matrix(self):
        g1 = m.sin(self.crystal.tb) * self.crystal.D[1] / 2 / self.crystal.l
        g2 = m.sin(self.crystal.tb) * self.crystal.D[1] / 2 / (
            self.crystal.l + self.crystal.D[1] * m.cos(self.crystal.tb))

        phi = self.detector.dim[0] / (self.crystal.l + self.leng + m.cos(self.crystal.tb) * self.crystal.D[0] / 2)
        n = 0
        for i in np.arange(-g1, g2, (g1 + g2) / self.source.number):
            if n % 10 == 0:
                print('{}/{}'.format(n, self.source.number))
            n += 1
            for j in np.arange(-phi, phi, 2 * phi / self.source.number):
                s = [j, i, 1]

                self.reflection_crystal(la.normalize(s))
        self.source.intensity_per_photon = self.source.intensity * (
            tl.solid_angle(self.crystal.D[0], self.crystal.l, self.crystal.tb) / (4 * m.pi)) / self.source.photons_total

    def shine_random(self):
        g1 = m.sin(self.crystal.tb) * self.crystal.D[1] / 2 / self.crystal.l
        g2 = m.sin(self.crystal.tb) * self.crystal.D[1] / 2 / (
            self.crystal.l + self.crystal.D[1] * m.cos(self.crystal.tb))

        phi = self.detector.dim[0] / (self.crystal.l + self.leng + m.cos(self.crystal.tb) * self.crystal.D[0] / 2)

        for i in range(self.source.number):
            if i % 100 == 0:
                print(i)

            done = False
            while not done:
                s = [phi * (random.random() - 0.5), random.random() * (g1 + g2) - g2, 1]
                done = self.reflection_crystal(la.normalize(s))

        # solid_angle = (m.cos(m.pi / 2 - g1) - m.cos(m.pi / 2 + g2)) * 2 * phi
        solid_angle = (m.fabs(g1) + m.fabs(g2)) * 2 * phi
        # print('sa: {}'.format(solid_angle))
        self.source.intensity_per_photon = self.source.intensity * (
            solid_angle / (4 * m.pi) / self.source.photons_total)

    def reflection_crystal(self, s):
        coeff = np.linalg.solve(np.array([s,
                                          self.crystal.ux,
                                          self.crystal.uy]).T,
                                [0, 0, self.crystal.l + m.cos(self.crystal.tb) * self.crystal.D[1] / 2])

        cp_loc = la.x(s, coeff[0])
        if cp_loc[0] < self.crystal.D[0] / 2:
            self.source.photons_total += 1

            if not self.reflection_point(cp_loc, s):
                return False

            return True
        return False

    def reflection_point(self, loc, ray: list):
        n = self.crystal.normal
        # out_intensity = self.curveSi.curve(self.crystal.rc, m.pi / 2 - m.acos(la.cos(ray, la.i(n))))
        out_intensity = 0

        for w in self.source.wl:
            out_intensity += Polycurve.curve(self.crystal.bragg_angle(w, 1),
                                             m.pi / 2 - m.acos(la.cos(ray, la.i(n))))
        self.crystal.points.append(loc)

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
                # return True
        return False

    def mesh_to_image(self):
        pos = self.detector.detected_position
        inte = self.detector.detected_intensity
        print('count on det {}'.format(len(inte)))
        for i in range(len(inte)):
            self.detector.image[pos[i][0], pos[i][1]] += inte[i] * self.source.intensity_per_photon
        pos.clear()
        inte.clear()

    def graph(self):
        f = open('name_spec.txt', 'r+')
        number = int(f.read())
        f.close()

        pilimage = Image.fromarray(self.detector.image.T)
        pilimage.save('images_spec/{}.tiff'.format(number))
        # misc.imsave('images/{}.png'.format(number), self.detector.image.T)

        text_info = '''
###################################
{}.txt
---------------------------------
c = Crystal(d={}, D={}, l={}, tb={})
d = Detector(dim={}, loc={}, res={})


        '''.format(
            number,
            self.crystal.d,
            self.crystal.D,
            self.crystal.l,
            self.crystal.tb,
            self.detector.dim,
            self.detector.loc,
            self.detector.res * 1e4,

        )
        text_info += '''
Source(loc={}, wavelength={}, intensity={}, number={})
        
        '''.format(
            self.source.loc,
            self.source.wl,
            self.source.intensity,
            self.source.number
        )
        text_info += '''
---------------------------------
Detector intensity: {}
Photon fraction on detector: {}
Photons on crystal {}
Photons reflected {}
        '''.format(
            self.detector.integral,
            self.detector.integral / self.source.intensity,
            self.source.photons_reached,
            self.crystal.count_reflected
        )

        info = open('info_spec.txt', 'a')
        info.write(text_info)
        info.close()

        f = open('name_spec.txt', 'w')
        f.write(str(number + 1))
        f.close()

    def statistics(self):
        print('--------------------\nStatisitcs')
        integral = self.detector.integral
        print('Intensity on detector: {}'.format(integral))
        print('Photon fraction on detector: {}'.format(integral / self.source.intensity))

    def work(self):
        t = time.time()
        # self.shine_matrix()

        self.shine_random()
        self.mesh_to_image()

        print('Elapsed time: {}s'.format(time.time() - t))

        self.graph()
        self.statistics()
