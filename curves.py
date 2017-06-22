from tools import Tools as tl
import math as m


class Curves:
    def __init__(self):
        self.files = [['reflectivity.csv', 79.1920], ['reflectivity_spec.csv', 3.4858226524]]
        self.data = list()
        self.bragg = list()
        self.bragg_lim = list()
        for f in self.files:
            table = [[tl.rad_from_deg(l[0] / 3600 + f[1]), l[1]] for l in self.list_from_csv(f[0])]
            self.data.append(table)
            bragg = table[0]

            for t in table:
                if t[1] > bragg[1]:
                    bragg = t

            self.bragg.append(bragg[0])
            self.bragg_lim.append([table[0][0], table[-1][0]])

    def list_from_csv(self, file):
        file = open(file, 'r')
        data = file.read()
        dlist = [d.split(';') for d in data.split('\n')]
        final = [[float(i) for i in d] for d in dlist]

        return final

    def curve(self, i, x):
        if not self.bragg_lim[i][0] < x < self.bragg_lim[i][1]:
            return 0
        lfc = self.data[i]

        for j in range(len(lfc)):
            if lfc[j][0] > x:
                return (lfc[j][1] - lfc[j - 1][1]) * ((x - lfc[j - 1][0]) / (lfc[j][0] - lfc[j - 1][0])) + lfc[j - 1][1]


class Polycurve:
    @staticmethod
    def lorentz(x, o, gamma):
        return 1 / m.pi * (0.5 * gamma) / ((x - o) ** 2 + (0.5 * gamma) ** 2) / (2 / gamma)

    @staticmethod
    def fwhm(theta):
        return tl.rad_from_deg(0.00655 / m.cos(theta) + 0.00866 * m.tan(theta))

    @staticmethod
    def curve(theta, x):
        # print(tl.deg_from_rad(theta))
        # print(tl.deg_from_rad(x))
        fwhm = Polycurve.fwhm(theta)
        if m.fabs(x - theta) > 5 * fwhm:
            return 0
        return Polycurve.lorentz(x, theta, fwhm)
