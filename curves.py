import numpy as np
import math as m
from tools import Tools as tl


class Curves:
    def __init__(self, file):
        self.lfc = [[tl.rad_from_deg(l[0] / 3600 + 79.1920), l[1]] for l in
                    self.list_from_csv(file)]
        bragg = self.lfc[0]

        for l in self.lfc:
            if l[1] > bragg[1]:
                bragg = l

        self.bragg = bragg[0]

    def list_from_csv(self, file):
        file = open(file, 'r')
        data = file.read()
        dlist = [d.split('  ') for d in data.split('\n')]
        final = [[float(i) for i in d] for d in dlist]

        return final

    def curve(self, x):
        if not self.lfc[0][0] < x < self.lfc[-1][0]:
            return 0
        for i in range(len(self.lfc)):
            if self.lfc[i][0] > x:
                return (self.lfc[i][1] - self.lfc[i - 1][1]) * (
                    (x - self.lfc[i - 1][0]) / (self.lfc[i][0] - self.lfc[i - 1][0])) + self.lfc[i - 1][1]

        # def rock_curve(x):
        #     return tl.gauss(x, mi=self.bragg, s=0.0014544410433286077 / 3)