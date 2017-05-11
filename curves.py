import numpy as np
import math as m


class Curves:
    def __init__(self, file):
        self.lfc = [[l[0] / 3600 + 79.1920, l[1]] for l in
                    self.list_from_csv(file)]

    def list_from_csv(self, file):
        file = open(file, 'r')
        data = file.read()
        dlist = [d.split('  ') for d in data.split('\n')]
        final = [[float(i) for i in d] for d in dlist]

        return final

    def curve(self, x):
        if not self.lfc[0][0] < x < self.lfce[-1][0]:
            return 0
        for i in range(len(self.lfc)):
            if self.lfc[i][0] > x:
                return (self.lfc[i][1] - self.lfc[i - 1][1]) * (
                (x - self.lfc[i - 1][0]) / (self.lfc[i][0] - self.lfc[i - 1][0])) + self.lfc[i - 1][1]
