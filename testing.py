from matplotlib.font_manager import FontProperties
import IPython
import matplotlib.pyplot as plt
import numpy as np
import math as m

import time
import timeit

from tracing import *
import winsound


def display(setup: SetUp):
    c = setup.crystal
    s = setup.source
    d = setup.detector

    x = [p.loc[0] for p in c.points]
    y = [p.loc[2] for p in c.points]

    x.append(0)
    y.append(0)

    for r in d.mesh:
        for p in r:
            x.append(p.loc[0])
            y.append(p.loc[2])

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
    ax.scatter(x, y, linewidth=2, color='green')
    ax.grid(True)
    fig.savefig('display/graph_xz.png', dpi=200, bbox_inches='tight')

    x = [p.loc[0] for p in c.points]
    y = [p.loc[1] for p in c.points]
    x.append(0)
    y.append(0)
    for r in d.mesh:
        for p in r:
            x.append(p.loc[0])
            y.append(p.loc[1])
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
    ax.scatter(x, y, linewidth=2, color='green')
    ax.grid(True)
    fig.savefig('display/graph_xy.png', dpi=200, bbox_inches='tight')

    x = [p.loc[1] for p in c.points]
    y = [p.loc[2] for p in c.points]
    x.append(0)
    y.append(0)
    for r in d.mesh:
        for p in r:
            x.append(p.loc[1])
            y.append(p.loc[2])
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
    ax.scatter(x, y, linewidth=2, color='green')
    ax.grid(True)
    fig.savefig('display/graph_yz.png', dpi=200, bbox_inches='tight')


def display_ref(c: Crystal):
    x0 = [p.loc[0] for p in c.points]
    y0 = [p.loc[1] for p in c.points]

    x = list()
    y = list()
    for p in c.points:
        if len(p.ray_out) != 0:
            x.append(p.loc[0])
            y.append(p.loc[1])

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
    ax.scatter(x0, y0, linewidth=1, color='green')
    ax.scatter(x, y, linewidth=1, color='red')
    ax.set_xlim([c.loc[0] - c.D, c.loc[0] + c.D])
    ax.set_ylim([c.loc[1] - c.D, c.loc[1] + c.D])
    ax.grid(True)
    fig.savefig('display/graph_crystal.png', dpi=200, bbox_inches='tight')


# c = Crystal(d=3.84031/2, D=5, r=30, loc=[0, 5, 80])
# s = Source(wavelength=1.91, intensity=1000, number=20)
# d = Detector(dim=[2, 10], loc=[0, 8.2, 80 - 20], res=200)
# d.rotate([0, 60, 0], 'd')

c = Crystal(d=2.3, D=5, r=30, loc=[0, 5, 80])
s = Source(wavelength=2.290, intensity=1000, number=20000)
d = Detector(dim=[1, 1], loc=[0, 5, 80 - 20], res=20)
print('detector prepared')
setup = SetUp(source=s, crystal=c, detector=d)

setup.work()

display_ref(c)
display(setup)

# winsound.Beep(440, 1000)
