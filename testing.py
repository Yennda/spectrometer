from matplotlib.font_manager import FontProperties
import IPython
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import numpy as np
import math as m

import time
import timeit

from tracing import *


def display(setup: SetUp):
    c = setup.crystal
    s = setup.source
    d = setup.detector

    x = [p.loc[0] for p in c.points]
    y = [p.loc[2] for p in c.points]

    x.append(s.loc[0])
    y.append(s.loc[2])

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
    x.append(s.loc[0])
    y.append(s.loc[1])
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
    x.append(s.loc[1])
    y.append(s.loc[2])
    for r in d.mesh:
        for p in r:
            x.append(p.loc[1])
            y.append(p.loc[2])
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
    ax.scatter(x, y, linewidth=2, color='green')
    ax.grid(True)
    fig.savefig('display/graph_yz.png', dpi=200, bbox_inches='tight')


def display_source(source: list):
    x0 = [s.loc[0] for s in source]
    y0 = [s.loc[1] for s in source]

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])

    ax.scatter(x0, y0, linewidth=1, color='green')
    # ax.set_xlim([c.loc[0] - c.D, c.loc[0] + c.D])
    # ax.set_ylim([c.loc[1] - c.D, c.loc[1] + c.D])
    ax.grid(True)
    fig.savefig('display/graph_source.png', dpi=200, bbox_inches='tight')


def display_crystal(c: Crystal):
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

    ax.add_patch(
        patches.Circle(
            (c.loc[0], c.loc[1]),
            c.D / 2,
            fill=False
        )
    )

    ax.scatter(x0, y0, linewidth=1, color='green')
    ax.scatter(x, y, linewidth=1, color='red')
    ax.set_xlim([c.loc[0] - c.D, c.loc[0] + c.D])
    ax.set_ylim([c.loc[1] - c.D, c.loc[1] + c.D])
    ax.grid(True)
    fig.savefig('display/graph_crystal.png', dpi=200, bbox_inches='tight')


def simple_point_source(c: Crystal, s: Source, d: Detector, z_dist):
    curveSi = Curves('reflectivity.csv')
    bragg = curveSi.bragg

    loc_cr = la.plus(s.loc, [z_dist / m.tan(bragg), 0, z_dist])

    det_dist = tl.imaging_equation(f=c.r / 2, a=la.norm(la.minus(loc_cr, s.loc)))

    loc_det = la.plus(loc_cr, la.x(la.normalize([z_dist / m.tan(bragg), 0, -z_dist]), det_dist))
    print(det_dist)
    print(la.norm(la.minus(loc_cr, loc_det)))

    print('Crystal: {}'.format(loc_cr))
    print('Detector: {}'.format(loc_det))


def generate_eliptical_source(dim: list, n: int):
    s = []
    for i in range(n):
        for j in range(n):
            if ((i - n / 2) ** 2 + (j - n / 2) ** 2) < ((n-1) / 2) ** 2:
                s.append(
                    Source(loc=[i / n * dim[0], j / n * dim[1], 0], wavelength=0.377207, intensity=1000, number=20))
    return s


t = time.time()

c = Crystal(d=0.384031, D=5, r=30, loc=[7.639207036831479, 0, 40])
s = generate_eliptical_source([0.07, 0.07], 30)
d = Detector(dim=[1, 1], loc=[12.09391264740393, 0.0, 15.82], res=5)


display_source(s)

# s = [Source(loc=[0, 0, 0], wavelength=0.377207, intensity=1000, number=200)]

setup = SetUp(source=s, crystal=c, detector=d)
setup.work()


# simple_point_source(s=s, c=c, d=d, z_dist=40)
# display_crystal(c)
# display(setup)

print('Elapsed time: {}s'.format(int(time.time() - t)))
