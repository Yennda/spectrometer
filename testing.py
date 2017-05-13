from matplotlib.font_manager import FontProperties
import IPython
import matplotlib.pyplot as plt
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


def simple_point_source(c: Crystal, s: Source, d: Detector, z_dist):
    setup = SetUp(source=s, crystal=c, detector=d)
    bragg = setup.curveSi.bragg
    # print(bragg)

    c.loc = la.plus(s.loc, [z_dist / m.tan(bragg), 0, z_dist])
    c.loc_centre = la.minus(c.loc, [0, 0, c.r])

    det_dist = setup.imaging_equation()
    print(det_dist)
    d.loc = la.plus(c.loc, la.x(la.normalize([z_dist / m.tan(bragg), 0, -z_dist]), det_dist))
    print(la.norm(la.minus(c.loc, d.loc)))

    print('Crystal: {}'.format(c.loc))
    print('Detector: {}'.format(d.loc))
    return SetUp(source=s, crystal=c, detector=d)


t = time.time()

# c = Crystal(d=3.84031 / 2, D=5, r=30, loc=[0, 15.4, 80])
# s = Source(loc=[0, 0, 0], wavelength=1.91, intensity=1000, number=1000)
# d = Detector(dim=[5, 5], loc=[0, 18.9, 80 - 18.72], res=50)
#
# setup = simple_point_source(c, s, d, 30)
# setup.work()

c = Crystal(d=1.920155, D=5, r=30, loc=[5.729405277623609, 0, 30])
s = Source(loc=[0, 0, 0], wavelength=1.91, intensity=1000, number=1000)
# d = Detector(dim=[1, 1], loc=[11.258935829952932, 0.0, 1.0465731552011377], res=20)
d = Detector(dim=[1, 2], loc=[11.258935829952932, 0.0, 2], res=20)
setup = SetUp(source=s, crystal=c, detector=d)

setup.work()
# for i in range(-5,20):
#     d = Detector(dim=[1, 1], loc=[11.258935829952932, 0.0, i/10], res=20)
#     setup = SetUp(source=s, crystal=c, detector=d)
#     setup.work()


# display_ref(c)
# display(setup)

print('Elapsed time: {}s'.format(int(time.time()-t)))
