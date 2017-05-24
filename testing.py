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

    # x = [p.loc[0] for p in c.points]
    # y = [p.loc[1] for p in c.points]
    # x.append(s.loc[0])
    # y.append(s.loc[1])
    # for r in d.mesh:
    #     for p in r:
    #         x.append(p.loc[0])
    #         y.append(p.loc[1])
    # fig = plt.figure(figsize=(6, 6))
    # ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
    # ax.scatter(x, y, linewidth=2, color='green')
    # ax.grid(True)
    # fig.savefig('display/graph_xy.png', dpi=200, bbox_inches='tight')
    #
    # x = [p.loc[1] for p in c.points]
    # y = [p.loc[2] for p in c.points]
    # x.append(s.loc[1])
    # y.append(s.loc[2])
    # for r in d.mesh:
    #     for p in r:
    #         x.append(p.loc[1])
    #         y.append(p.loc[2])
    # fig = plt.figure(figsize=(6, 6))
    # ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
    # ax.scatter(x, y, linewidth=2, color='green')
    # ax.grid(True)
    # fig.savefig('display/graph_yz.png', dpi=200, bbox_inches='tight')


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
    x0 = [p[0] for p in c.points]
    y0 = [p[1] for p in c.points]

    x = list()
    y = list()
    # for p in c.points:
    #     if len(p.ray_out) != 0:
    #         x.append(p.loc[0])
    #         y.append(p.loc[1])

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
    # ax.scatter(x, y, linewidth=1, color='red')
    ax.set_xlim([c.loc[0] - c.D, c.loc[0] + c.D])
    ax.set_ylim([c.loc[1] - c.D, c.loc[1] + c.D])
    ax.grid(True)
    fig.savefig('display/graph_crystal.png', dpi=200, bbox_inches='tight')


def adjustment(c: Crystal, z_dist, oq=0):
    curveSi = Curves('reflectivity.csv')
    bragg = curveSi.bragg

    loc_cr = [z_dist * m.cos(bragg), 0, m.sin(bragg) * z_dist]
    q = tl.imaging_equation(f=c.r / 2, a=la.norm(loc_cr))
    qm = z_dist * c.r * m.sin(bragg) / (2 * z_dist - c.r * m.sin(bragg))
    qs = z_dist * c.r / (2 * z_dist * m.sin(bragg) - c.r)

    loc_det = la.plus(loc_cr, la.x(la.normalize([z_dist / m.tan(bragg), 0, -z_dist]), q))
    loc_det_m = la.plus(loc_cr, la.x(la.normalize([z_dist / m.tan(bragg), 0, -z_dist]), qm))
    loc_det_s = la.plus(loc_cr, la.x(la.normalize([z_dist / m.tan(bragg), 0, -z_dist]), qs))
    loc_det_o = la.plus(loc_cr, la.x(la.normalize([z_dist / m.tan(bragg), 0, -z_dist]), oq))

    print(la.norm(la.minus(loc_cr, loc_det)))

    print('q={}'.format(q))
    print('qm={}'.format(qm))
    print('qs={}'.format(qs))

    print('Crystal: {}'.format(loc_cr))
    print('Detector: {}'.format(loc_det))
    print('Detector_m: {}'.format(loc_det_m))
    print('Detector_s: {}'.format(loc_det_s))
    print('Detector_o: {}'.format(loc_det_o))


def generate_eliptical_source(dim: list, n: int):
    s = []
    for i in range(n):
        for j in range(n):
            if ((i - n / 2) ** 2 + (j - n / 2) ** 2) < ((n - 1) / 2) ** 2:
                s.append(
                    Source(loc=[i / n * dim[0], j / n * dim[1], 0], wavelength=0.377207, intensity=1000, number=200))
    print(len(s))
    return s


def generate_square_source(dim: list, n: int):
    s = []
    for i in range(n):
        for j in range(n):
            s.append(
                Source(loc=[i / n * dim[0], j / n * dim[1], 0], wavelength=0.377207, intensity=1000, number=200))
    return s


t = time.time()

# x = [[m.cos(f / 10000 * m.pi * 2), m.sin(f / 10000 * m.pi * 2)] for f in range(10000)]
# x0 = [s[0] for s in x]
# y0 = [s[1] for s in x]
#
# fig = plt.figure(figsize=(6, 6))
# ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
#
# ax.plot(x0, y0, linewidth=1, color='green')
# # ax.set_xlim([c.loc[0] - c.D, c.loc[0] + c.D])
# # ax.set_ylim([c.loc[1] - c.D, c.loc[1] + c.D])
# ax.grid(True)
# fig.savefig('display/graph.png', dpi=200, bbox_inches='tight')


# print('Shoelace: {}'.format(tl.shoelace(x)))
c = Crystal(d=0.384031, D=2.4, r=38, loc=[5.627693772801094, 0, 29.46742375572686])

# s = generate_eliptical_source([0.05, 0.05], 100)
s = [Source(loc=[0, 0, 0], wavelength=0.377207, intensity=1000, number=10000)]

d = Detector(dim=[1, 1], loc=[15.367355795595522, 0.0, -23.3], res=5)
dm = Detector(dim=[1, 1], loc=[14.891626637551013, 0.0, -19.039877172681177], res=5)
ds = Detector(dim=[1, 1], loc=[15.842832108508592, 0.0, -24.020527475108953], res=5)

display_source(s)

setup = SetUp(source=s, crystal=c, detector=d)
setup.work()
setup = SetUp(source=s, crystal=c, detector=dm)
setup.work()
setup = SetUp(source=s, crystal=c, detector=ds)
setup.work()


display_crystal(c)
adjustment(c, z_dist=30, oq=51.92)
print('Final elapsed time: {}s'.format(int(time.time() - t)))
