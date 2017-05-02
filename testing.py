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


def display_ref(c):
    x0 = [p.loc[0] for p in c.points]
    y0 = [p.loc[1] for p in c.points]

    x = list()
    y = list()
    for p in c.points:
        if len(p.out) != 0:
            x.append(p.loc[0])
            y.append(p.loc[1])

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
    ax.scatter(x0, y0, linewidth=1, color='green')
    ax.scatter(x, y, linewidth=1, color='red')

    ax.grid(True)
    fig.savefig('display/graph_source.png', dpi=200, bbox_inches='tight')


print('Beginning')


#reference for editing
# c = Crystal(d=2.3, D=5, R=30, n=200)
# s = Source(loc=[0, 3, -30], wavelength=2.295)
# d = Detector(dim=[20, 20], loc=[0, 0, -30], res=1000)

c = Crystal(d=2.3, D=5, R=30, n=2000)
s = Source(loc=[0, 3, -30], wavelength=2.295)
d = Detector(dim=[10, 10], loc=[0, 0, 0], res=500)

setup = SetUp(source=s, crystal=c, detector=d)

# setup.compute_reflected()
setup.do()
# setup.compute_reflected()

display_ref(c)
display(setup)




# winsound.Beep(440, 1000)
