from matplotlib.font_manager import FontProperties
import IPython
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from tracing import *

from tools_setup import SetupTools as st
import numpy as np
import math as m

import time
import timeit

t = time.time()

c = Crystal(d=9.406767935906574, D=[3, 3], l=7, tb=tl.rad_from_deg(3.4858), rc=1)
s = Source(loc=[0, 0, 0], wavelength=[1.143892288781112], intensity=1e9, number=1000)
d = Detector(dim=[0.08, 2.4], loc=st.distances(D=3, d=2.4, l=7, tb=tl.rad_from_deg(3.4858)), res=800)
d.rotate([tl.rad_from_deg(2 * 3.4858), 0, 0])


print(st.distances(D=3, d=2.4, l=7, tb=tl.rad_from_deg(3.4858)))

setup = SetUp(source=s, crystal=c, detector=d)
setup.work()

print('count {}'.format(s.photons_total))

x0 = [p[1] for p in c.points]
y0 = [p[2] for p in c.points]

fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])

ax.scatter(x0, y0, linewidth=1, color='green')

ax.grid(True)
fig.savefig('display/graph_sp_crystal.png', dpi=200, bbox_inches='tight')

print('Final elapsed time: {}s'.format(int(time.time() - t)))
3.148781422795631e-09