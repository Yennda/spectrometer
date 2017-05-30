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

c = Crystal(d=0.384031, D=[3,3], l=5, tb=tl.rad_from_deg(-3.4858))
s = Source(loc=[0, 0, 0], wavelength=0.377207, intensity=1000, number=100)
d = Detector(dim=[30, 30], loc=[0,8.12, 75.509], res=50)
d.rotate([0,tl.rad_from_deg(2*3.4858),0])

# st.display_source(s)
setup = SetUp(source=s, crystal=c, detector=d)
setup.work()


x0 = [p[0] for p in c.points]
y0 = [p[1] for p in c.points]

fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])

ax.scatter(x0, y0, linewidth=1, color='green')

ax.grid(True)
fig.savefig('display/graph_sp_crystal.png', dpi=200, bbox_inches='tight')

# adjustment(c, z_dist=30, oq=51.92)
print('Final elapsed time: {}s'.format(int(time.time() - t)))
