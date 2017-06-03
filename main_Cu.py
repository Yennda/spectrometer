from tracing import *

from tools_setup import SetupTools as st
from tools import Tools as tl
import numpy as np
import math as m

import time

t = time.time()

c = Crystal(d=1.5414, D=2.4, r=38, loc=[5.627693772801094, 0, 29.46742375572686], rc=2)
s = Source(loc=[0, 0, 0], wavelength=tl.ang_from_kev(8.0478), intensity=1000, number=10000)
d = Detector(dim=[1, 1], loc=[15.367229373029803, 0.0, -21.530202323895065], res=5)
print(c.d)
print(s.wl)
print(tl.deg_from_rad(c.bragg_angle(s.wl, 2)))

# st.display_source(s)

# SETUPS
# setup = SetUp(source=s, crystal=c, detector=d)
# setup.work()

# st.display_crystal(c)

# st.adjustment(c, z_dist=30, oq=51.92)
print('Final elapsed time: {}s'.format(int(time.time() - t)))
