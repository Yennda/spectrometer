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

# CRYSTAL
c = Crystal(d=0.384031, D=2.4, r=38, loc=[5.627693772801094, 0, 29.46742375572686])


# SOURCE
# s = st.generate_eliptical_source([0.02, 0.02], 70)
# s = [Source(loc=[0, 0, 0], wavelength=0.377207, intensity=1000, number=10000)]
s = [Source(loc=[0, 0, 0], wavelength=0.377207, intensity=1000, number=10000)]


# DETECTORS
d = Detector(dim=[1, 1], loc=[15.367355795595522+0.1, 0.2, -21.530864290851103], res=5)
# dm = Detector(dim=[1, 1], loc=[14.891626637551013, 0.0, -19.039877172681177], res=5)
# ds = Detector(dim=[1, 1], loc=[15.842832108508592, 0.0, -24.020527475108953], res=5)

st.display_source(s)

# SETUPS
setup = SetUp(source=s, crystal=c, detector=d)
setup.work()
# setup = SetUp(source=s, crystal=c, detector=dm)
# setup.work()
# setup = SetUp(source=s, crystal=c, detector=ds)
# setup.work()



st.display_crystal(c)

# adjustment(c, z_dist=30, oq=51.92)
print('Final elapsed time: {}s'.format(int(time.time() - t)))
