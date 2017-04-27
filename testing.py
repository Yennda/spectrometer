from matplotlib.font_manager import FontProperties
import IPython
import matplotlib.pyplot as plt
import numpy as np
import math as m

import time
import timeit

from tracing import *
import winsound

t = time.time()
print('Beginning')

c = Crystal(d=2.3, D=5, R=30, n=40)
s = Source(loc=[0, 5, -50], wavelength=2.289)
d = Detector(dim=[80, 80], loc=[20, 10, -10], res=10000)

setup = SetUp(source=s, crystal=c, detector=d)
print('reflection')
setup.compute_reflected()
print('detector')
setup.intensity_for_detector()
setup.graph()

print(time.time() - t)
# winsound.Beep(440, 1000)
