from matplotlib.font_manager import FontProperties
import IPython
import matplotlib.pyplot as plt

import numpy as np
import math as m
import copy

import sys
import time

sys.path.append('C:/Users/jabuk/PycharmProjects')
from textable.datalist import DataList
from textable.main import TexTable

from tracing import *

t = time.time()

c = Crystal(d=2.3, D=5, R=30, n=2000)
s = Source(loc=[0, 5, -50], wavelength=2.289)


d = Detector(dim=[0.3, 0.3], loc=[0, 0, 0], res=1000)
d.generate_mesh()
print(d.mesh[1][0].loc)
print('-------------------------')
d.translate(np.matrix([1,0,0]))
print('-------------------------')

print(d.mesh[1][0].loc)
print('########')
print(np.matrix([1,0,0])+np.matrix([1,0,0]))

# d.rotate([14, 24, 45], u='d')


# setup = SetUp(source=s, crystal=c, detector=d)
# setup.compute_reflected()

print(time.time() - t)
