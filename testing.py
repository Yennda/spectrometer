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

c = Crystal(d=2.3, D=5, R=30, n=200)
s = Source(loc=[0, 5, -50], wavelength=2.289)
# d=Detector(dim=[])
setup = SetUp(source=s, crystal=c)
setup.compute_reflected()




print(time.time() - t)
