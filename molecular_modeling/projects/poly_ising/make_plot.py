import numpy as np
from matplotlib import pylab as py

A = np.loadtxt('energia')
kt = np.linspace(1, 8, 20)

sr = np.apply_along_axis( np.mean, 1, A )
odch = np.apply_along_axis( np.std, 1, A)

py.plot( kt, sr, 'b.')
py.show()
py.plot( kt, odch, 'b.')
py.show()
