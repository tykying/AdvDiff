import sys

data_folder = './unittest/q64/TTG_sinusoidal/'
# sys.path.insert(0, '/home/s1046972/opt/qgm2/python')

import numpy
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib
matplotlib.use('TKAgg')

import field_io

## Four panels
fig, axarr = plt.subplots(1, 2, figsize=(12,6))
def animate_full(ts):
  # Upper left panel: PDE solution
  ts_in = str(ts)
  field_in = 'q'
  field_data = data_folder + field_in + '_' + ts_in

  name, tt, arr, glayer, type_id = field_io.read_field(field_data)
  nx, ny, = arr.shape

  x_logical = numpy.linspace(1., nx+1, nx+1)
  y_logical = numpy.linspace(1., ny+1, ny+1)

  X, Y = numpy.meshgrid(x_logical, y_logical)

  axarr[0].pcolormesh(X, Y, numpy.transpose(arr), cmap='Reds')
  axarr[0].set_xlabel('Integral = '+str(numpy.sum(arr[:])))
  axarr[0].set_title('t = ' + str(tt))
  axarr[0].set_xlim(x_logical[0], x_logical[-1])
  axarr[0].set_ylim(y_logical[0], y_logical[-1])
  axarr[0].set_aspect('equal')

  # Lower left panel: negative part of PDE solution
  arrm = arr
  arrm[arrm > 0.0] = 0.0
  axarr[1].set_xlabel('Integral of -ve parts = '+str(numpy.sum(arrm[:])))
  #arrm[arrm < -1e-10] = -1.0
  axarr[1].pcolormesh(X, Y, -numpy.transpose(arrm), cmap='Blues')
  axarr[1].set_xlim(x_logical[0], x_logical[-1])
  axarr[1].set_ylim(y_logical[0], y_logical[-1])
  axarr[1].set_aspect('equal')

ani = FuncAnimation(fig, animate_full, frames=range(0, 150, 10)) #, blit=True)
#animate_full(0)
plt.show()



## Debugging
#ts_in = '24'
#field_data = data_folder + field_in + '_' + ts_in
#name, t, arr = field_io.read_field(field_data)

#numpy.set_printoptions(precision=1)
#print(arr)
