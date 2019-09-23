import sys

data_folder = './unittest/q64/TTG_sinusoidal/'
data_folder = './unittest/q64/QGM2_L1/'
data_folder = './unittest/cellular/'
data_folder = './unittest/cellular/'
data_folder = './unittest/validation/N4096_D256_I256/QGM2_L1_NPART2704/h32d/steps/'
#data_folder = './unittest/meanflow/'
# sys.path.insert(0, '/home/s1046972/opt/qgm2/python')

import numpy
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import BoundaryNorm
import matplotlib
matplotlib.use('TKAgg')

import field_io

## Four panels
fig, axarr = plt.subplots(1, 2, figsize=(12,6))
def animate_full(ts):
  # Upper left panel: PDE solution
  ts_in = str(ts)
  field_in = 'q146'
  field_in = 'rev_q146'
  field_in = 'q'
  #field_in = 'psi_test'
  #field_in = 'psi_filtered'
  field_data = data_folder + field_in + '_' + ts_in


  name, tt, arr, glayer, type_id = field_io.read_field(field_data)
  nx, ny, = arr.shape

  q_max = numpy.amax(numpy.abs(arr))
  levels   = numpy.linspace(-q_max, q_max, 200, endpoint = True)
  cmap = plt.get_cmap('bwr')
  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

  x_logical = numpy.linspace(1., nx+1, nx+1)
  y_logical = numpy.linspace(1., ny+1, ny+1)

  X, Y = numpy.meshgrid(x_logical, y_logical)

  tt_day = tt/(3600*24)
  axarr[0].pcolormesh(X, Y, numpy.transpose(arr), cmap='Reds')
  #axarr[0].contourf(numpy.transpose(arr))
  axarr[0].set_xlabel('Integral = '+str(numpy.sum(arr[:])))
  axarr[0].set_title('t = ' + str(tt_day) + 'days')
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

#ani = FuncAnimation(fig, animate_full, frames=range(0, 256+1, 8)) #, blit=True)
ani = FuncAnimation(fig, animate_full, frames=range(0, 2400+1, 1)) #, blit=True)
#animate_full(0)
plt.show()



## Debugging
#ts_in = '24'
#field_data = data_folder + field_in + '_' + ts_in
#name, t, arr = field_io.read_field(field_data)

#numpy.set_printoptions(precision=1)
#print(arr)
