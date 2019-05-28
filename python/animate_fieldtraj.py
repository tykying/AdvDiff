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
fig, axarr = plt.subplots(2, 2)
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

  axarr[0,0].pcolormesh(X, Y, numpy.transpose(arr), cmap='Reds')
  axarr[0,0].set_title('t = ' + str(tt))
  axarr[0,0].set_xlim(x_logical[0], x_logical[-1])
  axarr[0,0].set_ylim(y_logical[0], y_logical[-1])
  axarr[0,0].set_aspect('equal')

  # Lower left panel: negative part of PDE solution
  arrm = arr
  arrm[arrm > 0.0] = 0.0
  axarr[1,0].set_xlabel(str(numpy.sum(arrm[:])))
  axarr[1,0].pcolormesh(X, Y, -numpy.transpose(arrm), cmap='Blues')
  axarr[1,0].set_xlim(x_logical[0], x_logical[-1])
  axarr[1,0].set_ylim(y_logical[0], y_logical[-1])
  axarr[1,0].set_aspect('equal')

  # Upper right panel: particle positions
  traj_data = data_folder + 'traj' + '_' + ts_in
  name, t_text, x, y, t = field_io.read_traj(traj_data)

  axarr[0,1].clear()
  axarr[0,1].scatter(x[0:-1:100,0],y[0:-1:100,0], marker='.')
  axarr[0,1].set_title('t = ' + str(t[0,0]))
  axarr[0,1].set_xlim(x_logical[0], x_logical[-1])
  axarr[0,1].set_ylim(y_logical[0], y_logical[-1])
  axarr[0,1].set_aspect('equal')

  # Lower right
  x_edges = numpy.linspace(1., nx+1, nx+1)  # Same resolution as bins
  y_edges = numpy.linspace(1., ny+1, ny+1)  # Same resolution as bins
  #x_edges = numpy.linspace(1., nx+1, (nx+1)/2+1)  # Same resolution as bins
  #y_edges = numpy.linspace(1., ny+1, (ny+1)/2+1)  # Same resolution as bins

  axarr[1,1].hist2d(x[:,0],y[:,0], bins=[x_edges, y_edges], cmap='Reds')
  axarr[1,1].set_xlim(x_logical[0], x_logical[-1])
  axarr[1,1].set_ylim(y_logical[0], y_logical[-1])
  axarr[1,1].set_aspect('equal')

ani = FuncAnimation(fig, animate_full, frames=range(0, 100, 1)) #, blit=True)
#animate_full(0)
plt.show()



## Debugging
#ts_in = '24'
#field_data = data_folder + field_in + '_' + ts_in
#name, t, arr = field_io.read_field(field_data)

#numpy.set_printoptions(precision=1)
#print(arr)
