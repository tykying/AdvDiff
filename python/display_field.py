import sys

data_folder = './output/'
data_folder = './unittest/'
# sys.path.insert(0, '/home/s1046972/opt/qgm2/python')

from matplotlib import pyplot as plt
import numpy

import field_io

field_in = 'q'
for ts in range(0, 2, 1):
  ts_in = str(ts)

  field_data = data_folder + field_in + '_' + ts_in
  
  print(field_data)
  name, t, arr, glayer = field_io.read_field(field_data)

  nx, ny, = arr.shape

  x_logical = numpy.linspace(0., 1., nx)
  y_logical = numpy.linspace(0., 1., ny) 
  X, Y = numpy.meshgrid(x_logical, y_logical)
  
  # TODO: HARDCODED size of domain
  assert(nx == ny)
  ngrid = nx
  L = 3840.0E5
  nd_scale = L/ngrid
  
  x = x_logical*nd_scale
  y = y_logical*nd_scale
  t_physical = t*nd_scale

  X, Y = numpy.meshgrid(x, y)

  plt.pcolormesh(X*L/1E5, Y*L/1E5, numpy.transpose(arr))
  plt.title('time (in years) = ' + str(t_physical/(3600.0*24.0*365.25)))
  plt.xlabel('x (in km)')
  plt.ylabel('y (in km) ')
    
  fig_path = '../QGFIG_OUTPUT/' + field_in + '_' + ts_in + '.png'
  
  # plt.savefig(fig_path, format='png') 

  plt.show()
