import sys

data_folder = './unittest/q64/TTG_sinusoidal/'
data_folder = './unittest/convergence/'
data_folder = './output/N4096_D64_I64/TTG_sinusoidal/h60d/'
# sys.path.insert(0, '/home/s1046972/opt/qgm2/python')

import matplotlib
#matplotlib.use('TKAgg')
import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable

import field_io



## Four panels
fig, axarr = plt.subplots(1, 2, figsize=(12,6))
div = make_axes_locatable(axarr[0])
cax = div.append_axes('right', '5%', '5%')

frames = []
t_frames = []
X_frames = []
Y_frames = []
factor = 1
nframes = 50
i_final = 2**12
i_final = 1024


#for ts in range(1, 3200*factor+1, 3200*factor/nframes):
for ts in range(i_final*factor, i_final*factor+1, 1):
    ts_in = str(ts)
    field_in = 'q_final'
    field_data = data_folder + field_in + '_' + ts_in
    #field_in = 'q_initial_0'
    #field_data = data_folder + field_in

    name, tt, arr, glayer, type_id = field_io.read_field(field_data)
    nx, ny, = arr.shape
    frames.append(numpy.transpose(arr))
    t_frames.append(tt)

    x_l = numpy.linspace(0., 1., nx+1)
    y_l = numpy.linspace(0., 1., ny+1)

    X, Y = numpy.meshgrid(x_l, y_l)
    X_frames.append(X)
    Y_frames.append(Y)

cv0 = frames[0]
cv0 = frames[-1]
vmax     = numpy.max(cv0)
vmin     = -numpy.max(cv0)
levels   = numpy.linspace(vmin, vmax, 200, endpoint = True)
cmap = plt.get_cmap('bwr')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
cf = axarr[0].pcolormesh(X_frames[0], Y_frames[0], cv0, norm=norm)
cb = fig.colorbar(cf, cax=cax)
tx = axarr[0].set_title('Time 0')
tx_r = axarr[1].set_title('Error 0')

def animate(i):
    arr = frames[i]
    t = t_frames[i]
    X = X_frames[i]
    Y = Y_frames[i]
    levels   = numpy.linspace(vmin, vmax, 200, endpoint = True)
    cf = axarr[0].pcolormesh(X, Y, arr, norm=norm)
    cax.cla()
    fig.colorbar(cf, cax=cax)
    tx.set_text('Time {0}'.format(t))

    nx, ny, = arr.shape
    x_e = numpy.linspace(0., X[0,-1], nx+1)
    y_e = numpy.linspace(0., Y[-1,0], ny+1)
    x_e = 0.5*(x_e[0:-1]+x_e[1:])
    y_e = 0.5*(y_e[0:-1]+y_e[1:])
    X_e, Y_e = numpy.meshgrid(x_e, y_e)

    q_e = 1.0 +  numpy.pi**2 * (4.0*X_e*(1.0-Y_e))**2 * (4.0*Y_e*(1.0-X_e))**3 * numpy.cos(numpy.pi*t)
    #q_e = (X_e*(1.0-Y_e))**2 * (Y_e*(1.0-X_e))**3
    #q_e = (numpy.sin(2.0*numpy.pi*X_e) * numpy.sin(numpy.pi*Y_e))**2

    cf_r = axarr[1].pcolormesh(X, Y, q_e, norm=norm)
    #cf_r = axarr[1].pcolormesh(X, Y, q_e-arr)

    dx = x_e[1]-x_e[0]
    dy = y_e[1]-y_e[0]
    err = numpy.sqrt(dx * dy * numpy.sum((arr-q_e)**2))
    tx_r.set_text('L2 Error {0}'.format(err))



ani = FuncAnimation(fig, animate, frames=nframes)

plt.show()
