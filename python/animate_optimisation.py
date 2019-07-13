import sys

show_id = 20
#data_folder = './output/output16/K_TG_iso/'
#data_folder = './output/output16_ptwiseprop/K_TG_iso/'
data_folder = './output/output16/K_sinusoidal/'
data_folder = './output/output16_cornerK/K_sinusoidal/'
data_folder = './output/output16_ptwiseprop/K_sinusoidal/'
#data_folder = './output/output32/K_sinusoidal/'
#data_folder = './output/output32_ptwiseprop/K_sinusoidal/'
#data_folder = './output/output16/TTG_const/'
#data_folder = './output/output16/K_const/'
#data_folder = './output/output16ref/K_const/'
#data_folder = './output/output16_isoK/K_const/'
#data_folder = './output/output32/K_const/'
#data_folder = './output/output32_cornerK/K_const/'
#data_folder = './output/output16/K_sinusoidal/'
data_folder = './output/output16/TTG_sinusoidal/'
#data_folder = './output/output16/K_sinusoidal/'
#data_folder = './output/output16/TTG_sinusoidal_h6h/'
#data_folder = './output/output16/K_sinusoidal/'
#data_folder = './output/output16/K_sinusoidal_wor/'
#data_folder = './output/output16/K_sinusoidal_wM/'
#data_folder = './output/output32/QGM2_L1_NPART676/'
#data_folder = './output/output16/QGM2_L1_NPART676_iso/'
#data_folder = './output/output16/QGM2_L1_NPART676/'
#data_folder = './output/output16/QGM2_L1_NPART2704_iso/h64d/'
#data_folder = './output/output48/QGM2_L1_NPART2704_iso/h64d/'
data_folder = './output/N4096_D256_I256/QGM2_L1_NPART2704/h64d/'
# sys.path.insert(0, '/home/s1046972/opt/qgm2/python')

import numpy
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from numpy import pi, sin, cos, sqrt, power, exp

matplotlib.use('TKAgg')

import field_io

def visualise_field(ax, field_in, ts, max_level, resc_unit):
  pcolor1contour0 = 1
  levels = MaxNLocator(nbins=15).tick_values(-1.0*max_level/resc_unit, 1.0*max_level/resc_unit)
  #cmap = plt.get_cmap('Reds')
  cmap = plt.get_cmap('bwr')
  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

  ts_in = str(ts)

  field_data = data_folder + field_in + '_' + ts_in

  name, tt, arr, glayer, type_id = field_io.read_field(field_data)
  nx, ny, = arr.shape

  xl = numpy.linspace(0., 1., nx+1)
  yl = numpy.linspace(0., 1., ny+1)

  if (type_id == 1):
    # Values defined at corners
    xl_odd = numpy.linspace(0.-0.5/(nx-1), 1.+0.5/(nx-1), nx+1)
    yl_odd = numpy.linspace(0.-0.5/(ny-1), 1.+0.5/(ny-1), ny+1)
    if (pcolor1contour0 == 1):
      X, Y = numpy.meshgrid(xl_odd, yl_odd)
      #pcm = ax.pcolormesh(X, Y, (numpy.transpose(arr)-resc_unit)/resc_unit, cmap=cmap, norm=norm)
      pcm = ax.pcolormesh(X, Y, (numpy.transpose(arr))/resc_unit, cmap=cmap, norm=norm)
    else:
      X, Y = numpy.meshgrid(0.5*(xl_odd[0:-1]+xl_odd[1:]), 0.5*(yl_odd[0:-1]+yl_odd[1:]))
      #pcm = ax.contourf(X, Y, (numpy.transpose(arr)-resc_unit)/resc_unit, cmap=cmap, levels=levels)
      pcm = ax.contourf(X, Y, (numpy.transpose(arr))/resc_unit, cmap=cmap, levels=levels)
  else:
    if (pcolor1contour0 == 1):
      X, Y = numpy.meshgrid(xl, yl)
      #pcm = ax.pcolormesh(X, Y, (numpy.transpose(arr)-resc_unit)/resc_unit, cmap=cmap, norm=norm)
      pcm = ax.pcolormesh(X, Y, (numpy.transpose(arr))/resc_unit, cmap=cmap, norm=norm)
    else:
      X, Y = numpy.meshgrid(0.5*(xl[0:-1]+xl[1:]), 0.5*(yl[0:-1]+yl[1:]))
      #pcm = ax.contourf(X, Y, (numpy.transpose(arr)-resc_unit)/resc_unit, cmap=cmap, levels=levels)
      pcm = ax.contourf(X, Y, (numpy.transpose(arr))/resc_unit, cmap=cmap, levels=levels)

  ax.set_title(field_in + ': logPost = ' + str(tt))
  #ax.set_title(r'$\'+field_in+'$')
  ax.set_title(field_in)
  ax.set_xlim(xl[0], xl[-1])
  ax.set_ylim(yl[0], yl[-1])
  ax.set_aspect('equal')
  plt.axis('off')
  plt.colorbar(pcm, ax=ax)
  tight_layout()

  # theo_logLik = -386479.35214807594
  # init_logLik = -395073.95702572557
  # print('initial = ', init_logLik, 'Exact = ', theo_logLik, 'MAP = ', tt)

def animate_full(ts):
  # Upper left panel: Exact solution
  Pi = pi
  L = 3840.0*1000
  kappa_scale = 10000
  #kappa_scale = 5000
  psi_scale = 0.05*L
  kappa_unit = 10000
  psi_unit = 1

  xl = numpy.linspace(0., 1., 64+1)
  yl = numpy.linspace(0., 1., 64+1)

  X, Y = numpy.meshgrid(xl, yl)

  psi = 0*X

  if (data_folder.find('const') > -1):
    sigma1 = 0.*X + kappa_scale
    sigma2 = 0.*Y + kappa_scale
    phi = 0.*Y + 0.


  if (data_folder.find('TG_iso') > -1):
    Lx1 = pi; Ly1 = pi;
    Lx2 = pi; Ly2 = pi;

    sigma1 = sqrt(sin(Lx1*X) * sin(Ly1*Y));
    sigma2 = sqrt(sin(Lx2*X) * sin(Ly2*Y));
    phi = 0.*Y + 0.

  if (data_folder.find('sinusoidal') > -1):
    Lx1 = pi/2; Ly1 = -pi;
    Lx2 = pi/3; Ly2 = pi/2;
    Lxp = pi/3; Lyp = pi/4;

    sigma1 = sin(Lx1*X + Ly1*Y)
    sigma2 = cos(Lx2*X - Ly2*Y)
    phi = (pi/2)*sin(Lxp*X)*sin(Lyp*Y)

  if (data_folder.find('TTG') > -1):
    a = 3; b = -1;
    Lx = pi; Ly = 2*pi;
    k = -0.5;
    #psi_scale = 0.0521
    psi = psi_scale*exp(k*(a*X+b*Y)/(a+b))*sin(Lx*X) * sin(Ly*Y);

  if (data_folder.find('QGM2') > -1):
      sigma1 = 0.*X + kappa_scale
      sigma2 = 0.*Y + kappa_scale
      phi = 0.*Y + 0.

  Kxx = kappa_scale * (power(cos(phi) * sigma1, 2) + power(sin(phi) * sigma2, 2))
  Kyy = kappa_scale * (power(sin(phi) * sigma1, 2) + power(cos(phi) * sigma2, 2))
  Kxy = kappa_scale * (cos(phi) * sin(phi) * (power(sigma1,2) - power(sigma2, 2)) )

  Klevels = MaxNLocator(nbins=15).tick_values(-1.25*kappa_scale/kappa_unit, 1.25*kappa_scale/kappa_unit)
  psilevels = MaxNLocator(nbins=15).tick_values(-1.25*psi_scale/psi_unit, 1.25*psi_scale/psi_unit)
  #cmap = plt.get_cmap('Reds')
  cmap = plt.get_cmap('bwr')
  Knorm = BoundaryNorm(Klevels, ncolors=cmap.N, clip=True)
  psinorm = BoundaryNorm(psilevels, ncolors=cmap.N, clip=True)

  ## Two panels
  fig, axarr = plt.subplots(2, 4,figsize=(16,8))

  ax = axarr[0,0]
  # cf = ax.contourf(X, Y, Kxx/kappa_unit, levels=Klevels, cmap=cmap)
  cf = ax.pcolormesh(X, Y, Kxx/kappa_unit, cmap=cmap, norm=Knorm)
  ax.set_title('Exact')
  ax.set_xlim(0, 1)
  ax.set_ylim(0, 1)
  ax.set_aspect('equal')
  plt.colorbar(cf, ax=ax)

  ax = axarr[0,1]
  # cf = ax.contourf(X, Y, Kyy/kappa_unit, levels=Klevels, cmap=cmap)
  cf = ax.pcolormesh(X, Y, Kyy/kappa_unit, cmap=cmap, norm=Knorm)
  ax.set_title('Exact')
  ax.set_xlim(0, 1)
  ax.set_ylim(0, 1)
  ax.set_aspect('equal')
  plt.colorbar(cf, ax=ax)

  ax = axarr[0,2]
  # cf = ax.contourf(X, Y, Kxy/kappa_unit, levels=Klevels, cmap=cmap)
  cf = ax.pcolormesh(X, Y, Kxy/kappa_unit, cmap=cmap, norm=Knorm)
  ax.set_title('Exact')
  ax.set_xlim(0, 1)
  ax.set_ylim(0, 1)
  ax.set_aspect('equal')
  plt.colorbar(cf, ax=ax)

  ax = axarr[0,3]
  # cf = ax.contourf(X, Y, psi/psi_unit, levels=psilevels, cmap=cmap)
  cf = ax.pcolormesh(X, Y, psi/psi_unit, cmap=cmap, norm=psinorm)
  ax.set_title('Exact')
  ax.set_xlim(0, 1)
  ax.set_ylim(0, 1)
  ax.set_aspect('equal')
  plt.colorbar(cf, ax=ax)

  # Bottom panels
  visualise_field(axarr[1,0], 'K11', ts, kappa_scale, kappa_unit)
  visualise_field(axarr[1,1], 'K22', ts, kappa_scale, kappa_unit)
  visualise_field(axarr[1,2], 'K12', ts, kappa_scale, kappa_unit)
  visualise_field(axarr[1,3], 'psi', ts, psi_scale, psi_unit)

  if (data_folder.find('QGM2') > -1):
      fig, axarr = plt.subplots(1, 4,figsize=(16,4))
      # Bottom panels
      visualise_field(axarr[0], 'K11', ts, kappa_scale, kappa_unit)
      visualise_field(axarr[1], 'K22', ts, kappa_scale, kappa_unit)
      visualise_field(axarr[2], 'K12', ts, kappa_scale, kappa_unit)
      visualise_field(axarr[3], 'psi', ts, psi_scale, psi_unit)


#ani = FuncAnimation(fig, animate_full, frames=range(0, 15+1, 1)) #, blit=True)
animate_full(show_id)
plt.show()
