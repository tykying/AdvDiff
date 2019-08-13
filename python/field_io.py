#!/usr/bin/env python2
# -*- coding: UTF-8 -*-

from __future__ import division, print_function, unicode_literals

import array
import bz2
import gzip
import os
import pdb

import numpy
try:
  import lzma
  _lzma = True
except ImportError:
  _lzma = False

__all__ = [
    "IOException",
    "InvalidArgumentException",
    "open_compressed",
    "read_field"
  ]

class IOException(IOError):
  pass

class InvalidArgumentException(Exception):
  pass

def open_compressed(filename, binary = True):
  """
  Open a (possibly compressed) file for reading. Searches for [filename].gz,
  [filename].bz2, [filename].xz and [filename]. Raises an IOException if none or
  more than one such file is found.
  """

  if binary:
    mode = "rb"
  else:
    mode = "r"
  files = []
  try:
    files.append(gzip.GzipFile("%s.gz" % filename, mode))
  except IOError:
    pass
  try:
    files.append(bz2.BZ2File("%s.bz2" % filename, mode))
  except IOError:
    pass
  if _lzma:
    try:
      files.append(lzma.LZMAFile("%s.xz" % filename, mode))
    except IOError:
      pass
  else:
    if os.path.exists("%s.xz" % filename):
      raise IOException("xz not available")
  try:
    files.append(open(filename, mode))
  except IOError:
    pass

  if not len(files) == 1:
    raise IOException("Failed to open file '%s': %i found" % (filename, len(files)))

  return files[0]

def read_traj(filename):
  """
  Read a raw field
  """

  handle = open("%s.hdr" % filename, "r")
  line = handle.readline()
  if not line.strip() == "serial":
    raise IOException("Invalid header")
  name = handle.readline().strip()
  m, n = map(int, handle.readline().split())
  dim, = map(int, handle.readline().split())
  t_text, = map(float, handle.readline().split())
  if not len(handle.readline()) == 0:
    raise IOException("Invalid header")
  handle.close()

  handle = open_compressed("%s_x.dat" % filename, binary = True)
  x = array.array(b"d")
  x.fromstring(handle.read())
  x = numpy.array(x)
  x.shape = (n, m)
  handle.close()

  handle = open_compressed("%s_y.dat" % filename, binary = True)
  y = array.array(b"d")
  y.fromstring(handle.read())
  y = numpy.array(y)
  y.shape = (n, m)
  handle.close()

  handle = open_compressed("%s_t.dat" % filename, binary = True)
  t = array.array(b"d")
  t.fromstring(handle.read())
  t = numpy.array(t)
  t.shape = (n, m)
  handle.close()

  return name, t_text, x.transpose([1, 0]), y.transpose([1, 0]), t.transpose([1, 0])

def read_field(filename):
  """
  Read a raw field
  """

  handle = open("%s.hdr" % filename, "r")
  line = handle.readline()
  if not line.strip() == "serial":
    raise IOException("Invalid header")
  name = handle.readline().strip()
  m, n = map(int, handle.readline().split())
  glayer, = map(int, handle.readline().split())
  type_id, = map(int, handle.readline().split())
  t, = map(float, handle.readline().split())
  if not len(handle.readline()) == 0:
    raise IOException("Invalid header")
  handle.close()
  
  if (type_id == 2):
    type_id = 0

  handle = open_compressed("%s.dat" % filename, binary = True)
  arr = array.array(b"d")
  arr.fromstring(handle.read())
  arr = numpy.array(arr)
  arr.shape = (n+type_id+2*glayer, m+type_id+2*glayer)
  arr = arr[glayer:n+type_id+glayer, glayer:m+type_id+glayer]
  handle.close()

  return name, t, arr.transpose([1, 0]), glayer, type_id

def read_energy(filename):
  """
  Parse a log file for energy diagnostics.

  Arguments:
    filename: Log filename

  Returns:
    A tuple of arrays (t, ke, pe)
  """

  f = open_compressed(filename, binary = False)

  t, ke, pe = [], [], []
  layers = None
  line = f.readline()
  while len(line) > 0:
    if line.startswith("time (days) = "):
      lt = float(line.split()[-1])
      line = f.readline()

      lke = []
      l = 1
      while line.startswith("Layer %i KE = " % l):
        lke.append(float(line.split()[-1]))
        line = f.readline()
        l += 1
      if layers is None:
        layers = l - 1
      elif not layers == l - 1:
        break

      lpe = []
      l = 1
      while line.startswith("Layer %i PE = " % l):
        lpe.append(float(line.split()[-1]))
        line = f.readline()
        l += 1
      if not layers == l - 1:
        break

      t.append(lt)
      ke.append(lke)
      pe.append(lpe)

    line = f.readline()

  f.close()

  return (numpy.array(t, dtype = numpy.float64),
          numpy.array(ke, dtype = numpy.float64),
          numpy.array(pe, dtype = numpy.float64))
