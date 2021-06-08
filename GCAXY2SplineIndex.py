#!/usr/bin/python3
###!
# \file         GCAXY2SplineIndex.py
# \author       Bill Hill
# \date         June 2021
# \version      $Id$
# \par
# (C) University of Edinburgh, Edinburgh, UK
# (C) Heriot-Watt University, Edinburgh, UK
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be
# useful but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
# \brief
# Creates a single Woolz object in which the values are the spline
# index values corresponding to the x,y values of the coordinate value
# objects.
###

from __future__ import print_function
import sys
import re
import math as m
import json as jsn
import argparse as ap
import ctypes as c
from scipy.spatial import cKDTree
import Wlz as w

libc = c.CDLL("libc.so.6")

libc.fopen.restype = c.POINTER(w.FILE)

args = None
script_version = '0.0.1'
wlz_min_version = '1.8.0'
err_num = c.c_int(w.WLZ_ERR_NONE)

def isPython3():
  return (sys.version_info[0] == 3)
#end

def errExit(msg):
  global err_num
  e = re.sub(r'[^A-Z_0-9]', '', str(w.WlzStringFromErrorNum(err_num, None)))
  print(sys.argv[0] + ': Error - ' + msg + ' (' + e + ').', file=sys.stderr)
  sys.exit(1)
#end

def vrbMsg(msg):
  if bool(args) and args.verbose:
    print(sys.argv[0] + ': ' + msg, file=sys.stderr)
  #end
#end

def wlzVersionOK(min_version):
  wv = [int(x) for x in re.sub(r'[^0-9.]','', str(w.WlzVersion())).split('.')]
  mv = [int(x) for x in re.sub(r'[^0-9.]','', str(min_version)).split('.')]
  ok = ((wv[0] > mv[0]) or
       ((wv[0] == mv[0]) and ((wv[1] > mv[1]) or
       ((wv[1] == mv[1]) and (wv[2] >= mv[2])))))
  return(ok)
#end

def readWlzObj(fn):
  global err_num
  vrbMsg('Reading Woolz object from file ' + fn)
  obj = None
  fn = fn.encode('ascii') if (isPython3) else fn
  err_num = c.c_int(w.WLZ_ERR_FILE_OPEN)
  fp = libc.fopen(fn, 'rb')
  if(bool(fp)):
    obj = w.WlzReadObj(fp, c.byref(err_num))
    libc.fclose(fp)
  #end
  vrbMsg('obj = ' + str(obj) + ', err_num = ' +
      str(w.WlzStringFromErrorNum(err_num, None)) + ')')
  return obj
#end

def readJson(fn):
  global err_num
  vrbMsg('Reading JSON object from file ' + fn)
  obj = None
  fn = fn.encode('ascii') if (isPython3) else fn
  err_num = c.c_int(w.WLZ_ERR_FILE_OPEN)
  with open(fn, 'r') as f:
    obj = jsn.load(f)
  f.close()
  err_num = c.c_int(w.WLZ_ERR_NONE)
  return(obj)
#end

def writeWlzObj(fn, obj):
  global err_num
  vrbMsg('Writing Woolz object to file' + fn)
  fn = fn.encode('ascii') if (isPython3) else fn
  err_num = c.c_int(w.WLZ_ERR_FILE_OPEN)
  fp = libc.fopen(fn, 'wb')
  if(bool(fp)):
    err_num = w.WlzWriteObj(fp, obj)
    libc.fclose(fp)
  #end
#end

def parseArgs():
  global args
  parser = ap.ArgumentParser(
      description=
      'Creates a single Woolz object in which the values are the spline ' +
      'index values corresponding to the x,y values of the coordinate value ' +
      'objects.')
  parser.add_argument('-x', '--x-coords', type=str, required=True,
      help='Column (x) coordinate object.')
  parser.add_argument('-y', '--y-coords', type=str, required=True,
      help='Line (y) coordinate object.')
  parser.add_argument('-o', '--output', type=str, required=True,
      help='Output filename.')
  parser.add_argument('-v', '--verbose', action='store_true', default=False,
      help='Verbose output (possibly useful for debugging).')
  parser.add_argument('spline_tbl', default='',
      help='Input table of spline x,y values.')
  args = parser.parse_args()
  #end
#end

def main():
  global args
  global err_num
  gws = [None, None, None]
  iws = [None, None, None]
  obj = [None, None, None]
  # Read input x and y Woolz objects followed by the spline json object
  obj[0] = readWlzObj(args.x_coords)
  if not bool(err_num):
    obj[1] = readWlzObj(args.y_coords)
  #end
  if not bool(err_num):
    obj[2] = w.WlzCopyObject(obj[0], c.byref(err_num))
  #end
  if not bool(err_num):
    s_tbl = readJson(args.spline_tbl) 
  #end
  if not bool(err_num):
    # Build KDTree for point location queries into spline points
    kdt = cKDTree(s_tbl['points'])
    # For all coordinate pairs find nearest spline point in KDTree
    for i in range(0, len(obj)):
      gws[i] = w.WlzGreyWSpace()
      iws[i] = w.WlzIntervalWSpace()
      err_num = w.WlzInitGreyScan(obj[i], c.byref(iws[i]), c.byref(gws[i]))
      if bool(err_num):
        break
      #end
    #end
    while not bool(w.WlzNextGreyInterval(c.byref(iws[0]))):
      w.WlzNextGreyInterval(c.byref(iws[1]))
      w.WlzNextGreyInterval(c.byref(iws[2]))
      x = gws[0].u_grintptr.inp.contents.value
      y = gws[1].u_grintptr.inp.contents.value
      d, i = kdt.query([x, y], 1)
      gws[2].u_grintptr.inp.contents = c.c_int(i)
    #end
    for i in range(0, len(obj)):
      w.WlzEndGreyScan(c.byref(iws[i]), c.byref(gws[i]))
    #end
  #end
  # Write spline index object
  writeWlzObj(args.output, obj[2])
  for i in range(0, len(obj)):
    w.WlzFreeObj(obj[i])
  #end
#end

if __name__ == '__main__':
  status = 0
  # Proccess command line
  parseArgs()
  vrbMsg(sys.argv[0] + ' version: ' + script_version)
  vrbMsg('Woolz version: ' + re.sub(r'^b\'(.*)\'', r'\1', str(w.WlzVersion())))
  if(not wlzVersionOK(wlz_min_version)):
    err_num = c.c_int(w.WLZ_ERR_UNIMPLEMENTED)
    errExit('Requires PyWoolz which has been built with Woolz version >= ' +
            wlz_min_version)
  #end
  # Do the work
  main()
  sys.exit(0)
#end



