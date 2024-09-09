#!/usr/bin/python3
# -*- coding: utf-8 -*-

from __future__ import print_function

"""
GCAMakeMapTable.py

Creates a JSON encoded table of midline spline indices for closest
midline position at each vertex of the anatomy surface.
"""

__date__    = "August 2021"
__author__      = "Bill Hill"
__license__     = """
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be
  useful but WITHOUT ANY WARRANTY; without even the implied
  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public
  License along with this program; if not, write to the Free
  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
  Boston, MA  02110-1301, USA.
"""
__copyright__   = """
  (C) University of Edinburgh, Edinburgh, UK
  (C) Heriot-Watt University, Edinburgh, UK
"""

import sys
import re
import math as m
import argparse as ap
import ctypes as c
import json as jsn
from scipy.spatial import cKDTree
import Wlz as w

libc = c.CDLL('libc.so.6')
libc.fopen.restype = c.POINTER(w.FILE)

class WlzError(Exception):
  pass


line = 0
args = None
script_version = '0.0.1'
err_num = c.c_int(w.WLZ_ERR_NONE)

def isPython3():
  return (sys.version_info[0] == 3)
#end

def errExit(msg):
  e = re.sub(r'[^A-Z_0-9]', '', str(w.WlzStringFromErrorNum(err_num, None)))
  print(sys.argv[0] + ': Error - ' + msg + ' (' + e + ').', file=sys.stderr)
  sys.exit(1)
#end

def vrbMsg(msg):
  if bool(args) and args.verbose:
    print(sys.argv[0] + ': ' + msg, file=sys.stderr)
  #end
#end

def parseArgs():
  global args
  parser = ap.ArgumentParser(
      description=
      'Creates a JSON encoded table of midline spline indices for closest ' +
      'midline position at each vertex of the anatomy surface.');
  parser.add_argument('-o', '--output', type=str, default='-',
      help='Output file name.')
  parser.add_argument('-b', '--b-spline', type=str, default='-',
      help='B-Spline (JSON file format).')
  parser.add_argument('-s', '--surface', type=str, default='-',
      help='Surface (VTK ASCII file format).')
  parser.add_argument('-v', '--verbose', action='store_true', default=False,
      help='Verbose output (possibly useful for debugging).')
  parser.add_argument('nearest', type=str,
      help='Input domain object with nearest spline positions.')
  args = parser.parse_args()
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
  if bool(err_num):
    raise WlzError()
  #end
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
  return obj
#end

def getSurfVertexArray(surf_file):
  global err_num
  global line
  nv = -1
  vtx = []
  line = 0
  vrbMsg('Reading VTK file ' + surf_file + ' and parsing vertices into array')
  err_num = c.c_int(w.WLZ_ERR_FILE_OPEN)
  with open(surf_file, 'r') as f:
    for rec in f:    
      line = line + 1
      if(line == 0):
        res = re.match(r'#\s*vtk\s+datafile', rec, re.I)
        if(not bool(res)):
          err_num = c.c_int(w.WLZ_ERR_FILE_FORMAT)
          raise WlzError()
        #end
      else:
        if(nv < 0):
          res = re.match(r'points\s(\d+)\s+float', rec, re.I)
          if(bool(res)):
            nv = int(res.group(1))
          #end
        else:
          res = re.match(r'\s*(\d*\.?\d+)\s+(\d*\.?\d+)\s+(\d*\.?\d+)', 
              rec, re.I)
          if((not bool(res)) or (len(res.groups()) != 3)):
            err_num = c.c_int(w.WLZ_ERR_FILE_FORMAT)
            raise WlzError()
          #end
          vtx.append([float(res.group(1)), float(res.group(2)),
              float(res.group(3))])
          nv = nv - 1
          if(nv == 0):
            break
          #end
        #end
      #end
    #end
    err_num = c.c_int(w.WLZ_ERR_NONE)
  #end
  return vtx
#end

# a bodge because we can't use WlzMakeMain in PyWoolz
def wlzMakeMain(t, d):
  global err_num
  obj = w.WlzMakeEmpty(c.byref(err_num))
  if(bool(err_num)):
    raise WlzError()
  #end
  obj.contents.type = t
  obj.contents.domain = d
  d.core.contents.linkcount = d.core.contents.linkcount + 1
  return obj
#end

def vtxFromObjDom(obj):
  global err_num
  vtx = []
  vrbMsg('Extracting voxel coordinates from object.\n')
  pd = obj.contents.domain.p.contents
  for pl in range(pd.plane1, pd.lastpl + 1):
    idx = pl - pd.plane1
    dom2 = pd.domains[idx]
    if(bool(dom2)):
      obj2 = wlzMakeMain(w.WLZ_2D_DOMAINOBJ, dom2)
      itw = w.WlzIntervalWSpace() 
      err_num = c.c_int(
          w.WlzInitRasterScan(obj2, c.byref(itw), w.WLZ_RASTERDIR_ILIC))
      while(err_num.value == w.WLZ_ERR_NONE):
        err_num = c.c_int(w.WlzNextInterval(c.byref(itw)))
        if(not bool(err_num)):
          ln = itw.linpos
          for kl in range(itw.lftpos, itw.rgtpos + 1):
            vtx.append([kl, ln, pl])
          #end
        #end
      #end
    #end
    if(err_num.value == w.WLZ_ERR_EOO):
      err_num = c.c_int(w.WLZ_ERR_NONE)
    #end
    if(bool(err_num)):
      break
    #end
  #end
  if(bool(err_num)):
    raise WlzError()
  #end
  return vtx
#end

def getSplineIdxArray(srf_vtx, near_file, b_spl_file):
  global err_num
  spl_idx = []
  near_gvw = [None, None, None]
  vrbMsg('Creating array of closest spline vertex indices from the\n' +
         'surface vertex array, with data read from files:\n' +
         '  ' + near_file + ' and\n' +
         '  ' + b_spl_file + '\n')
  # Read nearest spline position and spline objects
  near_obj = readWlzObj(near_file)
  spl_dat = readJson(b_spl_file);
  near_cpd = c.cast(near_obj, c.POINTER(w.WlzCompoundArray))
  if((near_cpd.contents.type != w.WLZ_COMPOUND_ARR_1) or
     (near_cpd.contents.otype != w.WLZ_3D_DOMAINOBJ)):
    err_num = c.c_int(w.WLZ_ERR_OBJECT_TYPE)
    raise WlzError()
  elif(near_cpd.contents.n != 4):
    err_num = c.c_int(w.WLZ_ERR_OBJECT_DATA)
    raise WlzError()
  #end
  # Build KDTree for domain point location queries. The surface vertices
  # are often just outside of the woolz domain so use a KDTree populated
  # with domain boundary vertices to find the nearest domain vertex.
  near_gvw = [None, None, None]
  s_obj = w.WlzMakeSphereObject(c.c_int(w.WLZ_3D_DOMAINOBJ), 1.0,
      0.0, 0.0, 0.0, c.byref(err_num))
  if(bool(err_num)):
    raise WlzError()
  #end
  erd_obj = w.WlzStructErosion(near_cpd.contents.o[0], s_obj, c.byref(err_num))
  if(bool(err_num)):
    raise WlzError()
  #end
  df_obj = w.WlzDiffDomain(near_cpd.contents.o[0], erd_obj, c.byref(err_num))
  if(bool(err_num)):
    raise WlzError()
  #end
  dom_vtx = vtxFromObjDom(df_obj)
  dom_kdt = cKDTree(dom_vtx)
  # Build KDTree for point location queries into spline points
  spl_kdt = cKDTree(spl_dat['points'])
  # Create grey value workspaces for nearest spline coordinate random lookup
  for i in range(0,3):
    err_num = c.c_int(w.WLZ_ERR_NONE)
    near_gvw[i] = w.WlzGreyValueMakeWSp(near_cpd.contents.o[i],
        c.byref(err_num))
    if(bool(err_num)):
      raise WlzError()
    #end
  #end
  # For all surface vertices find nearest domain voxel position using
  # dom_kdt and then from it the index of nearest spline coordinate by
  # lookup and then spl_kdt
  for i in range(0, len(srf_vtx)):
    sv = srf_vtx[i]
    dst, idx = dom_kdt.query([sv[0], sv[1], sv[2]], k=1)
    dv = dom_vtx[idx]
    nv = [0, 0, 0]
    for j in range(0, 3):
      # Note parameter order WlzGreyValueGet(gv, pl, ln, kl)
      w.WlzGreyValueGet(near_gvw[j], dv[2], dv[1], dv[0]) 
      nv[j] = near_gvw[j].contents.gVal[0].inv
    #end
    vrbMsg(str(i) + ' ' + str(sv) + ' ' + str(dv) + ' ' + str(nv))
    dst, idx = spl_kdt.query([nv[0], nv[1], nv[2]], k=1)
    spl_idx.append(idx)
  #end
  return spl_idx
#end

def writeIdxData(out_file, spl_idx):
  global err_num
  vrbMsg('Writing B-spline index data to file ' + out_file)
  with open(out_file, 'w') as f:
    f.write('[\n')
    for i in range(0, len(spl_idx) - 1):
      v = spl_idx[i]
      f.write('  ' + str(v) + ',\n')
    #end
    v = spl_idx[-1]
    f.write('  ' + str(v) + '\n')
    f.write(']\n')
  #end
#end

if __name__ == '__main__':
  # Proccess command line
  parseArgs()
  vrbMsg(sys.argv[0] + ' version: ' + script_version)
  vrbMsg('Woolz version: ' + re.sub(r'^b\'(.*)\'', r'\1', str(w.WlzVersion())))
  try:
    surf_vtx = getSurfVertexArray(args.surface)
  except:
    errExit('Failed to read vertices from surface file ' +
         args.surface + ' (line ' + str(line) + ')')
  #end
  try:
    srf_idx = getSplineIdxArray(surf_vtx, args.nearest, args.b_spline)
  except:
    err_str = c.POINTER(c.c_char)()
    errExit('Failed to compute spline indices')
  writeIdxData(args.output, srf_idx)
#end

