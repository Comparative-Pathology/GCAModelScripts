#!/usr/bin/python3
###
# \file         GCAMapHuBMAPSampleDomains.py
# \author       Bill Hill
# \date         February 2022
# \version      $Id$
# \par
# Address:
#               Heriot-Watt University,
#               Edinburgh, Scotland, EH14 4AS, UK
# \par
# Copyright (C), [2022],
# Heriot-Watt University, Edinburgh, UK.
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
# 
# Computes the 3D spatial mapping of a HuBMAP male/female large/small Intestine
# sample to a corresponding HuBMAP v1.1 model.
###

from __future__ import print_function
import os
import re
import sys
import math as m
import json as jsn
import numpy as np
import argparse as ap
import ctypes as c
from datetime import datetime

libc = c.CDLL('libc.so.6')

haveWoolz = False
try:
  import Wlz as w
  haveWoolz = True
  libc.fopen.restype = c.POINTER(w.FILE)
except ImportError:
  pass
#end

args = None
SCRIPT_VERSION = '0.0.4'

valid_output_types = ['surface', 'domain', 'interval']

# Transform from HuBMAP registration interface coordinates to HuBMAP male v1.1
tr_rtomm_11 = np.array(
    [[0.00103615,  -3.98108e-5,   5.7803e-5,    -0.119858],
     [9.59889e-5,  -4.30641e-5,  -0.000967827,   0.0744525],
     [-2.21515e-5,  0.00100179,   1.95871e-6,   -0.0359116],
     [ 0.0,         0.0,          0.0,           1.0]])
# Transform from HuBMAP registration interface coordinates to HuBMAP model
tr_rtom_11 = tr_rtomm_11

def isPython3():
  """
  Is this python version 3?
  """
  return (sys.version_info[0] == 3)
#end

def fmtstr(s):
  """
  Simple string formating for *bold* and _underline_
  """
  c = {'b': '\033[1m', 'u': '\033[4m', 'e': '\033[0m'}
  r = re.sub(r'\*(\w*)\*', c['b'] + '\\1' + c['e'], 
             re.sub(r'_(\w*)_', c['u'] + '\\1' + c['e'], s))
  return r
#end

def degtorad(d):
  """
  Degrees to Radians
  """
  return d * m.pi / 180.
#end

def vrbMsg(msg):
  """
  Verbose messages
  Parameters
  ----------
  msg : string
        Verbose message for output
  """
  global args
  if bool(args) and bool(args.verbose):
    print(sys.argv[0] + ': ' + msg, file=sys.stderr)
  #end
#end

def errMsg(msg):
  """
  Error message then exit
  Parameters
  ----------
  msg : string
        Error message for output
  """
  global args
  if bool(args):
    print(sys.argv[0] + ': Error - ' + msg)
  #end
  sys.exit(1)
#end

def parseArgs():
  """
  Parse command line arguments
  """
  global args
  global haveWoolz
  global tr_rtom_11
  parser = ap.ArgumentParser(
      description=
      'Extracts the sample domain from a HuBMAP sample file. ' +
      'Version: ' + SCRIPT_VERSION)
  parser.add_argument('-c', '--config', type=str, default=None,
      help='Corresponding (JSON format) GCA configuration file for the ' +
           'HuBMAP model which contains the landmark positions. ' +
           'This is only required for computing midline intervals as the ' +
           'proportional distance between landmarks. ' +
           '(default: %(default)s)')
  parser.add_argument('-k', '--check-sample-file', action='store_true',
                      default=False,
      help='Check the sample file is valid before processing it')
  parser.add_argument('-p', '--path', type=str, default=None,
      metavar='PATH',
      help='Corresponding midline path file in JSON format. ' +
           'Only required for computing midline index interval. ' +
           '(default: %(default)s)')
  parser.add_argument('-o', '--output', type=str, default='-',
      help='Output filename for registered sample block. ' +
           '(default: %(default)s)')
  parser.add_argument('-r', '--transform', type=str, default=None,
      help='Woolz affine transform from the HuBMAP registration interface ' +
           'to the HuBMAP v1.1 model space ' + 
           '(default: built-in male colon transform ' +
           str(tr_rtom_11) +
           ')')
  parser.add_argument('-s', '--scale', type=float, default=1000.0,
      help='Scale factor, eg 1.0 for HuBMAP v1.1 space using metres and ' +
           '1000.0 for HuBMAP v1.1 space using millimetres. ' +
           '(default: %(default)s)')
  parser.add_argument('-T', '--transobj', action='store_true', default=False,
      help='Output a Woolz transobj if outputting a Woolz object. ' +
           '(default: %(default)s)')
  parser.add_argument('-t', '--output-type', type=str,
      default=valid_output_types[0],
      choices=valid_output_types,
      help=('Output file type for registered sample block. ' +
            'Where ' +
            fmtstr(
                '*surface* is encoded using the VTK ASCII polydata format, ' +
                '*domain* is a Woolz 3D domain object, ' +
                '*interval* may be either a simple midline index range or ' +
                'when a GCA config file is given proportional distances ' +
                'between landmarks')) + '.\n' +
            '(default: %(default)s)')
  parser.add_argument('-v', '--verbose', action='store_true', default=False,
      help='Verbose output for viewing parameter values and debugging. ' +
           '(default: %(default)s)')
  parser.add_argument('sam_file', default='',
      help='URL of input HuBMAP JSON encoded sample file (required).')
  args = parser.parse_args()
  vrbMsg('Parameter values:\n' + str(args))
  if bool(args.transform):
    if not haveWoolz:
      errMsg('Woolz is required to read affine transforms.')
    #end
    err_num = w.WlzErrorNum(w.WLZ_ERR_FILE_OPEN)
    fn = args.transform
    fn = fn.encode('ascii') if (isPython3) else fn
    fp = libc.fopen(fn, 'rb')
    if bool(fp):
      obj = w.WlzReadObj(fp, c.byref(err_num))
      libc.fclose(fp)
    #end
    if (not bool(err_num)) and (not obj.contents.type == w.WLZ_AFFINE_TRANS):
        err_num == w.WLZ_ERR_OBJECT_TYPE
    #end
    if not bool(err_num):
      mat = obj.contents.domain.t.contents.mat
      tr_rtom_11 = np.zeros((4,4))
      for i in range(0,4):
        for j in range(0,4):
          tr_rtom_11[i][j] = mat[i][j]
        #end
      #end
      vrbMsg('tr_rtom_11 = ' + str(tr_rtom_11))
    else:
      errMsg('Failed to read a Woolz affine transform from ' + args.transform)
    #end
  #end
#end

def readJson(fn):
  """
  Read JSON object from specified file
  Parameters
  ----------
  fn :  string
        Complete path of specified file
  """
  vrbMsg('Reading JSON object from file ' + fn)
  obj = None
  fn = fn.encode('ascii') if (isPython3) else fn
  with open(fn, 'r') as f:
    obj = jsn.load(f)
  f.close()
  return(obj)
#end

def checkSamFileValid(sam):
  """
  Check for valid HuBMAP sample object
  Parameters
  ----------
  sam : string
        Sample object
  """
  valid = False
  # Check contains valid anatomy
  valid_anat = set(['colon', 'large intestine',
                    'appendix', 'cecum', 'caecum',
                    'ascending colon', 'hepatic flexure',
                    'transverse colon' 'splenic flexure',
                    'descending colon', 'sigmoid colon',
                    'rectum', 'anus',
                    'small intestine', 'small bowel',
                    'duodenum', 'jejunum', 'ileum', 'terminal ileum'])
  if bool(sam):
    for i in range(0, 4):
      key = 'anatomy_' + str(i)
      if (key in sam):
        ana = sam[key]
        if(bool(valid_anat.intersection(set(ana)))):
          valid = True
          break
      #end
    #end
  #end
  # Check has a registered location
  if(valid):
    valid = ('rui_location' in sam) and ('placement' in sam['rui_location'])
  #end
  return(valid)
#end

def checkSamPlacement(loc, plc):
  """
  Check for valid HuBMAP placement in sample object
  Parameters
  ----------
  loc : object
        HuBMAP sample location object
  plc : object
        HuBMAP sample placement object
  """
  valid = ((loc['dimension_units'] == 'millimeter' or
            loc['dimension_units'] == 'millimetre') and
           (plc['rotation_units'] == 'degree') and
           (plc['translation_units'] == 'millimeter' or
            plc['translation_units'] == 'millimetre') and
           (plc['scaling_units'] == 'ratio') and
           (plc['rotation_order'] == 'XYZ'))
  return(valid)
#end

def makeBlock(block_rng):
  """
  Compute unaligned registration sample block centred on (0, 0, 0)
  Parameters
  ----------
  block_rng : array
              Minimum and maximum unaligned block coordinate values
  Returns
  -------
  array
            Vertices of the unalligned block
  """
  block = np.zeros((4,8))
  for i in range(0,8):
    block[0][i] = block_rng[0][(i & 1)]
    block[1][i] = block_rng[1][(i & 2) // 2]
    block[2][i] = block_rng[2][(i & 4) // 4]
    block[3][i] = 1.
  #end
  return block
#end

def outputSurface(block_rng, tr, block_name):
  """
  Output the aligned block as a VTK surface
  Parameters
  ----------
  block_rng :  array
               Minimum and maximum unaligned block coordinate values
  tr :         array
               Affine transform to apply to unaligned block
  block_name : string
               HuBMAP ID
  """
  global args
  # Compute unaligned registration sample block. This is centred on (0, 0, 0)
  block = makeBlock(block_rng)
  vrbMsg('block = \n' + str(block))
  # Transform block
  trblk = np.matmul(tr, block)
  vrbMsg('trblk = \n' + str(trblk))
  # Write out transformed block as a VTK surface.
  surf = ('# vtk DataFile Version 1.0\n' +
      'HuBMAPSampleCoords transformed block for ' + block_name + '\n' +
      'ASCII\n' +
      'DATASET POLYDATA\n' +
      'POINTS 8 float\n')
  for i in range(0,8):
    surf = surf + '{x:g} {y:g} {z:g}\n'.format(x=trblk[0][i],
        y = trblk[1][i], z = trblk[2][i]);
  #end
  surf = (surf + 'POLYGONS 12 48\n' +
          '3 4 0 5\n3 1 5 0\n3 7 6 4\n3 4 5 7\n3 3 2 6\n3 6 7 3\n' +
          '3 1 0 2\n3 2 3 1\n3 4 6 0\n3 0 6 2\n3 5 1 7\n3 3 7 1\n')
  f = None
  if args.output == '-':
    f = sys.stdout
  else:
    f = open(args.output, 'w')
  #end
  print(surf, file=f)
  if not (f == None or args.output == '-'):
    f.close()
  #end
#end

def outputDomain(block_rng, tr, block_name):
  """
  Output the aligned block as a Woolz domain object
  Parameters
  ----------
  block_rng :  array
               Minimum and maximum unaligned block coordinate values
  tr :         array
               Affine transform to apply to unaligned block
  block_name : string
               HuBMAP ID
  """
  global args
  global haveWoolz
  if not haveWoolz:
    errMsg('Woolz is required for domain output.')
  #end
  tfm = None
  blk_obj = None
  out_obj = None,
  rad = np.zeros(3)
  cen = np.zeros(3)
  err_num = c.c_int(w.WLZ_ERR_NONE)
  for i in range(0,3):
    rad[i] = (block_rng[i][1] - block_rng[i][0]) / 2.
    cen[i] = (block_rng[i][1] + block_rng[i][0]) / 2.
  #end
  vrbMsg('Domain centre = ' + str(cen) + ', radius = ' + str(rad))
  blk_obj = w.WlzAssignObject(
      w.WlzMakeCuboidObject(w.WLZ_3D_DOMAINOBJ,
          rad[0], rad[1], rad[2], cen[0], cen[1], cen[2],
          c.byref(err_num)), None)
  tfm = w.WlzAssignAffineTransform(
      w.WlzMakeAffineTransform(w.WLZ_TRANSFORM_3D_AFFINE, c.byref(err_num)),
      None)
  mat = tfm.contents.mat
  for i in range(0,4):
    for j in range(0,4):
      mat[i][j] = tr[i][j]
    #end
  #end
  if args.transobj:
    out_obj = w.WlzMakeEmpty(c.byref(err_num))
    out_obj.contents.type = w.WLZ_TRANS_OBJ
    out_obj.contents.domain.t = w.WlzAssignAffineTransform(tfm, None)
    out_obj.contents.values.obj = w.WlzAssignObject(blk_obj, None)
  else:
    out_obj = w.WlzAffineTransformObj(blk_obj, tfm, 
        w.WLZ_INTERPOLATION_NEAREST, c.byref(err_num))
  #end
  vrbMsg('Writing Woolz object to file ' + args.output)
  fp = None
  if args.output == '-':
    errMsg('Writting Woolz objects to stdout has not been implemented.')
  else:
    fn = args.output.encode('ascii') if (isPython3) else args.output
    fp = libc.fopen(fn, 'wb')
  #end
  if(bool(fp)):
    err_num = w.WlzWriteObj(fp, out_obj)
    if not (args.output == '-'):
      libc.fclose(fp)
    #end
  #end
  w.WlzFreeObj(blk_obj)
  w.WlzFreeAffineTransform(tfm)
  w.WlzFreeObj(out_obj)
#end

def closestPoint(points, q, i0, i1):
  """
  Finds the closest point on the midline to the given point over the
  given midline index range
  Parameters
  ----------
  points :  array
            Points on the midline
  q :       array
            Given point
  i0:       int
            Index of first point on midline at which to start search
  i1:       int
            Index of last point on midline at which to stop search
  """
  n = len(points)
  i0 = min(max(0, i0), n - 1)
  i1 = max(min(n - 1, i1), 0)
  i_min = i0
  p = points[i_min]
  vrbMsg('p = ' + str(p))
  d2_min = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2
  for i in range(i0 + 1, i1 + 1):
    p = points[i]
    d2 = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2
    if d2 < d2_min:
      d2_min = d2
      i_min = i
    #end
  #end
  return i_min
#end

def outputMidlineInterval(block_rng, tr, block_name, plc):
  """
  Output the midline index interval range
  Parameters
  ----------
  block_rng :  array
               Minimum and maximum unaligned block coordinate values
  tr :         array
               Affine transform to apply to unaligned block
  block_name : string
               HuBMAP ID
  plc :        HuBMAP sample placement
  """
  global args
  intervals = None
  positions = None
  if not bool(args.path):
    errMsg('Path file is required for computing midline path range.')
  #end
  # Read path file with points and tangents
  path = readJson(args.path)
  # Compute unregistered block
  block = makeBlock(block_rng)
  # Compute centre of unregistered sample block
  cen = np.zeros(4)
  for i in range(0,3):
    cen[i] = (block_rng[i][1] + block_rng[i][0]) / 2.
  #end
  cen[3] = 1.
  # Transform centre of block
  trcen = np.matmul(tr, cen)
  vrbMsg('trcen = \n' + str(trcen))
  # Find index of point closest to the centre of the registered sample block
  i_blk = [0, 0, 0, 0, 0, 0, 0, 0, 0]
  i_blk[8] = closestPoint(path['points'], trcen, 0, len(path['points']) - 1)
  vrbMsg('Closest midline point to centre [' +
         str(trcen[0])  + ', ' +
         str(trcen[1])  + ', ' +
         str(trcen[2]) + ']' +
         ' = points[' + str(i_blk[8]) + '] = [' +
         str(path['points'][i_blk[8]][0]) + ', ' +
         str(path['points'][i_blk[8]][1]) + ', ' +
         str(path['points'][i_blk[8]][2]) + ']')
  # Compute vertices of the registered sample block
  trblk = np.matmul(tr, block)
  vrbMsg('trblk = \n' + str(trblk))
  # Find indices on path matching registered sample block vertices
  if bool(args.config):
    vrbMsg('Computing proportional distance midline path interval.')
  else:
    vrbMsg('Computing simple midline path index interval.')
  #end
  maxrng = m.ceil(np.linalg.norm(trblk[:,0] - trblk[:,7])) * 2
  rng0 = i_blk[8] - maxrng
  rng1 = i_blk[8] + maxrng
  vrbMsg('[rng0, rng1] = [' + str(rng0) + ', ' + str(rng1) + ']')
  for i in range(0, 8):
    i_blk[i] = closestPoint(path['points'], trblk[:,i], rng0, rng1)
  #end
  vrbMsg('i_blk = ' + str(i_blk))
  interval = [min(i_blk), max(i_blk)]
  vrbMsg('interval = ' + str(interval))
  landmarks = []
  if bool(args.config):
    # Read the config file
    cfg = readJson(args.config)
    # Find which path this is
    path_id = findPathID(cfg, args.path)
    if not bool(path_id):
      errMsg('The corresponding path was not found in the given configuration '
             'file.')
    #end
    landmarks = findSortedLandmarks(cfg, path_id)
    if not bool(landmarks):
      errMsg('Unable to generate sorted list of landmarks for path ' + path_id)
    #end
    vrbMsg('landmarks = ' + str(landmarks))
    itv_lmk = [['', ''], ['', '']]
    itv_pos = [0.0, 0.0]
    tmp = indexToLandmarkPD(landmarks, interval[0])
    itv_lmk[0] = [tmp[0], tmp[1]]
    itv_pos[0] = tmp[2]
    tmp = indexToLandmarkPD(landmarks, interval[1])
    itv_lmk[1] = [tmp[0], tmp[1]]
    itv_pos[1] = tmp[2]
  #end
  f = None
  if args.output == '-':
    f = sys.stdout
  else:
    f = open(args.output, 'w')
  #end
  if bool(args.config):
    now = datetime.now()
    splc = str(plc).replace('\'', '"')
    print('{', file=f)
    print('"ident": "GCA LBPD Intervals 0.0.1",', file=f)
    print('"description": "GCA Mapped Spatial Data",', file=f)
    print('"created_by": "GCAMapHuBMAPSampleDomains.py ' +
        SCRIPT_VERSION + '",', file=f)
    print('"creation_date": "' + now.strftime('%d-%m-%Y %H:%M:%S') + '",',
        file=f)
    print('"provenance": "HuBMAP",', file=f)
    print('"external_sample_id": "' + block_name + '",', file=f)
    print('"external_placement": ' + splc + ',', file=f)
    print('"gca_sample_id": "' + block_name + '",', file=f)
    print('"gca_model_id": "' + cfg['id'] + '",', file=f)
    print('"transform": ' + transformToStr(tr) + ',', file=f)
    print('"intervals": [', file=f)
    print('  [{', file=f)
    print('    "landmarks": ["' + itv_lmk[0][0] + '",', file=f, end='')
    print(' "' + itv_lmk[0][1] + '"],', file=f)
    print('    "position":' + str(itv_pos[0]), file=f)
    print('  }, {', file=f)
    print('    "landmarks": ["' + itv_lmk[1][0] + '",', file=f, end='')
    print(' "' + itv_lmk[1][1] + '"],', file=f)
    print('    "position":' + str(itv_pos[1]), file=f)
    print('  }]', file=f)
    print('  ]', file=f)
    print('}', file=f)
  else:
    print(block_name + ' ' + str(interval[0]) + ' ' + str(interval[1]), file=f)
  #end
  if not (f == None or args.output == '-'):
    f.close()
  #end
#end

def transformToStr(t):
  """
  Creates a simple string representation of an 4x4 affine transfom matrix
  Parameters
  ----------
  t :          array
               Affine transform 
  Returns
  -------
  s :          string
               String representation of the affine transform
  """
  s = ('[[' + str(t[0][0]) + ', ' + str(t[0][1]) + ', ' +
              str(t[0][2]) + ', ' + str(t[0][3]) + '], ' +
       '[' +  str(t[1][0]) + ', ' + str(t[1][1]) + ', ' +
              str(t[1][2]) + ', ' + str(t[1][3]) + '], ' +
       '[' +  str(t[2][0]) + ', ' + str(t[2][1]) + ', ' +
              str(t[2][2]) + ', ' + str(t[2][3]) + '], ' +
       '[0.0, 0.0, 0.0, 1.0]]')
  return(s)
#end

def indexToLandmarkPD(landmarks, idx):
  """
  Given a spline index and the pre-sorted list of landmarks computes the
  landmark bounded proportional distance
  Parameters
  ----------
  landmarks :   array
                Pre-sorted list of landmarks each of which has a position and a
                GCA id. The landmarks must be pre-sorted by distance.
  idx:          The spline index.
  Returns
  -------
  lpd :         tuple
                Tuple of lower landmark, upper landmark and proportional
                distance from the lower to the upper landmark in terms of
                the spline index
  """
  li = 0
  lmk = [None, None]
  for i in range(0, len(landmarks)):
    l = landmarks[i]
    lp = l['position']
    if (lp <= idx):
      li = i
    else:
      break
    #end
  #end
  lmk[0] = landmarks[li]
  if li + 1 < len(landmarks):
    lmk[1] = landmarks[li + 1]
  else:
    lmk[1] = lmk[0]
  #end
  pd = 0.0
  if lmk[1]['position'] > lmk[0]['position']:
    pd = (idx - lmk[0]['position']) / (lmk[1]['position'] - lmk[0]['position'])
  #end
  lpd = (lmk[0]['id'], lmk[1]['id'], pd)
  return lpd
#end

def findSortedLandmarks(cfg, pth_id):
  """
  Finds and sorts (by position) the landmarks in the given GCA configuration
  with the given GCA path ID.
  Parameters
  ----------
  cfg :        GCA configuration
  pth_id :     GCA path ID
  Returns
  -------
  array of dictionary
               Sorted array of {'id', 'position'}
  """
  lmks = cfg['landmarks']
  sla = []
  for i in range(0, len(lmks)):
    lmk = lmks[i]
    path_ids = lmk['paths']
    for j in range(0, len(path_ids)):
      if pth_id == path_ids[j]:
        sla.append({'position': lmk['position'][j],
                    'id': lmk['anatomy'][j]['id']})
        break
      #end
    #end
  #end
  if len(sla) > 0:
    sla.sort(key= lambda item: item.get('position'))
  else:
    sla = None
  #end
  return sla
#end

def findPathID(cfg, pth):
  """
  Finds and returns the GCA path ID in the given GCA configuration which
  has the same filename as the given path file
  Parameters
  ----------
  cfg :        GCA configuration
  pth :        Given path file
  Returns
  -------
  string
              GCA path ID
  """
  pth_id = None
  paths = cfg['paths']
  # Strip directories from path
  p0 = os.path.basename(pth)
  for i in range(0, len(paths)):
    p = paths[i]
    if os.path.basename(p['spline_filename']) == p0:
      pth_id = p['id']
    #end 
  #end
  return pth_id
#end

def computeTransform(plc):
  """
  HuBMAP registration space is not the same as the v1.1 model space
  so use the transform tr_rtom_11. Also note that the registration uses
  mm but the v1.1 model uses metres again this is taken care of by the
  transform tr_rtom_11.
  Parameters
  ----------
  plc : object
        HuBMAP placement object
  """
  global tr_rtom_11
  # Compute translation transform tr_tln
  tr_tln = np.array(
      [[ 1,  0,  0,  plc['x_translation']],
       [ 0,  1,  0,  plc['y_translation']],
       [ 0,  0,  1,  plc['z_translation']],
       [ 0,  0,  0,  1]])
  vrbMsg('tr_tln = \n' + str(tr_tln))
  # Compute rotation transform tr_rot where
  # tr_rot = tr_rot_x * tr_rot_y * tr_rot_z for order XYZ
  r = degtorad(plc['x_rotation'])
  cr = m.cos(r)
  sr = m.sin(r)
  tr_rot_x = np.zeros((4,4))
  tr_rot_x[0][0] =  1.0
  tr_rot_x[1][1] =  cr
  tr_rot_x[1][2] = -sr
  tr_rot_x[2][1] =  sr
  tr_rot_x[2][2] =  cr
  tr_rot_x[3][3] =  1.0
  vrbMsg('tr_rot_x = \n' + str(tr_rot_x))
  tr_rot_y = np.zeros((4,4))
  r = degtorad(plc['y_rotation'])
  cr = m.cos(r)
  sr = m.sin(r)
  tr_rot_y[0][0] =  cr
  tr_rot_y[0][2] =  sr
  tr_rot_y[1][1] =  1.0
  tr_rot_y[2][0] = -sr
  tr_rot_y[2][2] =  cr
  tr_rot_y[3][3] =  1.0
  vrbMsg('tr_rot_y = \n' + str(tr_rot_y))
  tr_rot_z = np.zeros((4,4))
  r = degtorad(plc['z_rotation'])
  cr = m.cos(r)
  sr = m.sin(r)
  tr_rot_z[0][0] =  cr
  tr_rot_z[0][1] = -sr
  tr_rot_z[1][0] =  sr
  tr_rot_z[1][1] =  cr
  tr_rot_z[2][2] =  1.0
  tr_rot_z[3][3] =  1.0
  vrbMsg('tr_rot_z = \n' + str(tr_rot_z))
  tr_rot = np.matmul(tr_rot_x, np.matmul(tr_rot_y, tr_rot_z))
  vrbMsg('tr_rot = \n' + str(tr_rot))
  # Compute scaling transform tr_scl. This assumes scale is ratio
  tr_scl = np.zeros((4,4))
  tr_scl[0][0] = plc['x_scaling']
  tr_scl[1][1] = plc['y_scaling']
  tr_scl[2][2] = plc['z_scaling']
  tr_scl[3][3] = 1.0
  vrbMsg('tr_scl = \n' + str(tr_scl))
  # Commad line scale factor
  tr_csf = np.zeros((4,4))
  tr_csf[0][0] = args.scale
  tr_csf[1][1] = args.scale
  tr_csf[2][2] = args.scale
  tr_csf[3][3] = 1.0
  vrbMsg('tr_csf = \n' + str(tr_csf))
  # Compute overall transform from sample to the model
  tr = np.matmul(tr_csf,
       np.matmul(tr_rtom_11,
       np.matmul(tr_tln,
       np.matmul(tr_rot, tr_scl))))
  vrbMsg('tr = \n' + str(tr))
  return tr
#end

def main():
  global args
  sam = readJson(args.sam_file)
  if args.check_sample_file and (not checkSamFileValid(sam)):
    errMsg('Sample file does not have a relevant registered location.')
  #end
  loc = jsn.loads(sam['rui_location'])
  plc = loc['placement']
  if not checkSamPlacement(loc, plc):
    errMsg('Sample file does not include valid registration transform.')
  #end
  tr = computeTransform(plc)

  
  # Compute sample block range centred on (0, 0, 0)
  block_rng = np.zeros((4,2))
  for i in range(0,2):
    block_rng[0][i] = ((i & 1) - 0.5) * loc['x_dimension']
    block_rng[1][i] = ((i & 1) - 0.5) * loc['y_dimension']
    block_rng[2][i] = ((i & 1) - 0.5) * loc['z_dimension']
    block_rng[3][i] = 1.
  #end
  vrbMsg('block_rng = \n' + str(block_rng))
  # Compute and output required object
  if args.output_type == 'surface':
    outputSurface(block_rng, tr, sam['hubmap_id'])
  elif args.output_type == 'domain':
    outputDomain(block_rng, tr, sam['hubmap_id'])
  elif args.output_type == 'interval':
    outputMidlineInterval(block_rng, tr, sam['hubmap_id'], plc)
  #end


#end

if __name__ == '__main__':
  # Proccess command line
  parseArgs()
  vrbMsg(sys.argv[0] + ' version: ' + SCRIPT_VERSION)
  main()
  sys.exit(0)
#end





