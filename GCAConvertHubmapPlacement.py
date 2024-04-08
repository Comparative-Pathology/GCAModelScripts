#!/usr/bin/python3
###
# \file         GCAConvertHubmapPlacement.py
# \author       Bill Hill
# \date         March 2024
# \version      $Id$
# \par
# Address:
#               Heriot-Watt University,
#               Edinburgh, Scotland, EH14 4AS, UK
# \par
# Copyright (C), [2024],
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
# Converts sample placements downloaded from HuBMAP to Edinburgh Gut Cell Atlas
# landmark bounded proportional distance intervals.
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

args = None

VERBOSE_MAX = 3
SCRIPT_VERSION = '0.0.1'

file_req_text= '''
files:

Input sample files should be JSON format with HuBMAP samples of the form:
{
  {
    "samples": [
      {
        "sample_id": "STR",
        ...,
        "rui_location": {
          "rui_location": {
            placement": {
              "@type": "STR",
              "rotation_order": "XYZ",
              "rotation_units": "degree",
              "scaling_units": "ratio",
              "translation_units": "STR",
              "x_rotation": FPN,
              "x_scaling": FPN,
              "x_translation": FPN,
              "y_rotation": FPN,
              "y_scaling": FPN,
              "y_translation": FPN,
              "z_rotation": FPN,
              "z_scaling": FPN,
              "z_translation": FPN,
              ...
            },
            "dimension_units": "millimeter",
            "x_dimension": FPN,
            "y_dimension": FPN,
            "z_dimension": FPN,
            ...
          },
          ...
        }
      }
    ],
    ...
    [
      ...
    ]
  }
}

Paths and configuration files must be as used for the Edinburgh Gut Cell
Atlas models.

Output files will have the form:
{
  "ident": "GCA LBPD Intervals 0.0.1",
  "description": "GCA Mapped Spatial Data. Intervals are represented as
      ordered start then end positions, with each position being represented
      as a proportional distance from the first to the second landmark.
      All landmarks and positions directed from the anus towards the
      gastro-duodenum junction. A transform is from the unaligned source to
      the gca_model. ",
  "created_by": "GCAMapHuBMAPSampleDomains.py VSN",
  "creation_date": "STR
  "source": "HuBMAP",
  "source_sample_id": "STR
  "source_placement": "STR
  "gca_sample_id": "STR",
  "gca_model": "STR VSN",
  "gca_path": "STR VSN",
  "transform": [[FPN, FPN, FPN, FPN], [FPN, FPN, FPN, FPN],
      [FPN, FPN, FPN, FPN], [0.0, 0.0, 0.0, 1.0]],
  "intervals": [
    [{
      "landmarks": ["LMK", "LMK"],
      "position":FPN
    }, {
      "landmarks": ["LMK", "LMK"],
      "position":FPN
    }]
  ]
}

In all cases STR is a string, VSN a version number string and LMK Gut Cell
Atlas landmark string. Intervals are defined between the two landmark bounded
proportional positions, with the position the fractional distance along the
models midline from the first landmark to the second.
'''

rui_type_filters = {'e': 'SpatialEntity', 'p': 'SpatialPlacement'}

# Transform from HuBMAP registration interface coordinates to HuBMAP female
# large intestine v1.1
tr_htom11_fl = np.array(
    [[ 1.02093,      0.00443812,   0.0137555,    -129.303],
     [-0.078389,     0.0551648,   -1.0856,        154.624],
     [ 0.0138375,    0.995437,     0.0222402,    -16.9894],
     [ 0.0,          0.0,          0.0,           1.0]])

# Transform from HuBMAP registration interface coordinates to HuBMAP male 
# large intestine v1.1
tr_htom11_ml = np.array(
    [[ 0.00103615,  -3.98108e-5,   5.7803e-5,    -0.119858],
     [ 9.59889e-5,  -4.30641e-5,  -0.000967827,   0.0744525],
     [-2.21515e-5,   0.00100179,   1.95871e-6,   -0.0359116],
     [ 0.0,          0.0,          0.0,           1.0]])

# Transform from HuBMAP registration interface coordinates to HuBMAP female
# small intestine v1.1
tr_htom11_fs = np.array(
    [[ 0.994691,    -0.0276993,    0.0809571,    -88.7493],
     [-0.0155947,    0.0267578,   -1.12891,       83.952],
     [-0.0223554,    0.984934,     0.0377081,     79.7066],
     [ 0.0,          0.0,          0.0,           1.0]])

# Transform from HuBMAP registration interface coordinates to HuBMAP male
# small intestine v1.1
tr_htom11_ms = np.array(
    [[ 0.808391,     0.0486574,   -0.196351,    -43.2267],
     [ 0.0706348,   -0.0265005,   -0.790626,    -14.2504],
     [-0.264289,     1.06319,     -0.285726,   126.522],
     [ 0.0,          0.0,          0.0,          1.0]])


# Target model data
target_models = {'VHFColon':
                 {'sex': 'F', 'id': 'large intestine', 'version': '1.0',
                  'path file': None, 'path': None, 'config': None,
                  'transform': tr_htom11_fl},
                 'VHMColon':
                 {'sex': 'M', 'id': 'large intestine', 'version': '1.1',
                  'path file': None, 'path': None, 'config': None,
                  'transform': tr_htom11_ml},
                 'VHFLargeIntestine':
                 {'sex': 'F', 'id': 'large intestine', 'version': '1.1',
                  'path file': None, 'path': None, 'config': None,
                  'transform': tr_htom11_fl},
                 'VHMLargeIntestine':
                 {'sex': 'M', 'id': 'large intestine', 'version': '1.1',
                  'path file': None, 'path': None, 'config': None,
                  'transform': tr_htom11_ml},
                 'VHFSmallIntestine':
                 {'sex': 'F', 'id': 'small intestine', 'version': '1.1',
                  'path file': None, 'path': None, 'config': None,
                  'transform': tr_htom11_fs},
                 'VHMSmallIntestine':
                 {'sex': 'F', 'id': 'small intestine', 'version': '1.1',
                  'path file': None, 'path': None, 'config': None,
                  'transform': tr_htom11_ms}}



# Regular expressions
tgtmod_rx = re.compile('^ *http://purl.org/ccf/latest/ccf.owl#([A-Z][A-Za-z]+)')
mmetre_rx = re.compile('[Mm]illimet(?:re|er)')

sample_id_list = []

def isPython3():
  """
  Is this python version 3?
  """
  return (sys.version_info[0] == 3)
#end

def degToRad(d):
  """
  Degrees to Radians
  """
  return d * m.pi / 180.
#end

def vrbMsg(lvl, msg):
  """
  Verbose messages
  Parameters
  ----------
  lvl : int
        Verbosity level of the message
  msg : string
        Verbose message for output
  """
  global args
  if bool(args) and (lvl <= args.verbose):
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
  Parse commandline arguments
  """
  global args
  parser = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,
    description =
      '''
      Parses and converts HuBMAP RUI locations to Edinburgh GCA landmark
      bounded proportional distance intervals 
      ''' +
      'Version: ' + SCRIPT_VERSION,
    epilog=file_req_text)

  parser._action_groups.pop()
  required = parser.add_argument_group('required arguments')
  optional = parser.add_argument_group('optional arguments')
  required.add_argument('-c', '--conf-m', type=str, required=True,
      help='''
      Edinburgh GCA config file for the HuBMAP visible human male model.
      ''')
  required.add_argument('-d', '--conf-f', type=str, required=True,
      help='''
      Edinburgh GCA config file for the HuBMAP visible human female model.
      ''')
  optional.add_argument('-f', '--filter-types', type=str, default='e,p',
      help='''
      Filter input to only accept given rui_location types. Comma separated
      types are SpatialEntity (e) and SpatialPlacement (p)
      (default: %(default)s).
      ''')
  optional.add_argument('-i', '--id', type=str, default='',
      help='''
      Sample ID string. If given a single sample is processed else all
      samples in the file are processed.
      ''')
  required.add_argument('-P', '--path-ml', type=str, required=True,
      help='''
      Midline path file for the HuBMAP visible human male large intestine
      model.
      ''')
  required.add_argument('-Q', '--path-fl', type=str, required=True,
      help='''
      Midline path file for the HuBMAP visible human female large intestine
      model.
      ''')
  required.add_argument('-p', '--path-ms', type=str, required=True,
      help='''
      Midline path file for the HuBMAP visible human male small intestine
      model.
      ''')
  required.add_argument('-q', '--path-fs', type=str, required=True,
      help='''
      Midline path file for the HuBMAP visible human female small intestine
      model.
      ''')
  optional.add_argument('-o', '--output', type=str, default='-',
      help='Output filename for registered sample block. ' +
           '(default: %(default)s)')
  optional.add_argument('-v', '--verbose', type=int, default=0,
      help='Verbose output for viewing parameter values and debugging. ' +
           'Range [-1 - ' + str(VERBOSE_MAX) + ']; ' +
           '0 - just warnings, ' +
           str(VERBOSE_MAX) + ' - lots of verbose output. ' +
           'Use -1 for no warnings.'
           '(default: %(default)s)')
  parser.add_argument('sam_file', 
      help='URL of input HuBMAP JSON encoded sample file (required).')

  args = parser.parse_args()
  if args.verbose < -1:
    args.verbose = -1
  elif args.verbose > VERBOSE_MAX:
    args.verbose = VERBOSE_MAX
  #end
  vrbMsg(1, 'Parameter values:\n' + str(args))
#end

def readJson(fn):
  """
  Read JSON object from specified file
  Parameters
  ----------
  fn :  string
        Complete path of specified file
  """
  vrbMsg(1, 'Reading JSON object from file ' + fn)
  obj = None
  fn = fn.encode('ascii') if (isPython3) else fn
  with open(fn, 'r') as f:
    obj = jsn.load(f)
  f.close()
  return(obj)
#end

def findPathID(cfg, pth_file):
  """
  Finds and returns the GCA path ID in the given GCA configuration which
  has the same filename as the given path file
  Parameters
  ----------
  cfg :        GCA configuration
  pth_file :   Given path file
  Returns
  -------
  string
              GCA path ID
  """
  vrbMsg(3, 'Finding path id for path file ' + pth_file)
  pth_id = None
  paths = cfg['paths']
  vrbMsg(3, 'Config paths ' + str(paths))
  # Strip directories from path
  p0 = os.path.basename(pth_file)
  for i in range(0, len(paths)):
    p = paths[i]
    if os.path.basename(p['spline_filename']) == p0:
      pth_id = p['id']
    #end
  #end
  return pth_id
#end

def computeTransform(plm, mod):
  """
  HuBMAP registration space is not the same as the v1.1 model space
  so use the transform tr_htom11_ml or tr_htom11_fl. Also note that the
  registration uses mm but the v1.1 model uses metres again this is
  taken care of by the transform tr_htom11_ml.
  Parameters
  ----------
  plm : object
        HuBMAP rui placement object
  mod : object
        Appropriate model
  """
  # Compute translation transform tr_tln
  tr_tln = np.array(
      [[ 1,  0,  0,  plm['x_translation']],
       [ 0,  1,  0,  plm['y_translation']],
       [ 0,  0,  1,  plm['z_translation']],
       [ 0,  0,  0,  1]])
  vrbMsg(2, 'tr_tln = \n' + str(tr_tln))
  # Compute rotation transform tr_rot where
  # tr_rot = tr_rot_x * tr_rot_y * tr_rot_z for order XYZ
  r = degToRad(plm['x_rotation'])
  cr = m.cos(r)
  sr = m.sin(r)
  tr_rot_x = np.zeros((4,4))
  tr_rot_x[0][0] =  1.0
  tr_rot_x[1][1] =  cr
  tr_rot_x[1][2] = -sr
  tr_rot_x[2][1] =  sr
  tr_rot_x[2][2] =  cr
  tr_rot_x[3][3] =  1.0
  vrbMsg(2, 'tr_rot_x = \n' + str(tr_rot_x))
  tr_rot_y = np.zeros((4,4))
  r = degToRad(plm['y_rotation'])
  cr = m.cos(r)
  sr = m.sin(r)
  tr_rot_y[0][0] =  cr
  tr_rot_y[0][2] =  sr
  tr_rot_y[1][1] =  1.0
  tr_rot_y[2][0] = -sr
  tr_rot_y[2][2] =  cr
  tr_rot_y[3][3] =  1.0
  vrbMsg(2, 'tr_rot_y = \n' + str(tr_rot_y))
  tr_rot_z = np.zeros((4,4))
  r = degToRad(plm['z_rotation'])
  cr = m.cos(r)
  sr = m.sin(r)
  tr_rot_z[0][0] =  cr
  tr_rot_z[0][1] = -sr
  tr_rot_z[1][0] =  sr
  tr_rot_z[1][1] =  cr
  tr_rot_z[2][2] =  1.0
  tr_rot_z[3][3] =  1.0
  vrbMsg(2, 'tr_rot_z = \n' + str(tr_rot_z))
  tr_rot = np.matmul(tr_rot_x, np.matmul(tr_rot_y, tr_rot_z))
  vrbMsg(2, 'tr_rot = \n' + str(tr_rot))
  # Compute scaling transform tr_scl. This assumes scale is ratio
  tr_scl = np.zeros((4,4))
  tr_scl[0][0] = plm['x_scaling']
  tr_scl[1][1] = plm['y_scaling']
  tr_scl[2][2] = plm['z_scaling']
  tr_scl[3][3] = 1.0
  vrbMsg(2, 'tr_scl = \n' + str(tr_scl))
  # Scale factor (to metres)
  tr_csf = np.zeros((4,4))
  sc = 1000.0
  if ((mod['sex'] == 'F') or
      (('translation_units' in plm) and
       not bool(re.match(mmetre_rx, plm['translation_units'])))):
    sc = 1.0
  #end
  tr_csf[0][0] = sc
  tr_csf[1][1] = sc
  tr_csf[2][2] = sc
  tr_csf[3][3] = 1.0
  vrbMsg(2, 'tr_csf = \n' + str(tr_csf))
  # Compute overall transform from sample to the model
  tr_rto = mod['transform']
  tr = np.matmul(tr_csf,
       np.matmul(tr_rto,
       np.matmul(tr_tln,
       np.matmul(tr_rot, tr_scl))))
  vrbMsg(1, 'tr = \n' + str(tr))
  return tr
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
  vrbMsg(2, 'p = ' + str(p))
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

def outputMidlineInterval(block_rng, tr, block_name, plc, mod):
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
  plc :        HuBMAP rui sample placement
  mod          The appropriate HuBMAP model
  """
  global args
  intervals = None
  positions = None
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
  vrbMsg(2, 'trcen = \n' + str(trcen))
  # Find index of point closest to the centre of the registered sample block
  i_blk = [0, 0, 0, 0, 0, 0, 0, 0, 0]
  i_blk[8] = closestPoint(mod['path']['points'], trcen, 0,
      len(mod['path']['points']) - 1)
  vrbMsg(2, 'Closest midline point to centre [' +
         str(trcen[0])  + ', ' +
         str(trcen[1])  + ', ' +
         str(trcen[2]) + ']' +
         ' = points[' + str(i_blk[8]) + '] = [' +
         str(mod['path']['points'][i_blk[8]][0]) + ', ' +
         str(mod['path']['points'][i_blk[8]][1]) + ', ' +
         str(mod['path']['points'][i_blk[8]][2]) + ']')
  # Compute vertices of the registered sample block
  trblk = np.matmul(tr, block)
  vrbMsg(2, 'trblk = \n' + str(trblk))
  # Find indices on path matching registered sample block vertices
  vrbMsg(2, 'Computing proportional distance midline path interval.')
  maxrng = m.ceil(np.linalg.norm(trblk[:,0] - trblk[:,7])) * 2
  rng0 = i_blk[8] - maxrng
  rng1 = i_blk[8] + maxrng
  vrbMsg(2, '[rng0, rng1] = [' + str(rng0) + ', ' + str(rng1) + ']')
  for i in range(0, 8):
    i_blk[i] = closestPoint(mod['path']['points'], trblk[:,i], rng0, rng1)
  #end
  vrbMsg(2, 'i_blk = ' + str(i_blk))
  interval = [min(i_blk), max(i_blk)]
  vrbMsg(2, 'interval = ' + str(interval))
  landmarks = []
  # Find which path this is
  path_id = findPathID(mod['config'], mod['path file'])
  if not bool(path_id):
    errMsg('The corresponding path was not found in the given configuration '
           'file.')
  #end
  landmarks = findSortedLandmarks(mod['config'], path_id)
  if not bool(landmarks):
    errMsg('Unable to generate sorted list of landmarks for path ' + path_id)
  #end
  vrbMsg(2, 'landmarks = ' + str(landmarks))
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
    f = open(args.output + block_name + '.jsn', 'w')
  #end
  now = datetime.now()
  pv = ''
  mcv = ''
  if 'version' in mod['config']:
    mcv = str(mod['config']['version'])
  #end
  if 'sub_version' in mod['config']:
    mcv = mcv + '.' + str(mod['config']['sub_version'])
  #end
  if len(mcv) == 0:
    mcv = 'unkown'
  #end
  if 'version' in mod['path']:
    pv = mod['path']['version']
  else:
    pv = 'unknown'
  #end
  splc = str(plc).replace('\'', '"')
  print('{', file=f)
  print('"ident": "GCA LBPD Intervals 0.0.1",', file=f)
  print('"description": "GCA Mapped Spatial Data. ' +
        'Intervals are represented as ordered start then end positions, ' +
        'with each position being represented as a proportional distance ' +
        'from the first to the second landmark. All landmarks and ' +
        'positions directed from the anus towards the gastro-duodenum ' +
        'junction. ' +
        'A transform is from the unaligned source to the gca_model. ' +
        '",', file=f)
  print('"created_by": "GCAMapHuBMAPSampleDomains.py ' +
      SCRIPT_VERSION + '",', file=f)
  print('"creation_date": "' + now.strftime('%d-%m-%Y %H:%M:%S') + '",',
      file=f)
  print('"source": "HuBMAP",', file=f)
  print('"source_sample_id": "' + block_name + '",', file=f)
  print('"source_placement": ' + splc + ',', file=f)
  print('"gca_sample_id": "' + block_name + '",', file=f)
  print('"gca_model": "' + mod['config']['id'] + ' ' + mcv + '",', file=f)
  print('"gca_path": "' + path_id + ' ' + pv  + '",', file=f)
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

def main():
  global args
  global target_models
  global rui_type_filters
  global tgtmod_rx
  global sample_id_list
  sam_f = None
  finished = False
  rui_types = []
  for f in args.filter_types:
    if f in rui_type_filters:
      rui_types.append(rui_type_filters[f])
    #end
  #end
  vrbMsg(1, 'RUI types selected are ' + str(rui_types))
  # Read the config files for the Edinburgh GCA HuBMAP VHM and VHF models.
  cfg = readJson(args.conf_f)
  target_models['VHFColon']['config'] = cfg
  target_models['VHFLargeIntestine']['config'] = cfg
  target_models['VHFSmallIntestine']['config'] = cfg
  cfg = readJson(args.conf_m)
  target_models['VHMColon']['config'] = cfg
  target_models['VHMLargeIntestine']['config'] = cfg
  target_models['VHMSmallIntestine']['config'] = cfg
  # Read the path files for the Edinburgh GCA HuBMAP VHM and VHF large and
  # small intestine paths.
  target_models['VHFColon']['path file'] = args.path_fl
  target_models['VHFLargeIntestine']['path file'] = args.path_fl
  target_models['VHMColon']['path file'] = args.path_ml
  target_models['VHMLargeIntestine']['path file'] = args.path_ml
  target_models['VHFSmallIntestine']['path file'] = args.path_fs
  target_models['VHMSmallIntestine']['path file'] = args.path_ms
  pth = readJson(args.path_fl)
  target_models['VHFColon']['path'] = pth
  target_models['VHFLargeIntestine']['path'] = pth
  pth = readJson(args.path_ml)
  target_models['VHMColon']['path'] = pth
  target_models['VHMLargeIntestine']['path'] = pth
  pth = readJson(args.path_fs)
  target_models['VHFSmallIntestine']['path'] = pth
  pth = readJson(args.path_ms)
  target_models['VHMSmallIntestine']['path'] = pth
  # Read the sample file (JSON format)
  samples = None
  sam_o = readJson(args.sam_file)
  if 'samples' in sam_o:
    samples = sam_o['samples']
    if (not isinstance(samples, list)) or (len(samples) < 1):
      samples = None
    #end
  #end
  if not bool(samples):
    errMsg('Sample file ' + args.sam_file + ' has no samples.')
  else:
    for sam in samples:
      sid = None
      rui = None
      plm = None
      tgt = None
      mod = None
      if not 'sample_id' in sam:
        vrbMsg(0, 'Warning no sample id for sample.')
      else:
        sid = sam['sample_id']
        if (bool(sid) and ((len(args.id) == 0) or (args.id == sid))):
          vrbMsg(1, 'Sample ' + sid + ' matched.')
          if sid in sample_id_list:
            vrbMsg(0, 'Warning repeat of sample ' + sid)
          else:
            sample_id_list.append(sid)
          #end
        else:
          sid = None
        #end
      #end
      if bool(sid):
        vrbMsg(2, 'sid = ' + str(sid))
        if not 'rui_location' in sam:
          vrbMsg(0, 'Warning no rui location for sample ' + sid)
        else:
          rui = sam['rui_location']
        #end
      #end
      if bool(rui):
        vrbMsg(2, 'rui = ' + str(rui))
        if not 'placement' in rui:
          vrbMsg(0, 'Warning no placement for rui location in sample ' +
            sid)
        else:
          plm = rui['placement']
        #end
      #end
      if bool(plm):
        if not 'target' in plm:
          vrbMsg(0, 'Warning no target for placement in sample ' + sid)
        else:
          if '@type' in rui:
            if rui['@type'] in rui_types:
              t0 = re.findall(tgtmod_rx, plm['target'])
              tgt = t0[0]
              vrbMsg(2, 'rui type  ' + rui['@type'])
            else:
              vrbMsg(0, 'rui type  ' + rui['@type'] +
                  ' not in filter types for sample ' + sid)
            #end
          else:
            vrbMsg(0, 'Warning no type for sample ' + sid)
          #end
        #end
      #end
      if bool(tgt):
        vrbMsg(2, 'tgt = ' + str(tgt))
        if tgt in target_models:
          mod = target_models[tgt]
          vrbMsg(2, 'mod = ' + str(mod))
          tr = computeTransform(plm, mod)
          # Compute sample block range centred on (0, 0, 0)
          block_rng = np.zeros((4,2))
          for i in range(0,2):
            block_rng[0][i] = ((i & 1) - 0.5) * rui['x_dimension']
            block_rng[1][i] = ((i & 1) - 0.5) * rui['y_dimension']
            block_rng[2][i] = ((i & 1) - 0.5) * rui['z_dimension']
            block_rng[3][i] = 1.
          #end
          vrbMsg(1, 'block_rng = \n' + str(block_rng))
          outputMidlineInterval(block_rng, tr, sid, rui, mod)
        else:
          vrbMsg(0, 'Warning unknown rui target (' + tgt + ') for sample ' +
                 sid)
        #end
      #end
    #end
  #end
#end


if __name__ == '__main__':
  # Proccess command line
  parseArgs()
  vrbMsg(1, sys.argv[0] + ' version: ' + SCRIPT_VERSION)
  main()
  sys.exit(0)
#end

