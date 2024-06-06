#!/usr/bin/python3
###
# \file         mapSurfItvToHuBMAP.py
# \author       Bill Hill
# \date         May 2024
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
# Maps Edinburgh GCA landmark bounded proportional distance intervals
# to HuBMAP RUI placements
###

import sys
import re
import math as m
import json as jsn
import numpy as np
import trimesh as tm
import argparse as ap
import datetime as dt
from shapely import geometry as gm

args = None

VERBOSE_MAX = 3
SCRIPT_VERSION = '0.0.1'
DEFAULT_RADIUS = 3
CCF_OWL = 'http://purl.org/ccf/latest/ccf.owl'
CCF_ONTOLOGY='https://hubmapconsortium.github.io/ccf-ontology/ccf-context.json'

EPSILON = np.finfo(float).eps

HELP_TESTING = '''

This software uses Trimesh (https://trimesh.org/) to create and process
surfaces.

Set up environment, eg using micromamba

micromamba create -n trimesh -c conda-forge
micromamba activate trimesh
micromamba install -c conda-forge python=3.8 \\
    trimesh shapely meshpy py-triangle pyglet

Test using something like

GM="$GCA_MODEL_HOME"/models/HuBMAPVHM_3D_00070
GF="$GCA_MODEL_HOME"/models/HuBMAPVHF_3D_00080/

micromamba run -n trimesh python \\
./mapSurfItvToHuBMAP.py \\
    -v 1 -t m \\
    -c $GCA_MODEL_HOME/HuBMAPVHM_3D_00070_1_5.json \\
    -p $GM/paths/VH_M_Small_Intestine_v11_path-bs.jsn \\
    -P  $GM/paths/VHM_Colon_Low_path_11-bs.jsn \\
    00bf92f03e3eadf91552850d32eff9e4.jsn

micromamba run -n trimesh python \\
./mapSurfItvToHuBMAP.py \\
    -v 1 -t f \\
    -c $GCA_MODEL_HOME/HuBMAPVHF_3D_00080_1_6.json \\
    -p $GF/paths/VH_F_Small_Intestine_v11_path-bs.jsn \\
    -P $GF/paths/VH_F_colon_v11_path-bs.jsn \\
    00bf92f03e3eadf91552850d32eff9e4.jsn

'''


# Regular expressions
ws_rx = re.compile('[\\s\\n]+')
magic_rx = re.compile('#\\s+vtk\\s+datafile\\s+version\\s+1\\.0', re.IGNORECASE)
ascii_rx = re.compile('ascii', re.IGNORECASE)
dataset_rx = re.compile('dataset\\s+polydata', re.IGNORECASE)
points_rx = re.compile('points\\s+[1-9][0-9]*\\s+float', re.IGNORECASE)
polygons_rx = re.compile('POLYGONS\s+[1-9][0-9]*\\s+[1-9][0-9]*')

# Transform from HuBMAP registration interface coordinates to HuBMAP female
# large intestine v1.1
tr_htoe_fl = np.array(
    [[ 1.02093,      0.00443812,   0.0137555,    -129.303],
     [-0.078389,     0.0551648,   -1.0856,        154.624],
     [ 0.0138375,    0.995437,     0.0222402,    -16.9894],
     [ 0.0,          0.0,          0.0,           1.0]])

# Transform from HuBMAP registration interface coordinates to HuBMAP male
# large intestine v1.1
tr_htoe_ml = np.array(
    [[ 0.00103615,  -3.98108e-5,   5.7803e-5,    -0.119858],
     [ 9.59889e-5,  -4.30641e-5,  -0.000967827,   0.0744525],
     [-2.21515e-5,   0.00100179,   1.95871e-6,   -0.0359116],
     [ 0.0,          0.0,          0.0,           1.0]])

# Transform from HuBMAP registration interface coordinates to HuBMAP female
# small intestine v1.1
tr_htoe_fs = np.array(
    [[ 0.994691,    -0.0276993,    0.0809571,    -88.7493],
     [-0.0155947,    0.0267578,   -1.12891,       83.952],
     [-0.0223554,    0.984934,     0.0377081,     79.7066],
     [ 0.0,          0.0,          0.0,           1.0]])

# Transform from HuBMAP registration interface coordinates to HuBMAP male
# small intestine v1.1
tr_htoe_ms = np.array(
    [[ 0.808391,     0.0486574,   -0.196351,    -43.2267],
     [ 0.0706348,   -0.0265005,   -0.790626,    -14.2504],
     [-0.264289,     1.06319,     -0.285726,   126.522],
     [ 0.0,          0.0,          0.0,          1.0]])


model = {
    'target': None,    # Possible valid targets are (m|f)(l|s), eg 'fl', ...
    'config': None,    # Edinburgh configuration file
    'path_l': None,    # Edinburgh large intestine path file
    'path_s': None,    # Edinburgh samll intestine path file
    'tr_htoe_l': None, # Transform from HuBMAP to Edinburgh large int model
    'tr_htoe_s': None  # Transform from HuBMAP to Edinburgh small int model
}

def isPython3():
  """ 
  Is this python version 3?
  """
  return (sys.version_info[0] == 3)
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

def radToDeg(x):
  """
  Radians to degrees
  Parameters:
  ----------
    x : angle in radians
  Returns
  -------
    float : angle in degrees
  """
  return x * 180. / m.pi

def fmtstr(s):
  """
  Simple string formating for *bold* and _underline_
  Parameters:
  ----------
    s : string for formatting
  Returns
  ------
    Formatted string
  """
  c = {'b': '\033[1m', 'u': '\033[4m', 'e': '\033[0m'}
  r = re.sub(r'\*(\w*)\*', c['b'] + '\\1' + c['e'],
             re.sub(r'_(\w*)_', c['u'] + '\\1' + c['e'], s))
  return r
#end

def readJson(fn):
  """
  Read JSON object from specified file
  Parameters:
  ----------
  fn :  string
        Complete path of specified file
  Returns
  -------
  JSON object read from file
  
  """
  vrbMsg(1, 'Reading JSON object from file ' + fn)
  obj = None
  fn = fn.encode('ascii') if (isPython3) else fn
  with open(fn, 'r') as f:
    obj = jsn.load(f)
  f.close()
  return(obj)
#end


def pathIdxFromLBPD(mod, lbpd):
  """
  Computes the path index from a landmark bounded proportional distance
  Parameters
  ----------
  mod :  model dictionary
  lbpd: landmark bounded proportional distance
  Returns
  -------
  pnm, idx
  pnm: model path name ('path_l' or 'path_s') or None on error
  idx:   path index or -1 on error)
  """
  idx = -1
  pnm = None
  cfg = mod['config']
  cfg_paths = cfg['paths']
  cfg_lmks = cfg['landmarks']
  lbpd_lmks = lbpd['landmarks']
  pos = [-1, -1]
  pth = [None, None]
  # Find the landmarks of the landmark bounded proportional distance in the
  # config file
  for i in range(0, 2):
    for lmk in cfg_lmks:
      if lbpd_lmks[i] == lmk['anatomy'][0]['id']:
        pos[i] = lmk['position'][0]
        pth[i] = lmk['paths'][0]
      #end
    #end
  #end
  # For now, insist that an interval is on a single path (pth[0] == pth[1])
  if (pos[0] > -EPSILON) and (pos[1] > -EPSILON) and (pth[0] == pth[1]):
    idx = int(m.floor((pos[1] - pos[0]) * lbpd['position'] + pos[0]))
    if(idx < 0):
      idx = 0
    #end
  #end
  ubn_pth = None
  if   (pth[0] == cfg_paths[0]['id']):
    ubn_pth = cfg_paths[0]['anatomy']['external_id']
  elif (pth[0] == cfg_paths[1]['id']):
    ubn_pth = cfg_paths[1]['anatomy']['external_id']
  #end
  if ubn_pth == 'UBERON:0000059':
    pnm = 'path_l'
  elif ubn_pth == 'UBERON:0002108':
    pnm = 'path_s'
  #end
  if not bool(pnm):
    idx = -1
  #end
  return pnm, idx
#end

def surfFromItvls(model, path_n, itv):
  """
  Creates a spline extrusion following the midline between the
  interval ends.
  Parameters
  ----------
  model:  the model
  path_n: path name, either 'path_l' (large intestine) or
          'path_s' (small intestine)
  itv:    interval indices along the path
  Returns
  -------
  surf: an extruded cylinder along the midline path between covering the
        interval
  """

  lors = 'l' if path_n[0] == 'path_l' else 's'
  model['target'] = args.target + lors
  pth = model['path_' + lors]
  vrbMsg(2, 'lors = ' + lors)
  vrbMsg(2, 'itv = ' + str(itv))
  vrbMsg(2, 'points at ends if interval = ' +
      str(pth['points'][itv[0]]) + ', ' +
      str(pth['points'][itv[1]]))
  # Create a polygon to be extruded along the midline
  vec = np.array([0, 1]) * args.radius
  nseg = 64
  rotang = 2.0 * np.pi / nseg
  rotmat = np.array([
      [np.cos(rotang), -np.sin(rotang)],
      [np.sin(rotang), np.cos(rotang)]])
  perim = []
  for i in range(nseg):
    perim.append(vec)
    vec = np.dot(rotmat, vec)
  poly = gm.Polygon(perim)
  # Create midline interval array
  pnts = np.array(pth['points'][itv[0]:itv[1]])
  vrbMsg(3, 'pnts = ' + str(pnts))
  # Extrude polygon along path between intrvals to create surface 
  surf = tm.creation.sweep_polygon(poly, pnts, engine='triangle')
  #end
  # surf.show()
  surf.export('surf.stl')
  return surf
#end

def computeOBB(points):
  """
  Compute the bounding box extents and transform to the cartesian axes
  of an array of points. The extents and transform are ordered such
  that ext[0] >= ext[1] >= ext[2].
  Parameters
  ----------
  points: array of points
  Returns
  -------
  trx,ext
  trx: 4 x 4 numpy array, transform from OBB to AABB centred at the origin
  ext: 1 x 4 numpy array, extents of the AABB
  """
  vrbMsg(3, 'points[:3] = ' + str(points[:3]))
  v = points
  vrbMsg(3, 'v = ' + str(v))
  vc = np.mean(v, axis=0)
  vrbMsg(3, 'vc = ' + str(vc))
  vp = (v - vc).T
  vrbMsg(3, 'vp = ' + str(vp))
  cov = np.cov(vp)
  vrbMsg(3, 'cov = ' + str(cov))
  # covarience matrix is symmetric and +ve semi-definite so can use eigh
  evl, evc = np.linalg.eigh(cov)
  vrbMsg(2, 'evl, evc = \n' + str(evl) + ',\n' + str(evc))
  # Compute transform to the x-y-z axes
  evct = evc.T
  vpa = np.matmul(evct, vp)
  vca = np.matmul(evct, vc.T)
  vrbMsg(2, 'vpa = ' + str(vpa))
  # Transform eigen values and find extent
  ext = np.array([
      np.max(vpa[0, :]) - np.min(vpa[0, :]),
      np.max(vpa[1, :]) - np.min(vpa[1, :]),
      np.max(vpa[2, :]) - np.min(vpa[2, :])])
  # Function eigh returns vca st ext[2] >= ext[1] >= ext[0], but want the
  # ordering ext[0] >= ext[1] >= ext[2]. This is equivalent to a 90 degree
  # rotation about the y-axis and a shift along the z axis of 2 * ext[0], ie:
  # [[0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 2 * ext[0]], [0, 0, 0, 1]]
  # Now reorder the extents and the eigen vectors.
  prm_ext =  np.array([ext[2], ext[1], ext[0], 1])
  vrbMsg(2, 'prm_ext =\n' + str(prm_ext))
  # rot_trx rotates surface about centroid to AABB
  # off_trx takes surface centroid to origin
  off_trx = np.array([
      [1., 0., 0., -vc[0]],
      [0., 1., 0., -vc[1]],
      [0., 0., 1., -vc[2]],
      [0., 0., 0., 1.]])
  rot_trx = np.array([
      [evct[2, 0], evct[2, 1], evct[2, 2],    0.0],
      [evct[1, 0], evct[1, 1], evct[1, 2],    0.0],
      [-evct[0, 0], -evct[0, 1], -evct[0, 2], 0.0],
      [0., 0., 0., 1.]])
  prm_trx = np.matmul(rot_trx, off_trx)
  vrbMsg(2, 'prm_trx =\n' + str(prm_trx))
  return prm_trx, prm_ext
#end

def plmFromSurf(srf, pnm):
  """
  Computes a HuBMAP placement object for the given surface.
  Parameters
  ----------
  srf: trimesh surface in Edinburgh model space
  pnm: path name, either 'path_l' or 'path_s'
  Returns
  -------
  plm: placement object
  """

  plm = None
  vrbMsg(1, 'plmFromSurf pnm = ' + str(pnm))
  vrbMsg(1, 'srf.vertices[:3] =\n' + str(srf.vertices[:3]) + '\n...')
  tr_htoe = model['tr_htoe_l'] if pnm == 'path_l' else model['tr_htoe_s']
  # Transform surface to the HuBMAP space
  vrbMsg(2, 'tr_htoe = \n' + str(tr_htoe))
  tr_etoh = np.linalg.inv(tr_htoe)
  vrbMsg(2, 'tr_etoh = \n' + str(tr_etoh))
  if args.keep_edin_space:
    srf_h = srf
  else:
    srf_h = srf.apply_transform(tr_etoh)
    srf_h.export('surf_h.stl')
  #end
  # Compute oriented bounding box
  tr_sam, ext_sam = computeOBB(srf_h.vertices)
  vrbMsg(1, 'tr_sam = \n' + str(tr_sam))
  vrbMsg(1, 'ext_sam = ' + str(ext_sam))
  tr_isam = np.linalg.inv(tr_sam) # Sample AABB to OBB
  vrbMsg(1, 'tr_isam = \n' + str(tr_isam))
  # Compute transform angles going from AABB  to OBB
  # HuBMAP rotation order XYZ is Tx(Ty(Tz(x))) elsewhere this is ZYX
  rots = tm.transformations.euler_from_matrix(tr_isam, axes='szyx')
  vrbMsg(1, 'szyx rots = ' + str(rots))
  vrbMsg(1, '          = ' +
      str((radToDeg(rots[0]), radToDeg(rots[1]), radToDeg(rots[2]))))
  date = dt.datetime.now().strftime('%Y-%m-%d')
  target = (CCF_OWL + '#VH' +
      ('M' if 'm' in model['target'] else 'F') +
      ('Large' if 'l' in model['target'] else 'Small') +
      'Intestine')
  plm = {
      '@context': CCF_ONTOLOGY,
      '@id': 'TODO some-id',
      'placement': {
        '@id': 'TODO some-id',
        '@type': 'SpatialPlacement',
        'source': 'TODO source',
        'target': target,
        'x_scaling': 1,
        'y_scaling': 1,
        'z_scaling': 1,
        'scaling_units': 'ratio',
        'x_rotation': radToDeg(rots[2]),
        'y_rotation': radToDeg(rots[1]),
        'z_rotation': radToDeg(rots[0]),
        'rotation_order': 'XYZ',
        'rotation_units': 'degree',
        'x_translation': tr_isam[0,3],
        'y_translation': tr_isam[1,3],
        'z_translation': tr_isam[2,3],
        'translation_units': 'millimeter',
        'creation_date': date},
      '@type': 'SpatialEntity',
      'ccf_annotations': [
        'TODO'
        'http://purl.obolibrary.org/obo/UBERON_0001153',
        'http://purl.obolibrary.org/obo/UBERON_0001156'],
      'creator_first_name': 'TODO first-name',
      'creator_last_name': 'TODO last-name',
      'dimension_units': 'millimeter',
      'slice_count': 1,
      'slice_thickness': 1,
      'x_dimension': ext_sam[0],
      'y_dimension': ext_sam[1],
      'z_dimension': ext_sam[2],
      'creation_date': date,
      'creator': 'TODO first-name last-name'}
  return plm
#end

def writePlacement(plm):
  f = sys.stdout
  if args.output != '-':
    f = open(args.output, 'w') 
  #end
  jsn.dump(plm, f)
  if args.output != '-':
    f.close()
  #end
#end

def main():
  plm = None
  pnm = [None, None]
  itv = [-1, -1]
  vrbMsg(2, 'Reading Edinburgh HuBMAP ' + args.target +
      ' model configuration from file ' + args.conf)
  model['config'] = readJson(args.conf)
  vrbMsg(2, 'Reading Edinburgh HuBMAP ' + args.target + 
      ' large intestine path from file  ' + args.path_l)
  model['path_l'] = readJson(args.path_l)
  vrbMsg(2, 'Reading Edinburgh HuBMAP ' + args.target + 
      ' small intestine path from file  ' + args.path_s)
  model['path_s'] = readJson(args.path_s)
  vrbMsg(2, 'Reading interval from file ' + args.file)
  lbi = readJson(args.file)
  vrbMsg(2, 'Computing absolute path indices from landmark bounded ' +
      'proportional distances.')
  for i in range(0,2):
    pnm[i], itv[i] = pathIdxFromLBPD(model, lbi['intervals'][0][i])
  if (itv[0] < 0) or not bool(pnm[0]) or (itv[1] < 0) or not bool(pnm[1]):
    errMsg('Failed to find absolute path index interval on model.')
  #end
  vrbMsg(3, 'pnm = ' + str(pnm) + ', itv = ' + str(itv))

  vrbMsg(2, 'Computing surface segment bounded by the interval.')
  surf = surfFromItvls(model, pnm, itv)
  if   args.output_type == 's':
    surf.export(args.output)
  elif args.output_type == 'b':
    f = sys.stdout
    if not args.output == '-':
      f = open(args.output, 'wt')
    #end
    print(str(surf.bounds), file=f)
    if not args.output == '-':
      f.close()
    #end
  elif args.output_type == 'p':
    vrbMsg(2, 'Computing placement from surface segment.')
    plm = plmFromSurf(surf, pnm[0])
    writePlacement(plm)
  else:
    errMsg('Invalid output_type ' + args.output_type +
        ' specified, -h or --help for usage.')
  #end
#end

def parseArgs():
  global args
  parser = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,
    description = ('Maps Edinburgh GCA landmark bounded proportional ' +
        'distance intervals to HuBMAP RUI placements. Version: ' +
        SCRIPT_VERSION),
    epilog= fmtstr('*Instalation* and *Testing*') + HELP_TESTING)

  parser._action_groups.pop()
  required = parser.add_argument_group('required arguments')
  optional = parser.add_argument_group('optional arguments')
  required.add_argument('-c', '--conf', type=str, required=True,
      help='''
      Edinburgh GCA config file for the HuBMAP visible human male / female
      model.
      ''')
  required.add_argument('-p', '--path-s', type=str, required=True,
      help='''
      Midline path file for the HuBMAP visible human male / female
      small intestine model.
      ''')
  required.add_argument('-P', '--path-l', type=str, required=True,
      help='''
      Midline path file for the HuBMAP visible human male / female
      large intestine model. 
      ''')
  optional.add_argument('-e', '--keep-edin-space', action='store_true',
      default=False,
      help='''
      Keep in Edinburgh model space rather than transforming to the HuBMAP
      model space (may be useful for debugging, default: %(default)s).
      ''')
  optional.add_argument('-r', '--radius', type=str, default=DEFAULT_RADIUS,
    help='Cylinder radius, (default: %(default)s))')
  optional.add_argument('-t', '--target', type=str, default='m',
    help='Target HuBMAP model: m(ale) or (f)emale, (default: %(default)s))')
  optional.add_argument('-o', '--output', type=str, default='-',
    help='Output filename. ' +
         '(default: %(default)s)')
  optional.add_argument('-y', '--output-type', type=str, default='p',
    help='Output type: (b)ounding box, (p)lacement, ' +
        '(s)urface, (default: %(default)s))')
  optional.add_argument('-v', '--verbose', type=int, default=0,
      help='Verbose output for viewing parameter values and debugging. ' +
           'Range [-1 - ' + str(VERBOSE_MAX) + ']; ' +
           '0 - just warnings, ' +
           str(VERBOSE_MAX) + ' - lots of verbose output. ' +
           'Use -1 for no warnings.'
           '(default: %(default)s)')
  parser.add_argument('file',
      help='Edinburgh GCA landmark bounded proportional distance interval ' +
           '(JSON) file.')

  args = parser.parse_args()
  args.output_type = args.output_type.lower()
  if args.output_type not in ['b', 'p', 's']:
    errMsg('Invalid output type given, see usage (-h or --help)')
  #end
  args.target = args.target.lower()
  if args.target not in ['m', 'f']:
    errMsg('Invalid target given, see usage (-h or --help)')
  else:
    if args.target == 'm':
      model['tr_htoe_s'] = tr_htoe_ms;
      model['tr_htoe_l'] = tr_htoe_ml;
    else:
      model['tr_htoe_s'] = tr_htoe_fs;
      model['tr_htoe_l'] = tr_htoe_fl;
    #end
  #end
  if args.verbose < -1:
    args.verbose = -1
  elif args.verbose > VERBOSE_MAX:
    args.verbose = VERBOSE_MAX
  #end
  vrbMsg(1, 'Parameter values:\n' + str(args))
#end    

if __name__ == '__main__':
  if not isPython3():
    errMsg('Python 3 is required for this script.')
  #end
  parseArgs()
  vrbMsg(1, sys.argv[0] + ' version: ' + SCRIPT_VERSION)
  main()
  sys.exit(0)
#end

