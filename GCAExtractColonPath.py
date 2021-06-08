#!/usr/bin/python3
###!
# \file         GCAExtractColonPath.py
# \author       Bill Hill
# \date         September 2020
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
# Extracts the colon path from a Woolz segmented colon domain.
###

from __future__ import print_function
import sys
import re
import math as m
import argparse as ap
import ctypes as c
import Wlz as w

libc = c.CDLL("libc.so.6")
libc.fopen.restype = c.POINTER(w.FILE)

args = None
script_version = '0.0.3'
bspline_order = 3
err_num = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)

def isPython3():
  return (sys.version_info[0] == 3)

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

def vrbMsgEndpoints():
  if bool(args) and args.verbose:
    p0 = args.end_points[0]
    p1 = args.end_points[1]
    print(sys.argv[0] + ': Endpoints [(' +
      str(p0.vtX) + ',' + str(p0.vtY) + ',' + str(p0.vtZ) + '),('+
      str(p1.vtX) + ',' + str(p1.vtY) + ',' + str(p1.vtZ) + ')]',
      file=sys.stderr)
  #end
#end

def parseArgs():
  global args
  parser = ap.ArgumentParser(
      description=
      'Extracts the colon path from a Woolz segmented colon domain.')
  parser.add_argument('-e', '--end-points', type=str, default=None,
      help='End points for the path.')
  parser.add_argument('-o', '--output', type=str, default='',
      help='Output filename base.')
  parser.add_argument('-s', '--sample-scale', type=float, default=0.25,
      help='Sub-sampling scale.')
  parser.add_argument('-v', '--verbose', action='store_true', default=False,
      help='Verbose output (possibly useful for debugging).')
  parser.add_argument('in_domain', default='', 
      help='Input segmented colon domain.')
  args = parser.parse_args()
  if bool(args.end_points):
    try:
      s0 = args.end_points.split(',')
      e0 = [w.WlzIVertex3(int(s0[0]), int(s0[1]), int(s0[2])),
            w.WlzIVertex3(int(s0[3]), int(s0[4]), int(s0[5]))]
      args.end_points = e0
    except:
      args.end_points = None
      raise ap.ArgumentTypeError('Endpoints must be of the form ' +
        'x0,y0,z0,x1,y1,z1.') 
    #end
  #end
#end

def makeWlzMain(otype, dom, val, dst_err_num):
  # This function is a hack as 't use WlzMakeMain(). This is because
  # Python ctypes can't pass unions by value.
  obj = w.WlzMakeEmpty(dst_err_num)
  if bool(obj):
    obj.type = otype
    obj.contents.type = otype
    obj.contents.domain = dom
    obj.contents.value = val
    if bool(dom.core):
      dom.core.contents.linkcount = dom.core.contents.linkcount + 1
    #end
    if bool(val.core):
      val.core.contents.linkcount = val.core.contents.linkcount + 1
    #end
  #end
  return obj
#end

def freeWlzDomain(dom):
  # This function is a hack as 't use WlzFreeDomain. This is because
  # Python ctypes can't pass unions by value.
  if bool(dom.core):
    if dom.core.contents.linkcount >= 0:
      dom.core.contents.linkcount = dom.core.contents.linkcount - 1
      if dom.core.contents.linkcount <= 0:
        # Most (but not all) domains are free'd like this
        f = dom.core.contents.freeptr
        if bool(f):
          w.AlcFreeStackFree(f)
        #end
        w.AlcFree(dom.core)
      #end
    #end
  #end
#end

def readWlzObj(fn):
  global err_num
  vrbMsg('Reading Woolz object from file ' + fn)
  obj = None
  fn = fn.encode('ascii') if (isPython3) else fn
  err_num = w.enum__WlzErrorNum(w.WLZ_ERR_FILE_OPEN)
  fp = libc.fopen(fn, 'rb')
  if bool(fp):
    obj = w.WlzReadObj(fp, c.byref(err_num))
    libc.fclose(fp)
  #end
  vrbMsg('obj = ' + str(obj) + ', err_num = ' +
      str(w.WlzStringFromErrorNum(err_num, None)) + ')')
  return obj
#end

def subsampleObj(obj):
  global args, err_num
  vrbMsg('Subsampling object.') ###
  s = args.sample_scale
  st = w.WlzAffineTransformFromScale(w.WLZ_TRANSFORM_3D_AFFINE,
      s, s, s, c.byref(err_num))
  if not bool(err_num):
    obj_ss = w.WlzAffineTransformObj(obj, st, w.WLZ_INTERPOLATION_NEAREST,
        c.byref(err_num))
    
    w.WlzFreeAffineTransform(st)
  #end
  return obj_ss
#end

def findFarPoints(obj):
  global err_num
  vrbMsg('Finding farthest points in object.')
  p = [w.WlzDVertex3(0.0, 0.0, 0.0), w.WlzDVertex3(0.0, 0.0, 0.0)]
  p_obj = None
  dst_obj = None
  vrbMsg('Finding point farthest from boundary.')
  bnd_obj = w.WlzBoundaryDomain(obj, c.byref(err_num))
  if not bool(err_num):
    dst_obj = w.WlzDistanceTransform(obj, bnd_obj, w.WLZ_OCTAGONAL_DISTANCE,
        0.0, 0.0, c.byref(err_num))
  #end
  if not bool(err_num):
    cen = w.WlzGreyExtremumPos(dst_obj, 1, None, c.byref(err_num))
  #end
  w.WlzFreeObj(bnd_obj)
  bnd_obj = None
  w.WlzFreeObj(dst_obj)
  dst_obj = None
  if not bool(err_num):
    vrbMsg('Finding 1st farthest point.')
    p_obj = w.WlzMakeSinglePixelObject(w.WLZ_3D_DOMAINOBJ,
        cen.vtX, cen.vtY, cen.vtZ, c.byref(err_num))
  #end
  if not bool(err_num):
    dst_obj = w.WlzDistanceTransform(obj, p_obj, w.WLZ_OCTAGONAL_DISTANCE,
        0.0, 0.0, c.byref(err_num))
  #end
  if not bool(err_num):
    p[0] = w.WlzGreyExtremumPos(dst_obj, 1, None, c.byref(err_num))
  #end
  w.WlzFreeObj(dst_obj)
  dst_obj = None
  w.WlzFreeObj(p_obj)
  p_obj = None
  if not bool(err_num):
    vrbMsg('Finding 2nd farthest point.')
    p_obj = w.WlzMakeSinglePixelObject(w.WLZ_3D_DOMAINOBJ,
        p[0].vtX, p[0].vtY, p[0].vtZ, c.byref(err_num))
  #end
  if not bool(err_num):
    dst_obj = w.WlzDistanceTransform(obj, p_obj, w.WLZ_OCTAGONAL_DISTANCE,
        0.0, 0.0, c.byref(err_num))
  #end
  if not bool(err_num):
    p[1] = w.WlzGreyExtremumPos(dst_obj, 1, None, c.byref(err_num))
  #end
  w.WlzFreeObj(dst_obj)
  w.WlzFreeObj(p_obj)           
  if not bool(err_num):
    vrbMsg('Found farthest points [(' +
      str(p[0].vtX) + ',' + str(p[0].vtY) + ',' + str(p[0].vtZ) + '), (' +
      str(p[1].vtX) + ',' + str(p[1].vtY) + ',' + str(p[1].vtZ) + ')].')
  #end
  return p
#end

def fitBSpline(path_pts):
  global args, bspline_order, err_num
  vrbMsg('Fitting B-spline to path.')
  st = None
  path_bs = None
  fac = 1.0;
  tr_type = w.WLZ_TRANSFORM_2D_AFFINE
  dom = w.WlzDomain(None)
  val = w.WlzValues(None)
  bs_smoothing = path_pts.contents.domain.pts.contents.nPoints
  if ((path_pts.contents.domain.pts.contents.type == w.WLZ_POINTS_3I) or
      (path_pts.contents.domain.pts.contents.type == w.WLZ_POINTS_3D)):
    fac = 2.0
    w.WLZ_TRANSFORM_3D_AFFINE
  #end
  bs_smoothing = ((fac + bs_smoothing) / fac)
  vrbMsg('Smoothing factor = ' + str(bs_smoothing))
  dom.bs = w.WlzBSplineFromObj(path_pts, bspline_order, 0, bs_smoothing,
      c.byref(err_num))
  if not bool(err_num):
    path_bs = makeWlzMain(w.WLZ_SPLINE, dom, val, c.byref(err_num))
  #end
  if bool(err_num):
    w.WlzFreeDomain(dom)
    dom = w.WlzDomain(None)
  else:
    vrbMsg('Upscaling B-spline.') 
    s = 1.0 / args.sample_scale
    st = w.WlzAffineTransformFromScale(tr_type,
        s, s, s, c.byref(err_num))
    if not bool(err_num):
      t_bs = w.WlzAffineTransformObj(path_bs, st, w.WLZ_INTERPOLATION_NEAREST,
          c.byref(err_num))
      w.WlzFreeAffineTransform(st)
      w.WlzFreeObj(path_bs)
      path_bs = t_bs
    #end
  #end
  return path_bs
#end

def writeWlzObject(fn, obj):
  global err_num
  vrbMsg('Writing Woolz object to file ' + fn)
  fn = fn.encode('ascii') if (isPython3) else fn
  err_num = w.enum__WlzErrorNum(w.WLZ_ERR_FILE_OPEN)
  fp = libc.fopen(fn, 'wb')
  if(bool(fp)):
    err_num = w.WlzWriteObj(fp, obj)
    libc.fclose(fp)
  #end
  vrbMsg('err_num = ' + str(w.WlzStringFromErrorNum(err_num, None)))
#end

def writeJsnTnt2D(f, ps, tv, pe):
  l = m.sqrt(tv.vtX * tv.vtX + tv.vtY * tv.vtY)
  if l < 0.000001:
    l = 1.0
  #end
  nv = w.WlzDVertex2()
  nv.vtX = tv.vtX / l
  nv.vtY = tv.vtY / l
  writeJsnPnt2D(f, '{:1.8f}', ps, nv, pe)
#end

def writeJsnTnt3D(f, ps, tv, pe):
  l = m.sqrt(tv.vtX * tv.vtX + tv.vtY * tv.vtY + tv.vtZ * tv.vtZ)
  if l < 0.000001:
    l = 1.0
  #end
  nv = w.WlzDVertex3()
  nv.vtX = tv.vtX / l
  nv.vtY = tv.vtY / l
  nv.vtZ = tv.vtZ / l
  writeJsnPnt3D(f, '{:1.8f}', ps, nv, pe)
#end

def writeJsnPnt2D(f, fmt, ps, pv, pe):
  f.write(ps +
      '[' +
      fmt.format(pv.vtX) + ', ' +
      fmt.format(pv.vtY) +
      ']' +
      pe)
#end

def writeJsnPnt3D(f, fmt, ps, pv, pe):
  f.write(ps +
      '[' +
      fmt.format(pv.vtX) + ', ' +
      fmt.format(pv.vtY) + ', ' +
      fmt.format(pv.vtZ) +
      ']' +
      pe)
#end

def writeSplineData(fn, pts, tnt):
  global err_num
  vrbMsg('Writing B-spline path points and tangents to file.')
  n = pts.contents.nPoints
  dim3 = (pts.contents.type == w.WLZ_POINTS_3D)
  with open(fn, 'w') as f:
    f.write('{\n')
    f.write('  "n": ' + str(n) + ',\n')
    f.write('  "points": [\n')
    if dim3:
      pd = pts.contents.points.d3
      for i in range(0, n - 1):
        writeJsnPnt3D(f, '{:g}', '    ', pd[i], ',\n')
      #end
      writeJsnPnt3D(f, '{:g}', '    ', pd[n - 1], '\n')
    else:
      pd = pts.contents.points.d2
      for i in range(0, n - 1):
        writeJsnPnt2D(f, '{:g}', '    ', pd[i], ',\n')
      #end
      writeJsnPnt2D(f, '{:g}', '    ', pd[n - 1], '\n')
    #end
    f.write('   ],\n')
    f.write('  "tangents": [\n')
    if dim3:
      pd = tnt.contents.points.d3
      for i in range(0, n - 1):
        writeJsnTnt3D(f, '    ', pd[i], ',\n')
      #end
      writeJsnTnt3D(f, '    ', pd[n - 1], '\n')
    else:
      pd = pts.contents.points.d2
      for i in range(0, n - 1):
        writeJsnPnt2D(f, '{:g}', '    ', pd[i], ',\n')
      #end
      writeJsnPnt2D(f, '{:g}', '    ', pd[n - 1], '\n')
    #end
    f.write('   ]\n')
    f.write('}\n')
  #end
#end

def main():
  global args, err_num
  vrbMsg('Entering main() args are: ' + str(args)) 
  # Read input domain ###
  obj_in = readWlzObj(args.in_domain)
  if bool(err_num):
    errExit('Failed to read input domain from ' + args.in_domain)
  #end
  ### Sub-sample iput domain ###
  obj_ss = subsampleObj(obj_in)
  if bool(err_num):
    errExit('Failed to subsample input domain.')
  #end
  if(bool(args.end_points) and (len(args.end_points) == 2)):
    for i in range(0,2):
      s = args.sample_scale
      p = args.end_points[i]
      v = w.WlzIVertex3(int(s * p.vtX), int(s * p.vtY), int(s * p.vtZ))
      args.end_points[i] = v
    #end
  else:
    ### Finding endpoints within sub-sampled domain ###
    end_points = findFarPoints(obj_ss)
    if bool(err_num):
      errExit('Failed to find furthest points.')
    else:
      args.end_points = end_points
    #end
  #end
  vrbMsgEndpoints()
  # Compute minimum distance path between the endpoints
  vrbMsg('Computing minimum distance path.')
  path_pts = w.WlzLineSkeletonSegment(obj_ss, w.WLZ_POINTS,
          args.end_points[0], args.end_points[1], c.byref(err_num))
  if(bool(err_num)):
    errExit('Failed to compute line skeleton segment.')
  #end
  # Fit B-spline to the path
  path_bs = fitBSpline(path_pts)
  if bool(err_num):
    errExit('Failed to compute B-spline path.')
  #end
  # Free input object and the subsampled object
  w.WlzFreeObj(obj_in)
  w.WlzFreeObj(obj_ss)
  # Output the B-spline path object
  fn = args.output + '-bs.wlz'
  writeWlzObject(fn, path_bs)
  if bool(err_num):
    errExit('Failed to write output path B-spline to file ' + fn)
  #end
  # Compute the legth of the path
  vrbMsg('Compute length of the B-spline path.')
  err_num = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
  path_len = w.WlzBSplineLength(path_bs.contents.domain.bs, 0.0, 1.0,
      c.byref(err_num))
  if bool(err_num):
    errExit('Failed to compute the legth of the path.')
  #end
  path_len = int(path_len)
  vrbMsg('B-spline path length = ' + str(path_len))
  # Output an evaluation of the B-spline path and it's derivatives
  vrbMsg('Compute B-spline path as points.')
  path_bs_pts = w.WlzBSplineEvalPoints(path_bs.contents.domain.bs, path_len, 0,
      c.byref(err_num))
  if bool(err_num):
    errExit('Failed to compute B-spline path as points.')
  #end
  vrbMsg('Compute B-spline path tangents.')
  path_bs_tnt = w.WlzBSplineEvalPoints(path_bs.contents.domain.bs, path_len, 1,
      c.byref(err_num))
  if bool(err_num):
    errExit('Failed to compute B-spline path tangents.')
  #end
  writeSplineData(args.output + '-bs.jsn', path_bs_pts, path_bs_tnt)
  if bool(err_num):
    errExit('Failed to write B-spline path points and tangents to file.')
  #end
  dom = w.WlzDomain(None)
  dom.pts = path_bs_pts
  freeWlzDomain(dom)
  dom.pts = path_bs_tnt
  freeWlzDomain(dom)
  w.WlzFreeObj(path_bs)
#end

if __name__ == '__main__':
  status = 0
  # Proccess command line
  parseArgs()
  vrbMsg(sys.argv[0] + ' version: ' + script_version)
  vrbMsg('Woolz version: ' + re.sub(r'^b\'(.*)\'', r'\1', str(w.WlzVersion())))
  # Check version of Woolz is >= 1.8.0
  wv = [int(x) for x in re.sub(r'[^0-9.]','', str(w.WlzVersion())).split('.')]
  if (wv[0] < 1) or ((wv[0] < 2) and (wv[1] < 8)):
    err_num = w.enum__WlzErrorNum(w.WLZ_ERR_UNIMPLEMENTED)
    errExit('Requires PyWoolz which has been built with Woolz version >= 1.8.0')
  #end
  # Do the work
  main()
  sys.exit(0)
#end
