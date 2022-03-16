#!/usr/bin/python3
# -*- coding: utf-8 -*-

from __future__ import print_function

"""
GCAExtractColonPath.py

Extracts the colon path from a Woolz segmented colon domain.
"""

__date__        = "September 2020"
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
import Wlz as w

libc = c.CDLL("libc.so.6")
libc.fopen.restype = c.POINTER(w.FILE)

args = None
dim3 = None
script_version = '0.0.4'
bspline_order = 3
err_num = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)

def isPython3():
  """
  Is this Python version 3.X.X

  Returns:
    True or False
  """

  return (sys.version_info[0] == 3)
#end

def errExit(msg):
  """
  Prints error message including the Woolz error string to stderr and
  then exits.

  Parameters:
    msg (string): Error message
  """

  e = re.sub(r'[^A-Z_0-9]', '', str(w.WlzStringFromErrorNum(err_num, None)))
  print(sys.argv[0] + ': Error - ' + msg + ' (' + e + ').', file=sys.stderr)
  sys.exit(1)
#end

def vrbMsg(msg):
  """
  If verbose flag is set then prints verbose message to stderr.

  Parameters:
    msg (string): Verbose message
  """

  if bool(args) and args.verbose:
    print(sys.argv[0] + ': ' + msg, file=sys.stderr)
  #end
#end

def vrbMsgEndpoints():
  """
  If verbose flag is set then prints path end points to stderr.
  """

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
  """
  Parses the command line arguments.
  """

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
  """
  This function is a hack as we can't use WlzMakeMain(). This is because
  Python ctypes can't pass unions by value.
    Parameters:
      otype (c_int): Woolz object type
      dom (WlzDomain *): The Woolz domain pointer
      val (WlzValue *): The Woolz values pointer
      dst_err_num (c_int *): Destination error pointer
    Returns
      obj (WlzObject *): New Woolz object
  """

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
  """
  This function is a hack as we can't use WlzFreeDomain. This is because
  Python ctypes can't pass unions by value.
    Parameters:
      dom (WlzDomain *): The Woolz domain pointer
  """

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
  """
  Reads a Woolz object fr4om the given file
    Parameters:
      fn (string): file name
    Returns:
      obj (WlzObject *): object read from file
  """

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
  """
  Subsamples the given object using the sampling scale
    Parameters:
      obj (WlzObject *): Given object
    Returns:
      obj_ss (WlzObject *): Subsampled Woolz object
  """

  global args, err_num, dim3
  vrbMsg('Subsampling object.') ###
  s = args.sample_scale
  if dim3:
    tt = w.WLZ_TRANSFORM_3D_AFFINE 
  else:
    tt = w.WLZ_TRANSFORM_2D_AFFINE
  #end
  st = w.WlzAffineTransformFromScale(tt, s, s, s, c.byref(err_num))
  if not bool(err_num):
    obj_ss = w.WlzAffineTransformObj(obj, st, w.WLZ_INTERPOLATION_NEAREST,
        c.byref(err_num))
    
    w.WlzFreeAffineTransform(st)
  #end
  return obj_ss
#end

def findFarPoints(obj):
  """
  Finds furthest seperated points in the given object's domain
    Parameters:
      obj (WlzObject *): Woolz object with domain
    Returns:
      p ([WlzDVertex3, WlzDVertex3]): Furthest points
  """

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
  """
  Fits a B-spline to the given ordered points along the path
    Parameters:
      path_pts (WlzObject *): Given points object with ordered points
    Returns:
      path_bs (WlzObject *): B-spline object
  """

  global args, bspline_order, err_num, dim3
  vrbMsg('Fitting B-spline to path.')
  st = None
  path_bs = None
  fac = 1.0;
  tr_type = w.WLZ_TRANSFORM_2D_AFFINE
  if dim3:
    fac = 2.0
    tr_type = w.WLZ_TRANSFORM_3D_AFFINE
  #end
  dom = w.WlzDomain(None)
  val = w.WlzValues(None)
  bs_smoothing = path_pts.contents.domain.pts.contents.nPoints
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

def vecCrossProduct3D(a, b):
  """
  Vector cross product c = a x b
    Parameters:
      a (WlzDVertex3): Vector a
      b (WlzDVertex3): Vector b
    Returns:
    c (WlzDVertex3): Vector c
  """

  c = w.WlzDVertex3()
  c.vtX = (a.vtY * b.vtZ) - (b.vtY * a.vtZ)
  c.vtY = (a.vtZ * b.vtX) - (b.vtZ * a.vtX)
  c.vtZ = (a.vtX * b.vtY) - (b.vtX * a.vtY)
  return c
#end

def vecDotProduct2D(a, b):
  """
  Vector cross product c = a . b
    Parameters:
      a (WlzDVertex2): Vector a
      b (WlzDVertex2): Vector b
    Returns:
    c (WlzDVertex2): Vector c
  """

  c = (a.vtX * b.vtX) + (a.vtY * b.vtY)
  return c
#end

def vecDotProduct3D(a, b):
  """
  Vector cross product c = a . b
    Parameters:
      a (WlzDVertex3): Vector a
      b (WlzDVertex3): Vector b
    Returns:
    c (WlzDVertex3): Vector c
  """

  c = (a.vtX * b.vtX) + (a.vtY * b.vtY) + (a.vtZ * b.vtZ)
  return c
#end

def vecSub2D(a, b):
  """
  Subtracts 2nd vector from the 1st c = a - b
    Parameters:
      a (WlzDVertex2): Vector a
      b (WlzDVertex2): Vector b
    Returns:
    c (WlzDVertex2): Vector c
  """

  c = w.WlzDVertex2()
  c.vtX = a.vtX - b.vtX
  c.vtY = a.vtY - b.vtY
  return c
#end

def vecSub3D(a, b):
  """
  Subtracts 2nd vector from the 1st c = a - b
    Parameters:
      a (WlzDVertex3): Vector a
      b (WlzDVertex3): Vector b
    Returns:
    c (WlzDVertex3): Vector c
  """

  c = w.WlzDVertex3()
  c.vtX = a.vtX - b.vtX
  c.vtY = a.vtY - b.vtY
  c.vtZ = a.vtZ - b.vtZ
  return c
#end

def vecScale2D(a, v):
  """
  Scales the given vector by a scalar c = a v
    Parameters:
      a (float): Scalar a
      v (WlzDVertex2): Vector v
    Returns:
    c (WlzDVertex2): Vector c
  """

  u  = w.WlzDVertex2()
  u.vtX = a * v.vtX
  u.vtY = a * v.vtY
  return u
#end

def vecScale3D(a, v):
  """
  Scales the given vector by a scalar c = a v
    Parameters:
      a (float): Scalar a
      v (WlzDVertex3): Vector v
    Returns:
    c (WlzDVertex3): Vector c
  """

  u  = w.WlzDVertex3()
  u.vtX = a * v.vtX
  u.vtY = a * v.vtY
  u.vtZ = a * v.vtZ
  return u
#end

def vecNormalise2D(v):
  """
  Normalises a vector c = v / ||v||
    Parameters:
      v (WlzDVertex2): Vector v
    Returns:
      c (WlzDVertex2): Vector c
  """

  l = m.sqrt(vecDotProduct2D(v, v))
  if l < 0.000001:
    l = 0.0
  else:
    l = 1.0 / l
  #end
  v.vtX *= l
  v.vtY *= l
  return v
#end

def vecNormalise3D(v):
  """
  Normalises a vector c = v / ||v||
    Parameters:
      v (WlzDVertex3): Vector v
    Returns:
      c (WlzDVertex3): Vector c
  """

  l = m.sqrt(vecDotProduct3D(v, v))
  if l < 0.000001:
    l = 0.0
  else:
    l = 1.0 / l
  #end
  v.vtX *= l
  v.vtY *= l
  v.vtZ *= l
  return v
#end

def initialNormal(t):
  """
  Computes the initial normal of the path
    Parameters:
      t (WlzDVertex3): First tangent vector
    Returns:
      r (WlzDVertex3): First normal vector
  """
  # Create v1 not parallel to t
  v = w.WlzDVertex3()
  v.vtX = 1.0
  v.vtY = 1.0
  v.vtZ = 1.0
  if t.vtX > t.vtY:
    if t.vtX > t.vtZ:
      v.vtX = 0.0
    else:
      v.vtZ = 0.0
    #end
  else:
    if t.vtY > t.vtZ:
      v.vtY = 0.0
    else:
      v.vtZ = 0.0
    #end
  #end
  # Compute normal
  s = vecCrossProduct3D(t, v)
  r = vecCrossProduct3D(s, t)
  vecNormalise3D(r)
  return r
#end

def doubleReflectionNormal(p0, p1, r0, t0, t1):
  """
  Computes a normal (after the first) using double reflection
    Parameters:
      p0 (WlzDVertex3): Previous point
      p1 (WlzDVertex3): Current point
      r0 (WlzDVertex3): Previous normal vector
      t0 (WlzDVertex3): Previous tangent vector
      t1 (WlzDVertex3): Current tangent vector
    Returns:
      r1 (WlzDVertex3): Current normal vector
  """

  v1 = vecSub3D(p1, p0)
  c1 = vecDotProduct3D(v1, v1)
  rL = vecSub3D(r0, vecScale3D((2.0 / c1) * vecDotProduct3D(v1, r0), v1))
  tL = vecSub3D(t0, vecScale3D((2.0 / c1) * vecDotProduct3D(v1, t0), v1))
  v2 = vecSub3D(t1, tL)
  c2 = vecDotProduct3D(v2, v2)
  r1 = vecSub3D(rL, vecScale3D((2.0 / c2) * vecDotProduct3D(v2, rL), v2))
  return r1
#end

def  createRotMinNormals(pts, tnt):
  """
  If dim3 is True then computes normals so as to minimise their rotation along
  their path, otherwise if dim3 is False just returns None 
    Parameters:
      pts (WlzPoints *): Points domain with path point positions
      tnt (WlzPoints *): Points domain with path tangent vectors
    Returns:
      nd (WlzDVertex3 []): array of normals (dim3 True) or None (dim3 False)
  """

  global err_num, dim3
  if dim3:
    n = pts.contents.nPoints
    pd = pts.contents.points.d3
    td = tnt.contents.points.d3
    nd = []
    nd.append(initialNormal(td[0]))
    for i in range(1, n):
      j = i - 1
      nrm = doubleReflectionNormal(pd[j], pd[i], nd[j], td[j], td[i])
      nd.append(nrm)
    #end
  else:
    nd = None
  #end
  return nd
#end


def writeWlzObject(fn, obj):
  """
  Writes the given Woolz object to a file with the given path
    Parameters:
      fn (string): Complete file path
      obj (WlzObject *): Woolz object
  """

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

def writeJsnUnitVec2D(f, ps, tv, pe):
  """
  Writes a 2D unit vector as ASCII
    Parameters:
      f (File): File opened for writting
      ps (string): Start string
      tv (WlzDVertex2): Vector
      pe (string): End string
  """

  nv = vecNormalise2D(tv)
  writeJsnVec2D(f, '{:1.8f}', ps, nv, pe)
#end

def writeJsnUnitVec3D(f, ps, tv, pe):
  """
  Writes a 3D unit vector as ASCII
    Parameters:
      f (File): File opened for writting
      ps (string): Start string
      tv (WlzDVertex3): Vector
      pe (string): End string
  """

  nv = vecNormalise3D(tv)
  writeJsnVec3D(f, '{:1.8f}', ps, nv, pe)
#end

def writeJsnVec2D(f, fmt, ps, pv, pe):
  """
  Writes a 2D vector as ASCII
    Parameters:
      f (File): File opened for writting
      ps (string): Start string
      tv (WlzDVertex2): Vector
      pe (string): End string
  """
  f.write(ps +
      '[' +
      fmt.format(pv.vtX) + ', ' +
      fmt.format(pv.vtY) +
      ']' +
      pe)
#end

def writeJsnVec3D(f, fmt, ps, pv, pe):
  """
  Writes a 3D vector as ASCII
    Parameters:
      f (File): File opened for writting
      ps (string): Start string
      tv (WlzDVertex3): Vector
      pe (string): End string
  """

  f.write(ps +
      '[' +
      fmt.format(pv.vtX) + ', ' +
      fmt.format(pv.vtY) + ', ' +
      fmt.format(pv.vtZ) +
      ']' +
      pe)
#end

def writeSplineData(fn, pts, tnt, nrm):
  """
  Writes B-spline path points, tangents and normals to a JSON file
    Parameters:
      fn (string): JSON file path
      pts (WlzPoints *): Woolz points domain with point positions
      tnt (WlzPoints *): Woolz points domain with tangent vectors
      nrm (WlzDVertex3 [] or None): Array of WlzDVertex3 normal vectors
          (dim3 True) or None (dim3 False)
  """

  global err_num, dim3
  vrbMsg('Writing B-spline path points and tangents to file.')
  n = pts.contents.nPoints
  with open(fn, 'w') as f:
    f.write('{\n')
    f.write('  "n": ' + str(n) + ',\n')
    f.write('  "points": [\n')
    if dim3:
      pd = pts.contents.points.d3
      for i in range(0, n - 1):
        writeJsnVec3D(f, '{:g}', '    ', pd[i], ',\n')
      #end
      writeJsnVec3D(f, '{:g}', '    ', pd[n - 1], '\n')
    else:
      pd = pts.contents.points.d2
      for i in range(0, n - 1):
        writeJsnVec2D(f, '{:g}', '    ', pd[i], ',\n')
      #end
      writeJsnVec2D(f, '{:g}', '    ', pd[n - 1], '\n')
    #end
    f.write('   ],\n')
    f.write('  "tangents": [\n')
    if dim3:
      pd = tnt.contents.points.d3
      for i in range(0, n - 1):
        writeJsnUnitVec3D(f, '    ', pd[i], ',\n')
      #end
      writeJsnUnitVec3D(f, '    ', pd[n - 1], '\n')
    else:
      pd = pts.contents.points.d2
      for i in range(0, n - 1):
        writeJsnVec2D(f, '{:g}', '    ', pd[i], ',\n')
      #end
      writeJsnVec2D(f, '{:g}', '    ', pd[n - 1], '\n')
    #end
    if dim3:
      f.write('   ],\n')
      f.write('  "normals": [\n')
    else:
      f.write('   ]\n')
    #end
    if dim3:
      pd = nrm
      for i in range(0, n - 1):
        writeJsnUnitVec3D(f, '    ', pd[i], ',\n')
      #end
      writeJsnUnitVec3D(f, '    ', pd[n - 1], '\n')
      f.write('   ]\n')
    #end
    f.write('}\n')
  #end
#end

def main():
  """
  Computes the B-spline path through a Woolz domain object then writes the
  path both as a Woolz object and as a JSON encoded file
  """

  global args, err_num, dim3
  vrbMsg('Entering main() args are: ' + str(args)) 
  # Read input domain ###
  obj_in = readWlzObj(args.in_domain)
  if bool(err_num):
    errExit('Failed to read input domain from ' + args.in_domain)
  #end
  if obj_in.contents.type == w.WLZ_2D_DOMAINOBJ:
    dim3 = False
  elif obj_in.contents.type == w.WLZ_3D_DOMAINOBJ:
    dim3 = True
  else:
    errExit('Input object must be either a 2 or 3D spatial domain object.')
  #end
  ### Sub-sample input domain ###
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
  vrbMsg('Compute B-spline path normals.')
  path_bs_nrm = createRotMinNormals(path_bs_pts, path_bs_tnt)
  writeSplineData(args.output + '-bs.jsn', path_bs_pts, path_bs_tnt,
      path_bs_nrm)
  if bool(err_num):
    errExit('Failed to write B-spline path points, tangents and normals ' +
        'to file.')
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
