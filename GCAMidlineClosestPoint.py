#!/usr/bin/python3
###
# \file         GCAMidlineClosestPoint.py
# \author       Bill Hill
# \date         August 2024
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
# Simple script to find closest point in a JSON encoded mid-line
# to a given point.
###

from __future__ import print_function
import os
import re
import sys
import math as m
import json as jsn
import numpy as np
import argparse as ap
from datetime import datetime

args = None

SCRIPT_VERSION = '0.0.1'

file_req_text= '''
files:

Input sample files should be JSON encoded mil-line of the form:
{
  "n": 2870,
  "points": [
     [189.445, 240.376, 203.334],
     [189.639, 239.15, 203.18],
     ...
   ]
}
'''

def isPython3():
  """
  Is this python version 3?
  """
  return (sys.version_info[0] == 3)
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
  if args.verbose:
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
      Finds the closest point in a JSON encoded mid-line to a given point
      and then prints it using the format:
      <index of closest point> <decimal index of closest point> <closest point>
      ''' +
      'Version: ' + SCRIPT_VERSION,
    epilog=file_req_text)

  parser._action_groups.pop()
  required = parser.add_argument_group('required arguments')
  optional = parser.add_argument_group('optional arguments')
  required.add_argument('-p', '--point', type=str, required=True,
      nargs=3, metavar=('x', 'y', 'z'),
      help='''
      Given point of the form <x> <y> <z> (eg 1 23.4 5.6e+02)
      ''')
  ap._StoreAction(option_strings=['-p'], dest='point', nargs=3,
      const=None, default=None, type=float, choices=None,
      metavar=('x', 'y', 'z'))
  optional.add_argument('-v', '--verbose', action='store_true',
      default=False,
      help='Verbose output (mainly useful for debugging).')
  parser.add_argument('midline_file', 
      help='Input JSON encoded mid-line file (required).')
  args = parser.parse_args()
  args.point = [float(args.point[0]), float(args.point[1]),
                 float(args.point[2])]
  vrbMsg('Parameter values:\n' + str(args))
  return(args)
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


def main():
  global args

  args = parseArgs()
  vrbMsg(str(args))
  midline = readJson(args.midline_file)
  n1 = len(midline['points']) - 1
  i = closestPoint(midline['points'], args.point, 0, n1)
  vrbMsg(str(midline))
  print(str(i) + ' ' + str(float(i) / n1) + ' ' + str(midline['points'][i]))
#end


if __name__ == '__main__':
  if not isPython3():
    errMsg('Python 3 is required for this script.')
  #end
  # Proccess command line
  parseArgs()
  vrbMsg(sys.argv[0] + ' version: ' + SCRIPT_VERSION)
  main()
  sys.exit(0)
#end

