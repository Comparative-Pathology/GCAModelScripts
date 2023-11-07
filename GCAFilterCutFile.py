#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
GCAFilterCutFile.py

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
from datetime import date
import json

args = None
script_version = '0.0.1'

def errExit(msg):
  """
  Prints error message to stderr and then exits.

  Parameters:
    msg (string): Error message
  """

  print(sys.argv[0] + ': Error - ' + msg + '.', file=sys.stderr)
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

def parseAndSetRange(obj, p):
  """
  Parses floating point ranges from a pair of comma separated values
  (either of which may be missing) in a string of the given attribute
  of the given object. The floating point range array replaces the
  objects string attribute.
  """
  s = getattr(obj, p).split(',')
  if len(s) == 1:
    s = [s[0], '']
  #end
  r = [float('nan'), float('nan')]
  for i in range(0, 2):
    try:
      r[i] = float(s[i])
    except:
      pass
    #end
  #end
  setattr(obj, p, r)
#end

def parseArgs():
  """
  Parses the command line arguments.
  """

  global args
  parser = ap.ArgumentParser(
      description=
      '''
      Filters a given cut record file by outputting only either matched
      or non-matched tiles. The tiles are matched against a series of
      plane parameter ranges and expressions.
      Ranges are specified by a pair of optional comma separated values,
      eg "1", "1,", "1,3.2" and ",3.14"; giving ranges
      "[1,NaN]", "[1,NaN]", "[1,3.2]" and "[NaN, 3.14]",
      where "NaN" is used to indicate that there is no limit.
      Distance expressions are matched if True eg "(d-4)%5==0" would
      match all tiles with distance d for which (d - 4) mod 5 = 0. 
      ''')
  parser.add_argument('-d', '--distance', type=str, default=',',
      help='Distance range.');
  parser.add_argument('-D', '--dist-expr', type=str, default='',
      help='Distance expression.');
  parser.add_argument('-i', '--inverse', action='store_true', default=False,
      help='Inverse filtering, ie output tiles that don\'t match.')
  parser.add_argument('-o', '--output', type=str, default='-',
      help='Filtered output cut record file.')
  parser.add_argument('-p', '--phi', type=str, default=',',
      help='Angle phi range.');
  parser.add_argument('-t', '--theta', type=str, default=',',
      help='Angle theta range.');
  parser.add_argument('-z', '--zeta', type=str, default=',',
      help='Angle zeta range.');
  parser.add_argument('-v', '--verbose', action='store_true', default=False,
      help='Verbose output (possibly useful for debugging).')
  parser.add_argument('-x', '--index', type=str, default=',',
      help='Index range.');
  parser.add_argument('input', default='', 
      help='Input cut record file.')
  args = parser.parse_args()
  for p in ['distance', 'index', 'phi', 'theta', 'zeta']:
    parseAndSetRange(args, p)
  #end
  vrbMsg('Parsed arguments are: ' +
      ' distance = ' + str(args.distance) +
      ' dist_expr = ' + str(args.dist_expr) +
      ' inverse = ' + str(args.inverse) +
      ' output = ' + str(args.output) +
      ' phi = ' + str(args.phi) +
      ' theta = ' + str(args.theta) +
      ' zeta = ' + str(args.zeta) +
      ' verbose = ' + str(args.verbose) +
      ' index = ' + str(args.index) +
      ' input = ' + str(args.input))
#end

def main():
  global args
  eps = 7./3 - 4./3 - 1 # float epsilon for the current machine

  # Read and parse input JSON file
  with open(args.input, 'r') as fin:
    cr = json.load(fin, strict=False)
  #end
  fout = sys.stdout
  if not args.output == '-':
    fout = open(args.output, 'w')
  #end
  # Modify the object parsed from the JSON file
  print('{', file=fout)
  print('"commandline": ' + cr['commandline'] + '; ' + ' '.join(sys.argv) +
      ',', file=fout)
  print('"creation_date": ' + '"' + date.today().ctime() + '",', file=fout)
  print('"tiles": [', file=fout)
  # Filter and output the cut tiles
  sep = '' # output plane separator to avoid comma after last
  for t in cr['tiles']:
    idx = t['index']
    pln = t['plane']
    dist = pln['dist']
    phi = pln['phi']
    theta = pln['theta']
    zeta = pln['zeta']
    d = dist
    match = (
        (m.isnan(args.index[0]) or (args.index[0] < idx + eps)) and
        (m.isnan(args.index[1]) or (args.index[1] > idx - eps)) and
        ((len(args.dist_expr) == 0) or bool(eval(args.dist_expr))) and
        (m.isnan(args.distance[0]) or (args.distance[0] < dist + eps)) and
        (m.isnan(args.distance[1]) or (args.distance[1] > dist - eps)) and
        (m.isnan(args.phi[0]) or (args.phi[0] < phi + eps)) and
        (m.isnan(args.phi[1]) or (args.phi[1] > phi - eps)) and
        (m.isnan(args.theta[0]) or (args.theta[0] < theta + eps)) and
        (m.isnan(args.theta[1]) or (args.theta[1] > theta - eps)) and
        (m.isnan(args.zeta[0]) or (args.zeta[0] < zeta + eps)) and
        (m.isnan(args.zeta[1]) or (args.zeta[1] > zeta - eps)))
    vrbMsg('tile index = ' + str(idx) + ' match = ' + str(match))
    if (match and (not args.inverse)) or (args.inverse and (not match)):
      print(sep + json.dumps(t), end='', file=fout)
      if len(sep) == 0:
        sep = ',\n'
      #end
   #end
  #end
  # Tidy up and close the cut record file
  print(']', file=fout)
  print('}', file=fout)
  if not args.output == '-':
    fout.close()
  #end
#end


if __name__ == '__main__':
  status = 0
  # Proccess command line
  parseArgs()
  vrbMsg(sys.argv[0] + ' version: ' + script_version)
  # Do the work
  main()
  sys.exit(0)
#end
