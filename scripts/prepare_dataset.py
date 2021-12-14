from common import *
import sys

fp = open(sys.argv[1], 'r')

for l in fp:
  es = l.strip().split()
  if len(es) == 1:
      conf = es[0]
      ixyz = int_from_xyz(conf)
  else:
      seq = es[0]
      eng = int(es[1])
      print seq, ixyz, eng
