import sys
from common import *

##### main #########
fn_gen = sys.argv[1]
shift = int(sys.argv[2])

lines = open(fn_gen, 'r').readlines()

if sys.argv[3] == "int":
    for l in lines[shift:]:
        s = l.strip()
        print cal_Rg2(real_conf_from_int(s))
else:
    for l in lines[shift:]:
        s = l.strip()
        print cal_Rg2(real_conf_from_xyz(s))
