import sys
import numpy as np
from common import *

##### main #########
fn_gen = sys.argv[1]
shift = int(sys.argv[2])

lines = open(fn_gen, 'r').readlines()

rg2_lst = []
if sys.argv[3] == "int":
    for l in lines[shift:]:
        s = l.strip()
        rg2_lst.append( cal_Rg2(real_conf_from_int(s)) )
else:
    for l in lines[shift:]:
        s = l.strip()
        rg2 = cal_Rg2(real_conf_from_xyz(s))
        rg2_lst.append( rg2 )

lines = open(sys.argv[4], 'r').readlines()
n_states = int(lines[0].split()[0])
prob_RGs = np.zeros([n_states])
for i in range(1, n_states+1):
    es = lines[i].strip().split()
    prob_RGs[i-1] = float(es[1]) + 1e-10
    #assert( int(float(es[0])*2) == i-1 )

pv = np.ones(shape=prob_RGs.shape)
for rg2 in rg2_lst:
    #i = int(rg2*2)
    i = int(rg2)
    if i>=len(pv): i=len(pv)-1
    pv[i] += 1
pv /= np.sum(pv)
KL = np.sum(prob_RGs * np.log(prob_RGs/pv))

print KL, np.mean(rg2_lst), np.mean([x**2 for x in rg2_lst]), np.mean([x**3 for x in rg2_lst])
