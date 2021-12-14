import sys
lines = open(sys.argv[1], 'r').readlines()
n = len(lines)
sum_N = 0
for l in lines:
    sum_N += int(l.split()[2])
print( "%d %d" % (n, sum_N))
for l in lines:
    es = l.split()
    print("%s %f" % (es[0], float(es[2])/sum_N))
