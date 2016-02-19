#!  /bin/python

import sys;
#import regex;
from collections import defaultdict

invert = 0

vertices = input();
edges    = input();

for line in sys.stdin:
    edge = line.split(",")
    v1 = int(edge[0])
    v2 = int(edge[1])
    print str(v1) + "\t" + str(v2)
