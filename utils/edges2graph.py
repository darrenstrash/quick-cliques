#!  /bin/python

import sys;
#import regex;
from collections import defaultdict

invert = 0

vertices = input();
edges    = input();

edges = edges/2

if invert:
    edges = (vertices-1)*vertices/2 - edges

print vertices, edges

neighbors = defaultdict(list)

for line in sys.stdin:
    edge = line.split(",")
    v1 = int(edge[0])
    v2 = int(edge[1])
    neighbors[v1].append(v2)

for i in range(0,vertices):
    neighbors[i].sort()
    nonneighbor=i+1
    for j in range(0, len(neighbors[i])):
        sys.stdout.write(str(neighbors[i][j] + 1));
        if j != len(neighbors[i]) - 1:
            sys.stdout.write(" ")
    sys.stdout.write('\n')
