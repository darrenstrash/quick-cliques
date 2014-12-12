#!  /bin/python

import sys;
#import regex;
from collections import defaultdict

invert = 1

for line in sys.stdin:
    if (line.startswith("c")):
        print line
        continue

    if line.startswith("p"):
        splitline = line.split()
        startIndex = len(splitline) - 2
        vertices  = int(splitline[startIndex])
        edges     = int(splitline[startIndex+1])
        break; 

print "c"

if invert:
    edges = (vertices-1)*vertices/2 - edges

print "p edge", vertices, edges

neighbors = defaultdict(list)

for line in sys.stdin:
    if (line.startswith("c")):
        continue
    edge = line.split()
    v1 = int(edge[1])
    v2 = int(edge[2])
    if (v1 < v2):
        neighbors[v1].append(v2)

if invert == 1:
    for i in range(1,vertices):
        neighbors[i].sort()
        nonneighbor=i+1
        for neighbor in neighbors[i]:
            while nonneighbor < neighbor:
                print "e", i, nonneighbor
                nonneighbor = nonneighbor + 1
            nonneighbor = nonneighbor + 1

        while nonneighbor <= vertices:
            print "e", i, nonneighbor
            nonneighbor = nonneighbor + 1
