#!  /bin/python

import sys;
#import regex;
from collections import defaultdict

invert = 0

if len(sys.argv) > 1:
    if sys.argv[1] == "--invert":
#        print "Inverting..."
        invert = 1

vertices = 0
edges = 0

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

if invert:
    edges = (vertices-1)*vertices/2 - edges

edges = 2*edges

print vertices
print edges

#edge = input()
#while len(edge) != 0:
#    if (edge[0] < edge[1]):
#        print "e", edge[0]+1, edge[1]+1
#    edge = input()

neighbors = defaultdict(list)

for line in sys.stdin:
    edge = line.split()
    v1 = int(edge[1])
    v2 = int(edge[2])
    if invert == 1:
        neighbors[v1].append(v2)
        neighbors[v2].append(v1)
    else:
        print str(v1-1) + "," + str(v2-1)
        print str(v2-1) + "," + str(v1-1)

if invert == 1:
    for i in range(0,vertices-1):
        neighbors[i].sort()
        nonneighbor=0
        for neighbor in neighbors[i]:
            while nonneighbor < neighbor:
                print "e", i+1, nonneighbor + 1
                nonneighbor = nonneighbor + 1
            nonneighbor = nonneighbor + 1

        while nonneighbor < vertices:
            print "e", i+1, nonneighbor + 1
            nonneighbor = nonneighbor + 1
