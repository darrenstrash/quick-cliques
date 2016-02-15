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
        #print line
        continue

    if line.startswith("p"):
        splitline = line.split()
        startIndex = len(splitline) - 2
        vertices  = int(splitline[startIndex])
        edges     = int(splitline[startIndex+1])
        break; 

if invert == 1:
    edges = ((vertices-1)*vertices/2 - edges - vertices)*2
else:
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

    if v1 == v2:
        print "ERROR: found a loop in dimacs file."

    if invert == 1:
        neighbors[v1-1].append(v2-1)
        neighbors[v2-1].append(v1-1)
    else:
        print str(v1-1) + "," + str(v2-1)
        print str(v2-1) + "," + str(v1-1)

if invert == 1:
    for i in range(0,vertices):
        neighbors[i].sort()
        nonneighbor=0
        for neighbor in neighbors[i]:
            while nonneighbor < neighbor:
                if i != nonneighbor:
                    print str(i) + "," + str(nonneighbor)
                nonneighbor = nonneighbor + 1
            nonneighbor = nonneighbor + 1

        while nonneighbor < vertices:
            if i != nonneighbor:
                print str(i) + "," + str(nonneighbor)
            nonneighbor = nonneighbor + 1
