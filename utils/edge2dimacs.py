#!  /bin/python

import sys;
#import regex;
from collections import defaultdict

invert = 0

if len(sys.argv) > 1:
    if sys.argv[1] == "--invert":
#        print "Inverting..."
        invert = 1

vertices = input();
edges    = input();

print "c"

edges = edges/2

if invert:
    edges = (vertices-1)*vertices/2 - edges

print "p edge", vertices, edges

#edge = input()
#while len(edge) != 0:
#    if (edge[0] < edge[1]):
#        print "e", edge[0]+1, edge[1]+1
#    edge = input()

neighbors = defaultdict(list)

for line in sys.stdin:
    edge = line.split(",")
    v1 = int(edge[0])
    v2 = int(edge[1])
    if (v1 < v2):
        if invert == 1:
            neighbors[v1].append(v2)
        else:
            print "e", v1+1, v2+1

if invert == 1:
    for i in range(0,vertices):
        neighbors[i].sort()
        nonneighbor=i+1
        for neighbor in neighbors[i]:
            while nonneighbor < neighbor:
                print "e", i+1, nonneighbor+1
                nonneighbor = nonneighbor + 1
            nonneighbor = nonneighbor + 1

        while nonneighbor < vertices:
            print "e", i+1, nonneighbor+1
            nonneighbor = nonneighbor + 1
