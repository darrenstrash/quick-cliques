#include "PartitionTools.h"

#include <cassert>
#include <cstdio>
#include <ctime>
#include <fstream> // ifstream
#include <iostream>
#include <cmath>

using namespace std;

void PartitionTools::ReadPartitionFile(
    string              const &sFilename,
    vector<int>               &vVertexToPartition,
    vector<vector<int>>       &vPartitions)
{

    vVertexToPartition.clear();
    vPartitions.clear();

    ifstream instream(sFilename.c_str());
    int partition(-1); // endvertices, to read edges.
    int numPartitions(0);
    int intsRead(0);
    while (instream >> partition) {
        numPartitions = max(partition+1,numPartitions);
        intsRead++;

        vVertexToPartition.push_back(partition);
    }

    vPartitions.resize(numPartitions);
    size_t largestPartition(0);
    for (size_t vertex=0; vertex < intsRead; ++vertex) {
        vPartitions[vVertexToPartition[vertex]].push_back(vertex);
        largestPartition = max(largestPartition, vPartitions[vVertexToPartition[vertex]].size());
    }

    cout << "Read        " << intsRead         << " integers." << endl << flush;
    cout << "And         " << numPartitions    << " partitions." << endl << flush;
    cout << "Largest has " << largestPartition << " vertices." << endl << flush;
}
