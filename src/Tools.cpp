/* 
    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by 
    the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version. 
 
    This program is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
    GNU General Public License for more details. 
 
    You should have received a copy of the GNU General Public License 
    along with this program.  If not, see <http://www.gnu.org/licenses/> 
*/

#include <cassert>
#include <cstdio>
#include <ctime>
#include <fstream> // ifstream
//#include <csys/resource.h>

#include "Tools.h"
#include <list>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "MemoryManager.h"
#include "Algorithm.h"

using namespace std;

/*! \file Tools.cpp

    \brief A collection of useful comparators and print functions

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

/*! \brief compare integers return -1,0,1 for <,=,>

    \param node1 an integer

    \param node2 an integer

    \return -1 if <, 0 if =, and 1 if >.
*/

int nodeComparator(int node1, int node2)
{
    if (node1 < node2)
        return -1;
    if(node1 > node2)
        return 1;

    return 0;
}

/*! \brief compare integer pointers; return -1,0,1 for <,=,>;
           used for calling sort().

    \param node1 a pointer to an integer

    \param node2 a pointer to an integer

    \return -1 if <, 0 if =, and 1 if >.
*/

int sortComparator(int node1, int node2)
{
    if (node1 < node2)
        return -1;
    if(node1 > node2)
        return 1;

    return 0;
}

/*! \brief print an array of integers to standard out.

    \param array the array to print

    \param size the length of the array
*/

void printArray(int* array, int size)
{
    int i = 0;
    while(i<size)
        printf("%d ", array[i++]);
    printf("\n");
}

void printArrayWithIndexArrows(int* array, int size, int index1, int index2, int index3)
{
    printArray(array, size);
    int i = 0;
    while (i++ < index1)
        printf(" ");
    printf("^");

    while (i++ < index2)
        printf(" ");
    printf("^");

    while (i++ < index3)
        printf(" ");
    printf("^");

    printf("\n");
}

/*! \brief print an abbreviated version of an adjacency list

    \param listOfLists the adjacency list

    \param size the number of vertices in the graph
*/

void printArrayOfLinkedLists(vector<list<int>> const &arrayOfLists, int size)
{
    // list graph contents

    int i = 0;
    while (i < arrayOfLists.size())
    {
        if (!arrayOfLists[i].empty())
        {
            printf("%d:", i);
            printListAbbv(arrayOfLists[i], &Tools::printInt);
        }
        i++;
    }
}

/*! \brief print a clique, that is formatted as an integer
           array ending with -1.

    \param clique the clique.
*/

void printClique(int* clique)
{
    int i = 0;
    while(clique[i]!=-1)
    {
        printf("%d", clique[i]);
        if(clique[i+1]!=-1)
            printf(" ");
        i++;
    }
    printf("\n");
}

/*! \brief print an integer 

    \param integer an integer cast to a void*
*/

void Tools::printInt(int integer)
{
    printf("%d", integer);
}

/*! \brief destroy a linked list of integer arrays that have
           -1 in the last cell, have have been allocated by
           the user.

    \param cliques the linked list of integer arrays
*/

void destroyCliqueResults(list<list<int>> &cliques)
{
    cliques.clear();
}

/*! \brief read in a graph from stdin and return an 
           adjacency list, as an array of linked lists
           of integers.

    \param n this will be the number of vertices in the
             graph when this function returns.

    \param m this will be 2x the number of edges in the
             graph when this function returns.

    \return an array of linked lists of integers (adjacency list) 
            representation of the graph
*/

vector<list<int>> readInGraphAdjList(int* n, int* m)
{
    int u, v; // endvertices, to read edges.

    if(scanf("%d", n)!=1)
    {
        fprintf(stderr, "problem with line 1 in input file\n");
        exit(1);
    }

    if(scanf("%d", m)!=1)
    {
        fprintf(stderr, "problem with line 2 in input file\n");
        exit(1);
    }

#ifdef DEBUG
    printf("Number of vertices: %d\n", *n);
    printf("Number of edges: %d\n", *m);
#endif
    
    vector<list<int>> adjList(*n);

    int i = 0;
    while(i < *m)
    {
        if(scanf("%d,%d", &u, &v)!=2)
        {
            printf("problem with line %d in input file\n", i+2);
            exit(1);
        }
        assert(u < *n && u > -1);
        assert(v < *n && v > -1);
        if(u==v)
            printf("%d=%d\n", u, v);
        assert(u != v);

        adjList[u].push_back(v);

        i++;
    }

#ifdef DEBUG
    printArrayOfLinkedLists(adjList, *n);
#endif

    return adjList;
}

vector<list<int>> readInGraphAdjListEdgesPerLine(int &n, int &m, string const &fileName)
{
    ifstream instream(fileName.c_str());

    if (instream.good() && !instream.eof()) {
        string line;
        std::getline(instream, line);
////        cout << "Read Line: " << line << endl << flush;
        while((line.empty() || line[0] == '%') && instream.good() && !instream.eof()) {
            std::getline(instream, line);
        }
        stringstream strm(line);
        strm >> n >> m;
    } else {
        fprintf(stderr, "ERROR: Problem reading number of vertices and edges in file %s\n", fileName.c_str());
        exit(1);
    }

#ifdef DEBUG
    printf("Number of vertices: %d\n", n);
    printf("Number of edges: %d\n", m);
#endif
    
    vector<list<int>> adjList(n);

    int u, v; // endvertices, to read edges.
    int i = 0;
    while (i < n) {
        if (!instream.good()  || instream.eof()) {
            fprintf(stderr, "ERROR: Problem reading line %d in file %s\n", i+1, fileName.c_str());
            exit(1);
        }

        string line;
        std::getline(instream, line);
        u = i; // TODO/DS: remove.
        stringstream strm(line);
////        bool debug(true); ////u == 40656 || u == 40653);
////if (debug)        cout << (u+1) << " : " << endl << flush;
////if (debug)        cout << "Read     Line: " << line << endl << flush;
////if (debug)        cout << "Actually Read: ";
        while (!line.empty() && strm.good() && !strm.eof()) {
            strm >> v;
            ////if (!strm.good()) break;
////if (debug)            cout << v << " ";
            v--;

            assert(u < n && u > -1);
            assert(v < n && v > -1);
            if (u==v) {
                fprintf(stderr, "ERROR: Detected loop %d->%d\n", u + 1, v + 1);
                exit(1);
            }

            adjList[u].push_back(v);
        }
////if (debug)        cout << endl << flush;

        i++;
    }

#ifdef DEBUG
    printArrayOfLinkedLists(adjList, n);
#endif

    return adjList;
}


vector<list<int>> readInGraphAdjList(int &n, int &m, string const &fileName)
{

    ifstream instream(fileName.c_str());

    if (instream.good() && !instream.eof())
        instream >> n;
    else {
        fprintf(stderr, "problem with line 1 in input file\n");
        exit(1);
    }


    if (instream.good() && !instream.eof())
        instream >> m;
    else {

        fprintf(stderr, "problem with line 2 in input file\n");
        exit(1);
    }

#ifdef DEBUG
    printf("Number of vertices: %d\n", n);
    printf("Number of edges: %d\n", m);
#endif
    
    vector<list<int>> adjList(n);

    int u, v; // endvertices, to read edges.
    int i = 0;
    while(i < m)
    {
        char comma;
        if (instream.good() && !instream.eof()) {
            instream >> u >> comma >> v;
        } else {
            fprintf(stderr, "problem with line %d in input file\n", i+2);
            exit(1);
        }
        assert(u < n && u > -1);
        assert(v < n && v > -1);
        if(u==v)
            fprintf(stderr, "Detected loop %d->%d\n", u, v);
        assert(u != v);

        adjList[u].push_back(v);

        i++;
    }

#ifdef DEBUG
    printArrayOfLinkedLists(adjList, n);
#endif

    return adjList;
}

#if 0
vector<list<int>> readInGraphAdjListDimacs(int &n, int &m, string const &fileName)
{

    std::getline(instream, line);
    ifstream instream(fileName.c_str());

    if (instream.good() && !instream.eof())
        instream >> m;
    else {

        fprintf(stderr, "problem with line 2 in input file\n");
        exit(1);
    }

#ifdef DEBUG
    printf("Number of vertices: %d\n", n);
    printf("Number of edges: %d\n", m);
#endif
    
    vector<list<int>> adjList(n);

    int u, v; // endvertices, to read edges.
    int i = 0;
    while(i < m)
    {
        char comma;
        if (instream.good() && !instream.eof()) {
            instream >> u >> comma >> v;
        } else {
            fprintf(stderr, "problem with line %d in input file\n", i+2);
            exit(1);
        }
        assert(u < n && u > -1);
        assert(v < n && v > -1);
        if(u==v)
            fprintf(stderr, "Detected loop %d->%d\n", u, v);
        assert(u != v);

        adjList[u].push_back(v);

        i++;
    }

#ifdef DEBUG
    printArrayOfLinkedLists(adjList, n);
#endif

    return adjList;
}
#endif



/*! \brief execute an maximal clique listing algorithm
           that takes an adjacency matrix as input, time it,
           and print the number of maximal cliques found
           along with time information.

    \param function a function that computes all maximal cliques
                    and returns the number of maximal cliques found

    \param algName a zero-terminated character string that will
                   be printed with the statistics of the algorithm
                   run.

    \param adjMatrix the input graph in the adjacency matrix format.

    \param n the number of vertices in the input graph

*/

void runAndPrintStatsMatrix(long (*function)(char**,
                                             int),
                            const char* algName,
                            char** adjMatrix,
                            int n )
{
    fprintf(stderr, "%s: ", algName);
    fflush(stderr);

    clock_t start = clock();

    long cliqueCount = function(adjMatrix, n);

    clock_t end = clock();

    fprintf(stderr, "%ld maximal cliques, ", cliqueCount);
    fprintf(stderr, "in %f seconds\n", (double)(end-start)/(double)(CLOCKS_PER_SEC));
    fflush(stderr);
}

void RunAndPrintStats(Algorithm *pAlgorithm, list<list<int>> &cliques, bool const outputLatex)
{
    fprintf(stderr, "%s: ", pAlgorithm->GetName().c_str());
    fflush(stderr);

    clock_t start = clock();

    long const cliqueCount = pAlgorithm->Run(cliques);

    clock_t end = clock();

    if (!outputLatex) {
        fprintf(stderr, "%ld maximal cliques, ", cliqueCount);
        fprintf(stderr, "in %f seconds\n", (double)(end-start)/(double)(CLOCKS_PER_SEC));
    } else {
        printf("%.2f", (double)(end-start)/(double)(CLOCKS_PER_SEC));
    }
    fflush(stderr);
}

/*! \brief Print the items in the linked list.

    \param linkedList A linked list.

    \param printFunc A function to print the data elements in
                     the linked list.
*/

void Tools::printList(list<int> const &linkedList, void (*printFunc)(int))
{
#ifdef DEBUG
    printf("printList...\n");
#endif
    int count = 0;
    for (int const value : linkedList) {
        printFunc(value);
        if (count != linkedList.size()) {
            printf(" ");
        }
    }

    printf("\n");
   
}

/*! \brief Print the first 10 items in the linked list

    \param linkedList A linked list.

    \param printFunc A function to print the data elements in
                     the linked list.
*/

void printListAbbv(list<int> const &linkedList, void (*printFunc)(int))
{
#ifdef DEBUG
    printf("printListAbbv...\n");
#endif
    int count = 0;

    for (list<int>::const_iterator cit = linkedList.begin();
         cit != linkedList.end() && count != 10; ++cit)
    {
        count++;
        printFunc(*cit);
        if(count != linkedList.size())
        {
            printf(" ");
        }
    }

    if(count != linkedList.size())
    {
        printf("... plus %lu more", linkedList.size()-10);
    }

    printf("\n");
}

void DescribeVertex(int const lineNumber, int *vertexSets, int *vertexLookup, int const size, int const vertex, int const beginX, int const beginD, int const beginP, int const beginR)
{
    int const vertexLocation(vertexLookup[vertex]);

    cout << lineNumber << ": vertex " << vertex << " is in position " << vertexLocation << (vertexSets[vertexLocation] == vertex ? "(consistent)" : "(inconsistent: " + to_string(vertexSets[vertexLocation]) + " is there)" ) << " in set ";

    if (vertexLocation < beginX) {
        cout << "(before X)" << endl;
    }
    if (vertexLocation >= beginX && vertexLocation < beginD) {
        cout << "X" << endl;
    }
    if (vertexLocation >= beginD && vertexLocation < beginP) {
        cout << "D" << endl;
    }
    if (vertexLocation >= beginP && vertexLocation < beginR) {
        cout << "P" << endl;
    }
    if (vertexLocation >= beginR) {
        cout << "R" << endl;
    }
}

void DescribeSet(string const &setName, int const begin, int const end)
{
    cout << " " << setName << "=[" << begin << "->" << end << "]";
}

void DescribeState(int const lineNumber, int *vertexSets, int *vertexLookup, int const size, int const beginX, int const beginD, int const beginP, int const beginR)
{
    cout << lineNumber << ": Size " << size;
    DescribeSet("X", beginX, beginD-1);
    DescribeSet("D", beginD, beginP-1);
    DescribeSet("P", beginP, beginR-1);
    DescribeSet("R", beginR, size-1);
    cout << endl << flush;
}

void CheckConsistency(int const lineNumber, size_t const recursionNumber, int *vertexSets, int *vertexLookup, int const size)
{
    //if (recursionNumber > 2) return;
    //cout << recursionNumber << "1 is in position " << ver
    for (int i=0; i < size; ++i) {
        if (vertexSets[vertexLookup[i]] != i) {
            cout << recursionNumber << "(line " << lineNumber << ") : inconsistency -- vertex " << i  << " is supposed to be in position " << vertexLookup[i] << " but vertex " <<  vertexSets[vertexLookup[i]] << " is there." << endl;
        }
    }
}

void CheckReverseConsistency(int const lineNumber, size_t const recursionNumber, int *vertexSets, int *vertexLookup, int const size)
{
    //if (recursionNumber > 2) return;
    //cout << recursionNumber << "1 is in position " << ver
    for (int i=0; i < size; ++i) {
        if (vertexLookup[vertexSets[i]] != i) {
            cout << recursionNumber << "(line " << lineNumber << ") : inconsistency -- vertex " << vertexSets[i]  << " is supposed to be in position " << vertexLookup[vertexSets[i]] << " but it is in position " << i << "." << endl;
        }
    }
}

////bool IsMaximalClique(std::list<int> const &clique, std::vector<std::vector<int>> const &adjacencyList)
////{
////    set<int>
////    // is clique
////    int index;
////    for (int const vertex : partialClique) {
////    }
////
////    // is maximal
////}


void InvertGraph(vector<list<int>> const &adjList)
{
    int const n(adjList.size());
    cout << n << endl;
    size_t numEdgesInInverse(0);
    for (list<int> const &neighbors : adjList) {
        numEdgesInInverse += n - neighbors.size() - 1; // all non-edges except loops
    }

    cout << numEdgesInInverse << endl;

    for (int i = 0; i < adjList.size(); ++i) {
        set<int> setNeighbors;
        setNeighbors.insert(adjList[i].begin(), adjList[i].end());
        for (int neighbor=0; neighbor < adjList.size(); neighbor++) {
            if (setNeighbors.find(neighbor) == setNeighbors.end() && neighbor != i) {
                cout << "(" << i << "," << neighbor << i << ")" << endl;
            }
        }
    }
}

string Tools::GetTimeInSeconds(clock_t delta, bool const brackets) {
    stringstream strm;

    strm.precision(2);
    strm.setf(std::ios::fixed, std::ios::floatfield);
    if (brackets) {
        strm << "[" << (double)(delta)/(double)(CLOCKS_PER_SEC) << "s]";
    } else {
        strm << (double)(delta)/(double)(CLOCKS_PER_SEC) << "s";
    }
    return strm.str();
}

vector<int> Tools::ReadMetisOrdering(string const &fileName)
{
    ifstream instream(fileName.c_str());

    vector<int> ordering;

    int u, v; // endvertices, to read edges.

    if (!instream.good()  || instream.eof()) {
        fprintf(stderr, "ERROR: Problem reading line 1 in file %s\n", fileName.c_str());
        exit(1);
    }

    while (instream.good() && !instream.eof()) {

        string line;
        std::getline(instream, line);
        stringstream strm(line);
////        bool debug(true); ////u == 40656 || u == 40653);
////if (debug)        cout << (u+1) << " : " << endl << flush;
////if (debug)        cout << "Read     Line: " << line << endl << flush;
////if (debug)        cout << "Actually Read: ";
        if (!line.empty() && strm.good() && !strm.eof()) {
            strm >> v;
            ////if (!strm.good()) break;
////if (debug)            cout << v << " ";
            cout << "read: " << v << endl;

            assert(v > -1);

            ordering.push_back(v);
        }
////if (debug)        cout << endl << flush;
    }

#ifdef DEBUG
    printArrayOfLinkedLists(adjList, n);
#endif

    return ordering;
}

