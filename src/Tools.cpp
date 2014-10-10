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
//#include <csys/resource.h>

#include "Tools.h"
#include <list>
#include <vector>
#include "MemoryManager.h"

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
            printListAbbv(arrayOfLists[i], &printInt);
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

void printInt(int integer)
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

    \param cliques A linked list of cliques to return. <b>(only available when compiled 
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \param n the number of vertices in the input graph

*/

void runAndPrintStatsMatrix(long (*function)(char**,
                                             #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                             list<list<int>> &,
                                             #endif
                                             int),
                            const char* algName,
                            char** adjMatrix,
                            #ifdef RETURN_CLIQUES_ONE_BY_ONE
                            list<list<int>> const &cliques,
                            #endif
                            int n )
{
    fprintf(stderr, "%s: ", algName);
    fflush(stderr);

    clock_t start = clock();

    long cliqueCount = function(adjMatrix,
                                #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                cliques,
                                #endif
                                n);

    clock_t end = clock();

    fprintf(stderr, "%ld maximal cliques, ", cliqueCount);
    fprintf(stderr, "in %f seconds\n", (double)(end-start)/(double)(CLOCKS_PER_SEC));
    fflush(stderr);
}

/*! \brief execute an maximal clique listing algorithm
           that takes an adjacency list as input, time it,
           and print the number of maximal cliques found
           along with time information.

    \param function a function that computes all maximal cliques
                    and returns the number of maximal cliques found

    \param algName a zero-terminated character string that will
                   be printed with the statistics of the algorithm
                   run.

    \param adjListLinked the input graph in the array of linked lists
                         (adjacency list) format.

    \param adjListArray the input graph in the array of arrays
                         (adjacency list) format.

    \param cliques A linked list of cliques to return. <b>(only available when compiled 
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \param degree An array, indexed by vertex, containing the degree of that vertex.

    \param n the number of vertices in the input graph

*/

void runAndPrintStatsListList( long (*function)(vector<list<int>> const &, 
                                                int**, 
                                                #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                                list<list<int>> &,
                                                #endif
                                                int*, int),
                               const char* algName,
                               vector<list<int>> const &adjListLinked,
                               int** adjListArray,
                               #ifdef RETURN_CLIQUES_ONE_BY_ONE
                               list<list<int>> &cliques,
                               #endif
                               int* degree,
                               int n )
{
    fprintf(stderr, "%s: ", algName);
    fflush(stderr);

    clock_t const start = clock();

    long const cliqueCount = function(adjListLinked, 
                                      adjListArray, 
                                      #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                      cliques,
                                      #endif
                                      degree, n);

    clock_t const end = clock();

    fprintf(stderr, "%ld maximal cliques, ", cliqueCount);
    fprintf(stderr, "in %f seconds\n", (double)(end-start)/(double)(CLOCKS_PER_SEC));
    fflush(stderr);
}

/*! \brief Print the items in the linked list.

    \param linkedList A linked list.

    \param printFunc A function to print the data elements in
                     the linked list.
*/

void printList(list<int> const &linkedList, void (*printFunc)(int))
{
#ifdef DEBUG
    printf("printList...\n");
#endif
    int count = 0;
    for (int const value : linkedList)
    {
        printFunc(value);
        if(count != linkedList.size())
        {
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
        printf("... plus %d more", linkedList.size()-10);
    }

    printf("\n");
}

