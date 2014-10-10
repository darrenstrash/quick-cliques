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

#include<assert.h>
#include<stdio.h>
#include<time.h>
#include<sys/resource.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"

/*! \file misc.c

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

int nodeComparator(void* node1, void* node2)
{
    if ((int)node1 < (int)node2)
        return -1;
    if((int)node1 > (int)node2)
        return 1;

    return 0;
}

/*! \brief compare integer pointers; return -1,0,1 for <,=,>;
           used for calling sort().

    \param node1 a pointer to an integer

    \param node2 a pointer to an integer

    \return -1 if <, 0 if =, and 1 if >.
*/

int sortComparator(void* node1, void* node2)
{
    if (*(int*)node1 < *(int*)node2)
        return -1;
    if(*(int*)node1 > *(int*)node2)
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

void printArrayOfLinkedLists(LinkedList** listOfLists, int size)
{
    // list graph contents

    int i=0;

    while(i<size)
    {
        if(!isEmpty(listOfLists[i]))
        {
            printf("%d:", i);
            printListAbbv(listOfLists[i], &printInt);
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

void printInt(void* integer)
{
    printf("%d", (int)integer);
}

/*! \brief destroy a linked list of integer arrays that have
           -1 in the last cell, have have been allocated by
           the user.

    \param cliques the linked list of integer arrays
*/

void destroyCliqueResults(LinkedList* cliques)
{
    Link* curr = cliques->head->next;
    while(!isTail(curr))
    {
        int* clique = (int*)curr->data;

#ifdef DEBUG
        int i=0;
        while(clique[i] != -1)
        {
            printf("%d", clique[i]);
            if(clique[i+1] != -1)
                printf(" ");
            i++;
        }
        printf("\n");
#endif
        Free(clique);
        curr = curr->next;
    } 

    destroyLinkedList(cliques); 
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

LinkedList** readInGraphAdjList(int* n, int* m)
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
    
    LinkedList** adjList = Calloc(*n, sizeof(LinkedList*));

    int i = 0;
    while(i < *n)
    {
        adjList[i] = createLinkedList();
        i++;
    }

    i = 0;

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

        addLast(adjList[u], (void*)v);

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
                                             LinkedList*,
                                             #endif
                                             int),
                            const char* algName,
                            char** adjMatrix,
                            #ifdef RETURN_CLIQUES_ONE_BY_ONE
                            LinkedList* cliques,
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

void runAndPrintStatsListList( long (*function)(LinkedList**, 
                                                int**, 
                                                #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                                LinkedList*,
                                                #endif
                                                int*, int),
                               const char* algName,
                               LinkedList** adjListLinked,
                               int** adjListArray,
                               #ifdef RETURN_CLIQUES_ONE_BY_ONE
                               LinkedList* cliques,
                               #endif
                               int* degree,
                               int n )
{
    fprintf(stderr, "%s: ", algName);
    fflush(stderr);

    clock_t start = clock();

    long cliqueCount = function(adjListLinked, 
                                adjListArray, 
                                #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                cliques,
                                #endif
                                degree, n);

    clock_t end = clock();

    fprintf(stderr, "%ld maximal cliques, ", cliqueCount);
    fprintf(stderr, "in %f seconds\n", (double)(end-start)/(double)(CLOCKS_PER_SEC));
    fflush(stderr);
}

/*! \brief process a clique, which may include printing it in
           one of several formats and/or adding the 
           clique to a linked list.

    \param cliques A linked list of cliques to return. <b>(only available when compiled 
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \param clique the clique to add to the linked list

*/

inline void processClique(
                          #ifdef RETURN_CLIQUES_ONE_BY_ONE
                          LinkedList *cliques,
                          #endif
                          LinkedList *clique)
{
    #ifdef PRINT_CLIQUES_TOMITA_STYLE
    printf("c ");
    #endif

    #ifdef PRINT_CLIQUES_ONE_BY_ONE
    printList(clique, &printInt);
    #endif

    #ifdef RETURN_CLIQUES_ONE_BY_ONE
    {
    int* cliqueArray = Calloc(length(clique) + 1, sizeof(int));
    int i = 0;

    Link* curr = clique->head->next;
    while(!isTail(curr))
    {
        cliqueArray[i] = (int)curr->data;
        i++;
        curr = curr->next;
    }

    cliqueArray[i] = -1;

    addLast(cliques, (void*)cliqueArray);
    }
    #endif
}
