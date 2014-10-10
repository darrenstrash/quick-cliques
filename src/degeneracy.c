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
#include<stdlib.h>
#include<time.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include"degeneracy_algorithm.h"

/*! \file degeneracy.c

    \brief Execute the algorithm in degeneracy_algorithm.c
           and print the number of cliques found and wall clock
           execution time.

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

int main()
{

#ifdef MEMORY_DEBUG
    fprintf(stderr, "WARNING: MEMORY_DEBUG is defined, timing will be off.\n");
#endif

    int n; // number of vertices
    int m; // 2x number of edges

    LinkedList** adjacencyList = readInGraphAdjList(&n,&m);

    int i;

    int** adjList; // = Calloc(n, sizeof(int*));
    int* degree; // = Calloc(n, sizeof(int));

    #ifdef RETURN_CLIQUES_ONE_BY_ONE
    LinkedList* cliques = createLinkedList();
    #endif

    runAndPrintStatsListList( &listAllMaximalCliquesDegeneracy,
                              "degeneracy",
                              adjacencyList, adjList, 
                              #ifdef RETURN_CLIQUES_ONE_BY_ONE
                              cliques,
                              #endif
                              degree, n );

    #ifdef RETURN_CLIQUES_ONE_BY_ONE
    destroyCliqueResults(cliques);
    #endif

    // Free up memory from adjacency list.

    i = 0;
    while(i<n)
    {
        destroyLinkedList(adjacencyList[i]);
        i++;
    }

    Free(adjacencyList); 

    return 0;
}

/*! \mainpage Quick Cliques: A package to efficiently compute all maximal cliques in sparse graphs

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly

    \section intro_sec Introduction

    This package contains the code that was used to generate the results of 
    Eppstein and Strash (2011). (See \ref recommended_reading for references I am using.)
 
    This package contains code to generate all maximal cliques of graphs. It contains implementations
    of the following algorithms: 

    -# <b>(tomita)</b> The algorithm of Tomita et al. (2006), which is known to work well in practice,
    but consumes much memory because it uses an adjacency matrix representation of the input graph.
    For <b> small </b> sparse graphs, this algorithm is an excellent choice, and it is worst-case
    optimal. This algorithm also works very fast in practice for dense graphs.
    \see tomita.c tomita_algorithm.c <br><br>

    -# <b>(tomita-adjacency-list)</b> The algorithm of Tomita et al. (2006), modified to use adjacency 
    lists, therefore only requiring space linear in the size of the graph. This version is explained in the paper
    by Eppstein and Strash (2011), and is referred to in their paper as the "max-degree" algorithm. For sparse graphs,
    this algorithm can be quite fast on some inputs, but very slow on others. In particular, there are some sparse graphs 
    that cause this algorithm to run for over 5 hours, when the <b>hybrid</b> and <b>degeneracy</b> algorithms described
    next takes 10 and 4 minutes, respectively. This behavior is likely because there are some high-degree vertices whose 
    neighbor lists are iterated througha many times. Therefore, <i>use caution</i> when running this algorithm. 
    \see adjlist.c adjlist_algorithm.c <br><br>

    -# <b>(hybrid)</b> An implementation of the algorithm of Eppstein, Löffler, and Strash (2010), which 
    computes a degeneracy ordering of the input graph, and uses this ordering to speed up computation without any 
    fancy data structuring. This algorithm is called "hybrid" in Eppstein and Strash (2011). For sparse graphs, this 
    algorithm is consistently fast in practice, it is especially fast when the degeneracy of the graph is <b>very small</b>
    (for example, less than 6), or with graphs whose vertices, on average, have very few later neighbors in the degeneracy
    ordering. However, with larger degeneracy graphs, and graphs whose vertices have, on average, many more later neighbors in
    the degeneracy ordering, using <b>degeneracy</b>, the algorithm described next, is typically faster. Why? It uses a more 
    complex data structure that lets us gets us better speed in this case.
    \see hybrid.c hybrid_algorithm.c <br><br>

    -# <b>(degeneracy)</b> The algorithm of Eppstein, Löffler, and Strash (2010), which computes a degeneracy ordering of the
    input graph, and builds a data structure to speed up computation. This algorithm is called "degeneracy" in Eppstein
    and Strash (2011). For sparse graphs, this algorithm is consistently fast
    in practice, it is particular fast when the degeneracy of the graph is higher than 6, or with graphs whose vertices,
    on average, have many later neighbors in the degeneracy ordering.
    \see degeneracy.c degeneracy_algorithm.c <br><br>

    \section install_sec Compilation and Execution

        Under Linux systems, or *nix, executing "make" or "make all" from the top-level directory should compile all the source
        files in the <i>src</i> directory, place all object files in the <i>obj</i> directory, and place all
        executables in the <i>bin</i> directory. The makefile is not very sophisticated, so you may want to
        always run "make clean" first. If you have gcc and gnu make, this should just work. This code has not been tested
        in other environments.

        Executing the bash shell script "test.sh" in the top-level directory will execute every executable 
        on the data sets in <i>data</i> and print statistics about the data sets.

        \subsection Defines

        You can modify the defines in the Makefile to alter the behavior of the executables.

        - <b>DEBUG</b> - print verbose debug information; may be annoying
        - <b>MEMORY_DEBUG</b> - check for malloc or calloc returning NULL; if so, exit gracefully (causes slightly slower code)
        - <b>RETURN_CLIQUES_ONE_BY_ONE</b> - add cliques to a linked list that is returned by the clique listing algorithms
        - <b>PRINT_CLIQUES_ONE_BY_ONE</b> - print cliques to standard output one per line
            - (don't use with define PRINT_CLIQUES_TOMITA_STYLE)
        - <b>PRINT_CLIQUES_TOMITA_STYLE</b> - print cliques as described by Tomita et al. (2006)
            - prints to standard output: print the vertex when it is added to the partial clique, 
              print "b" when a vertex is removed, and print "c" when
              a maximal clique is found. (don't use with define PRINT_CLIQUES_ONE_BY_ONE)
        - <b>ALLOW_ALLOC_ZERO_BYTES</b> - some systems behave strangely when you allocate 0 bytes with either
          malloc or calloc, if this is not defined, then we always allocated at least one byte.

        The Defines used to generate the results of Eppstein and Strash (2011) are -DALLOW_ALLOC_ZERO_BYTES and -DPRINT_CLIQUES_TOMITA_STYLE
        with redirecting the standard output to /dev/null.

        \subsection Executables

        Each executable reads in the input graph from standard input, and writes the main product to 
        standard output, any statistical information or by-products are printed to standard error.

        - Statistics:
            - <i>bin/printnm</i> - print the number of vertices and edges in the data set
            - <i>bin/compdegen</i> - compute and print the degeneracy of the data set

        - Maximal clique listing: (print number of maximal cliques, and algorithm running time)
            - <i>bin/tomita</i> - execute algorithm <b>tomita</b> on the data set
            - <i>bin/adjlist</i> - execute algorithm <b>tomita-adjacency-list</b> on the data set
            - <i>bin/hybrid</i> - execute algorithm <b>hybrid</b> on the data set
            - <i>bin/degeneracy</i> - execute algorithm <b>degeneracy</b> on the data set

    \section data_sets Data Sets
    
    Several data sets are in the directory <i>data</i>, in a subdirectory named according to their origin.

        \subsection Sources

        I have included several data sets from BioGRID version 3.0.65 (http://thebiogrid.org/) in subdirectory <i>biogrid</i>.

See Eppstein and Strash (2011) for more sources

        \subsection Format

        Each file contains the following:
        -# The number of vertices in the graph
        -# Twice the number of edges in the graph.
        -# a list of edges in x,y format, where 0<=x,y<number of vertices.
            - the graph is symmetric, so x,y and y,x are in the list.

    \section recommended_reading Recommended Reading

    - "The Worst-Case Time Complexity for Generating All Maximal Cliques" by Esuji Tomita, Akira Tanaka, and Haruhisha Takahashi, <i>Theoretical Computer Science</i>, 363(1), 2006. (http://dx.doi.org/10.1016/j.tcs.2006.06.015)

    - "Listing All Maximal Cliques in Near-Optimal Time" by David Eppstein, Maarten Löffler, and Darren Strash, <i>ISAAC 2010</i>, LNCS volume 6506, pp. 403-416, 2010. (http://dx.doi.org/10.1007/978-3-642-17517-6_36 or http://arxiv.org/abs/1006.5440)

    - "Listing All Maximal Cliques in Large Sparse Real-World Graphs" by David Eppstein and Darren Strash, <i>SEA 2011</i>, LNCS volume 6630, pp. 364-375, 2011. (http://dx.doi.org/10.1007/978-3-642-20662-7_31 or http://arxiv.org/abs/1103.0318)
 
*/
