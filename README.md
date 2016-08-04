# **Quick Cliques**: Quickly compute all maximal cliques in sparse graphs

[![license](https://img.shields.io/badge/license-GPL%20v3.0-blue.svg)](http://www.gnu.org/licenses/)

The original intent of this software was to provide exact reproducibility of experimental results from two papers:

*Listing All Maximal Cliques in Large Sparse Real-World Graphs in Near-Optimal Time*,  
**D. Eppstein, M. Löffler, and D. Strash**,  
Journal of Experimental Algorithmics, 18 (3): 3.1, 2013,   
[doi:10.1145/2543629](https://doi.org/10.1145/2543629)  

and

*Listing All Maximal Cliques in Large Sparse Real-World Graphs*,   
**D. Eppstein and D. Strash**  
Proceedings of the 10th International Conference on Experimental Algorithms (SEA 2011), LNCS vol. 6630, pp. 403-414.  
[doi:10.1007/978-3-642-20662-7_31](https://doi.org/10.1007/978-3-642-20662-7_31)  
[arXiv:1103.0318](https://arxiv.org/abs/1103.0318)  

However, in the original implementation (written in C) was not designed to be integrated with other code. It has since been upgraded to C++11, is becoming more hardened, and has better support for integration into other software packages.

Want the exact code from the original experiments? See release v1.0.

###This package includes:

 - C++ code for four implementations for enumerating all maximal cliques of a graph
 - A small set of data used for testing (see ./data) [ For a larger set of data (too large for GitHub), you can download the following: http://www.ics.uci.edu/~dstrash/data.tar.gz ]
 - The test script used to build and run all algorithms on all data sets (./test.sh)

Please feel free to contact me with any questions!

### Version
2.0beta

### Building

```sh
$ git clone https://github.com/darrenstrash/quick-cliques.git
$ cd quick-cliques
$ make
```

### Running
```sh
$ ./bin/qc --input-file=<input graph> --algorithm=<tomita|adjlist|degeneracy|hybrid>
```

or

```sh
$ ./test.sh
```

to run all algorithms on all data sets in ./data directory.

### Printing Cliques

Cliques can pe printed in two formats:

 - One clique per line, as space separated integers.

   This is activated by the `PRINT_CLIQUES_ONE_BY_ONE` define in `makefile`
 - The recursion tree format used by [Tomita et al. (2006)](http://dx.doi.org/10.1016/j.tcs.2006.06.015). That is:
   - Print a number when search evaluates a new vertex
   - Print a `c` when a new maximal clique is found
   - Print a `b` when search backtracks.

   This is activated by the `PRINT_CLIQUES_TOMITA_STYLE` define in `makefile`

### Graph Format

First, note that loops and directed edges are not supported, and could lead to errors.

Currently, two formats are supported:

 - **The unweighted METIS format**: Which consists of

   `<# vertices> <# edges> 1`

   followed by `<# vertices>` lines of space-separated vertices,  where the `i`-th line consists of 
   all neighbors of `i`. All vertices range from `1` to `<# vertices>`

 - **A custom format**: Which consists of

   `<# vertices> 2*<# edges>`

   followed by `2*<# edges>` lines of comma-separated vertices of the form 
 
   `u,v`
 
   where `u` and `v` are vertices represented by integers from `0` to `<# vertices>-1`. If `u,v` is in the list, then `v,u` must be in the list also.

There are plans to phase out the custom format and put in support for *snap* and *dimacs* files.

Copyright
----

Copyright (c) 2011-2016 Darren Strash.


License
----

This code is released under the GNU Public License (GPL) 3.0.

To read the GPL 3.0, read the file COPYING in this directory.

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

Contact
----

**Darren Strash** (first name DOT last name AT gmail DOT com)
