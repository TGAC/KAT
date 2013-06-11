#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>


#include <zlib/zlib.h>

#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/compacted_hash.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/fstream_default.hpp>

#include <kseq/kseq.h>

#include <gnuplot/gnuplot_i.hpp>

#include <matrix/sparse_matrix.hpp>

#include "comp_args.hpp"
#include "comp_main.hpp"

using std::vector;
using std::string;
using std::cout;
using std::cerr;



// Start point
int compStart(int argc, char *argv[])
{
    // Parse args
    CompArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
    {
        args.print();
    }

    // Create handles to the memory mapped hash files
    mapped_file dbf1(args.db1_arg);
    mapped_file dbf2(args.db2_arg);

    // Advise kernel on how we will use this memory (i.e. random access and we'll
    // need it available soon.)
    dbf1.random().will_need();
    dbf2.random().will_need();

    // Get jellyfish has type
    char type1[8];
    char type2[8];
    memcpy(type1, dbf1.base(), sizeof(type1));
    memcpy(type2, dbf2.base(), sizeof(type2));

    // Create sparse matrix
    SparseMatrix<long> matrix(1000, 1000);

    matrix(100,200) = 5L;

    cout << matrix(100,200) << "\n";

    return 0;
}
