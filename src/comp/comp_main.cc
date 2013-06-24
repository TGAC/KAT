#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/compacted_hash.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/fstream_default.hpp>
#include <jellyfish/jellyfish_helper.hpp>

#include "comp.hpp"
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
        args.print();

    // Create the sequence coverage object
    Comp<hash_query_t> comp(args.db1_arg, args.db2_arg, 0, 0, args.threads_arg, args.xscale_arg, args.yscale_arg);

    // Output comp parameters to stderr if requested
    if (args.verbose)
        comp.printVars(cerr);

    // Do the work
    comp.do_it();

    // Send matrix to output file
    ofstream_default matrix_out(args.output_arg, cout);
    comp.printMatrix(matrix_out);
    matrix_out.close();

    // Send kmer statistics to stdout
    comp.printCounters(cout);


    return 0;
}
