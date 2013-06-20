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

#include <kseq/kseq.h>

#include <gnuplot/gnuplot_i.hpp>

#include "plot_args.hpp"
#include "plot_main.hpp"

using std::vector;
using std::string;
using std::cout;
using std::cerr;



void flamePlot()
{

}

void asmPlot()
{

}




// Start point
int plotStart(int argc, char *argv[])
{
    // Parse args
    PlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    // Pass remaining args to relevant child tool
    if (args.mode_arg.compare("flame") == 0)
    {
        flamePlot();
    }
    else if (args.mode_arg.compare("asm") == 0)
    {
        asmPlot();
    }

    /*// Create the sequence coverage object
    Comp<hash_query_t> comp(&hash1, &hash2, 0, 0, args.threads_arg, args.xscale_arg, args.yscale_arg);

    // Output comp parameters to stderr if requested
    if (args.verbose)
        comp.printVars(cerr);

    // Do the work
    comp.do_it();

    // Send matrix to output file
    ofstream_default matrix_out(args.output_arg, std::cout);
    comp.printMatrix(matrix_out);
    matrix_out.close();

    // Send kmer statistics to stdout
    comp.printCounters(cout);
   */

    return 0;
}
