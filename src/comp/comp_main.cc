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
    Comp<hash_query_t> comp(args.db1_arg, args.db2_arg, args.db3_arg, args.threads_arg,
                            args.d1_scale_arg, args.d2_scale_arg,
                            args.d1_bins, args.d2_bins,
                            args.verbose);

    // Output comp parameters to stderr if requested
    if (args.verbose)
        comp.printVars(cerr);

    // Do the work
    comp.do_it();

    // Send main matrix to output file
    std::ostringstream main_mx_out_path;
    main_mx_out_path << args.output_prefix_arg << "_main.mx";
    ofstream_default main_mx_out_stream(main_mx_out_path.str().c_str(), cout);
    comp.printMainMatrix(main_mx_out_stream);
    main_mx_out_stream.close();

    // Output ends matricies if required
    if (args.db3_arg)
    {
        // Ends matrix
        std::ostringstream ends_mx_out_path;
        ends_mx_out_path << args.output_prefix_arg << "_ends.mx";
        ofstream_default ends_mx_out_stream(ends_mx_out_path.str().c_str(), cout);
        comp.printEndsMatrix(ends_mx_out_stream);
        ends_mx_out_stream.close();

        // Middle matrix
        std::ostringstream middle_mx_out_path;
        middle_mx_out_path << args.output_prefix_arg << "_middle.mx";
        ofstream_default middle_mx_out_stream(middle_mx_out_path.str().c_str(), cout);
        comp.printMiddleMatrix(middle_mx_out_stream);
        middle_mx_out_stream.close();

        // Mixed matrix
        std::ostringstream mixed_mx_out_path;
        mixed_mx_out_path << args.output_prefix_arg << "_mixed.mx";
        ofstream_default mixed_mx_out_stream(mixed_mx_out_path.str().c_str(), cout);
        comp.printMixedMatrix(mixed_mx_out_stream);
        mixed_mx_out_stream.close();
    }

    // Send kmer statistics to file
    std::ostringstream stats_out_path;
    stats_out_path << args.output_prefix_arg << ".stats";
    ofstream_default stats_out_stream(stats_out_path.str().c_str(), cout);
    comp.printCounters(stats_out_stream);
    stats_out_stream.close();

    // Send kmer statistics to stdout as well
    comp.printCounters(cout);


    return 0;
}
