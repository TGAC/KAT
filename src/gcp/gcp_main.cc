#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdio.h>

#include <jellyfish/mer_counting.hpp>

#include "gcp.hpp"
#include "gcp_args.hpp"
#include "gcp_main.hpp"

//using std::cerr;

// Start point
int gcpStart(int argc, char *argv[])
{
    // Parse args
    GcpArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    // Create the sequence coverage object
    Gcp<hash_query_t> gcp(&args);

    // Output seqcvg parameters to stderr if requested
    if (args.verbose)
        gcp.printVars(cerr);

    // Do the work (outputs data to files as it goes)
    gcp.do_it();

    // Send main matrix to output file
    std::ostringstream main_mx_out_path;
    main_mx_out_path << args.output_prefix << ".mx";
    ofstream_default main_mx_out_stream(main_mx_out_path.str().c_str(), cout);
    gcp.printMainMatrix(main_mx_out_stream);
    main_mx_out_stream.close();

    return 0;
}
