//  ********************************************************************
//  This file is part of KAT - the Kmer Analysis Toolkit.
//
//  KAT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  KAT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with KAT.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

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

using kat::GcpArgs;
using kat::Gcp;

// Start point
int kat::gcpStart(int argc, char *argv[])
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
