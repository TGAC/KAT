//  ********************************************************************
//  This file is part of KAT - the K-mer Analysis Toolkit.
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

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

namespace bfs = boost::filesystem;
using boost::shared_ptr;

#include "comp.hpp"
#include "comp_args.hpp"
#include "comp_main.hpp"

using std::vector;
using std::string;
using std::cout;
using std::cerr;

using kat::CompArgs;
using kat::Comp;

// Start point

int kat::compStart(int argc, char *argv[]) {
    
    // Parse args
    CompArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    // Check input file exists
    if (!bfs::exists(args.db1_path) && !bfs::symbolic_link_exists(args.db1_path)) {
        cerr << endl << "Could not find first jellyfish hash file at: " << args.db1_path << "; please check the path and try again." << endl << endl;
        return 1;
    }

    // Check input file exists
    if (!bfs::exists(args.db2_path) && !bfs::symbolic_link_exists(args.db2_path)) {
        cerr << endl << "Could not find second jellyfish hash file at: " << args.db2_path << "; please check the path and try again." << endl << endl;
        return 1;
    }

    // Check input file exists
    if (!args.db3_path.empty() && !bfs::exists(args.db3_path) && !bfs::symbolic_link_exists(args.db3_path)) {
        cerr << endl << "Could not find third jellyfish hash file at: " << args.db3_path << "; please check the path and try again." << endl << endl;
        return 1;
    }

    // Create the sequence coverage object
    Comp comp(args);

    // Output comp parameters to stderr if requested
    if (args.verbose)
        comp.printVars(cerr);

    // Do the work
    comp.execute();

    // Send main matrix to output file
    std::ostringstream main_mx_out_path;
    main_mx_out_path << args.output_prefix << "_main.mx";
    ofstream_default main_mx_out_stream(main_mx_out_path.str().c_str(), cout);
    comp.printMainMatrix(main_mx_out_stream);
    main_mx_out_stream.close();

    // Output ends matricies if required
    if (!(args.db3_path.empty())) {
        // Ends matrix
        std::ostringstream ends_mx_out_path;
        ends_mx_out_path << args.output_prefix << "_ends.mx";
        ofstream_default ends_mx_out_stream(ends_mx_out_path.str().c_str(), cout);
        comp.printEndsMatrix(ends_mx_out_stream);
        ends_mx_out_stream.close();

        // Middle matrix
        std::ostringstream middle_mx_out_path;
        middle_mx_out_path << args.output_prefix << "_middle.mx";
        ofstream_default middle_mx_out_stream(middle_mx_out_path.str().c_str(), cout);
        comp.printMiddleMatrix(middle_mx_out_stream);
        middle_mx_out_stream.close();

        // Mixed matrix
        std::ostringstream mixed_mx_out_path;
        mixed_mx_out_path << args.output_prefix << "_mixed.mx";
        ofstream_default mixed_mx_out_stream(mixed_mx_out_path.str().c_str(), cout);
        comp.printMixedMatrix(mixed_mx_out_stream);
        mixed_mx_out_stream.close();
    }

    // Send K-mer statistics to file
    std::ostringstream stats_out_path;
    stats_out_path << args.output_prefix << ".stats";
    ofstream_default stats_out_stream(stats_out_path.str().c_str(), cout);
    comp.printCounters(stats_out_stream);
    stats_out_stream.close();

    // Send K-mer statistics to stdout as well
    comp.printCounters(cout);


    return 0;
}
