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

#include <iostream>

#include <jellyfish/mer_counting.hpp>

#include <file_utils.hpp>

#include "sect.hpp"
#include "sect_args.hpp"
#include "sect_main.hpp"

using std::cerr;

using kat::SectArgs;
using kat::Sect;

// Start point
int kat::sectStart(int argc, char *argv[])
{
    // Parse args
    SectArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    // Check input files exist
    if (!fileExists(args.jellyfish_hash))
    {
        cerr << endl << "Could not find jellyfish hash file at: " << args.jellyfish_hash << "; please check the path and try again." << endl << endl;
        return 2;
    }

    if (!fileExists(args.seq_file))
    {
        cerr << endl << "Could not find sequence file at: " << args.seq_file << "; please check the path and try again." << endl << endl;
        return 3;
    }

    // Create the sequence coverage object
    Sect<hash_query_t> sect(&args);

    // Output seqcvg parameters to stderr if requested
    if (args.verbose)
        sect.printVars(cerr);

    // Do the work (outputs data to files as it goes)
    sect.do_it();

    // Return the Sect result code... hopefully should be 0
    return sect.getResultCode();
}
