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

#include "sect.hpp"
#include "sect_args.hpp"
#include "sect_main.hpp"

using std::vector;
using std::string;
using std::cout;
using std::cerr;


// Start point
int sectStart(int argc, char *argv[])
{
    // Parse args
    SectArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    // Create the sequence coverage object
    Sect<hash_query_t> sect(&args);

    // Output seqcvg parameters to stderr if requested
    if (args.verbose)
        sect.printVars(cerr);

    // Do the work (outputs data to files as it goes)
    sect.do_it();

    return 0;
}
