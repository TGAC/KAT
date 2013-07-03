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

#include <jellyfish/mer_counting.hpp>

#include <kseq/kseq_helper.h>

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

    // Load enitre fasta file into memory
    vector<string> names;
    vector<string> seqs;
    readFasta(args.fasta_arg, names, seqs, args.verbose);

    // Create the sequence coverage object
    Sect<hash_query_t> sect(args.db_arg, &names, &seqs, args.gc_bins, args.cvg_bins, args.cvg_logscale, args.threads_arg, args.verbose);

    // Output seqcvg parameters to stderr if requested
    if (args.verbose)
        sect.printVars(cerr);

    // Do the work
    sect.do_it();

    // Send sequence kmer counts to file
    std::ostringstream count_path;
    count_path << args.output_prefix << "_counts.cvg";
    ofstream_default count_path_stream(count_path.str().c_str(), cout);
    sect.printCounts(count_path_stream);
    count_path_stream.close();


    // Send average sequence coverage and GC% scores to file
    std::ostringstream cvg_gc_path;
    cvg_gc_path << args.output_prefix << "_cvg-gc.csv";
    ofstream_default cvg_gc_stream(cvg_gc_path.str().c_str(), cout);
    sect.printStatTable(cvg_gc_stream);
    cvg_gc_stream.close();

    // Send contamination matrix to file
    std::ostringstream contamination_mx_path;
    contamination_mx_path << args.output_prefix << "_contamination.mx";
    ofstream_default contamination_mx_stream(contamination_mx_path.str().c_str(), cout);
    sect.printContaminationMatrix(contamination_mx_stream);
    contamination_mx_stream.close();

    return 0;
}
