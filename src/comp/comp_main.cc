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
    {
        args.print();
    }

    // Create handles to the memory mapped hash files
    mapped_file dbf1(args.db1_arg);
    mapped_file dbf2(args.db2_arg);

    // Advise kernel on how we will use this memory.  For hash1 treat it as sequential as we will iterate through it
    // hash2 should be random access as we will look up kmers from hash1.  It therefore makes sense for hash1 to be
    // the bigger hash, for example if comparing reads to an assembly.
    dbf1.sequential().will_need();
    dbf2.random().will_need();

    // Get jellyfish has type
    char type1[8];
    char type2[8];
    memcpy(type1, dbf1.base(), sizeof(type1));
    memcpy(type2, dbf2.base(), sizeof(type2));


    if(!strncmp(type1, jellyfish::compacted_hash::file_type, sizeof(type1)) &&
       !strncmp(type2, jellyfish::compacted_hash::file_type, sizeof(type2)))
    {
            if (args.verbose)
                cerr << "Compacted hashes detected.  Setting up query structure.\n\n";

            // Load the jellyfish hashes
            hash_query_t hash1(dbf1);
            hash_query_t hash2(dbf2);

            // Output jellyfish has details if requested
            if (args.verbose)
            {
                cerr << "hash1 mer length  = " << hash1.get_mer_len() << "\n"
                          << "hash1 hash size   = " << hash1.get_size() << "\n"
                          << "hash1 max reprobe = " << hash1.get_max_reprobe() << "\n"
                          << "hash1 matrix      = " << hash1.get_hash_matrix().xor_sum() << "\n"
                          << "hash1 inv_matrix  = " << hash1.get_hash_inverse_matrix().xor_sum() << "\n\n";

                cerr << "hash2 mer length  = " << hash2.get_mer_len() << "\n"
                          << "hash2 hash size   = " << hash2.get_size() << "\n"
                          << "hash2 max reprobe = " << hash2.get_max_reprobe() << "\n"
                          << "hash2 matrix      = " << hash2.get_hash_matrix().xor_sum() << "\n"
                          << "hash2 inv_matrix  = " << hash2.get_hash_inverse_matrix().xor_sum() << "\n\n";
            }

            // Check kmer lengths are the same for both hashes.  We can't continue if they are not.
            if (hash1.get_mer_len() != hash2.get_mer_len())
            {
                cerr << "Cannot process hashes that were created with different kmer lengths.\n";
                return 1;
            }

            // Create the sequence coverage object
            Comp<hash_query_t> comp(&hash1, &hash2, 0, 0, args.threads_arg, args.xscale_arg, args.yscale_arg, args.noindex);

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
    }
    else
    {
        cerr << "Can't process jellyfish hashes.  Wrong type.  Can only process compacted hashes.\n";
    }

    return 0;
}
