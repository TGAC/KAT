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
#include "sect_args.hpp"
#include "sect.hpp"
#include "sect_main.hpp"

using std::vector;
using std::string;
using std::cout;
using std::cerr;

KSEQ_INIT(gzFile, gzread)


// Loads Fasta file into memory.  Two vectors hold names and sequences respectively.  Assumes no reordering will
// take place
void readFasta(const char *fastaPath, vector<string>& fastaNames, vector<string>& fastaSeqs, bool verbose)
{
    if (verbose)
    {
        cerr << "Fasta file load: " << fastaPath << "\n";
    }

    gzFile fp;
    fp = gzopen(fastaPath, "r"); // STEP 2: open the file handler
    kseq_t *seq = kseq_init(fp); // STEP 3: initialize seq
    int l;
    while ((l = kseq_read(seq)) >= 0)   // STEP 4: read sequence
    {
        fastaNames.push_back(seq->name.s);
        fastaSeqs.push_back(seq->seq.s);
    }
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp); // STEP 6: close the file handler
}

// Debugging routine.. not required for normal use
void printFastaData(vector<string>& names, vector<string>& seqs)
{
    if (names.size() != seqs.size())
    {
        cerr << "Somehow the fasta names and sequence vector went out of sync!  Names size: " << names.size() << "; Seqs size: " << seqs.size() << "\n";
        return;
    }

    cerr << "Printing fasta data\n";

    for (int i = 0; i < names.size(); i++)
    {
        cerr << "name: " << names[i].c_str() << "\n";
        cerr << "seq:  " << seqs[i].c_str() << "\n";
    }
}

// Start point
int sectStart(int argc, char *argv[])
{
    // Parse args
    SectArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
    {
        args.print();
    }

    // Create handle to the memory mapped hash file into memory
    mapped_file dbf(args.db_arg);

    // Advise kernel on how we will use this memory (i.e. random access and we'll
    // need it available soon.)
    dbf.random().will_need();

    // Load enitre fasta file into memory
    vector<string> names;
    vector<string> seqs;
    readFasta(args.fasta_arg, names, seqs, args.verbose);

    // Get jellyfish has type
    char type[8];
    memcpy(type, dbf.base(), sizeof(type));

    // Process data differently depending on jellyfish hash type
    if(!strncmp(type, jellyfish::raw_hash::file_type, sizeof(type)))
    {
        cerr << "Raw hash detected but not supported yet\n";
    }
    else if(!strncmp(type, jellyfish::compacted_hash::file_type, sizeof(type)))
    {
        if (args.verbose)
        {
            cerr << "Compacted hash detected.  Setting up query structure.\n";
        }

        // Load the jellyfish hash object
        hash_query_t hash(dbf);

        // Output jellyfish has details if requested
        if (args.verbose)
        {
            cerr << "mer length  = " << hash.get_mer_len() << "\n"
                      << "hash size   = " << hash.get_size() << "\n"
                      << "max reprobe = " << hash.get_max_reprobe() << "\n"
                      << "matrix      = " << hash.get_hash_matrix().xor_sum() << "\n"
                      << "inv_matrix  = " << hash.get_hash_inverse_matrix().xor_sum() << "\n";
        }

        // Create the sequence coverage object
        Sect<hash_query_t> sect(&hash, &names, &seqs, args.kmer_arg, args.threads_arg);

        // Output seqcvg parameters to stderr if requested
        if (args.verbose)
        {
            sect.printVars();
        }

        // Do the work
        sect.do_it();

        // Send sequence kmer counts to file if requested
        if (args.outputGiven())
        {
            ofstream_default count_out(args.outputGiven() ? args.output_arg : 0, std::cout);
            sect.printCounts(count_out);
            count_out.close();
        }

        // Send sequence coverage scores to standard out
        sect.printCoverages(std::cout);
    }
    else
    {
        die << "Invalid jellyfish file type '" << err::substr(type, sizeof(type)) << "'.";
    }

    return 0;
}
