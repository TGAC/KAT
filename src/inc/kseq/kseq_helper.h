#ifndef __KSEQ_HELPER_HPP__
#define __KSEQ_HELPER_HPP__

#include <string.h>
#include <iostream>
#include <vector>

#include <zlib/zlib.h>

#include <kseq/kseq.h>

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
        cerr << "Fasta file load: " << fastaPath << "\n";

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

#endif
