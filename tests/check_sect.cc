#include <iostream>
#include <vector>

#include <kseq/kseq_helper.h>

#include "../src/sect/sect.hpp"


#define TEST_JF_HASH_1 "data/jf_hash"
#define TEST_FASTA_1 "data/test1.fa"


int testFastaLoad()
{
    vector<string> names;
    vector<string> seqs;
    readFasta(TEST_FASTA_1, names, seqs, false);

    if (names.size() != 6)
        return 1;

    if (seqs.size() != 6)
        return 1;

    return 0;
}

int main (void)
{
    int fastaLoad = testFastaLoad();

    if (fastaLoad != 0)
    {
        std::cerr << "Fasta Load failed" << std::endl;
        return 1;
    }

    return 0;
}

