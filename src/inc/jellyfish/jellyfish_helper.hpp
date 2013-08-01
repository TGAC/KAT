#pragma once

#include <string.h>
#include <iostream>

#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/compacted_hash.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/fstream_default.hpp>
#include <jellyfish/jellyfish_helper.hpp>

using std::string;
using std::cerr;
using std::endl;

class JellyfishHelper
{
private:
    const string jfHashPath;
    mapped_file* dbf;
    hash_query_t* hash;

public:

    JellyfishHelper(const string _jfHashPath) :
        jfHashPath(_jfHashPath)
    {
        dbf = new mapped_file(_jfHashPath.c_str());
    }

    ~JellyfishHelper()
    {
        if (dbf)
            delete dbf;

        dbf = NULL;

        if (hash)
            delete hash;

        hash = NULL;
    }


    hash_query_t* loadHash(const bool sequential, std::ostream* out)
    {

        // Advise kernel on how we will use this memory.
        if (sequential)
            dbf->sequential();
        else
            dbf->random();

        dbf->will_need();

        // Get jellyfish has type
        char type[8];
        memcpy(type, dbf->base(), sizeof(type));

        if(!strncmp(type, jellyfish::compacted_hash::file_type, sizeof(type)))
        {
            if (out)
                *out << endl << "Compacted hash detected.  Setting up query structure." << endl;

            // Load the jellyfish hashes
            hash = new hash_query_t(*dbf);

            // Output jellyfish has details if requested
            if (out)
            {
                *out << "Jellyfish hash: " << jfHashPath.c_str() << endl
                    << " - mer length  = " << hash->get_mer_len() << endl
                    << " - hash size   = " << hash->get_size() << endl
                    << " - max reprobe = " << hash->get_max_reprobe() << endl
                    << " - matrix      = " << hash->get_hash_matrix().xor_sum() << endl
                    << " - inv_matrix  = " << hash->get_hash_inverse_matrix().xor_sum() << endl;
            }

            return hash;
        }
        else
        {
            cerr << "Can't process jellyfish hash.  Wrong hash type.  Can only process compacted hash.\n";
            throw;
        }
    }

    hash_query_t* getHash()
    {
        return hash;
    }


};
