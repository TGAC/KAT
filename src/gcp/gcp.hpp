#pragma once

#include <stdint.h>
#include <iostream>

#include <jellyfish/hash.hpp>
#include <jellyfish/counter.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/jellyfish_helper.hpp>

#include <matrix/matrix_metadata_extractor.hpp>
#include <matrix/threaded_sparse_matrix.hpp>

#include "gcp_args.hpp"

using std::ostream;

#define MAX_COUNT 1000

template<typename hash_t>
class Gcp : public thread_exec
{
private:

    // Input args
    const GcpArgs                   *args;
    const uint64_t                  nb_slices;
    counter_t                       slice_id;
    uint8_t                         kmer_length;


    // Variables that live for the lifetime of this object
    JellyfishHelper                 *jfh;
    hash_t                          *hash;		// Jellyfish hash

    // Stores results
    ThreadedSparseMatrix<uint64_t>  *gcp_mx;  // Stores cumulative base count for each sequence where GC and CVG are binned

public:
    Gcp(GcpArgs *_args) :
        args(_args), nb_slices(args->threads_arg * 100)
    {
        // Setup handle to jellyfish hash
        jfh = new JellyfishHelper(args->db_arg);

        // No output as yet
        gcp_mx = NULL;
    }

    ~Gcp()
    {
        if (jfh)
            delete jfh;

        jfh = NULL;


        if (gcp_mx)
            delete gcp_mx;

        gcp_mx = NULL;
    }


    void do_it()
    {
        // Setup output stream for jellyfish initialisation
        std::ostream* out_stream = args->verbose ? &cerr : (std::ostream*)0;

        // Load the jellyfish hash for sequential access
        hash = jfh->loadHash(true, out_stream);

        // Create matrix of appropriate size
        gcp_mx = new ThreadedSparseMatrix<uint64_t>(hash->get_mer_len(), MAX_COUNT + 1, args->threads_arg);

        // Process batch with worker threads
        // Process each sequence is processed in a different thread.
        // In each thread lookup each kmer in the hash
        exec_join(args->threads_arg);

        // Merge the contamination matrix
        gcp_mx->mergeThreadedMatricies();
    }

    void start(int th_id)
    {
        SparseMatrix<uint64_t>* mx = gcp_mx->getThreadMatrix(th_id);

        for(size_t i = slice_id++; i < nb_slices; i = slice_id++)
        {
            typename hash_t::iterator it = hash->iterator_slice(i, nb_slices);
            while(it.next())
            {
                string kmer = it.get_dna_str();
                uint64_t kmer_count = it.get_val();

                uint16_t g_or_c = 0;

                for(uint16_t i = 0; i < kmer.length(); i++)
                {
                    char c = kmer[i];

                    if (c == 'G' || c == 'g' || c == 'C' || c == 'c')
                        g_or_c++;
                }


                if(it.get_val() > MAX_COUNT)
                    mx->inc(g_or_c, MAX_COUNT, 1);
                else
                    mx->inc(g_or_c, kmer_count, 1);
            }
        }
    }


    void printVars(ostream &out)
    {
        out << "GCP parameters:" << endl;
        out << " - Hash File Path: " << args->db_arg << endl;
        out << " - Hash: " << (hash ? "loaded" : "not loaded") << endl;
        out << " - Threads: " << args->threads_arg << endl;
    }

    // Print kmer comparison matrix
    void printMainMatrix(ostream &out)
    {
        SparseMatrix<uint64_t>* mx = gcp_mx->getFinalMatrix();

        out << mme::KEY_TITLE << "Kmer coverage vs GC count plot for: " << args->db_arg << endl;
        out << mme::KEY_X_LABEL << "Kmer multiplicity" << endl;
        out << mme::KEY_Y_LABEL << "GC count" << endl;
        out << mme::KEY_NB_COLUMNS << mx->height() << endl;
        out << mme::KEY_NB_ROWS << mx->width() << endl;
        out << mme::KEY_MAX_VAL << mx->getMaxVal() << endl;
        out << mme::MX_META_END << endl;

        mx->printMatrix(out);
    }

};
