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

#pragma once

#include <stdint.h>
#include <iostream>
#include <math.h>

#include <jellyfish/hash.hpp>
#include <jellyfish/counter.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/jellyfish_helper.hpp>

#include <matrix/matrix_metadata_extractor.hpp>
#include <matrix/threaded_sparse_matrix.hpp>

#include "gcp_args.hpp"

using std::ostream;

namespace kat
{
    template<typename hash_t>
    class Gcp : public thread_exec
    {
    private:

        // Input args
        const GcpArgs                   *args;
        const uint64_t                  nb_slices;
        counter_t                       slice_id;


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

            // Whether to treat this hash as double stranded or not.
            // Ideally it would be nice to determine this directly from the hash but I'm
            // not sure how to do that at the moment... it might not be possible
            if (args->both_strands)
                hash->set_canonical(true);

            // Create matrix of appropriate size (adds 1 to cvg bins to account for 0)
            gcp_mx = new ThreadedSparseMatrix<uint64_t>(hash->get_mer_len(), args->cvg_bins + 1, args->threads_arg);

            // Process batch with worker threads
            // Process each sequence is processed in a different thread.
            // In each thread lookup each K-mer in the hash
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

                    // Apply scaling factor
                    uint64_t cvg_pos = kmer_count == 0 ? 0 : ceil((double)kmer_count * args->cvg_scale);

                    if(cvg_pos > args->cvg_bins)
                        mx->inc(g_or_c, args->cvg_bins, 1);
                    else
                        mx->inc(g_or_c, cvg_pos, 1);
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

        // Print K-mer comparison matrix
        void printMainMatrix(ostream &out)
        {
            SparseMatrix<uint64_t>* mx = gcp_mx->getFinalMatrix();

            out << mme::KEY_TITLE << "K-mer coverage vs GC count plot for: " << args->db_arg << endl;
            out << mme::KEY_X_LABEL << "K-mer multiplicity" << endl;
            out << mme::KEY_Y_LABEL << "GC count" << endl;
            out << mme::KEY_Z_LABEL << "Distinct K-mers per bin" << endl;
            out << mme::KEY_NB_COLUMNS << mx->height() << endl;
            out << mme::KEY_NB_ROWS << mx->width() << endl;
            out << mme::KEY_MAX_VAL << mx->getMaxVal() << endl;
            out << mme::MX_META_END << endl;

            mx->printMatrix(out);
        }

    };
}
