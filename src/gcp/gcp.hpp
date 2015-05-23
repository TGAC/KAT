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
#include <memory>
#include <thread>
using std::shared_ptr;
using std::make_shared;
using std::thread;

#include <jellyfish/mer_dna.hpp>
#include <jellyfish_helper.hpp>

#include <matrix/matrix_metadata_extractor.hpp>
#include <matrix/threaded_sparse_matrix.hpp>

#include "gcp_args.hpp"

using std::ostream;

namespace kat {

    class Gcp {
    private:

        // Input args
        GcpArgs args;
        const uint64_t nb_slices;
        uint64_t slice_id;

        // Variables that live for the lifetime of this object
        shared_ptr<JellyfishHelper> jfh;
        
        // Stores results
        shared_ptr<ThreadedSparseMatrix> gcp_mx; // Stores cumulative base count for each sequence where GC and CVG are binned

    public:

        Gcp(GcpArgs& _args) :
        args(_args), nb_slices(args.threads_arg * 100) {
            // Setup handle to jellyfish hash
            jfh = make_shared<JellyfishHelper>(args.db_arg, AccessMethod::SEQUENTIAL);
        }

        ~Gcp() {            
        }

        void execute() {
            
            // Setup output stream for jellyfish initialisation
            std::ostream* out_stream = args.verbose ? &cerr : (std::ostream*)0;

            // Create matrix of appropriate size (adds 1 to cvg bins to account for 0)
            gcp_mx = make_shared<ThreadedSparseMatrix>(jfh->getKeyLen(), args.cvg_bins + 1, args.threads_arg);

            // Process batch with worker threads
            // Process each sequence is processed in a different thread.
            // In each thread lookup each K-mer in the hash
            startAndJoinThreads(args.threads_arg);

            // Merge the contamination matrix
            gcp_mx->mergeThreadedMatricies();
        }
        
        void printVars(ostream &out) {
            out << "GCP parameters:" << endl;
            out << " - Hash File Path: " << args.db_arg << endl;
            out << " - Hash: " << (jfh != nullptr ? "loaded" : "not loaded") << endl;
            out << " - Threads: " << args.threads_arg << endl;
        }

        // Print K-mer comparison matrix

        void printMainMatrix(ostream &out) {
            SM64 mx = gcp_mx->getFinalMatrix();

            out << mme::KEY_TITLE << "K-mer coverage vs GC count plot for: " << args.db_arg << endl;
            out << mme::KEY_X_LABEL << "K-mer multiplicity" << endl;
            out << mme::KEY_Y_LABEL << "GC count" << endl;
            out << mme::KEY_Z_LABEL << "Distinct K-mers per bin" << endl;
            out << mme::KEY_NB_COLUMNS << mx->height() << endl;
            out << mme::KEY_NB_ROWS << mx->width() << endl;
            out << mme::KEY_MAX_VAL << mx->getMaxVal() << endl;
            out << mme::KEY_TRANSPOSE << "0" << endl;
            out << mme::MX_META_END << endl;

            mx->printMatrix(out);
        }
        
    protected:
        
         void startAndJoinThreads(const uint16_t nbThreads) {
            
            thread t[nbThreads];
            
            for(int i = 0; i < nbThreads; i++) {
                t[i] = thread(&Gcp::start, this, i);
            }
            
            for(int i = 0; i < nbThreads; i++){
                t[i].join();
            }
        }

        void start(int th_id) {
            SM64 mx = gcp_mx->getThreadMatrix(th_id);

            for (size_t i = slice_id++; i < nb_slices; i = slice_id++) {
                lha::region_iterator it = jfh->getSlice(i, nb_slices);
                while (it.next()) {
                    string kmer = ((jellyfish::mer_dna)it.key()).to_str();
                    uint64_t kmer_count = it.val();

                    uint16_t g_or_c = 0;

                    for (uint16_t i = 0; i < kmer.length(); i++) {
                        char c = kmer[i];

                        if (c == 'G' || c == 'g' || c == 'C' || c == 'c')
                            g_or_c++;
                    }

                    // Apply scaling factor
                    uint64_t cvg_pos = kmer_count == 0 ? 0 : ceil((double) kmer_count * args.cvg_scale);

                    if (cvg_pos > args.cvg_bins)
                        mx->inc(g_or_c, args.cvg_bins, 1);
                    else
                        mx->inc(g_or_c, cvg_pos, 1);
                }
            }
        }

        

    };
}
