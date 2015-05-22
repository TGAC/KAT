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

#include <iostream>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <math.h>
#include <memory>
using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::ostream;
using std::shared_ptr;
using std::make_shared;

#include <matrix/matrix_metadata_extractor.hpp>
#include <matrix/sparse_matrix.hpp>
#include <matrix/threaded_sparse_matrix.hpp>

#include <jellyfish_helper.hpp>

#include "comp_args.hpp"


namespace kat {

    class CompCounters {
    public:
        uint64_t hash1_total;
        uint64_t hash2_total;
        uint64_t hash3_total;
        uint64_t hash1_distinct;
        uint64_t hash2_distinct;
        uint64_t hash3_distinct;
        uint64_t hash1_only_total;
        uint64_t hash2_only_total;
        uint64_t hash1_only_distinct;
        uint64_t hash2_only_distinct;
        uint64_t shared_hash1_total;
        uint64_t shared_hash2_total;
        uint64_t shared_distinct;

        const char* hash1_path;
        const char* hash2_path;
        const char* hash3_path;

        CompCounters(const char* _hash1_path, const char* _hash2_path, const char* _hash3_path) :
        hash1_path(_hash1_path), hash2_path(_hash2_path), hash3_path(_hash3_path) {

            hash1_total = 0;
            hash2_total = 0;
            hash3_total = 0;
            hash1_distinct = 0;
            hash2_distinct = 0;
            hash3_distinct = 0;
            hash1_only_total = 0;
            hash2_only_total = 0;
            hash1_only_distinct = 0;
            hash2_only_distinct = 0;
            shared_hash1_total = 0;
            shared_hash2_total = 0;
            shared_distinct = 0;
        }

        void updateHash1Counters(uint64_t hash1_count, uint64_t hash2_count) {
            hash1_total += hash1_count;
            hash1_distinct++;

            if (!hash2_count) {
                hash1_only_total += hash1_count;
                hash1_only_distinct++;
            }
        }

        void updateHash2Counters(uint64_t hash1_count, uint64_t hash2_count) {
            hash2_total += hash2_count;
            hash2_distinct++;

            if (!hash1_count) {
                hash2_only_total += hash2_count;
                hash2_only_distinct++;
            }
        }

        void updateHash3Counters(uint64_t hash3_count) {

            hash3_total += hash3_count;
            hash3_distinct++;
        }

        void updateSharedCounters(uint64_t hash1_count, uint64_t hash2_count) {

            if (hash1_count && hash2_count) {
                shared_hash1_total += hash1_count;
                shared_hash2_total += hash2_count;
                shared_distinct++;
            }
        }

        void printCounts(ostream &out) {

            out << "K-mer statistics for: " << endl;
            out << " - Hash 1: " << hash1_path << endl;
            out << " - Hash 2: " << hash2_path << endl;

            if (hash3_total > 0)
                out << " - Hash 3: " << hash3_path << endl;

            out << endl;

            out << "Total K-mers in: " << endl;
            out << " - Hash 1: " << hash1_total << endl;
            out << " - Hash 2: " << hash2_total << endl;

            if (hash3_total > 0)
                out << " - Hash 3: " << hash3_total << endl;

            out << endl;

            out << "Distinct K-mers in:" << endl;
            out << " - Hash 1: " << hash1_distinct << endl;
            out << " - Hash 2: " << hash2_distinct << endl;
            if (hash3_total > 0)
                out << " - Hash 3: " << hash3_distinct << endl;

            out << endl;

            out << "Total K-mers only found in:" << endl;
            out << " - Hash 1: " << hash1_only_total << endl;
            out << " - Hash 2: " << hash2_only_total << endl;
            out << endl;

            out << "Distinct K-mers only found in:" << endl;
            out << " - Hash 1: " << hash1_only_distinct << endl;
            out << " - Hash 2: " << hash2_only_distinct << endl << endl;

            out << "Shared K-mers:" << endl;
            out << " - Total shared found in hash 1: " << shared_hash1_total << endl;
            out << " - Total shared found in hash 2: " << shared_hash2_total << endl;
            out << " - Distinct shared K-mers: " << shared_distinct << endl << endl;
        }
    };

    class Comp {
    private:
        // Args passed in
        CompArgs args;

        // Jellyfish mapped file hash vars
        shared_ptr<JellyfishHelper> jfh1;
        shared_ptr<JellyfishHelper> jfh2;
        shared_ptr<JellyfishHelper> jfh3;

        // Threaded matrix data
        shared_ptr<ThreadedSparseMatrix> main_matrix;
        shared_ptr<ThreadedSparseMatrix> ends_matrix;
        shared_ptr<ThreadedSparseMatrix> middle_matrix;
        shared_ptr<ThreadedSparseMatrix> mixed_matrix;

        // Final data (created by merging thread results)
        shared_ptr<CompCounters> final_comp_counters;

        // Thread specific data
        shared_ptr<vector<shared_ptr<CompCounters> > > thread_comp_counters;


    public:

        Comp(CompArgs& _args) : args(_args) {
            if (args.verbose)
                cerr << "Setting up comp tool..." << endl;

            // Setup handles to load hashes
            jfh1 = make_shared<JellyfishHelper>(args.db1_path, AccessMethod::RANDOM);
            jfh2 = make_shared<JellyfishHelper>(args.db2_path, AccessMethod::RANDOM);
            jfh3 = !(args.db3_path.empty()) ?
                    make_shared<JellyfishHelper>(args.db3_path, AccessMethod::RANDOM) :
                    NULL;

            // Create the final K-mer counter matrices
            main_matrix = shared_ptr<ThreadedSparseMatrix>(
                    new ThreadedSparseMatrix(args.d1_bins, args.d2_bins, args.threads));

            // Initialise extra matrices for hash3 (only allocates space if required)
            if (!(args.db3_path.empty())) {
                if (args.verbose)
                    cerr << " - Setting up matrices for hash 3" << endl;

                ends_matrix = make_shared<ThreadedSparseMatrix>(args.d1_bins, args.d2_bins, args.threads);
                middle_matrix = make_shared<ThreadedSparseMatrix>(args.d1_bins, args.d2_bins, args.threads);
                mixed_matrix = make_shared<ThreadedSparseMatrix>(args.d1_bins, args.d2_bins, args.threads);
            }

            // Create the final comp counters
            final_comp_counters = make_shared<CompCounters>(
                    args.db1_path.c_str(), args.db2_path.c_str(), args.db3_path.c_str());

            // Create the comp counters for each thread
            thread_comp_counters = make_shared<vector<shared_ptr<CompCounters>>>(args.threads);
            for (int i = 0; i < args.threads; i++) {
                thread_comp_counters->push_back(make_shared<CompCounters>(
                        args.db1_path.c_str(), args.db2_path.c_str(), args.db3_path.c_str()));
            }

            if (args.verbose)
                cerr << "Comp tool setup successfully." << endl;
        }

        ~Comp() {
        }

        void execute() {

            if (args.verbose) {
                cerr << "Loading hashes..." << endl;
            }

            std::ostream* out_stream = args.verbose ? &cerr : (std::ostream*)0;

            // Load the hashes
            /*hash1 = jfh1.loadHash(true, out_stream);
            hash2 = jfh2.loadHash(false, out_stream);

            if (jfh3)
                hash3 = jfh3.loadHash(false, out_stream);*/

            // Check K-mer lengths are the same for both hashes.  We can't continue if they are not.
            if (jfh1->getKeyLen() != jfh2->getKeyLen()) {
                cerr << "Cannot process hashes that were created with different K-mer lengths.  "
                        << "Hash1: " << jfh1->getKeyLen() << ".  Hash2: " << jfh2->getKeyLen() << "." << endl;
                throw;
            } else if (jfh3 && (jfh3->getKeyLen() != jfh1->getKeyLen())) {
                cerr << "Cannot process hashes that were created with different K-mer lengths.  "
                        << "Hash1: " << jfh1->getKeyLen() << ".  Hash3: " << jfh3->getKeyLen() << "." << endl;
                throw;
            }

            if (args.verbose)
                cerr << endl
                    << "All hashes loaded successfully." << endl
                    << "Starting threads...";

            // Run the threads
            exec_join(args.threads);

            if (args.verbose)
                cerr << "done." << endl
                    << "Merging results...";

            // Merge results from the threads
            merge();

            if (args.verbose)
                cerr << "done." << endl;
        }

        void start(int th_id) {

            // Get handle to sparse matrix for this thread
            shared_ptr<SparseMatrix<uint64_t> > thread_matrix = main_matrix->getThreadMatrix(th_id);

            // Get handle on this thread's comp counter
            shared_ptr<CompCounters> comp_counters = (*thread_comp_counters)[th_id];

            // Setup iterator for this thread's chunk of hash1
            typename hash_t::iterator hash1Iterator = jfh1->getReader()->->pos() hash1->iterator_slice(th_id, args.threads);

            // Go through this thread's slice for hash1
            while (hash1Iterator.next()) {
                // Get the current K-mer count for hash1
                uint64_t hash1_count = hash1Iterator.get_val();

                // Get the count for this K-mer in hash2 (assuming it exists... 0 if not)
                uint64_t hash2_count = (*hash2)[hash1Iterator.get_key()];

                // Get the count for this K-mer in hash3 (assuming it exists... 0 if not)
                uint64_t hash3_count = hash3 ? (*hash3)[hash1Iterator.get_key()] : 0;

                // Increment hash1's unique counters
                comp_counters->updateHash1Counters(hash1_count, hash2_count);

                // Increment shared counters
                comp_counters->updateSharedCounters(hash1_count, hash2_count);

                // Scale counters to make the matrix look pretty
                uint64_t scaled_hash1_count = scaleCounter(hash1_count, args.d1_scale);
                uint64_t scaled_hash2_count = scaleCounter(hash2_count, args.d2_scale);
                uint64_t scaled_hash3_count = scaleCounter(hash3_count, args.d2_scale);

                // Modifies hash counts so that K-mer counts larger than MATRIX_SIZE are dumped in the last slot
                if (scaled_hash1_count >= args.d1_bins) scaled_hash1_count = args.d1_bins - 1;
                if (scaled_hash2_count >= args.d2_bins) scaled_hash2_count = args.d2_bins - 1;
                if (scaled_hash3_count >= args.d2_bins) scaled_hash3_count = args.d2_bins - 1;

                // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
                thread_matrix->inc(scaled_hash1_count, scaled_hash2_count, 1);

                // Update hash 3 related matricies if hash 3 was provided
                if (hash3) {
                    if (scaled_hash2_count == scaled_hash3_count)
                        ends_matrix->getThreadMatrix(th_id)->inc(scaled_hash1_count, scaled_hash3_count, 1);
                    else if (scaled_hash3_count > 0)
                        mixed_matrix->getThreadMatrix(th_id)->inc(scaled_hash1_count, scaled_hash3_count, 1);
                    else
                        middle_matrix->getThreadMatrix(th_id)->inc(scaled_hash1_count, scaled_hash3_count, 1);
                }
            }

            // Setup iterator for this thread's chunk of hash2
            // We setup hash2 for random access, so hopefully performance isn't too bad here...
            // Hash2 should be smaller than hash1 in most cases so hopefully we can get away with this.
            typename hash_t::iterator hash2Iterator = hash2->iterator_slice(th_id, args.threads);

            // Iterate through this thread's slice of hash2
            while (hash2Iterator.next()) {
                // Get the current K-mer count for hash2
                uint64_t hash2_count = hash2Iterator.get_val();

                // Get the count for this K-mer in hash1 (assuming it exists... 0 if not)
                uint64_t hash1_count = (*hash1)[hash2Iterator.get_key()];

                // Increment hash2's unique counters (don't bother with shared counters... we've already done this)
                comp_counters->updateHash2Counters(hash1_count, hash2_count);

                // Only bother updating thread matrix with K-mers not found in hash1 (we've already done the rest)
                if (!hash1_count) {
                    // Scale counters to make the matrix look pretty
                    uint64_t scaled_hash2_count = scaleCounter(hash2_count, args.d2_scale);

                    // Modifies hash counts so that K-mer counts larger than MATRIX_SIZE are dumped in the last slot
                    if (scaled_hash2_count >= args.d2_bins) scaled_hash2_count = args.d2_bins - 1;

                    // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
                    thread_matrix->inc(0, scaled_hash2_count, 1);
                }
            }

            // Only update hash3 counters if hash3 was provided
            if (hash3) {
                // Setup iterator for this thread's chunk of hash3
                typename hash_t::iterator hash3Iterator = hash3->iterator_slice(th_id, args->threads);

                // Iterate through this thread's slice of hash2
                while (hash3Iterator.next()) {
                    // Get the current K-mer count for hash2

                    uint64_t hash3_count = hash3Iterator.get_val();

                    // Increment hash3's unique counters (don't bother with shared counters... we've already done this)
                    comp_counters->updateHash3Counters(hash3_count);
                }
            }
        }



        // Threaded matrix data

        SM64 getMainMatrix() {

            return main_matrix ? main_matrix->getFinalMatrix() : NULL;
        }

        SM64 getEndsMatrix() {

            return ends_matrix ? ends_matrix->getFinalMatrix() : NULL;
        }

        SM64 getMiddleMatrix() {

            return middle_matrix ? middle_matrix->getFinalMatrix() : NULL;
        }

        SM64 getMixedMatrix() {

            return mixed_matrix ? mixed_matrix->getFinalMatrix() : NULL;
        }

        // Final data (created by merging thread results)

        shared_ptr<CompCounters> getCompCounters() {

            return final_comp_counters;
        }

        // Print Comp setup

        void printVars(ostream &out) {

            out << endl
                    << "Comp parameters:" << endl
                    << " - Hash 1: " << (jfh1 ? "mapped file configured" : "not specified") << endl
                    << " - Hash 2: " << (jfh2 ? "mapped file configured" : "not specified") << endl
                    << " - Hash 3: " << (jfh3 ? "mapped file configured" : "not specified") << endl
                    << " - Threads: " << args.threads << endl
                    << " - Dataset 1 scaling factor: " << args.d1_scale << endl
                    << " - Dataset 2 scaling factor: " << args.d2_scale << endl
                    << " - Dataset 1 bins: " << args.d1_bins << endl
                    << " - Dataset 2 bins: " << args.d2_bins << endl
                    << endl;
        }


        // Print K-mer comparison matrix

        void printMainMatrix(ostream &out) {

            SM64 mx = main_matrix->getFinalMatrix();

            out << mme::KEY_TITLE << "K-mer comparison plot" << endl
                    << mme::KEY_X_LABEL << "K-mer multiplicity for: " << args.db1_path << endl
                    << mme::KEY_Y_LABEL << "K-mer multiplicity for: " << args.db2_path << endl
                    << mme::KEY_Z_LABEL << "Distinct K-mers per bin" << endl
                    << mme::KEY_NB_COLUMNS << mx->height() << endl
                    << mme::KEY_NB_ROWS << mx->width() << endl
                    << mme::KEY_MAX_VAL << mx->getMaxVal() << endl
                    << mme::KEY_TRANSPOSE << "1" << endl
                    << mme::MX_META_END << endl;

            mx->printMatrix(out);
        }

        // Print K-mer comparison matrix

        void printEndsMatrix(ostream &out) {

            out << "# Each row represents K-mer multiplicity for: " << args.db1_path << endl;
            out << "# Each column represents K-mer multiplicity for sequence ends: " << args.db3_path << endl;

            ends_matrix->getFinalMatrix()->printMatrix(out);
        }

        // Print K-mer comparison matrix

        void printMiddleMatrix(ostream &out) {

            out << "# Each row represents K-mer multiplicity for: " << args.db1_path << endl;
            out << "# Each column represents K-mer multiplicity for sequence middles: " << args.db2_path << endl;

            middle_matrix->getFinalMatrix()->printMatrix(out);
        }

        // Print K-mer comparison matrix

        void printMixedMatrix(ostream &out) {

            out << "# Each row represents K-mer multiplicity for hash file 1: " << args.db1_path << endl;
            out << "# Each column represents K-mer multiplicity for mixed: " << args.db2_path << " and " << args.db3_path << endl;

            mixed_matrix->getFinalMatrix()->printMatrix(out);
        }

        // Print K-mer statistics

        void printCounters(ostream &out) {

            final_comp_counters->printCounts(out);
        }



    private:


        // Scale counters to make the matrix look pretty

        uint64_t scaleCounter(uint64_t count, double scale_factor) {

            return count == 0 ? 0 : (uint64_t) ceil((double) count * scale_factor);
        }

        // Combines each threads matrix into a single matrix

        void merge() {

            main_matrix->mergeThreadedMatricies();

            if (hash3) {
                ends_matrix->mergeThreadedMatricies();
                middle_matrix->mergeThreadedMatricies();
                mixed_matrix->mergeThreadedMatricies();
            }

            // Merge counters
            for (int k = 0; k < args.threads; k++) {
                shared_ptr<CompCounters> thread_comp_counter = (*thread_comp_counters)[k];

                final_comp_counters->hash1_total += thread_comp_counter->hash1_total;
                final_comp_counters->hash2_total += thread_comp_counter->hash2_total;
                final_comp_counters->hash3_total += thread_comp_counter->hash3_total;
                final_comp_counters->hash1_distinct += thread_comp_counter->hash1_distinct;
                final_comp_counters->hash2_distinct += thread_comp_counter->hash2_distinct;
                final_comp_counters->hash3_distinct += thread_comp_counter->hash3_distinct;
                final_comp_counters->hash1_only_total += thread_comp_counter->hash1_only_total;
                final_comp_counters->hash2_only_total += thread_comp_counter->hash2_only_total;
                final_comp_counters->hash1_only_distinct += thread_comp_counter->hash1_only_distinct;
                final_comp_counters->hash2_only_distinct += thread_comp_counter->hash2_only_distinct;
                final_comp_counters->shared_hash1_total += thread_comp_counter->shared_hash1_total;
                final_comp_counters->shared_hash2_total += thread_comp_counter->shared_hash2_total;
                final_comp_counters->shared_distinct += thread_comp_counter->shared_distinct;
            }
        }

    };
}
