#pragma once

#include <iostream>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <math.h>

#include <matrix/matrix_metadata_extractor.hpp>
#include <matrix/sparse_matrix.hpp>
#include <matrix/threaded_sparse_matrix.hpp>

#include <jellyfish/hash.hpp>
#include <jellyfish/counter.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/jellyfish_helper.hpp>

#include "comp_args.hpp"

using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::ostream;


class CompCounters{
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
        hash1_path(_hash1_path), hash2_path(_hash2_path), hash3_path(_hash3_path)
    {
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

    void updateHash1Counters(uint64_t hash1_count, uint64_t hash2_count)
    {
        hash1_total += hash1_count;
        hash1_distinct++;

        if (!hash2_count)
        {
            hash1_only_total += hash1_count;
            hash1_only_distinct++;
        }
    }

    void updateHash2Counters(uint64_t hash1_count, uint64_t hash2_count)
    {
        hash2_total += hash2_count;
        hash2_distinct++;

        if (!hash1_count)
        {
            hash2_only_total += hash2_count;
            hash2_only_distinct++;
        }
    }

    void updateHash3Counters(uint64_t hash3_count)
    {
        hash3_total += hash3_count;
        hash3_distinct++;
    }


    void updateSharedCounters(uint64_t hash1_count, uint64_t hash2_count)
    {
        if (hash1_count && hash2_count)
        {
            shared_hash1_total += hash1_count;
            shared_hash2_total += hash2_count;
            shared_distinct++;
        }
    }

    void printCounts(ostream &out)
    {
        out << "Kmer statistics for: " << endl;
        out << " - Hash 1: " << hash1_path << endl;
        out << " - Hash 2: " << hash2_path << endl;

        if (hash3_total > 0)
            out << " - Hash 3: " << hash3_path << endl;

        out << endl;

        out << "Total kmers in: " << endl;
        out << " - Hash 1: " << hash1_total << endl;
        out << " - Hash 2: " << hash2_total << endl;

        if (hash3_total > 0)
            out << " - Hash 3: " << hash3_total << endl;

        out << endl;

        out << "Distinct kmers in:" << endl;
        out << " - Hash 1: " << hash1_distinct << endl;
        out << " - Hash 2: " << hash2_distinct << endl;
        if (hash3_total > 0)
            out << " - Hash 3: " << hash3_distinct << endl;

        out << endl;

        out << "Total kmers only found in:" << endl;
        out << " - Hash 1: " << hash1_only_total << endl;
        out << " - Hash 2: " << hash2_only_total << endl;
        out << endl;

        out << "Distinct kmers only found in:" << endl;
        out << " - Hash 1: " << hash1_only_distinct << endl;
        out << " - Hash 2: " << hash2_only_distinct << endl << endl;

        out << "Shared kmers:" << endl;
        out << " - Total shared found in hash 1: " << shared_hash1_total << endl;
        out << " - Total shared found in hash 2: " << shared_hash2_total << endl;
        out << " - Distinct shared kmers: " << shared_distinct << endl << endl;
    }
};

template<typename hash_t>
class Comp : public thread_exec
{
private:
    // Args passed in
    const CompArgs          *args;

    // Jellyfish mapped file hash vars
    JellyfishHelper         *jfh1;
    JellyfishHelper         *jfh2;
    JellyfishHelper         *jfh3;
    hash_t                  *hash1;		// Jellyfish hash 1
    hash_t                  *hash2;     // Jellyfish hash 2
    hash_t                  *hash3;     // Jellyfish hash 3

    // Threaded matrix data
    ThreadedSparseMatrix<uint64_t>    *main_matrix;
    ThreadedSparseMatrix<uint64_t>    *ends_matrix;
    ThreadedSparseMatrix<uint64_t>    *middle_matrix;
    ThreadedSparseMatrix<uint64_t>    *mixed_matrix;

    // Final data (created by merging thread results)
    CompCounters            *final_comp_counters;

    // Thread specific data
    CompCounters            **thread_comp_counters;


public:
    Comp(CompArgs * _args) : args(_args)
    {
        if (args->verbose)
            cerr << "Setting up comp tool..." << endl;

        // Setup handles to load hashes
        jfh1 = new JellyfishHelper(args->db1_path);
        jfh2 = new JellyfishHelper(args->db2_path);
        jfh3 = args->db3_path ? new JellyfishHelper(args->db3_path) : NULL;

        // Ensure all hashes are null at this stage (we'll load them later)
        hash1 = NULL;
        hash2 = NULL;
        hash3 = NULL;

        // Create the final kmer counter matricies
        main_matrix = new ThreadedSparseMatrix<uint64_t>(args->d1_bins, args->d2_bins, args->threads);

        // Initialise extra matricies for hash3 (only allocates space if required)
        if (args->db3_path)
        {
            if (args->verbose)
                cerr << " - Setting up matricies for hash 3" << endl;

            ends_matrix = new ThreadedSparseMatrix<uint64_t>(args->d1_bins, args->d2_bins, args->threads);
            middle_matrix = new ThreadedSparseMatrix<uint64_t>(args->d1_bins, args->d2_bins, args->threads);
            mixed_matrix = new ThreadedSparseMatrix<uint64_t>(args->d1_bins, args->d2_bins, args->threads);
        }
        else
        {
            ends_matrix = NULL;
            middle_matrix = NULL;
            mixed_matrix = NULL;
        }

        // Create the final comp counters
        final_comp_counters = new CompCounters(args->db1_path, args->db2_path, args->db3_path);

        // Create the comp counters for each thread
        thread_comp_counters = new CompCounters*[args->threads];
        for(int i = 0; i < args->threads; i++) {
            thread_comp_counters[i] = new CompCounters(args->db1_path, args->db2_path, args->db3_path);
        }

        if (args->verbose)
            cerr << "Comp tool setup without error." << endl;
    }

    ~Comp()
    {
        destroyFinalMatricies();

        destroyCounters();

        destroyJellyfishHashes();
    }


    void do_it()
    {
        if (args->verbose)
        {
            cerr << "Loading hashes..." << endl;
        }

        std::ostream* out_stream = args->verbose ? &cerr : (std::ostream*)0;

        // Load the hashes
        hash1 = jfh1->loadHash(true, out_stream);
        hash2 = jfh2->loadHash(false, out_stream);

        if (jfh3)
            hash3 = jfh3->loadHash(false, out_stream);

        // Whether to treat this hash as double stranded or not.
        // Ideally it would be nice to determine this directly from the hash but I'm
        // not sure how to do that at the moment... it might not be possible
        if (args->both_strands)
        {
            hash1->set_canonical(true);
            hash2->set_canonical(true);
            if (hash3)
                hash3->set_canonical(true);
        }

        // Check kmer lengths are the same for both hashes.  We can't continue if they are not.
        if (hash1->get_mer_len() != hash2->get_mer_len())
        {
            cerr << "Cannot process hashes that were created with different kmer lengths.  "
                 << "Hash1: " << hash1->get_mer_len() << ".  Hash2: " << hash2->get_mer_len() << "." << endl;
            throw;
        }
        else if(hash3 && (hash3->get_mer_len() != hash1->get_mer_len()))
        {
            cerr << "Cannot process hashes that were created with different kmer lengths.  "
                 << "Hash1: " << hash1->get_mer_len() << ".  Hash3: " << hash3->get_mer_len() << "." << endl;
            throw;
        }

        if (args->verbose)
            cerr << endl
                 << "All hashes loaded successfully." << endl
                 << "Starting threads...";

        // Run the threads
        exec_join(args->threads);

        if (args->verbose)
            cerr << "done." << endl
                 << "Merging results...";

        // Merge results from the threads
        merge();

        if (args->verbose)
            cerr << "done." << endl;
    }

    void start(int th_id)
    {
        // Get handle to sparse matrix for this thread
        SparseMatrix<uint64_t>* thread_matrix = main_matrix->getThreadMatrix(th_id);

        // Get handle on this thread's comp counter
        CompCounters* comp_counters = thread_comp_counters[th_id];

        // Setup iterator for this thread's chunk of hash1
        typename hash_t::iterator hash1Iterator = hash1->iterator_slice(th_id, args->threads);

        // Go through this thread's slice for hash1
        while (hash1Iterator.next())
        {
            // Get the current kmer count for hash1
            uint64_t hash1_count = hash1Iterator.get_val();

            // Get the count for this kmer in hash2 (assuming it exists... 0 if not)
            uint64_t hash2_count = (*hash2)[hash1Iterator.get_key()];

            // Get the count for this kmer in hash3 (assuming it exists... 0 if not)
            uint64_t hash3_count = hash3 ? (*hash3)[hash1Iterator.get_key()] : 0;

            // Increment hash1's unique counters
            comp_counters->updateHash1Counters(hash1_count, hash2_count);

            // Increment shared counters
            comp_counters->updateSharedCounters(hash1_count, hash2_count);

            // Scale counters to make the matrix look pretty
            uint64_t scaled_hash1_count = scaleCounter(hash1_count, args->d1_scale);
            uint64_t scaled_hash2_count = scaleCounter(hash2_count, args->d2_scale);
            uint64_t scaled_hash3_count = scaleCounter(hash3_count, args->d2_scale);

            // Modifies hash counts so that kmer counts larger than MATRIX_SIZE are dumped in the last slot
            if (scaled_hash1_count >= args->d1_bins) scaled_hash1_count = args->d1_bins - 1;
            if (scaled_hash2_count >= args->d2_bins) scaled_hash2_count = args->d2_bins - 1;
            if (scaled_hash3_count >= args->d2_bins) scaled_hash3_count = args->d2_bins - 1;

            // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
            thread_matrix->inc(scaled_hash1_count, scaled_hash2_count, 1);

            // Update hash 3 related matricies if hash 3 was provided
            if (hash3)
            {
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
        typename hash_t::iterator hash2Iterator = hash2->iterator_slice(th_id, args->threads);

        // Iterate through this thread's slice of hash2
        while (hash2Iterator.next())
        {
            // Get the current kmer count for hash2
            uint64_t hash2_count = hash2Iterator.get_val();

            // Get the count for this kmer in hash1 (assuming it exists... 0 if not)
            uint64_t hash1_count = (*hash1)[hash2Iterator.get_key()];

            // Increment hash2's unique counters (don't bother with shared counters... we've already done this)
            comp_counters->updateHash2Counters(hash1_count, hash2_count);

            // Only bother updating thread matrix with kmers not found in hash1 (we've already done the rest)
            if (!hash1_count)
            {
                // Scale counters to make the matrix look pretty
                uint64_t scaled_hash2_count = scaleCounter(hash2_count, args->d2_scale);

                // Modifies hash counts so that kmer counts larger than MATRIX_SIZE are dumped in the last slot
                if (scaled_hash2_count >= args->d2_bins) scaled_hash2_count = args->d2_bins - 1;

                // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
                thread_matrix->inc(0, scaled_hash2_count, 1);
            }
        }

        // Only update hash3 counters if hash3 was provided
        if (hash3)
        {
            // Setup iterator for this thread's chunk of hash3
            typename hash_t::iterator hash3Iterator = hash3->iterator_slice(th_id, args->threads);

            // Iterate through this thread's slice of hash2
            while (hash3Iterator.next())
            {
                // Get the current kmer count for hash2
                uint64_t hash3_count = hash3Iterator.get_val();

                // Increment hash3's unique counters (don't bother with shared counters... we've already done this)
                comp_counters->updateHash3Counters(hash3_count);
            }
        }
    }



    // Print Comp setup
    void printVars(ostream &out)
    {
        out << endl;
        out << "Comp parameters:" << endl;
        out << " - Hash 1: " << (jfh1 ? "mapped file configured" : "not specified") << endl;
        out << " - Hash 2: " << (jfh2 ? "mapped file configured" : "not specified") << endl;
        out << " - Hash 3: " << (jfh3 ? "mapped file configured" : "not specified") << endl;
        out << " - Threads: " << args->threads << endl;
        out << " - Dataset 1 scaling factor: " << args->d1_scale << endl;
        out << " - Dataset 2 scaling factor: " << args->d2_scale << endl;
        out << " - Dataset 1 bins: " << args->d1_bins << endl;
        out << " - Dataset 2 bins: " << args->d2_bins << endl;
        out << endl;
   }


    // Print kmer comparison matrix
    void printMainMatrix(ostream &out)
    {
        SparseMatrix<uint64_t>* mx = main_matrix->getFinalMatrix();

        out << mme::KEY_TITLE << "Kmer comparison plot" << endl;
        out << mme::KEY_X_LABEL << "Kmer multiplicity for: " << args->db1_path << endl;
        out << mme::KEY_Y_LABEL << "Kmer multiplicity for: " << args->db2_path << endl;
        out << mme::KEY_Z_LABEL << "Distinct Kmers per bin" << endl;
        out << mme::KEY_NB_COLUMNS << mx->height() << endl;
        out << mme::KEY_NB_ROWS << mx->width() << endl;
        out << mme::KEY_MAX_VAL << mx->getMaxVal() << endl;
        out << mme::MX_META_END << endl;

        mx->printMatrix(out);
    }

    // Print kmer comparison matrix
    void printEndsMatrix(ostream &out)
    {
        out << "# Each row represents kmer multiplicity for: " << args->db1_path << endl;
        out << "# Each column represents kmer multiplicity for sequence ends: " << args->db3_path << endl;

        ends_matrix->getFinalMatrix()->printMatrix(out);
    }

    // Print kmer comparison matrix
    void printMiddleMatrix(ostream &out)
    {
        out << "# Each row represents kmer multiplicity for: " << args->db1_path << endl;
        out << "# Each column represents kmer multiplicity for sequence middles: " << args->db2_path << endl;

        middle_matrix->getFinalMatrix()->printMatrix(out);
    }

    // Print kmer comparison matrix
    void printMixedMatrix(ostream &out)
    {
        out << "# Each row represents kmer multiplicity for hash file 1: " << args->db1_path << endl;
        out << "# Each column represents kmer multiplicity for mixed: " << args->db2_path << " and " << args->db3_path << endl;

        mixed_matrix->getFinalMatrix()->printMatrix(out);
    }

    // Print kmer statistics
    void printCounters(ostream &out)
    {
       final_comp_counters->printCounts(out);
    }



private:


    // Scale counters to make the matrix look pretty
    uint64_t scaleCounter(uint64_t count, double scale_factor)
    {
        return count == 0 ? 0 : (uint64_t)ceil((double)count * scale_factor);
    }

    // Combines each threads matrix into a single matrix
    void merge() {

        main_matrix->mergeThreadedMatricies();

        if (hash3)
        {
            ends_matrix->mergeThreadedMatricies();
            middle_matrix->mergeThreadedMatricies();
            mixed_matrix->mergeThreadedMatricies();
        }

        // Merge counters
        for(int k = 0; k < args->threads; k++)
        {
            CompCounters* thread_comp_counter = thread_comp_counters[k];

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


    void destroyCounters()
    {
        if (final_comp_counters)
            delete final_comp_counters;

        final_comp_counters = NULL;

        if (thread_comp_counters)
        {
            for(int i = 0; i < args->threads; i++)
            {
                if (thread_comp_counters[i])
                {
                    delete thread_comp_counters[i];
                    thread_comp_counters[i] = NULL;
                }
            }

            delete thread_comp_counters;
            thread_comp_counters = NULL;
        }
    }

    void destroyJellyfishHashes()
    {
        if (jfh1)
            delete jfh1;

        if (jfh2)
            delete jfh2;

        if (jfh3)
            delete jfh3;

        jfh1 = NULL;
        jfh2 = NULL;
        jfh3 = NULL;
    }

    void destroyFinalMatricies()
    {
        if(main_matrix)
            delete main_matrix;

        if(ends_matrix)
            delete ends_matrix;

        if(middle_matrix)
            delete middle_matrix;

        if(mixed_matrix)
            delete mixed_matrix;

        main_matrix = NULL;
        ends_matrix = NULL;
        middle_matrix = NULL;
        mixed_matrix = NULL;
    }


};
