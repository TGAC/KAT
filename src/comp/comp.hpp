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
    const char              *jf_hash_path_1;
    const char              *jf_hash_path_2;
    const char              *jf_hash_path_3;
    const uint16_t          threads;	// Number of threads to use
    const float             d1_scale, d2_scale;  // Scaling factors to make the matrix look pretty.
    const uint16_t          d1_bins, d2_bins;  // Scaling factors to make the matrix look pretty.
    const bool              verbose;

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
    Comp(const char *_jfHashPath1, const char *_jfHashPath2, const char *_jfHashPath3,
        uint16_t _threads, float _xscale, float _yscale, uint16_t _d1_bins, uint16_t _d2_bins, bool _verbose) :
        jf_hash_path_1(_jfHashPath1), jf_hash_path_2(_jfHashPath2), jf_hash_path_3(_jfHashPath3),
        threads(_threads),
        d1_scale(_xscale), d2_scale(_yscale),
        d1_bins(_d1_bins), d2_bins(_d2_bins),
        verbose(_verbose)
    {
        if (verbose)
            cerr << "Setting up comp tool..." << endl;

        // Setup handles to load hashes
        jfh1 = new JellyfishHelper(jf_hash_path_1);
        jfh2 = new JellyfishHelper(jf_hash_path_2);
        jfh3 = jf_hash_path_3 ? new JellyfishHelper(jf_hash_path_3) : NULL;

        // Ensure all hashes are null at this stage (we'll load them later)
        hash1 = NULL;
        hash2 = NULL;
        hash3 = NULL;

        // Create the final kmer counter matricies
        main_matrix = new ThreadedSparseMatrix<uint64_t>(d1_bins, d2_bins, threads);

        // Initialise extra matricies for hash3 (only allocates space if required)
        if (jf_hash_path_3)
        {
            if (verbose)
                cerr << " - Setting up matricies for hash 3" << endl;

            ends_matrix = new ThreadedSparseMatrix<uint64_t>(d1_bins, d2_bins, threads);
            middle_matrix = new ThreadedSparseMatrix<uint64_t>(d1_bins, d2_bins, threads);
            mixed_matrix = new ThreadedSparseMatrix<uint64_t>(d1_bins, d2_bins, threads);
        }
        else
        {
            ends_matrix = NULL;
            middle_matrix = NULL;
            mixed_matrix = NULL;
        }

        // Create the final comp counters
        final_comp_counters = new CompCounters(jf_hash_path_1, jf_hash_path_2, jf_hash_path_3);

        // Create the comp counters for each thread
        thread_comp_counters = new CompCounters*[threads];
        for(int i = 0; i < threads; i++) {
            thread_comp_counters[i] = new CompCounters(jf_hash_path_1, jf_hash_path_2, jf_hash_path_3);
        }

        if (verbose)
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
        if (verbose)
        {
            cerr << "Loading hashes..." << endl;
        }

        std::ostream* out_stream = verbose ? &cerr : (std::ostream*)0;

        // Load the hashes
        hash1 = jfh1->loadHash(true, out_stream);
        hash2 = jfh2->loadHash(false, out_stream);

        if (jfh3)
            hash3 = jfh3->loadHash(false, out_stream);

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

        if (verbose)
            cerr << endl
                 << "All hashes loaded successfully." << endl
                 << "Starting threads...";

        // Run the threads
        exec_join(threads);

        if (verbose)
            cerr << "done." << endl
                 << "Merging results...";

        // Merge results from the threads
        merge();

        if (verbose)
            cerr << "done." << endl;
    }

    void start(int th_id)
    {
        // Get handle to sparse matrix for this thread
        SparseMatrix<uint64_t>* thread_matrix = main_matrix->getThreadMatrix(th_id);

        // Get handle on this thread's comp counter
        CompCounters* comp_counters = thread_comp_counters[th_id];

        // Setup iterator for this thread's chunk of hash1
        typename hash_t::iterator hash1Iterator = hash1->iterator_slice(th_id, threads);

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
            uint64_t scaled_hash1_count = scaleCounter(hash1_count, d1_scale);
            uint64_t scaled_hash2_count = scaleCounter(hash2_count, d2_scale);
            uint64_t scaled_hash3_count = scaleCounter(hash3_count, d2_scale);

            // Modifies hash counts so that kmer counts larger than MATRIX_SIZE are dumped in the last slot
            if (scaled_hash1_count >= d1_bins) scaled_hash1_count = d1_bins - 1;
            if (scaled_hash2_count >= d2_bins) scaled_hash2_count = d2_bins - 1;
            if (scaled_hash3_count >= d2_bins) scaled_hash3_count = d2_bins - 1;

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
        typename hash_t::iterator hash2Iterator = hash2->iterator_slice(th_id, threads);

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
                uint64_t scaled_hash2_count = scaleCounter(hash2_count, d2_scale);

                // Modifies hash counts so that kmer counts larger than MATRIX_SIZE are dumped in the last slot
                if (scaled_hash2_count >= d2_bins) scaled_hash2_count = d2_bins - 1;

                // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
                thread_matrix->inc(0, scaled_hash2_count, 1);
            }
        }

        // Only update hash3 counters if hash3 was provided
        if (hash3)
        {
            // Setup iterator for this thread's chunk of hash3
            typename hash_t::iterator hash3Iterator = hash3->iterator_slice(th_id, threads);

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
        out << " - Threads: " << threads << endl;
        out << " - Dataset 1 scaling factor: " << d1_scale << endl;
        out << " - Dataset 2 scaling factor: " << d2_scale << endl;
        out << " - Dataset 1 bins: " << d1_bins << endl;
        out << " - Dataset 2 bins: " << d2_bins << endl;
        out << endl;
   }


    // Print kmer comparison matrix
    void printMainMatrix(ostream &out)
    {
        SparseMatrix<uint64_t>* mx = main_matrix->getFinalMatrix();

        out << mme::KEY_TITLE << "Kmer comparison plot" << endl;
        out << mme::KEY_X_LABEL << "Kmer multiplicity for: " << jf_hash_path_1 << endl;
        out << mme::KEY_Y_LABEL << "Kmer multiplicity for: " << jf_hash_path_2 << endl;
        out << mme::KEY_NB_COLUMNS << mx->height() << endl;
        out << mme::KEY_NB_ROWS << mx->width() << endl;
        out << mme::KEY_MAX_VAL << mx->getMaxVal() << endl;
        out << mme::MX_META_END << endl;

        mx->printMatrix(out);
    }

    // Print kmer comparison matrix
    void printEndsMatrix(ostream &out)
    {
        out << "# Each row represents kmer multiplicity for: " << jf_hash_path_1 << endl;
        out << "# Each column represents kmer multiplicity for sequence ends: " << jf_hash_path_3 << endl;

        ends_matrix->getFinalMatrix()->printMatrix(out);
    }

    // Print kmer comparison matrix
    void printMiddleMatrix(ostream &out)
    {
        out << "# Each row represents kmer multiplicity for: " << jf_hash_path_1 << endl;
        out << "# Each column represents kmer multiplicity for sequence middles: " << jf_hash_path_2 << endl;

        middle_matrix->getFinalMatrix()->printMatrix(out);
    }

    // Print kmer comparison matrix
    void printMixedMatrix(ostream &out)
    {
        out << "# Each row represents kmer multiplicity for hash file 1: " << jf_hash_path_1 << endl;
        out << "# Each column represents kmer multiplicity for mixed: " << jf_hash_path_2 << " and " << jf_hash_path_3 << endl;

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
        for(int k = 0; k < threads; k++)
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
            for(int i = 0; i < threads; i++)
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
