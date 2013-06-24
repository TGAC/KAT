#if !defined(DEF_COMP_H)
#define DEF_COMP_H

#include <iostream>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <math.h>

#include <matrix/sparse_matrix.hpp>

#include <jellyfish/counter.hpp>
#include <jellyfish/thread_exec.hpp>

using std::vector;
using std::string;
using std::cerr;
using std::endl;

class CompCounters{
public:
    uint64_t hash1_total;
    uint64_t hash2_total;
    uint64_t hash1_distinct;
    uint64_t hash2_distinct;
    uint64_t hash1_only_total;
    uint64_t hash2_only_total;
    uint64_t hash1_only_distinct;
    uint64_t hash2_only_distinct;
    uint64_t shared_hash1_total;
    uint64_t shared_hash2_total;
    uint64_t shared_distinct;

    CompCounters()
    {
        hash1_total = 0;
        hash2_total = 0;
        hash1_distinct = 0;
        hash2_distinct = 0;
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

    void updateSharedCounters(uint64_t hash1_count, uint64_t hash2_count)
    {
        if (hash1_count && hash2_count)
        {
            shared_hash1_total += hash1_count;
            shared_hash2_total += hash2_count;
            shared_distinct++;
        }
    }

    void printCounts(std::ostream &out)
    {
        out << "Kmer statistics for hash1 and hash2:" << endl << endl;

        out << " - Total kmers in hash 1: " << hash1_total << endl;
        out << " - Total kmers in hash 2: " << hash2_total << endl << endl;

        out << " - Distinct kmers in hash 1: " << hash1_distinct << endl;
        out << " - Distinct kmers in hash 2: " << hash2_distinct << endl << endl;

        out << " - Total kmers only found in hash 1: " << hash1_only_total << endl;
        out << " - Total kmers only found in hash 2: " << hash2_only_total << endl << endl;

        out << " - Distinct kmers only found in hash 1: " << hash1_only_distinct << endl;
        out << " - Distinct kmers only found in hash 2: " << hash2_only_distinct << endl << endl;

        out << " - Total shared kmers found in hash 1: " << shared_hash1_total << endl;
        out << " - Total shared kmers found in hash 2: " << shared_hash2_total << endl;
        out << " - Distinct shared kmers: " << shared_distinct << endl << endl;
    }
};

#define KMER_LIMIT 1000
#define MATRIX_SIZE KMER_LIMIT+1        // Adds one to kmer limit so we can add absent kmers

template<typename hash_t>
class Comp : public thread_exec
{
    const hash_t            *hash1;		// Jellyfish hash 1
    const hash_t            *hash2;     // Jellyfish hash 2
    const vector<string>    *names;	    // Names of fasta sequences (in same order as seqs)
    const vector<string>    *seqs;	    // Sequences in fasta sequences (in same order as names)
    const uint_t            threads;	// Number of threads to use
    const float             xscale, yscale;  // Scaling factors to make the matrix look pretty.
    const bool              noindex;    // Whether or not to output an index value for the first column of the matrix.

    // Final data (created by merging thread results)
    SparseMatrix<uint_t>    *final_matrix;
    CompCounters            *final_comp_counters;

    // Thread specific data
    SparseMatrix<uint_t>    **thread_matricies;
    CompCounters            **thread_comp_counters;


public:
    Comp(const hash_t *_hash1, const hash_t *_hash2, const vector<string> *_names, const vector<string> *_seqs,
        uint_t _threads, float _xscale, float _yscale, bool _noindex) :
        hash1(_hash1), hash2(_hash2), names(_names), seqs(_seqs),
        threads(_threads),
        xscale(_xscale), yscale(_yscale), noindex(_noindex)
    {
        // Create the final kmer counter matrix
        final_matrix = new SparseMatrix<uint_t>(MATRIX_SIZE, MATRIX_SIZE);

        // Create kmer counter matricies for each thread
        thread_matricies = new SparseMatrix<uint_t>*[threads];
        for(int i = 0; i < threads; i++) {
            thread_matricies[i] = new SparseMatrix<uint_t>(MATRIX_SIZE, MATRIX_SIZE);
        }

        // Create the final comp counters
        final_comp_counters = new CompCounters();

        // Create the comp counters for each thread
        thread_comp_counters = new CompCounters*[threads];
        for(int i = 0; i < threads; i++) {
            thread_comp_counters[i] = new CompCounters();
        }
    }

    ~Comp()
    {
        if(final_matrix)
            delete final_matrix;

        if (thread_matricies)
        {
            for(int i = 0; i < threads; i++)
            {
                if (thread_matricies[i])
                    delete thread_matricies[i];
            }
        }

        if (final_comp_counters)
            delete final_comp_counters;

        if (thread_comp_counters)
        {
            for(int i = 0; i < threads; i++)
            {
                if (thread_comp_counters[i])
                    delete thread_comp_counters[i];
            }
        }
    }


    void do_it()
    {
        exec_join(threads);

        merge();
    }

    void start(int th_id)
    {
        // Get handle to sparse matrix for this thread
        SparseMatrix<uint_t>* thread_matrix = thread_matricies[th_id];

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

            // Increment hash1's unique counters
            comp_counters->updateHash1Counters(hash1_count, hash2_count);

            // Increment shared counters
            comp_counters->updateSharedCounters(hash1_count, hash2_count);

            // Scale counters to make the matrix look pretty
            uint64_t scaled_hash1_count = scaleCounter(hash1_count, xscale);
            uint64_t scaled_hash2_count = scaleCounter(hash2_count, yscale);

            // Modifies hash counts so that kmer counts larger than MATRIX_SIZE are dumped in the last slot
            if (scaled_hash1_count > KMER_LIMIT) scaled_hash1_count = KMER_LIMIT;
            if (scaled_hash2_count > KMER_LIMIT) scaled_hash2_count = KMER_LIMIT;

            // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
            thread_matrix->inc(scaled_hash1_count, scaled_hash2_count, 1);
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
                uint64_t scaled_hash2_count = scaleCounter(hash2_count, yscale);

                // Modifies hash counts so that kmer counts larger than MATRIX_SIZE are dumped in the last slot
                if (scaled_hash2_count > KMER_LIMIT) scaled_hash2_count = KMER_LIMIT;

                // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
                thread_matrix->inc(0, scaled_hash2_count, 1);
            }
        }
    }



    // Print Comp setup
    void printVars(std::ostream &out)
    {
        out << "Comp parameters:" << endl;
        out << " - Hash1: " << (hash1 ? "present" : "not specified") << endl;
        out << " - Hash2: " << (hash2 ? "present" : "not specified") << endl;
        out << " - Threads: " << threads << endl;
        out << " - X scaling factor: " << xscale << endl;
        out << " - Y scaling factor: " << yscale << endl;
        out << " - No index: " << noindex << endl;
   }


    // Print kmer comparison matrix
    void printMatrix(std::ostream &out)
    {
        uint_t iMax = noindex ? MATRIX_SIZE : MATRIX_SIZE + 1;

        for(int i = 0; i < iMax; i++)
        {
            if (!noindex)
                out << i << " ";

            for(int j = 0; j < MATRIX_SIZE; j++)
            {
                out << final_matrix->get(i, j);

                if (j != MATRIX_SIZE - 1)
                    out << " ";
            }
            out << endl;
        }
    }

    // Print kmer statistics
    void printCounters(std::ostream &out)
    {
       final_comp_counters->printCounts(out);
    }



private:


    // Scale counters to make the matrix look pretty
    uint64_t scaleCounter(uint64_t count, float scale_factor)
    {
        return count == 0 ? 0 : ceil(count * scale_factor);
    }

    // Combines each threads matrix into a single matrix
    void merge() {

        // Merge matrix
        for(int i = 0; i < KMER_LIMIT; i++)
        {
            for(int j = 0; j < KMER_LIMIT; j++)
            {
                for(int k = 0; k < threads; k++)
                {
                    final_matrix->inc(i, j, thread_matricies[k]->get(i, j));
                }
            }
        }

        // Merge counters
        for(int k = 0; k < threads; k++)
        {
            CompCounters* thread_comp_counter = thread_comp_counters[k];

            final_comp_counters->hash1_total += thread_comp_counter->hash1_total;
            final_comp_counters->hash2_total += thread_comp_counter->hash2_total;
            final_comp_counters->hash1_distinct += thread_comp_counter->hash1_distinct;
            final_comp_counters->hash2_distinct += thread_comp_counter->hash2_distinct;
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

#endif //DEF_COMP_H
