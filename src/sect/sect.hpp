#if !defined(DEF_SECT_H)
#define DEF_SECT_H

#include <iostream>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <math.h>

#include <jellyfish/counter.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/jellyfish_helper.hpp>

#include <matrix/threaded_sparse_matrix.hpp>

using std::vector;
using std::string;
using std::cerr;
using std::endl;

template<typename hash_t>
class Sect : public thread_exec
{
private:
    const char                      *jfHashPath;
    const vector<string>            *names;	    // Names of fasta sequences (in same order as seqs)
    const vector<string>            *seqs;	    // Sequences in fasta sequences (in same order as names)
    const uint16_t                  gc_bins;
    const uint16_t                  cvg_bins;
    const bool                      cvg_logscale;
    const uint16_t                  threads;	// Number of threads to use
    const bool                      verbose;
    const size_t                    bucket_size, remaining;	// Chunking vars

    JellyfishHelper                 *jfh;
    hash_t                          *hash;		// Jellyfish hash
    vector<vector<uint64_t>*>       *counts;    // Kmer counts for each kmer window in sequence (in same order as seqs and names; built by this class)
    vector<float>                   *coverages; // Overall coverage calculated for each sequence from the kmer windows.
    vector<float>                   *gcs;       // GC% for each sequence
    vector<uint32_t>                *lengths;       // GC% for each sequence
    ThreadedSparseMatrix<uint64_t>  *contamination_mx;  // Stores cumulative base count for each sequence where GC and CVG are binned

public:
    Sect(const char *_jfHashPath, const vector<string> *_names, const vector<string> *_seqs,
        uint16_t _gc_bins, uint16_t _cvg_bins, float _cvg_logscale, uint16_t _threads, bool _verbose) :
        jfHashPath(_jfHashPath), names(_names), seqs(_seqs),
        gc_bins(_gc_bins), cvg_bins(_cvg_bins), cvg_logscale(_cvg_logscale),
        threads(_threads), verbose(_verbose),
        bucket_size(seqs->size() / threads),
        remaining(seqs->size() % (bucket_size < 1 ? 1 : threads))
    {
        // Setup handle to jellyfish hash
        jfh = new JellyfishHelper(jfHashPath);

        // Setup space for storing output
        counts = new vector<vector<uint64_t>*>(seqs->size());
        coverages = new vector<float>(seqs->size());
        gcs = new vector<float>(seqs->size());
        lengths = new vector<uint32_t>(seqs->size());
        contamination_mx = new ThreadedSparseMatrix<uint64_t>(gc_bins, cvg_bins, threads);
    }

    ~Sect()
    {
        // Jellyfish helper takes care of freeing the hash
        //if (hash)
        //    delete hash;

        if (jfh)
            delete jfh;

        jfh = NULL;

        if(counts)
        {
            for(uint_t i = 0; i < counts->size(); i++)
            {
                delete (*counts)[i];
                (*counts)[i] = NULL;
            }
            delete counts;
            counts = NULL;
        }

        if (coverages)
            delete coverages;

        coverages = NULL;

        if (gcs)
            delete gcs;

        gcs = NULL;

        if (lengths)
            delete lengths;

        lengths = NULL;

        if (contamination_mx)
            delete contamination_mx;

        contamination_mx = NULL;

    }


    void do_it()
    {
        std::ostream* out_stream = verbose ? &cerr : (std::ostream*)0;

        // Load the jellyfish hash
        hash = jfh->loadHash(true, out_stream);

        // Process each fasta sequence in a different thread.
        // In each thread lookup each kmer in the hash
        exec_join(threads);

        contamination_mx->mergeThreadedMatricies();
    }

    void start(int th_id)
    {
        // Check to see if we have useful work to do for this thread, return if not
        if (bucket_size < 1 && th_id >= seqs->size())
        {
            return;
        }

        //processInBlocks(th_id);
        processInterlaced(th_id);
    }


    void printVars(std::ostream &out)
    {
        out << "SECT parameters:" << endl;
        out << " - Hash: " << (hash ? "present" : "not specified") << endl;
        out << " - Sequences to process: " << seqs->size() << endl;
        out << " - Threads: " << threads << endl;
        out << " - Bucket size: " << bucket_size << endl;
        out << " - Remaining: " << remaining << endl << endl;
    }


    void printCounts(std::ostream &out)
    {
        for(int i = 0; i < names->size(); i++)
        {
            out << ">" << (*names)[i].c_str() << endl;

            vector<uint64_t>* seqCounts = (*counts)[i];

            if (seqCounts != NULL && !seqCounts->empty())
            {
                out << (*seqCounts)[0];

                for(uint_t j = 1; j < seqCounts->size(); j++)
                {
                    out << " " << (*seqCounts)[j];
                }

                out << endl;
            }
            else
            {
                out << "0" << endl;
            }
        }
    }

    void printStatTable(std::ostream &out)
    {
        out << "seq_name coverage gc% seq_length" << endl;
        for(int i = 0; i < names->size(); i++)
        {
            out << (*names)[i] << " " << (*coverages)[i] << " " << (*gcs)[i] << " " << (*lengths)[i] << endl;
        }
    }

    // Print kmer comparison matrix
    void printContaminationMatrix(std::ostream &out)
    {
        out << "# Fasta file processed: " << "???" << endl;
        out << "# Jellyfish hash: " << jfHashPath << endl;
        out << "# Each column represents the GC%, with " << gc_bins << " bins. First bin is 0% and last bin in 100%" << endl;
        out << "# Each row represents the sequence coverage, with " << cvg_bins << " bins.  First bin in 0x last bin in " << endl;

        contamination_mx->getFinalMatrix()->printMatrix(out);
    }


private:

    // This method won't be optimal in most cases... Fasta files are normally sorted by length (largest first)
    // So first thread will be asked to do more work than the rest
    void processInBlocks(uint_t th_id)
    {
        size_t start = bucket_size < 1 ? th_id : th_id * bucket_size;
        size_t end = bucket_size < 1 ? th_id : start + bucket_size - 1;
        for(size_t i = start; i <= end; i++)
        {
            processSeq(i);
        }

        // Process a remainder if required
        if (th_id < remaining)
        {
            size_t rem_idx = (threads * bucket_size) + th_id;
            processSeq(rem_idx);
        }
    }

    // This method is probably makes more efficient use of multiple cores on a length sorted fasta file
    void processInterlaced(uint_t th_id)
    {
        size_t start = th_id;
        size_t end = seqs->size();
        for(size_t i = start; i < end; i+=threads)
        {
            processSeq(i, th_id);
        }
    }

    void processSeq(const size_t index, const uint_t th_id)
    {
        uint_t kmer = hash->get_mer_len();
        string seq = (*seqs)[index];
        uint_t nbCounts = seq.length() - kmer + 1;
        float mean_cvg = 0;

        if (seq.length() < kmer)
        {
            cerr << (*names)[index].c_str() << ": " << seq << "  is too short to compute coverage.  Setting sequence coverage to 0." << endl;
        }
        else
        {

            vector<uint_t>* seqCounts = new vector<uint_t>(nbCounts);

            uint64_t sum = 0;

            for(uint_t i = 0; i < nbCounts; i++)
            {
                string merstr = seq.substr(i, kmer);

                // Jellyfish compacted hash does not support Ns so if we find one set this mer count to 0
                if (merstr.find("N") != string::npos)
                {
                    (*seqCounts)[i] = 0;
                }
                else
                {
                    const char* mer = merstr.c_str();
                    uint_t count = (*hash)[mer];
                    sum += count;

                    (*seqCounts)[i] = count;
                }
            }

            (*counts)[index] = seqCounts;

            // Assumes simple mean calculation for sequence coverage for now... plug in Bernardo's method later.
            mean_cvg = (float)sum / (float)nbCounts;
            (*coverages)[index] = mean_cvg;
        }

        // Add length
        (*lengths)[index] = seq.length();

        // Calc GC%
        uint32_t gs = 0;
        uint32_t cs = 0;
        uint32_t ns = 0;

        for(uint32_t i = 0; i < seq.length(); i++)
        {
            char c = seq[i];

            if (c == 'G' || c == 'g')
                gs++;
            else if (c == 'C' || c == 'c')
                cs++;
            else if (c == 'N' || c == 'n')
                ns++;
        }

        float gc_perc = ((float)(gs + cs)) / ((float)(seq.length() - ns));
        (*gcs)[index] = gc_perc;

        float log_cvg = cvg_logscale ? log10(mean_cvg) : mean_cvg;

        // Assume log_cvg 5 is max value
        float compressed_cvg = cvg_logscale ? log_cvg * (cvg_bins / 5.0) : mean_cvg * 0.1;

        uint16_t x = gc_perc * gc_bins;  // Convert float to 1.dp
        uint16_t y = compressed_cvg >= cvg_bins ? cvg_bins - 1 : compressed_cvg;      // Simply cap the y value

        // Add bases to matrix
        contamination_mx->getThreadMatrix(th_id)->inc(x, y, seq.length());
    }

};

#endif //DEF_SECT_H
