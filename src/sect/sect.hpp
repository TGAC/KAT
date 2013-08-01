#pragma once

#include <iostream>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <math.h>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <jellyfish/hash.hpp>
#include <jellyfish/counter.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/jellyfish_helper.hpp>

#include <matrix/threaded_sparse_matrix.hpp>

#include "sect_args.hpp"

using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::stringstream;

using seqan::SequenceStream;
using seqan::StringSet;
using seqan::CharString;
using seqan::Dna5String;


#define BATCH_SIZE 1024

template<typename hash_t>
class Sect : public thread_exec
{
private:

    // Input args
    const SectArgs                  *args;
    const size_t                    bucket_size, remaining;	// Chunking vars

    // Variables that live for the lifetime of this object
    JellyfishHelper                 *jfh;
    hash_t                          *hash;		// Jellyfish hash
    ThreadedSparseMatrix<uint64_t>  *contamination_mx;  // Stores cumulative base count for each sequence where GC and CVG are binned
    uint32_t                        offset;
    uint16_t                        recordsInBatch;

    // Variables that are refreshed for each batch
    StringSet<CharString>           *names;
    StringSet<Dna5String>           *seqs;
    vector<vector<uint64_t>*>       *counts;    // Kmer counts for each kmer window in sequence (in same order as seqs and names; built by this class)
    vector<float>                   *coverages; // Overall coverage calculated for each sequence from the kmer windows.
    vector<float>                   *gcs;       // GC% for each sequence
    vector<uint32_t>                *lengths;   // Length in nucleotides for each sequence

public:
    Sect(SectArgs *_args) :
        args(_args),
        bucket_size(BATCH_SIZE / args->threads_arg),
        remaining(BATCH_SIZE % (bucket_size < 1 ? 1 : args->threads_arg))
    {
        // Setup handle to jellyfish hash
        jfh = new JellyfishHelper(args->db_arg);

        // Setup space for storing output
        offset = 0;
        recordsInBatch = 0;

        contamination_mx = new ThreadedSparseMatrix<uint64_t>(args->gc_bins, args->cvg_bins, args->threads_arg);
    }

    ~Sect()
    {
        if (jfh)
            delete jfh;

        jfh = NULL;


        if (contamination_mx)
            delete contamination_mx;

        contamination_mx = NULL;

        if(names)
            delete names;

        names = NULL;

        if(seqs)
            delete seqs;

        seqs = NULL;

        // Destroy any remaining batch variables
        destroyBatchVars();
    }


    void do_it()
    {
        // Setup output stream for jellyfish initialisation
        std::ostream* out_stream = args->verbose ? &cerr : (std::ostream*)0;

        // Load the jellyfish hash
        hash = jfh->loadHash(false, out_stream);

        // Setup stream to sequence file and check all is well
        SequenceStream seqStream(args->fasta_arg, SequenceStream::READ);
        if (!isGood(seqStream))
        {
            std::cerr << "ERROR: Could not open the sequence file: " << args->fasta_arg << endl;
            return;
        }


        // Setup output streams for files

        // Sequence kmer counts output stream
        std::ostringstream count_path;
        count_path << args->output_prefix << "_counts.cvg";
        ofstream_default count_path_stream(count_path.str().c_str(), cout);

        // Average sequence coverage and GC% scores output stream
        std::ostringstream cvg_gc_path;
        cvg_gc_path << args->output_prefix << "_stats.csv";
        ofstream_default cvg_gc_stream(cvg_gc_path.str().c_str(), cout);
        cvg_gc_stream << "seq_name coverage gc% seq_length" << endl;



        // Processes sequences in batches of records to reduce memory requirements
        while(!atEnd(seqStream))
        {
            // Create object to contain records to process for this batch
            names = new StringSet<CharString>();
            seqs = new StringSet<Dna5String>();

            // Read batch of records
            if (readBatch(*names, *seqs, seqStream, BATCH_SIZE) != 0)
            {
                std::cerr << "ERROR: Could not read batch of records! Offset: " << offset << endl;
                return;
            }

            // Record the number of records in this batch (may not be BATCH_SIZE) if we are at the end of the file
            recordsInBatch = length(*seqs);

            // Allocate memory for output produced by this batch
            createBatchVars(recordsInBatch);

            // Process batch with worker threads
            // Process each sequence is processed in a different thread.
            // In each thread lookup each kmer in the hash
            exec_join(args->threads_arg);

            // Output findings for this batch
            printCounts(count_path_stream);
            printStatTable(cvg_gc_stream);

            // Remove any batch specific variables from memory
            destroyBatchVars();

            // Increment batch management vars
            offset += recordsInBatch;
        }


        // Close output streams
        count_path_stream.close();
        cvg_gc_stream.close();

        // Merge the contamination matrix
        contamination_mx->mergeThreadedMatricies();

        // Send contamination matrix to file
        std::ostringstream contamination_mx_path;
        contamination_mx_path << args->output_prefix << "_contamination.mx";
        ofstream_default contamination_mx_stream(contamination_mx_path.str().c_str(), cout);
        printContaminationMatrix(contamination_mx_stream, args->fasta_arg);
        contamination_mx_stream.close();
    }

    void start(int th_id)
    {
        // Check to see if we have useful work to do for this thread, return if not
        if (bucket_size < 1 && th_id >= recordsInBatch)
        {
            return;
        }

        //processInBlocks(th_id);
        processInterlaced(th_id);
    }


    void printVars(std::ostream &out)
    {
        out << "SECT parameters:" << endl;
        out << " - Sequence File Path: " << args->fasta_arg << endl;
        out << " - Hash File Path: " << args->db_arg << endl;
        out << " - Hash: " << (hash ? "loaded" : "not loaded") << endl;
        out << " - Threads: " << args->threads_arg << endl;
        out << " - Bucket size: " << bucket_size << endl;
        out << " - Remaining: " << remaining << endl << endl;
    }





private:

    void destroyBatchVars()
    {
        if(names)
            delete names;

        names = NULL;

        if(seqs)
            delete seqs;

        seqs = NULL;


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
    }

    void createBatchVars(uint16_t batchSize)
    {
        counts = new vector<vector<uint64_t>*>(batchSize);
        coverages = new vector<float>(batchSize);
        gcs = new vector<float>(batchSize);
        lengths = new vector<uint32_t>(batchSize);
    }


    void printCounts(std::ostream &out)
    {
        for(int i = 0; i < recordsInBatch; i++)
        {
            out << ">" << toCString((*names)[i]) << endl;

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
        for(int i = 0; i < recordsInBatch; i++)
        {
            out << (*names)[i] << " " << (*coverages)[i] << " " << (*gcs)[i] << " " << (*lengths)[i] << endl;
        }
    }

    // Print kmer comparison matrix
    void printContaminationMatrix(std::ostream &out, const char* seqFile)
    {
        out << "# Sequence file processed: " << seqFile << endl;
        out << "# Jellyfish hash: " << args->db_arg << endl;
        out << "# Each column represents the GC%, with " << args->gc_bins << " bins. First bin is 0% and last bin in 100%" << endl;
        out << "# Each row represents the sequence coverage, with " << args->cvg_bins << " bins.  First bin in 0x last bin in " << endl;

        contamination_mx->getFinalMatrix()->printMatrix(out);
    }

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
            size_t rem_idx = (args->threads_arg * bucket_size) + th_id;
            processSeq(rem_idx);
        }
    }

    // This method is probably makes more efficient use of multiple cores on a length sorted fasta file
    void processInterlaced(uint_t th_id)
    {
        size_t start = th_id;
        size_t end = recordsInBatch;
        for(size_t i = start; i < end; i += args->threads_arg)
        {
            processSeq(i, th_id);
        }
    }

    void processSeq(const size_t index, const uint_t th_id)
    {
        uint_t kmer = hash->get_mer_len();

        // There's no substring functionality in SeqAn in this version (1.4.1).  So we'll just
        // use regular c++ string's for this bit.  The next version of SeqAn may offer substring
        // functionality, at which time I might change this code to make it run faster using
        // SeqAn's datastructures.
        stringstream ssSeq;
        ssSeq << (*seqs)[index];
        string seq = ssSeq.str();

        uint16_t seqLength = seq.length();
        uint_t nbCounts = seqLength - kmer + 1;
        float mean_cvg = 0;

        if (seqLength < kmer)
        {
            cerr << (*names)[index] << ": " << seq << " is too short to compute coverage.  Sequence length is "
                 << seqLength << " and kmer length is " << kmer << ". Setting sequence coverage to 0." << endl;
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
        (*lengths)[index] = seqLength;

        // Calc GC%
        uint32_t gs = 0;
        uint32_t cs = 0;
        uint32_t ns = 0;

        for(uint32_t i = 0; i < seqLength; i++)
        {
            char c = seq[i];

            if (c == 'G' || c == 'g')
                gs++;
            else if (c == 'C' || c == 'c')
                cs++;
            else if (c == 'N' || c == 'n')
                ns++;
        }

        float gc_perc = ((float)(gs + cs)) / ((float)(seqLength - ns));
        (*gcs)[index] = gc_perc;

        float log_cvg = args->cvg_logscale ? log10(mean_cvg) : mean_cvg;

        // Assume log_cvg 5 is max value
        float compressed_cvg = args->cvg_logscale ? log_cvg * (args->cvg_bins / 5.0) : mean_cvg * 0.1;

        uint16_t x = gc_perc * args->gc_bins;  // Convert float to 1.dp
        uint16_t y = compressed_cvg >= args->cvg_bins ? args->cvg_bins - 1 : compressed_cvg;      // Simply cap the y value

        // Add bases to matrix
        contamination_mx->getThreadMatrix(th_id)->inc(x, y, seqLength);

    }

};
