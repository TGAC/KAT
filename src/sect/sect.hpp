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
using std::stringstream;
using std::shared_ptr;
using std::make_shared;

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish_helper.hpp>

#include <matrix/matrix_metadata_extractor.hpp>
#include <matrix/threaded_sparse_matrix.hpp>

#include "sect_args.hpp"


namespace kat {
    const uint16_t BATCH_SIZE = 1024;

    class Sect {
    private:

        // Input args
        const SectArgs args;
        const size_t bucket_size, remaining; // Chunking vars

        // Variables that live for the lifetime of this object
        shared_ptr<JellyfishHelper> jfh;
        shared_ptr<ThreadedSparseMatrix> contamination_mx; // Stores cumulative base count for each sequence where GC and CVG are binned
        uint32_t offset;
        uint16_t recordsInBatch;

        // Variables that are refreshed for each batch
        seqan::StringSet<seqan::CharString> names;
        seqan::StringSet<seqan::Dna5String> seqs;
        vector<vector<uint64_t>*> *counts; // K-mer counts for each K-mer window in sequence (in same order as seqs and names; built by this class)
        vector<double> *coverages; // Overall coverage calculated for each sequence from the K-mer windows.
        vector<double> *gcs; // GC% for each sequence
        vector<uint32_t> *lengths; // Length in nucleotides for each sequence

        int resultCode;

    public:

        Sect(SectArgs& _args) :
        args(_args),
        bucket_size(BATCH_SIZE / args.threads_arg),
        remaining(BATCH_SIZE % (bucket_size < 1 ? 1 : args.threads_arg)) {

            // Setup handle to jellyfish hash
            jfh = make_shared<JellyfishHelper>(args.jellyfish_hash, AccessMethod::RANDOM);

            // Setup space for storing output
            offset = 0;
            recordsInBatch = 0;

            names = seqan::StringSet<seqan::CharString>();
            seqs = seqan::StringSet<seqan::Dna5String>();

            contamination_mx = make_shared<ThreadedSparseMatrix>(args.gc_bins, args.cvg_bins, args.threads_arg);

            resultCode = 0;
        }

        ~Sect() {
        }

        void execute() {
            // Setup output stream for jellyfish initialisation
            std::ostream* out_stream = args.verbose ? &cerr : (std::ostream*)0;

            // Open file, create RecordReader and check all is well
            std::fstream in(args.seq_file.c_str(), std::ios::in);
            seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);

            // Create the AutoSeqStreamFormat object and guess the file format.
            seqan::AutoSeqStreamFormat formatTag;
            if (!seqan::guessStreamFormat(reader, formatTag)) {
                std::cerr << "ERROR: Could not detect file format for: " << args.seq_file << endl;
                return;
            }

            // Setup output streams for files
            if (args.verbose)
                *out_stream << endl;

            // Sequence K-mer counts output stream
            ofstream_default* count_path_stream = NULL;
            if (!args.no_count_stats) {
                std::ostringstream count_path;
                count_path << args.output_prefix << "_counts.cvg";
                count_path_stream = new ofstream_default(count_path.str().c_str(), cout);
            }

            // Average sequence coverage and GC% scores output stream
            std::ostringstream cvg_gc_path;
            cvg_gc_path << args.output_prefix << "_stats.csv";
            ofstream_default cvg_gc_stream(cvg_gc_path.str().c_str(), cout);
            cvg_gc_stream << "seq_name coverage gc% seq_length" << endl;

            int res = 0;

            // Processes sequences in batches of records to reduce memory requirements
            while (!seqan::atEnd(reader) && res == 0) {
                if (args.verbose)
                    *out_stream << "Loading Batch of sequences... ";

                res = loadBatch(reader, formatTag, recordsInBatch);

                if (args.verbose)
                    *out_stream << "Loaded " << recordsInBatch << " records.  Processing batch... ";

                // Allocate memory for output produced by this batch
                createBatchVars(recordsInBatch);

                // Process batch with worker threads
                // Process each sequence is processed in a different thread.
                // In each thread lookup each K-mer in the hash
                //exec_join(args.threads_arg);

                // Output counts for this batch if (not not) requested
                if (!args.no_count_stats)
                    printCounts(*count_path_stream);

                // Output stats
                printStatTable(cvg_gc_stream);

                // Remove any batch specific variables from memory
                destroyBatchVars();

                // Increment batch management vars
                offset += recordsInBatch;

                if (args.verbose)
                    *out_stream << "done" << endl;
            }

            // Close output streams
            if (!args.no_count_stats) {
                count_path_stream->close();
                delete count_path_stream;
            }

            cvg_gc_stream.close();

            // Merge the contamination matrix
            contamination_mx->mergeThreadedMatricies();

            // Send contamination matrix to file
            std::ostringstream contamination_mx_path;
            contamination_mx_path << args.output_prefix << "_contamination.mx";
            ofstream_default contamination_mx_stream(contamination_mx_path.str().c_str(), cout);
            printContaminationMatrix(contamination_mx_stream, args.seq_file.c_str());
            contamination_mx_stream.close();

            // If there was a problem reading the data notify the user, otherwise output
            // the contamination matrix
            if (res != 0) {
                cerr << "ERROR: SECT could not analyse all sequences in the provided sequence file." << endl;
                resultCode = 1;
            }
        }

        void start(int th_id) {
            // Check to see if we have useful work to do for this thread, return if not
            if (bucket_size < 1 && th_id >= recordsInBatch) {
                return;
            }

            //processInBlocks(th_id);
            processInterlaced(th_id);
        }

        void printVars(std::ostream &out) {
            out << "SECT parameters:" << endl;
            out << " - Sequence File Path: " << args.seq_file << endl;
            out << " - Hash File Path: " << args.jellyfish_hash << endl;
            out << " - Threads: " << args.threads_arg << endl;
            out << " - Bucket size: " << bucket_size << endl;
            out << " - Remaining: " << remaining << endl << endl;
        }

        int getResultCode() {
            return resultCode;
        }



    private:

        void destroyBatchVars() {
            if (counts != NULL) {
                for (uint16_t i = 0; i < counts->size(); i++) {
                    vector<uint64_t>* ci = (*counts)[i];
                    if (ci != NULL)
                        delete ci;
                    ci = NULL;
                }
                delete counts;
                counts = NULL;
            }

            if (coverages != NULL)
                delete coverages;

            coverages = NULL;

            if (gcs != NULL)
                delete gcs;

            gcs = NULL;

            if (lengths != NULL)
                delete lengths;

            lengths = NULL;
        }

        void createBatchVars(uint16_t batchSize) {
            counts = new vector<vector<uint64_t>*>(batchSize);
            coverages = new vector<double>(batchSize);
            gcs = new vector<double>(batchSize);
            lengths = new vector<uint32_t>(batchSize);
        }

        uint16_t loadBatch(seqan::RecordReader<std::fstream, seqan::SinglePass<> >& reader, seqan::AutoSeqStreamFormat& formatTag, uint16_t& recordsRead) {
            // Create object to contain records to process for this batch
            seqan::clear(names);
            seqan::clear(seqs);

            int res = 0;
            uint32_t recordIndex = offset;

            // Read in a batch (Don't use readBatch... this seems to be broken!)
            // Room for improvement here... we could load the next batch while processing the previous one.
            // This would require double buffering.
            for (unsigned i = 0; (res == 0) && (i < BATCH_SIZE) && !seqan::atEnd(reader); ++i) {
                seqan::CharString id;
                seqan::CharString seq;

                res = seqan::readRecord(id, seq, reader, formatTag);
                if (res == 0) {
                    seqan::appendValue(names, id);

                    // Hopefully auto converts chars not in {A,T,G,C,N} to N
                    seqan::Dna5String dnaSeq = seq;

                    seqan::appendValue(seqs, dnaSeq);

                    recordIndex++;
                } else {
                    cerr << endl << "ERROR: SECT cannot finish processing all records in file.  SECT encountered an error reading file at record: " << recordIndex
                            << "; Error code: " << res << "; Last sequence ID: " << id << "; Continuing to process currently loaded records." << endl;
                }
            }

            // Record the number of records in this batch (may not be BATCH_SIZE) if we are at the end of the file
            recordsRead = seqan::length(names);

            return res;
        }

        void printCounts(std::ostream &out) {
            for (int i = 0; i < recordsInBatch; i++) {
                out << ">" << seqan::toCString(names[i]) << endl;

                vector<uint64_t>* seqCounts = (*counts)[i];

                if (seqCounts != NULL && !seqCounts->empty()) {
                    out << (*seqCounts)[0];

                    for (size_t j = 1; j < seqCounts->size(); j++) {
                        out << " " << (*seqCounts)[j];
                    }

                    out << endl;
                } else {
                    out << "0" << endl;
                }
            }
        }

        void printStatTable(std::ostream &out) {
            for (int i = 0; i < recordsInBatch; i++) {
                out << names[i] << " " << (*coverages)[i] << " " << (*gcs)[i] << " " << (*lengths)[i] << endl;
            }
        }

        // Print K-mer comparison matrix

        void printContaminationMatrix(std::ostream &out, const char* seqFile) {
            SM64 mx = contamination_mx->getFinalMatrix();

            out << mme::KEY_TITLE << "Contamination Plot for " << args.seq_file << " and " << args.jellyfish_hash << endl;
            out << mme::KEY_X_LABEL << "GC%" << endl;
            out << mme::KEY_Y_LABEL << "Average K-mer Coverage" << endl;
            out << mme::KEY_Z_LABEL << "Base Count per bin" << endl;
            out << mme::KEY_NB_COLUMNS << args.gc_bins << endl;
            out << mme::KEY_NB_ROWS << args.cvg_bins << endl;
            out << mme::KEY_MAX_VAL << mx->getMaxVal() << endl;
            out << mme::KEY_TRANSPOSE << "0" << endl;
            out << mme::MX_META_END << endl;

            contamination_mx->getFinalMatrix()->printMatrix(out);
        }

        // This method won't be optimal in most cases... Fasta files are normally sorted by length (largest first)
        // So first thread will be asked to do more work than the rest

        void processInBlocks(uint16_t th_id) {
            size_t start = bucket_size < 1 ? th_id : th_id * bucket_size;
            size_t end = bucket_size < 1 ? th_id : start + bucket_size - 1;
            for (size_t i = start; i <= end; i++) {
                processSeq(i, th_id);
            }

            // Process a remainder if required
            if (th_id < remaining) {
                size_t rem_idx = (args.threads_arg * bucket_size) + th_id;
                processSeq(rem_idx, th_id);
            }
        }

        // This method is probably makes more efficient use of multiple cores on a length sorted fasta file

        void processInterlaced(uint16_t th_id) {
            size_t start = th_id;
            size_t end = recordsInBatch;
            for (size_t i = start; i < end; i += args.threads_arg) {
                processSeq(i, th_id);
            }
        }

        void processSeq(const size_t index, const uint16_t th_id) {

            unsigned int kmer = jfh->getKeyLen();

            // There's no substring functionality in SeqAn in this version (1.4.1).  So we'll just
            // use regular c++ string's for this bit.  The next version of SeqAn may offer substring
            // functionality, at which time I might change this code to make it run faster using
            // SeqAn's datastructures.
            stringstream ssSeq;
            ssSeq << seqs[index];
            string seq = ssSeq.str();

            uint64_t seqLength = seq.length();
            uint64_t nbCounts = seqLength - kmer + 1;
            double average_cvg = 0.0;

            if (seqLength < kmer) {

                cerr << names[index] << ": " << seq << " is too short to compute coverage.  Sequence length is "
                        << seqLength << " and K-mer length is " << kmer << ". Setting sequence coverage to 0." << endl;
            } else {

                vector<uint64_t>* seqCounts = new vector<uint64_t>(nbCounts);

                uint64_t sum = 0;

                for (uint64_t i = 0; i < nbCounts; i++) {

                    string merstr = seq.substr(i, kmer);

                    // Jellyfish compacted hash does not support Ns so if we find one set this mer count to 0
                    if (merstr.find("N") != string::npos) {
                        (*seqCounts)[i] = 0;
                    } else {
                        mer_dna mer(merstr.c_str());
                        uint64_t count = jfh->getCount(mer);
                        sum += count;

                        (*seqCounts)[i] = count;
                    }
                }

                (*counts)[index] = seqCounts;

                if (args.median) {
                    
                    // Create a copy of the counts, and sort it first, then take median value
                    vector<uint64_t> sortedSeqCounts = *seqCounts;                    
                    std::sort(sortedSeqCounts.begin(), sortedSeqCounts.end());
                    average_cvg = (double)(sortedSeqCounts[sortedSeqCounts.size() / 2]);                    
                }
                else {
                    // Calculate the mean
                    average_cvg = (double)sum / (double)nbCounts;                    
                }
                
                (*coverages)[index] = average_cvg;

            }

            // Add length
            (*lengths)[index] = seqLength;

            // Calc GC%
            uint64_t gs = 0;
            uint64_t cs = 0;
            uint64_t ns = 0;

            for (uint64_t i = 0; i < seqLength; i++) {
                char c = seq[i];

                if (c == 'G' || c == 'g')
                    gs++;
                else if (c == 'C' || c == 'c')
                    cs++;
                else if (c == 'N' || c == 'n')
                    ns++;
            }

            double gc_perc = ((double) (gs + cs)) / ((double) (seqLength - ns));
            (*gcs)[index] = gc_perc;

            double log_cvg = args.cvg_logscale ? log10(average_cvg) : average_cvg;

            // Assume log_cvg 5 is max value
            double compressed_cvg = args.cvg_logscale ? log_cvg * (args.cvg_bins / 5.0) : average_cvg * 0.1;

            uint16_t x = gc_perc * args.gc_bins; // Convert double to 1.dp
            uint16_t y = compressed_cvg >= args.cvg_bins ? args.cvg_bins - 1 : compressed_cvg; // Simply cap the y value

            // Add bases to matrix
            contamination_mx->getThreadMatrix(th_id)->inc(x, y, seqLength);
        }
    };
}
