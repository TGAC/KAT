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
#include <thread>
using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::stringstream;
using std::shared_ptr;
using std::make_shared;
using std::thread;

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;
using boost::lexical_cast;

#include <jellyfish/mer_dna.hpp>
#include <jellyfish_helper.hpp>

#include <matrix/matrix_metadata_extractor.hpp>
#include <matrix/threaded_sparse_matrix.hpp>

typedef boost::error_info<struct SectError,string> SectErrorInfo;
struct SectException: virtual boost::exception, virtual std::exception { };

namespace kat {
    
    
    class Sect {
    private:

        static const uint16_t BATCH_SIZE = 1024;

        // Input args
        vector<path>    countsFiles;
        path            seqFile;
        path            outputPrefix;
        uint16_t        gcBins;
        uint16_t        cvgBins;
        bool            cvgLogscale;
        uint16_t        threads;
        bool            canonical;
        uint16_t        merLen;
        uint64_t        hashSize;
        bool            noCountStats;
        bool            median;
        bool            verbose;
            
        // Chunking vars
        size_t bucket_size, remaining; 

        // Variables that live for the lifetime of this object
        LargeHashArrayPtr hash;
        shared_ptr<ThreadedSparseMatrix> contamination_mx; // Stores cumulative base count for each sequence where GC and CVG are binned
        uint32_t offset;
        uint16_t recordsInBatch;
        path hashFile;

        // Variables that are refreshed for each batch
        seqan::StringSet<seqan::CharString> names;
        seqan::StringSet<seqan::Dna5String> seqs;
        shared_ptr<vector<shared_ptr<vector<uint64_t>>>> counts; // K-mer counts for each K-mer window in sequence (in same order as seqs and names; built by this class)
        shared_ptr<vector<double>> coverages; // Overall coverage calculated for each sequence from the K-mer windows.
        shared_ptr<vector<double>> gcs; // GC% for each sequence
        shared_ptr<vector<uint32_t>> lengths; // Length in nucleotides for each sequence

    public:

        Sect(const path _counts_file, const path _seq_file) : 
            Sect(vector<path>(), _seq_file) {            
            countsFiles.push_back(_counts_file);            
        }
        
        Sect(const vector<path> _counts_files, const path _seq_file) :
            countsFiles(_counts_files), seqFile(_seq_file) {
            outputPrefix = "kat-sect";
            gcBins = 1001;
            cvgBins = 1001;
            cvgLogscale = false;
            threads = 1;
            canonical = false;
            merLen = DEFAULT_MER_LEN;
            hashSize = DEFAULT_HASH_SIZE;
            noCountStats = false;
            median = true;
            verbose = false;
        }

        ~Sect() {
        }

        vector<path> getCountsFiles() const {
            return countsFiles;
        }

        void setCountsFiles(vector<path> countsFiles) {
            this->countsFiles = countsFiles;
        }

        
        path getOutputPrefix() const {
            return outputPrefix;
        }

        void setOutputPrefix(path outputPrefix) {
            this->outputPrefix = outputPrefix;
        }

        
        bool isCanonical() const {
            return canonical;
        }

        void setCanonical(bool canonical) {
            this->canonical = canonical;
        }

        uint16_t getCvgBins() const {
            return cvgBins;
        }

        void setCvgBins(uint16_t cvg_bins) {
            this->cvgBins = cvg_bins;
        }

        bool isCvgLogscale() const {
            return cvgLogscale;
        }

        void setCvgLogscale(bool cvg_logscale) {
            this->cvgLogscale = cvg_logscale;
        }

        uint16_t getGcBins() const {
            return gcBins;
        }

        void setGcBins(uint16_t gc_bins) {
            this->gcBins = gc_bins;
        }

        bool isMedian() const {
            return median;
        }

        void setMedian(bool median) {
            this->median = median;
        }

        bool isNoCountStats() const {
            return noCountStats;
        }

        void setNoCountStats(bool no_count_stats) {
            this->noCountStats = no_count_stats;
        }

        path getSeqFile() const {
            return seqFile;
        }

        void setSeqFile(path seq_file) {
            this->seqFile = seq_file;
        }

        uint16_t getThreads() const {
            return threads;
        }

        void setThreads(uint16_t threads) {
            this->threads = threads;
        }
        
        uint64_t getHashSize() const {
            return hashSize;
        }

        void setHashSize(uint64_t hashSize) {
            this->hashSize = hashSize;
        }

        uint16_t getMerLen() const {
            return merLen;
        }

        void setMerLen(uint16_t merLen) {
            this->merLen = merLen;
        }

        bool isVerbose() const {
            return verbose;
        }

        void setVerbose(bool verbose) {
            this->verbose = verbose;
        }

        
        void execute();


    private:

        void loadHashes();
        
        void startAndJoinThreads();
        
        void start(int th_id);
        
        void destroyBatchVars();

        void createBatchVars(uint16_t batchSize);

        void printCounts(std::ostream &out);

        void printStatTable(std::ostream &out);

        // Print K-mer comparison matrix

        void printContaminationMatrix(std::ostream &out, const path seqFile);

        // This method won't be optimal in most cases... Fasta files are normally sorted by length (largest first)
        // So first thread will be asked to do more work than the rest
        void processInBlocks(uint16_t th_id);
        
        // This method is probably makes more efficient use of multiple cores on a length sorted fasta file
        void processInterlaced(uint16_t th_id);

        void processSeq(const size_t index, const uint16_t th_id);
        
        static string helpMessage() {            
        
            return string(  "Usage: kat sect [options] <sequence_file> <counts_file>\n\n") +
                            "Estimates coverage levels for a collection of sequences using jellyfish K-mer counts.\n\n" \
                            "This tool will produce a fasta style file containing K-mer coverage counts mapped across each " \
                            "sequence.  In addition, a space separated table file containing the mean coverage score and GC " \
                            "of each sequence is produced.  The row order is identical to the original sequence file.\n\n" \
                            "Note: K-mers containing any Ns derived from sequences in the sequence file not be included.\n\n" \
                            "Options";

        }
        
    public:
        
        static int main(int argc, char *argv[]);
    };
}
