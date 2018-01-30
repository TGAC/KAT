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
#include <boost/exception/exception.hpp>
#include <boost/exception/info.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::path;
using boost::lexical_cast;

#include <jellyfish/mer_dna.hpp>

#include <kat/matrix_metadata_extractor.hpp>
#include <kat/jellyfish_helper.hpp>
#include <kat/input_handler.hpp>
#include <kat/sparse_matrix.hpp>
using kat::InputHandler;
using kat::ThreadedSparseMatrix;

typedef boost::error_info<struct SectError,string> SectErrorInfo;
struct SectException: virtual boost::exception, virtual std::exception { };

namespace kat {


    class Sect {
    private:

        static const uint16_t BATCH_SIZE = 1024;

        // Input args
        InputHandler    input;
        path            seqFile;
        path            outputPrefix;
        uint16_t        gcBins;
        uint16_t        cvgBins;
        bool            cvgLogscale;
        uint16_t        threads;
        bool            noCountStats;
        bool            outputGCStats;
        bool            extractNR;
        bool            extractR;
        uint32_t        minRepeat;
        uint32_t        maxRepeat;
        bool            verbose;

        // Chunking vars
        size_t bucket_size, remaining;

        // Variables that live for the lifetime of this object
        shared_ptr<ThreadedSparseMatrix> contamination_mx; // Stores cumulative base count for each sequence where GC and CVG are binned
        uint32_t offset;
        uint16_t recordsInBatch;
        path hashFile;

        // Variables that are refreshed for each batch
        seqan::StringSet<seqan::CharString> names;
        seqan::StringSet<seqan::CharString> seqs;
        shared_ptr<vector<shared_ptr<vector<uint64_t>>>> counts; // K-mer counts for each K-mer window in sequence (in same order as seqs and names; built by this class)
        shared_ptr<vector<shared_ptr<vector<int16_t>>>> gc_counts; // GC counts for each K-mer window in sequence (in same order as seqs and names; built by this class)
        shared_ptr<vector<uint32_t>> medians; // Overall coverage calculated for each sequence from the K-mer windows.
        shared_ptr<vector<double>> means; // Overall coverage calculated for each sequence from the K-mer windows.
        shared_ptr<vector<double>> gcs; // GC% for each sequence
        shared_ptr<vector<uint32_t>> lengths; // Length in nucleotides for each sequence
        shared_ptr<vector<uint32_t>> nonZero;
        shared_ptr<vector<double>> percentNonZero;
        shared_ptr<vector<uint32_t>> invalid;
        shared_ptr<vector<double>> percentInvalid;
        shared_ptr<vector<double>> percentNonZeroCorrected;


    public:

        Sect(const vector<path> _counts_files, const path _seq_file);

        virtual ~Sect() {
        }

        path getOutputPrefix() const {
            return outputPrefix;
        }

        void setOutputPrefix(path outputPrefix) {
            this->outputPrefix = outputPrefix;
        }

        void setTrim(const vector<uint16_t>& _5ptrim) {
            this->input.set5pTrim(_5ptrim);
        }

        bool isCanonical() const {
            return input.canonical;
        }

        void setCanonical(bool canonical) {
            this->input.canonical = canonical;
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

        bool isNoCountStats() const {
            return noCountStats;
        }

        void setNoCountStats(bool no_count_stats) {
            this->noCountStats = no_count_stats;
        }

        bool isOutputGCStats() const {
            return outputGCStats;
        }

        void setOutputGCStats(bool outputGCStats) {
            this->outputGCStats = outputGCStats;
        }

        bool isExtractNR() const {
            return extractNR;
        }

        void setExtractNR(bool extractNR) {
            this->extractNR = extractNR;
        }

        bool isExtractR() const {
            return extractR;
        }

        void setExtractR(bool extractR) {
            this->extractR = extractR;
        }

        uint32_t getMaxRepeat() const {
            return maxRepeat;
        }

        void setMaxRepeat(uint32_t maxRepeat) {
            this->maxRepeat = maxRepeat;
        }
        
        uint32_t getMinRepeat() const {
            return minRepeat;
        }

        void setMinRepeat(uint32_t minRepeat) {
            this->minRepeat = minRepeat;
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
            return input.hashSize;
        }

        void setHashSize(uint64_t hashSize) {
            this->input.hashSize = hashSize;
        }

        uint16_t getMerLen() const {
            return input.merLen;
        }

        void setMerLen(uint16_t merLen) {
            this->input.merLen = merLen;
        }

        bool isDumpHash() const {
            return input.dumpHash;
        }

        void setDumpHash(bool dumpHash) {
            this->input.dumpHash = dumpHash;
        }

        bool isVerbose() const {
            return verbose;
        }

        void setVerbose(bool verbose) {
            this->verbose = verbose;
        }


        void execute();

        void save();


    private:

        void processSeqFile();

        void analyseBatch();

        void analyseBatchSlice(int th_id);

        void merge();

        void destroyBatchVars();

        void createBatchVars(uint16_t batchSize);

        void printCounts(std::ostream &out);

        void printGCCounts(std::ostream &out);

        void printRegions(std::ostream &out, const uint32_t min_count, const uint32_t max_count);

        void printStatTable(std::ostream &out);

        // Print K-mer comparison matrix

        void printContaminationMatrix(std::ostream &out, const path seqFile);

        // This method won't be optimal in most cases... Fasta files are normally sorted by length (largest first)
        // So first thread will be asked to do more work than the rest
        void processInBlocks(uint16_t th_id);

        // This method is probably makes more efficient use of multiple cores on a length sorted fasta file
        void processInterlaced(uint16_t th_id);

        void processSeq(const size_t index, const uint16_t th_id);

        double gcCountToPercentage(int16_t count);

        static string helpMessage() {

            return string(  "Usage: kat sect [options] <sequence_file> (<input>)+\n\n") +
                            "Estimates coverage levels across sequences in the provided input sequence file.\n\n" \
                            "This tool will produce a fasta style representation of the input sequence file containing " \
                            "K-mer coverage counts mapped across each sequence.  K-mer coverage is determined from the " \
                            "provided counts input file, which can be either one jellyfish hash, or one or more FastA / " \
                            "FastQ files.  In addition, a space separated table file containing the mean coverage score and GC " \
                            "of each sequence is produced.  The row order is identical to the original sequence file.\n\n" \
                            "NOTE: K-mers containing any Ns derived from sequences in the sequence file not be included.\n\n" \
                            "WARNING: The <sequence_file> cannot be gzipped compressed.\n\n" \
                            "Options";

        }

    public:

        static int main(int argc, char *argv[]);
    };
}
