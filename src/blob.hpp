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


namespace kat {

    typedef boost::error_info<struct BlobError,string> BlobErrorInfo;
    struct BlobException: virtual boost::exception, virtual std::exception { };

    const string     DEFAULT_BLOB_PLOT_OUTPUT_TYPE     = "png";


    class Blob {
    private:

        static const uint16_t BATCH_SIZE = 1024;

        // Input args
        InputHandler    reads;
        InputHandler    assembly;
        path            outputPrefix;
        uint16_t        gcBins;
        uint16_t        cvgBins;
        uint16_t        threads;
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
        shared_ptr<vector<uint32_t>> medians; // Overall coverage calculated for each sequence from the K-mer windows.
        shared_ptr<vector<double>> means; // Overall coverage calculated for each sequence from the K-mer windows.
        shared_ptr<vector<uint32_t>> asmCns; // Overall coverage calculated for each sequence from the K-mer windows.
        shared_ptr<vector<double>> gcs; // GC% for each sequence
        shared_ptr<vector<uint32_t>> lengths; // Length in nucleotides for each sequence
        shared_ptr<vector<uint32_t>> nonZero;
        shared_ptr<vector<double>> percentNonZero;
        shared_ptr<vector<uint32_t>> invalid;
        shared_ptr<vector<double>> percentInvalid;
        shared_ptr<vector<double>> percentNonZeroCorrected;


    public:

        Blob(const vector<path> _reads_files, const path _asm_file);

        virtual ~Blob() {
        }

        path getOutputPrefix() const {
            return outputPrefix;
        }

        void setOutputPrefix(path outputPrefix) {
            this->outputPrefix = outputPrefix;
        }

        void setReadsTrim(const vector<uint16_t>& _5ptrim) {
            this->reads.set5pTrim(_5ptrim);
        }

        uint16_t getCvgBins() const {
            return cvgBins;
        }

        void setCvgBins(uint16_t cvg_bins) {
            this->cvgBins = cvg_bins;
        }

        uint16_t getGcBins() const {
            return gcBins;
        }

        void setGcBins(uint16_t gc_bins) {
            this->gcBins = gc_bins;
        }

        uint16_t getThreads() const {
            return threads;
        }

        void setThreads(uint16_t threads) {
            this->threads = threads;
        }

        uint64_t getHashSize() const {
            return reads.hashSize;
        }

        void setHashSize(uint64_t hashSize) {
            this->reads.hashSize = hashSize;
            this->assembly.hashSize = hashSize / 2;
        }

        uint16_t getMerLen() const {
            return reads.merLen;
        }

        void setMerLen(uint8_t merLen) {
            reads.merLen = merLen;
            assembly.merLen = merLen;
        }

        bool dumpHashes() const {
            return reads.dumpHash;
        }

        void setDumpHashes(bool dumpHashes) {
            this->reads.dumpHash = dumpHashes;
            this->assembly.dumpHash = dumpHashes;
        }

        bool hashGrowDisabled() const {
            return reads.disableHashGrow;
        }

        void setDisableHashGrow(bool disableHashGrow) {
            this->reads.disableHashGrow = disableHashGrow;
        }

        bool isVerbose() const {
            return verbose;
        }

        void setVerbose(bool verbose) {
            this->verbose = verbose;
        }


        void execute();

		void plot(const string& output_type);


    private:

        void processSeqFile();

        void analyseBatch();

        void analyseBatchSlice(int th_id);

        void destroyBatchVars();

        void createBatchVars(uint16_t batchSize);

        void printStatTable(std::ostream &out);


        // This method won't be optimal in most cases... Fasta files are normally sorted by length (largest first)
        // So first thread will be asked to do more work than the rest
        void processInBlocks(uint16_t th_id);

        // This method is probably makes more efficient use of multiple cores on a length sorted fasta file
        void processInterlaced(uint16_t th_id);

        void processSeq(const size_t index, const uint16_t th_id);

        double gcCountToPercentage(int16_t count);

        static string helpMessage() {

            return string(  "Usage: kat blob [options] <assembly> <reads>\n\n") +
                            "Calculated median read k-mer coverage, assembly k-mer coverage and GC% across each sequence in the provided assembly. " \
                            "The, assuming plotting is enabled, the results are converted into something something similar to a blobplot as " \
                            "would be produced by blobtools.  Each blob is colored according to a similar scheme used in spectra-cn plots.\n\n " \
                            "The <assembly> should be a fasta file that is NOT gzipped compressed.  The <reads> can be any number of <fasta/q> " \
                            "files, which CAN be gzipped compressed, or a pre-counted hash.\n\n" \
                            "Options";

        }

    public:

        static int main(int argc, char *argv[]);
    };
}
