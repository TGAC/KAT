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

#include <string.h>
#include <stdint.h>
#include <vector>
#include <memory>
#include <mutex>
using std::vector;
using std::string;
using std::shared_ptr;
using std::mutex;

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::path;

#include <jellyfish/large_hash_iterator.hpp>

#include <matrix/matrix_metadata_extractor.hpp>
#include <matrix/sparse_matrix.hpp>
#include <matrix/threaded_sparse_matrix.hpp>

#include "jellyfish_helper.hpp"
#include "input_handler.hpp"
using kat::JellyfishHelper;
using kat::InputHandler;

typedef boost::error_info<struct CompError,string> CompErrorInfo;
struct CompException: virtual boost::exception, virtual std::exception { };

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

        path hash1_path;
        path hash2_path;
        path hash3_path;

        CompCounters();
        
        CompCounters(const path& _hash1_path, const path& _hash2_path, const path& _hash3_path);
        
        CompCounters(const CompCounters& o);

        void updateHash1Counters(uint64_t hash1_count, uint64_t hash2_count);

        void updateHash2Counters(uint64_t hash1_count, uint64_t hash2_count);

        void updateHash3Counters(uint64_t hash3_count);

        void updateSharedCounters(uint64_t hash1_count, uint64_t hash2_count);

        void printCounts(ostream &out);
    };
    
    class ThreadedCompCounters {
    private: 
        uint16_t threads;

        shared_ptr<CompCounters> final_matrix;
        shared_ptr<vector<shared_ptr<CompCounters>>> threaded_counters;
        
    public:
        ThreadedCompCounters(const path& _hash1_path, const path& _hash2_path) :
           ThreadedCompCounters(_hash1_path, _hash2_path, path()) {}
           
        ThreadedCompCounters(const path& _hash1_path, const path& _hash2_path, const path& _hash3_path);                
        
        void printCounts(ostream &out);
        
        void add(shared_ptr<CompCounters> cc);
        
        size_t size() {
            return threaded_counters->size();
        }
        
        shared_ptr<CompCounters> getFinalMatrix() {
            return final_matrix;
        }
        
        shared_ptr<CompCounters> getThreadedMatrixAt(uint16_t index) {
            return threaded_counters->at(index);
        }
        
        void merge();
    
    };

    class Comp {
    private:
        
                
        // Args passed in
        vector<InputHandler> input;
        path outputPrefix;
        double d1Scale;
        double d2Scale;
        uint16_t d1Bins;
        uint16_t d2Bins;
        uint16_t threads;
        uint16_t merLen;
        bool parallelIO;
        bool verbose;

        // Threaded matrix data
        shared_ptr<ThreadedSparseMatrix> main_matrix;
        shared_ptr<ThreadedSparseMatrix> ends_matrix;
        shared_ptr<ThreadedSparseMatrix> middle_matrix;
        shared_ptr<ThreadedSparseMatrix> mixed_matrix;

        // Final data (created by merging thread results)
        shared_ptr<ThreadedCompCounters> comp_counters;
        
        std::mutex mu;


    public:

        Comp();
        
        Comp(const path& _input1, const path& _input2);
        
        Comp(const path& _input1, const path& _input2, const path& _input3);

        virtual ~Comp() {}
        
        bool doThirdHash() {
            return input.size() == 3;
        }
        
        bool isCanonical(uint16_t index) const {
            return input[index].canonical;
        }

        void setCanonical(uint16_t index, bool canonical) {
            this->input[index].canonical = canonical;
        }

        uint16_t getD1Bins() const {
            return d1Bins;
        }

        void setD1Bins(uint16_t d1Bins) {
            this->d1Bins = d1Bins;
        }

        double getD1Scale() const {
            return d1Scale;
        }

        void setD1Scale(double d1Scale) {
            this->d1Scale = d1Scale;
        }

        uint16_t getD2Bins() const {
            return d2Bins;
        }

        void setD2Bins(uint16_t d2Bins) {
            this->d2Bins = d2Bins;
        }

        double getD2Scale() const {
            return d2Scale;
        }

        void setD2Scale(double d2Scale) {
            this->d2Scale = d2Scale;
        }

        uint64_t getHashSize(uint16_t index) const {
            return input[index].hashSize;
        }

        void setHashSize(uint16_t index, uint64_t hashSize) {
            this->input[index].hashSize = hashSize;
        }

        path getInput(uint16_t index) const {
            return input[index].input[0];
        }

        void setInput(uint16_t index, path input) {
            this->input[index].input[0] = input;
        }

        uint8_t getMerLen() const {
            return merLen;
        }

        void setMerLen(uint8_t merLen) {
            this->merLen = merLen;
        }

        path getOutputPrefix() const {
            return outputPrefix;
        }

        void setOutputPrefix(path outputPrefix) {
            this->outputPrefix = outputPrefix;
        }

        uint16_t getThreads() const {
            return threads;
        }

        void setThreads(uint16_t threads) {
            this->threads = threads;
        }
        
        bool isParallelIO() const {
            return parallelIO;
        }

        void setParallelIO(bool parallelIO) {
            this->parallelIO = parallelIO;
        }

        bool isDumpHashes() const {
            return input[0].dumpHash;
        }

        void setDumpHashes(bool dumpHashes) {
            for(int i = 0; i < input.size(); i++) {
                this->input[i].dumpHash = dumpHashes;
            }
        }

        bool isVerbose() const {
            return verbose;
        }

        void setVerbose(bool verbose) {
            this->verbose = verbose;
        }

        
        void execute();
        

        
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

        
        // Print K-mer comparison matrix

        void printMainMatrix(ostream &out);

        // Print K-mer comparison matrix

        void printEndsMatrix(ostream &out);

        // Print K-mer comparison matrix

        void printMiddleMatrix(ostream &out);

        // Print K-mer comparison matrix

        void printMixedMatrix(ostream &out);

        // Print K-mer statistics

        void printCounters(ostream &out);


    private:

        void loadHashes();
        
        void compare();
        
        void compareSlice(int th_id);

        void merge();
        
        // Scale counters to make the matrix look pretty
        uint64_t scaleCounter(uint64_t count, double scale_factor) {

            return count == 0 ? 0 : (uint64_t) ceil((double) count * scale_factor);
        }
        
        static string helpMessage() {            
        
            return string(  "Usage: kat comp [options] <input_1> <input_2> [<input_3>]\n\n") +
                            "Compares jellyfish K-mer count hashes.\n\n" \
                            "The most common use case for this tool is to compare two (or three) K-mer hashes.  The typical use case for " \
                            "this tool is to compare K-mers from two K-mer hashes both representing K-mer counts for reads.  However, " \
                            "it is also common to compare K-mers generated from reads to those generated from an assembly.\n" \
                            "If comparing K-mers from reads to K-mers from an assembly, the larger (most likely the read) K-mer hash " \
                            "should be provided first, then the assembly K-mer hash second.\n" \
                            "The third optional jellyfish hash acts as a filter, restricting the analysis to the K-mers present on that " \
                            "set.  The manual contains more details on specific use cases.\n\n" \
                            "Options";

        }
        
    public:
        
        
        static int main(int argc, char *argv[]);

    };
}
