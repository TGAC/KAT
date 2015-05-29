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
#include <mutex>
using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::ostream;
using std::shared_ptr;
using std::make_shared;

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::path;

#include <jellyfish/large_hash_iterator.hpp>

#include <matrix/matrix_metadata_extractor.hpp>
#include <matrix/sparse_matrix.hpp>
#include <matrix/threaded_sparse_matrix.hpp>

#include <jellyfish_helper.hpp>
using kat::JellyfishHelper;


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
        path input1;
        path input2;
        path input3;
        path outputPrefix;
        double d1Scale;
        double d2Scale;
        uint16_t d1Bins;
        uint16_t d2Bins;
        uint16_t threads;
        uint16_t merLen;
        bool canonical1;
        bool canonical2;
        bool canonical3;
        uint64_t hashSize1;
        uint64_t hashSize2;
        uint64_t hashSize3;
        bool verbose;

        // Jellyfish mapped file hash vars
        shared_ptr<JellyfishHelper> jfh1;
        shared_ptr<JellyfishHelper> jfh2;
        shared_ptr<JellyfishHelper> jfh3;

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
        
        bool isCanonical1() const {
            return canonical1;
        }

        void setCanonical1(bool canonical1) {
            this->canonical1 = canonical1;
        }

        bool isCanonical2() const {
            return canonical2;
        }

        void setCanonical2(bool canonical2) {
            this->canonical2 = canonical2;
        }

        bool isCanonical3() const {
            return canonical3;
        }

        void setCanonical3(bool canonical3) {
            this->canonical3 = canonical3;
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

        uint64_t getHashSize1() const {
            return hashSize1;
        }

        void setHashSize1(uint64_t hashSize1) {
            this->hashSize1 = hashSize1;
        }

        uint64_t getHashSize2() const {
            return hashSize2;
        }

        void setHashSize2(uint64_t hashSize2) {
            this->hashSize2 = hashSize2;
        }

        uint64_t getHashSize3() const {
            return hashSize3;
        }

        void setHashSize3(uint64_t hashSize3) {
            this->hashSize3 = hashSize3;
        }

        path getInput1() const {
            return input1;
        }

        void setInput1(path input1) {
            this->input1 = input1;
        }

        path getInput2() const {
            return input2;
        }

        void setInput2(path input2) {
            this->input2 = input2;
        }

        path getInput3() const {
            return input3;
        }

        void setInput3(path input3) {
            this->input3 = input3;
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

        
        
        void startAndJoinThreads();
        
        void start(int th_id);

        // Scale counters to make the matrix look pretty

        uint64_t scaleCounter(uint64_t count, double scale_factor) {

            return count == 0 ? 0 : (uint64_t) ceil((double) count * scale_factor);
        }

        // Combines each threads matrix into a single matrix

        
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
