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
using std::ostream;
using std::shared_ptr;
using std::make_shared;
using std::thread;

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
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

        const char* hash1_path;
        const char* hash2_path;
        const char* hash3_path;

        CompCounters(const char* _hash1_path, const char* _hash2_path, const char* _hash3_path) :
        hash1_path(_hash1_path), hash2_path(_hash2_path), hash3_path(_hash3_path) {

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

        void updateHash1Counters(uint64_t hash1_count, uint64_t hash2_count) {
            hash1_total += hash1_count;
            hash1_distinct++;

            if (!hash2_count) {
                hash1_only_total += hash1_count;
                hash1_only_distinct++;
            }
        }

        void updateHash2Counters(uint64_t hash1_count, uint64_t hash2_count) {
            hash2_total += hash2_count;
            hash2_distinct++;

            if (!hash1_count) {
                hash2_only_total += hash2_count;
                hash2_only_distinct++;
            }
        }

        void updateHash3Counters(uint64_t hash3_count) {

            hash3_total += hash3_count;
            hash3_distinct++;
        }

        void updateSharedCounters(uint64_t hash1_count, uint64_t hash2_count) {

            if (hash1_count && hash2_count) {
                shared_hash1_total += hash1_count;
                shared_hash2_total += hash2_count;
                shared_distinct++;
            }
        }

        void printCounts(ostream &out) {

            out << "K-mer statistics for: " << endl;
            out << " - Hash 1: " << hash1_path << endl;
            out << " - Hash 2: " << hash2_path << endl;

            if (hash3_total > 0)
                out << " - Hash 3: " << hash3_path << endl;

            out << endl;

            out << "Total K-mers in: " << endl;
            out << " - Hash 1: " << hash1_total << endl;
            out << " - Hash 2: " << hash2_total << endl;

            if (hash3_total > 0)
                out << " - Hash 3: " << hash3_total << endl;

            out << endl;

            out << "Distinct K-mers in:" << endl;
            out << " - Hash 1: " << hash1_distinct << endl;
            out << " - Hash 2: " << hash2_distinct << endl;
            if (hash3_total > 0)
                out << " - Hash 3: " << hash3_distinct << endl;

            out << endl;

            out << "Total K-mers only found in:" << endl;
            out << " - Hash 1: " << hash1_only_total << endl;
            out << " - Hash 2: " << hash2_only_total << endl;
            out << endl;

            out << "Distinct K-mers only found in:" << endl;
            out << " - Hash 1: " << hash1_only_distinct << endl;
            out << " - Hash 2: " << hash2_only_distinct << endl << endl;

            out << "Shared K-mers:" << endl;
            out << " - Total shared found in hash 1: " << shared_hash1_total << endl;
            out << " - Total shared found in hash 2: " << shared_hash2_total << endl;
            out << " - Distinct shared K-mers: " << shared_distinct << endl << endl;
        }
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
        shared_ptr<CompCounters> final_comp_counters;

        // Thread specific data
        shared_ptr<vector<shared_ptr<CompCounters> > > thread_comp_counters;


    public:

        Comp(path _input1, path _input2) : 
            Comp(_input1, _input2, path()) {}
        
        Comp(path _input1, path _input2, path _input3) : 
            input1(_input1), input2(_input2), input3(_input3) {
            outputPrefix = "kat-comp";
            d1Scale = 1.0;
            d2Scale = 1.0;
            d1Bins = 1001;
            d2Bins = 1001;
            threads = 1;
            merLen = DEFAULT_MER_LEN;
            canonical1 = false;
            canonical2 = false;
            canonical3 = false;
            hashSize1 = DEFAULT_HASH_SIZE;
            hashSize2 = DEFAULT_HASH_SIZE;
            hashSize3 = DEFAULT_HASH_SIZE;
            verbose = false;
        }

        ~Comp() {
        }

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

        
        void execute() {

            // Check input file exists
            if (!bfs::exists(input1) && !bfs::symbolic_link_exists(input1)) {
                BOOST_THROW_EXCEPTION(CompException() << CompErrorInfo(string(
                        "Could not find first input file at: ") + input1.string() + "; please check the path and try again."));
            }

            // Check input file exists
            if (!bfs::exists(input2) && !bfs::symbolic_link_exists(input2)) {
                BOOST_THROW_EXCEPTION(CompException() << CompErrorInfo(string(
                        "Could not find second jellyfish hash file at: ") + input2.string() + "; please check the path and try again."));
            }

            // Check input file exists
            if (!input3.empty() && !bfs::exists(input3) && !bfs::symbolic_link_exists(input3)) {
                BOOST_THROW_EXCEPTION(CompException() << CompErrorInfo(string(
                        "Could not find third jellyfish hash file at: ") + input3.string() + "; please check the path and try again."));
            }
            
            // Create the final K-mer counter matrices
            main_matrix = shared_ptr<ThreadedSparseMatrix>(
                    new ThreadedSparseMatrix(d1Bins, d2Bins, threads));

            // Initialise extra matrices for hash3 (only allocates space if required)
            if (!input3.empty()) {
                
                ends_matrix = make_shared<ThreadedSparseMatrix>(d1Bins, d2Bins, threads);
                middle_matrix = make_shared<ThreadedSparseMatrix>(d1Bins, d2Bins, threads);
                mixed_matrix = make_shared<ThreadedSparseMatrix>(d1Bins, d2Bins, threads);
            }

            // Create the final comp counters
            final_comp_counters = make_shared<CompCounters>(
                    input1.c_str(), input2.c_str(), input3.c_str());

            // Create the comp counters for each thread
            thread_comp_counters = make_shared<vector<shared_ptr<CompCounters>>>(threads);
            for (int i = 0; i < threads; i++) {
                thread_comp_counters->push_back(make_shared<CompCounters>(
                        input1.c_str(), input2.c_str(), input3.c_str()));
            }

            std::ostream* out_stream = verbose ? &cerr : (std::ostream*)0;

            // Count kmers if necessary
            path i1Ext = input1.extension();
            path i2Ext = input2.extension();
            path i3Ext = !input3.empty() ? input3.extension() : path();
            
            path i1 = input1;
            path i2 = input2;
            path i3 = input3;
            
            string merLenStr = lexical_cast<string>(merLen);
            
            if (JellyfishHelper::isSequenceFile(i1Ext)) {
                
                cout << "Input 1 is a sequence file.  Executing jellyfish to count kmers." << endl;
                
                i1 = path(outputPrefix.string() + string(".jf") + merLenStr);
                
                JellyfishHelper::jellyfishCount(input1, i1, merLen, hashSize1, threads, canonical1, true);
            }
            
            if (JellyfishHelper::isSequenceFile(i2Ext)) {
                
                cout << "Input 2 is a sequence file.  Executing jellyfish to count kmers." << endl;
                
                i2 += path(outputPrefix.string() + string(".jf") + merLenStr);
                
                JellyfishHelper::jellyfishCount(input1, i2, merLen, hashSize2, threads, canonical2, true);
            }
            
            if (!input3.empty() && JellyfishHelper::isSequenceFile(i3Ext)) {
                
                cout << "Input 3 is a sequence file.  Executing jellyfish to count kmers." << endl;
                
                i3 += path(outputPrefix.string() + string(".jf") + merLenStr);
                
                JellyfishHelper::jellyfishCount(input1, i3, merLen, hashSize3, threads, canonical3, true);
            }
            
            
            // Load the hashes
            jfh1 = make_shared<JellyfishHelper>(input1, AccessMethod::RANDOM);
            jfh2 = make_shared<JellyfishHelper>(input2, AccessMethod::RANDOM);
            jfh3 = !(input3.empty()) ?
                    make_shared<JellyfishHelper>(input3, AccessMethod::RANDOM) :
                    nullptr;

            // Check K-mer lengths are the same for both hashes.  We can't continue if they are not.
            if (jfh1->getKeyLen() != jfh2->getKeyLen()) {
                BOOST_THROW_EXCEPTION(CompException() << CompErrorInfo(string(
                        "Cannot process hashes that were created with different K-mer lengths.  Hash1: ") +
                        lexical_cast<string>(jfh1->getKeyLen()) + 
                        ".  Hash2: " + lexical_cast<string>(jfh2->getKeyLen()) + "."));
            } else if (jfh3 && (jfh3->getKeyLen() != jfh1->getKeyLen())) {
                BOOST_THROW_EXCEPTION(CompException() << CompErrorInfo(string(
                        "Cannot process hashes that were created with different K-mer lengths.  Hash1: ") +
                        lexical_cast<string>(jfh1->getKeyLen()) + 
                        ".  Hash3: " + lexical_cast<string>(jfh3->getKeyLen()) + "."));
            }

            if (verbose)
                cerr << endl
                    << "All hashes loaded successfully." << endl
                    << "Starting threads...";

            // Run the threads
            startAndJoinThreads();

            if (verbose)
                cerr << "done." << endl
                    << "Merging results...";

            // Merge results from the threads
            merge();

            if (verbose)
                cerr << "done." << endl;
        }
        

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

        // Final data (created by merging thread results)

        shared_ptr<CompCounters> getCompCounters() {

            return final_comp_counters;
        }


        // Print K-mer comparison matrix

        void printMainMatrix(ostream &out) {

            SM64 mx = main_matrix->getFinalMatrix();

            out << mme::KEY_TITLE << "K-mer comparison plot" << endl
                    << mme::KEY_X_LABEL << "K-mer multiplicity for: " << input1 << endl
                    << mme::KEY_Y_LABEL << "K-mer multiplicity for: " << input2 << endl
                    << mme::KEY_Z_LABEL << "Distinct K-mers per bin" << endl
                    << mme::KEY_NB_COLUMNS << mx->height() << endl
                    << mme::KEY_NB_ROWS << mx->width() << endl
                    << mme::KEY_MAX_VAL << mx->getMaxVal() << endl
                    << mme::KEY_TRANSPOSE << "1" << endl
                    << mme::MX_META_END << endl;

            mx->printMatrix(out);
        }

        // Print K-mer comparison matrix

        void printEndsMatrix(ostream &out) {

            out << "# Each row represents K-mer multiplicity for: " << input1 << endl;
            out << "# Each column represents K-mer multiplicity for sequence ends: " << input3 << endl;

            ends_matrix->getFinalMatrix()->printMatrix(out);
        }

        // Print K-mer comparison matrix

        void printMiddleMatrix(ostream &out) {

            out << "# Each row represents K-mer multiplicity for: " << input1 << endl;
            out << "# Each column represents K-mer multiplicity for sequence middles: " << input2 << endl;

            middle_matrix->getFinalMatrix()->printMatrix(out);
        }

        // Print K-mer comparison matrix

        void printMixedMatrix(ostream &out) {

            out << "# Each row represents K-mer multiplicity for hash file 1: " << input1 << endl;
            out << "# Each column represents K-mer multiplicity for mixed: " << input2 << " and " << input3 << endl;

            mixed_matrix->getFinalMatrix()->printMatrix(out);
        }

        // Print K-mer statistics

        void printCounters(ostream &out) {

            final_comp_counters->printCounts(out);
        }



    private:

        
        
        void startAndJoinThreads() {
            
            thread t[threads];
            
            for(int i = 0; i < threads; i++) {
                t[i] = thread(&Comp::start, this, i);
            }
            
            for(int i = 0; i < threads; i++){
                t[i].join();
            }
        }
        
        void start(int th_id) {

            // Get handle to sparse matrix for this thread
            shared_ptr<SparseMatrix<uint64_t> > thread_matrix = main_matrix->getThreadMatrix(th_id);

            // Get handle on this thread's comp counter
            shared_ptr<CompCounters> comp_counters = (*thread_comp_counters)[th_id];

            // Setup iterator for this thread's chunk of hash1
            lha::region_iterator hash1Iterator = jfh1->getSlice(th_id, threads);

            // Go through this thread's slice for hash1
            while (hash1Iterator.next()) {
                // Get the current K-mer count for hash1
                uint64_t hash1_count = hash1Iterator.val();

                // Get the count for this K-mer in hash2 (assuming it exists... 0 if not)
                uint64_t hash2_count = jfh1->getCount(hash1Iterator.key());

                // Get the count for this K-mer in hash3 (assuming it exists... 0 if not)
                uint64_t hash3_count = jfh3 != nullptr ? jfh3->getCount(hash1Iterator.key()) : 0;

                // Increment hash1's unique counters
                comp_counters->updateHash1Counters(hash1_count, hash2_count);

                // Increment shared counters
                comp_counters->updateSharedCounters(hash1_count, hash2_count);

                // Scale counters to make the matrix look pretty
                uint64_t scaled_hash1_count = scaleCounter(hash1_count, d1Scale);
                uint64_t scaled_hash2_count = scaleCounter(hash2_count, d2Scale);
                uint64_t scaled_hash3_count = scaleCounter(hash3_count, d2Scale);

                // Modifies hash counts so that K-mer counts larger than MATRIX_SIZE are dumped in the last slot
                if (scaled_hash1_count >= d1Bins) scaled_hash1_count = d1Bins - 1;
                if (scaled_hash2_count >= d2Bins) scaled_hash2_count = d2Bins - 1;
                if (scaled_hash3_count >= d2Bins) scaled_hash3_count = d2Bins - 1;

                // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
                thread_matrix->inc(scaled_hash1_count, scaled_hash2_count, 1);

                // Update hash 3 related matricies if hash 3 was provided
                if (jfh3 != nullptr) {
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
            lha::region_iterator hash2Iterator = jfh2->getSlice(th_id, threads);

            // Iterate through this thread's slice of hash2
            while (hash2Iterator.next()) {
                // Get the current K-mer count for hash2
                uint64_t hash2_count = hash2Iterator.val();

                // Get the count for this K-mer in hash1 (assuming it exists... 0 if not)
                uint64_t hash1_count = jfh1->getCount(hash2Iterator.key());

                // Increment hash2's unique counters (don't bother with shared counters... we've already done this)
                comp_counters->updateHash2Counters(hash1_count, hash2_count);

                // Only bother updating thread matrix with K-mers not found in hash1 (we've already done the rest)
                if (!hash1_count) {
                    // Scale counters to make the matrix look pretty
                    uint64_t scaled_hash2_count = scaleCounter(hash2_count, d2Scale);

                    // Modifies hash counts so that K-mer counts larger than MATRIX_SIZE are dumped in the last slot
                    if (scaled_hash2_count >= d2Bins) scaled_hash2_count = d2Bins - 1;

                    // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
                    thread_matrix->inc(0, scaled_hash2_count, 1);
                }
            }

            // Only update hash3 counters if hash3 was provided
            if (jfh3 != nullptr) {
                // Setup iterator for this thread's chunk of hash3
                lha::region_iterator hash3Iterator = jfh3->getSlice(th_id, threads);

                // Iterate through this thread's slice of hash2
                while (hash3Iterator.next()) {
                    // Get the current K-mer count for hash2

                    uint64_t hash3_count = hash3Iterator.val();

                    // Increment hash3's unique counters (don't bother with shared counters... we've already done this)
                    comp_counters->updateHash3Counters(hash3_count);
                }
            }
        }

        // Scale counters to make the matrix look pretty

        uint64_t scaleCounter(uint64_t count, double scale_factor) {

            return count == 0 ? 0 : (uint64_t) ceil((double) count * scale_factor);
        }

        // Combines each threads matrix into a single matrix

        void merge() {

            main_matrix->mergeThreadedMatricies();

            if (jfh3 != nullptr) {
                ends_matrix->mergeThreadedMatricies();
                middle_matrix->mergeThreadedMatricies();
                mixed_matrix->mergeThreadedMatricies();
            }

            // Merge counters
            for (int k = 0; k < threads; k++) {
                shared_ptr<CompCounters> thread_comp_counter = (*thread_comp_counters)[k];

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
        
        
        static int main(int argc, char *argv[]) {
            
            path input1;
            path input2;
            path input3;
            path output_prefix;
            double d1_scale;
            double d2_scale;
            uint16_t d1_bins;
            uint16_t d2_bins;
            uint16_t threads;
            uint16_t mer_len;
            bool canonical_1;
            bool canonical_2;
            bool canonical_3;
            uint64_t hash_size_1;
            uint64_t hash_size_2;
            uint64_t hash_size_3;
            bool verbose;
            bool help;
        
            // Declare the supported options.
            po::options_description generic_options(Comp::helpMessage(), 100);
            generic_options.add_options()
                    ("output_prefix,o", po::value<path>(&output_prefix)->default_value("kat-comp"), 
                        "Path prefix for files generated by this program.")
                    ("threads,t", po::value<uint16_t>(&threads)->default_value(1),
                        "The number of threads to use")
                    ("d1_scale,x", po::value<double>(&d1_scale)->default_value(1.0),
                        "Scaling factor for the first dataset - float multiplier")
                    ("d2_scale,y", po::value<double>(&d2_scale)->default_value(1.0),
                        "Scaling factor for the second dataset - float multiplier")
                    ("d1_bins,i", po::value<uint16_t>(&d1_bins)->default_value(1001),
                        "Number of bins for the first dataset.  i.e. number of rows in the matrix")
                    ("d2_bins,j", po::value<uint16_t>(&d2_bins)->default_value(1001),
                        "Number of bins for the second dataset.  i.e. number of rows in the matrix")
                    ("canonical1,c1", po::bool_switch(&canonical_1)->default_value(false),
                        "Whether the jellyfish hash for input 1 contains K-mers produced for both strands.  If this is not set to the same value as was produced during jellyfish counting then output from sect will be unpredicatable.")
                    ("canonical2,c2", po::bool_switch(&canonical_2)->default_value(false),
                        "Whether the jellyfish hash for input 2 contains K-mers produced for both strands.  If this is not set to the same value as was produced during jellyfish counting then output from sect will be unpredicatable.")
                    ("canonical3,c3", po::bool_switch(&canonical_3)->default_value(false),
                        "Whether the jellyfish hash for input 3 contains K-mers produced for both strands.  If this is not set to the same value as was produced during jellyfish counting then output from sect will be unpredicatable.  Only applicable if you are using a third input.")
                    ("mer_len,m", po::value<uint16_t>(&mer_len)->default_value(DEFAULT_MER_LEN),
                        "The kmer length to use in the kmer hashes.  Larger values will provide more discriminating power between kmers but at the expense of additional memory and lower coverage.")
                    ("hash_size_1,h1", po::value<uint64_t>(&hash_size_1)->default_value(DEFAULT_HASH_SIZE),
                        "If kmer counting is required for input 1, then use this value as the hash size.  It is important this is larger than the number of distinct kmers in your set.  We do not try to merge kmer hashes in this version of KAT.")
                    ("hash_size_2,h2", po::value<uint64_t>(&hash_size_2)->default_value(DEFAULT_HASH_SIZE),
                        "If kmer counting is required for input 2, then use this value as the hash size.  It is important this is larger than the number of distinct kmers in your set.  We do not try to merge kmer hashes in this version of KAT.")
                    ("hash_size_3,h3", po::value<uint64_t>(&hash_size_3)->default_value(DEFAULT_HASH_SIZE),
                        "If kmer counting is required for input 3, then use this value as the hash size.  It is important this is larger than the number of distinct kmers in your set.  We do not try to merge kmer hashes in this version of KAT.")
                    ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                        "Print extra information.")
                    ("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
                    ;

            // Hidden options, will be allowed both on command line and
            // in config file, but will not be shown to the user.
            po::options_description hidden_options("Hidden options");
            hidden_options.add_options()
                    ("input_1,i1", po::value<path>(&input1), "Path to the first input file.  Can be either FastA, FastQ or a jellyfish hash (non bloom filtered)")
                    ("input_2,i2", po::value<path>(&input2), "Path to the second input file.  Can be either FastA, FastQ or a jellyfish hash (non bloom filtered)")
                    ("input_3,i3", po::value<path>(&input3), "Path to the third input file.  Can be either FastA, FastQ or a jellyfish hash (non bloom filtered)")
                    ;

            // Positional option for the input bam file
            po::positional_options_description p;
            p.add("input_1", 1);
            p.add("input_2", 1);
            p.add("input_3", 1);
            

            // Combine non-positional options
            po::options_description cmdline_options;
            cmdline_options.add(generic_options).add(hidden_options);

            // Parse command line
            po::variables_map vm;
            po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
            po::notify(vm);

            // Output help information the exit if requested
            if (help || argc <= 1) {
                cout << generic_options << endl;
                return 1;
            }
        
        

            auto_cpu_timer timer(1, "KAT COMP completed.\nTotal runtime: %ws\n\n");        

            cout << "Running KAT in COMP mode" << endl
                 << "------------------------" << endl << endl;
        

            // Create the sequence coverage object
            Comp comp(input1, input2, input3);
            comp.setOutputPrefix(output_prefix);
            comp.setD1Scale(d1_scale);
            comp.setD2Scale(d2_scale);
            comp.setD1Bins(d1_bins);
            comp.setD2Bins(d2_bins);
            comp.setThreads(threads);
            comp.setMerLen(mer_len);
            comp.setCanonical1(canonical_1);
            comp.setCanonical2(canonical_2);
            comp.setCanonical3(canonical_3);
            comp.setHashSize1(hash_size_1);
            comp.setHashSize2(hash_size_2);
            comp.setHashSize3(hash_size_3);
            comp.setVerbose(verbose);
            
            // Do the work
            comp.execute();

            // Send main matrix to output file
            std::ostringstream main_mx_out_path;
            main_mx_out_path << output_prefix << "_main.mx";
            ofstream_default main_mx_out_stream(main_mx_out_path.str().c_str(), cout);
            comp.printMainMatrix(main_mx_out_stream);
            main_mx_out_stream.close();

            // Output ends matricies if required
            if (!(input3.empty())) {
                // Ends matrix
                std::ostringstream ends_mx_out_path;
                ends_mx_out_path << output_prefix << "_ends.mx";
                ofstream_default ends_mx_out_stream(ends_mx_out_path.str().c_str(), cout);
                comp.printEndsMatrix(ends_mx_out_stream);
                ends_mx_out_stream.close();

                // Middle matrix
                std::ostringstream middle_mx_out_path;
                middle_mx_out_path << output_prefix << "_middle.mx";
                ofstream_default middle_mx_out_stream(middle_mx_out_path.str().c_str(), cout);
                comp.printMiddleMatrix(middle_mx_out_stream);
                middle_mx_out_stream.close();

                // Mixed matrix
                std::ostringstream mixed_mx_out_path;
                mixed_mx_out_path << output_prefix << "_mixed.mx";
                ofstream_default mixed_mx_out_stream(mixed_mx_out_path.str().c_str(), cout);
                comp.printMixedMatrix(mixed_mx_out_stream);
                mixed_mx_out_stream.close();
            }

            // Send K-mer statistics to file
            std::ostringstream stats_out_path;
            stats_out_path << output_prefix << ".stats";
            ofstream_default stats_out_stream(stats_out_path.str().c_str(), cout);
            comp.printCounters(stats_out_stream);
            stats_out_stream.close();

            // Send K-mer statistics to stdout as well
            comp.printCounters(cout);


            return 0;
        }

    };
}
