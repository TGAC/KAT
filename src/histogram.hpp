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

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <thread>
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

#include <matrix/matrix_metadata_extractor.hpp>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish_helper.hpp>

typedef boost::error_info<struct HistogramError,string> HistogramErrorInfo;
struct HistogramException: virtual boost::exception, virtual std::exception { };

namespace kat {

    class Histogram {
    private:
        
        // Arguments from user
        vector<path>    inputs;
        path            outputPrefix;
        uint16_t        threads;
        uint64_t        low;
        uint64_t        high;
        bool            canonical;
        uint16_t        merLen;
        uint64_t        hashSize;            
        bool            verbose;

        // Jellyfish mapped file hash vars
        shared_ptr<JellyfishHelper> jfh;
        path hashFile;
        
        // Internal vars
        uint64_t base, ceil, inc, nb_buckets, nb_slices;
        vector<uint64_t> data;
        vector<shared_ptr<vector<uint64_t>>> threadedData;
        uint64_t slice_id;

    public:

        Histogram(vector<path> _inputs, uint64_t _low, uint64_t _high, uint64_t _inc);

        virtual ~Histogram() {
        }
        
        bool isCanonical() const {
            return canonical;
        }

        void setCanonical(bool canonical) {
            this->canonical = canonical;
        }
        
        uint64_t getHashSize() const {
            return hashSize;
        }

        void setHashSize(uint64_t hash_size) {
            this->hashSize = hash_size;
        }
        
        uint16_t getMerLen() const {
            return merLen;
        }

        void setMerLen(uint16_t merLen) {
            this->merLen = merLen;
        }

        path getOutputPrefix() const {
            return outputPrefix;
        }

        void setOutputPrefix(path outputPrefix) {
            this->outputPrefix = outputPrefix;
        }

        uint64_t getHigh() const {
            return high;
        }

        void setHigh(uint64_t high) {
            this->high = high;
        }

        uint64_t getInc() const {
            return inc;
        }

        void setInc(uint64_t inc) {
            this->inc = inc;
        }

        vector<path> getInputs() const {
            return inputs;
        }

        void setInputs(vector<path> inputs) {
            this->inputs = inputs;
        }

        uint64_t getLow() const {
            return low;
        }

        void setLow(uint64_t low) {
            this->low = low;
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
        
        void print(std::ostream &out);
        
    protected:
        
        uint64_t calcBase() {
            return low > 1 ? (1 >= low ? 1 : low - 1) : 1;
        }

        uint64_t calcCeil() {
            return high + 1;
        }
        
        void merge();
        
        void loadHashes();
        
        void startAndJoinThreads();
         
        void start(int th_id);
        
        static string helpMessage(){
            
            return string("Usage: kat hist [options] <jellyfish_hash>\n\n") +
                            "Create an histogram of k-mer occurrences in a sequence file.\n\n" +
                            "Create an histogram with the number of k-mers having a given count. In bucket 'i' are tallied the k-mers " \
                            "which have a count 'c' satisfying 'low+i*inc <= c < low+(i+1)'. Buckets in the output are labeled by the " \
                            "low end point (low+i).\n" \
                            "The last bucket in the output behaves as a catchall: it tallies all k-mers with a count greater or equal to " \
                            "the low end point of this bucket.\n" \
                            "This tool is very similar to the \"histo\" tool in jellyfish itself.  The primary difference being that the " \
                            "output contains metadata that make the histogram easier for the user to plot.\n\n" \
                            "Options";

        }
      
    public:
        
        static int main(int argc, char *argv[]);
    };
}
