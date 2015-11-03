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

#include <memory>
using std::unique_ptr;

#include "input_handler.hpp"
using kat::InputHandler;
 

typedef boost::error_info<struct FilterKmerError,string> FilterKmerErrorInfo;
struct FilterKmerException: virtual boost::exception, virtual std::exception { };

namespace kat { namespace filter {
    
const uint64_t  DEFAULT_FILT_KMER_LOW_COUNT    = 0;
const uint64_t  DEFAULT_FILT_KMER_HIGH_COUNT   = 10000;
const uint16_t  DEFAULT_FILT_KMER_LOW_GC       = 0;
const uint16_t  DEFAULT_FILT_KMER_HIGH_GC      = 31;
const bool      DEFAULT_FILT_KMER_INVERT       = false;
const bool      DEFAULT_FILT_KMER_SEPARATE     = false;
    

class Counter {
public:
    uint64_t distinct;
    uint64_t total;

    Counter() : distinct(0), total(0)
    {}

    ~Counter() {}
    
    void increment(const uint64_t total_inc) {
        distinct++;
        total += total_inc;
    }
    
    string toString() const;
};

class ThreadedCounter {
private:
    vector<Counter> counter;

public:
    
    ThreadedCounter() : ThreadedCounter(1) {};
    ThreadedCounter(const uint16_t threads) {
        counter.resize(threads);
    }
    
    ~ThreadedCounter() {}

    void increment(const uint16_t th_id, const uint64_t total_inc);

    unique_ptr<Counter> merge();
    
    void resize(uint16_t threads) {
        counter.resize(threads);
    }
};

class FilterKmer
{
private:
    // Args
    InputHandler    input;
    path            output_prefix;

    uint64_t    low_count;
    uint64_t    high_count;
    uint32_t    low_gc;
    uint32_t    high_gc;
    bool        invert;
    bool        separate;
    uint16_t    threads;
    uint16_t    merLen;
    bool        verbose;

    ThreadedCounter all;
    ThreadedCounter in;
    ThreadedCounter out;

    void init(const vector<path>& _input);

public:
    FilterKmer(const path& _input);

    FilterKmer(const vector<path>& _input);

    ~FilterKmer() {}

    uint64_t getHigh_count() const {
        return high_count;
    }

    void setHigh_count(uint64_t high_count) {
        this->high_count = high_count;
    }

    uint32_t getHigh_gc() const {
        return high_gc;
    }

    void setHigh_gc(uint32_t high_gc) {
        this->high_gc = high_gc;
    }

    bool isInvert() const {
        return invert;
    }

    void setInvert(bool invert) {
        this->invert = invert;
    }

    uint64_t getLow_count() const {
        return low_count;
    }

    void setLow_count(uint64_t low_count) {
        this->low_count = low_count;
    }

    uint32_t getLow_gc() const {
        return low_gc;
    }

    void setLow_gc(uint32_t low_gc) {
        this->low_gc = low_gc;
    }

    uint16_t getMerLen() const {
        return merLen;
    }

    void setMerLen(uint16_t merLen) {
        this->merLen = merLen;
    }

    path getOutputPrefix() const {
        return output_prefix;
    }

    void setOutputPrefix(path output) {
        this->output_prefix = output;
    }

    bool isSeparate() const {
        return separate;
    }

    void setSeparate(bool separate) {
        this->separate = separate;
    }
    
    bool isCanonical() const {
        return input.canonical;
    }

    void setCanonical(bool canonical) {
        this->input.canonical = canonical;
    }

    path getOutput_prefix() const {
        return output_prefix;
    }

    void setOutput_prefix(path output_prefix) {
        this->output_prefix = output_prefix;
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
            
    bool isVerbose() const {
        return verbose;
    }

    void setVerbose(bool verbose) {
        this->verbose = verbose;
    }


    void execute();
    
    void plot();


protected:

    void filter(HashCounter& inCounter, HashCounter& outCounter);

    void filterSlice(int th_id, HashCounter& inCounter, HashCounter& outCounter);

    bool inBounds(const string& kmer_seq, const uint64_t& kmer_count);

    void dump(path& out_path, HashCounter* hash, file_header& header);

    void merge();
    
    static string helpMessage() {            
        
        return string(  "Usage: kat filter kmer [options] <input>\n\n") +
                        "Filter kmers to those within defined bounds and those outside.\n\n" \
                        "The user can produce K-mer hashes, within and outside user defined GC and k-mer coverage bounds.\n" \
                        "This is useful for isolating k-mers that could be attributable to contamination, or for contamination\n" \
                        "removal.  Normally, the user would identify such regions using plots from the GCP tool.\n\n" \
                        "Options";

    }

public:

    static int main(int argc, char *argv[]);

};
}}

