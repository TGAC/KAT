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
#include <memory>
using std::unique_ptr;
using std::stringstream;

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "input_handler.hpp"
using kat::InputHandler;
 

typedef boost::error_info<struct FilterSeqError,string> FilterSeqErrorInfo;
struct FilterSeqException: virtual boost::exception, virtual std::exception { };

namespace kat { namespace filter {
    
const string    DEFAULT_FILT_SEQ_OUTPUT_PREFIX  = "kat.filter.seq";
const double    DEFAULT_FILT_SEQ_THRESHOLD      = 0.1;
const bool      DEFAULT_FILT_SEQ_INVERT         = false;
const bool      DEFAULT_FILT_SEQ_SEPARATE       = false;

struct SeqStats {
    int64_t index;
    uint64_t matches;
    uint64_t nb_kmers;
    
    SeqStats() : SeqStats(-1, 0, 0) {}
    SeqStats(const int32_t index, const uint64_t matches, const uint64_t nb_kmers) : 
        index(index), matches(matches), nb_kmers(nb_kmers) {}
    
    string toString() const {
        stringstream ss;
        ss << index << "\t" << nb_kmers << "\t" << matches;
        return ss.str();
    }
    
    double calcRatio() {
        return (double)matches / (double)nb_kmers;
    }
};

struct SeqStatsComparator {
    inline bool operator() (shared_ptr<SeqStats> k1, shared_ptr<SeqStats> k2) {
        return (k1->index < k2->index);
    }
};

class SeqFilterCounter {
    
public:
        
    vector<shared_ptr<SeqStats> > seq_stats;
    
    SeqFilterCounter() {}
    
    void add(shared_ptr<SeqStats> k) {
        seq_stats.push_back(std::move(k));
    }
    
    uint32_t calcKeepers(const double threshold) {
        uint32_t nb_found = 0;
        for(const auto& k : seq_stats) {
            if (k->calcRatio() >= threshold) {
                nb_found++;
            }
        }
        
        return nb_found;
    }
    
    size_t size() const {
        return seq_stats.size();
    }
    
    void sort() {
        // Sort keeper indices
        std::sort(seq_stats.begin(), seq_stats.end(), SeqStatsComparator()); 
    }
    
};

class ThreadedSeqStatsCounters {
private:
    vector<SeqFilterCounter> counters;

public:
    
    ThreadedSeqStatsCounters() : ThreadedSeqStatsCounters(1) {};
    ThreadedSeqStatsCounters(const uint16_t threads) {
        counters.resize(threads);
    }
    
    ~ThreadedSeqStatsCounters() {}

    void add(const uint32_t th_id, unique_ptr<SeqStats> keeper) {
        counters[th_id].add(std::move(keeper));            
    };

    unique_ptr<SeqFilterCounter> merge();
    
    void resize(uint16_t threads) {
        counters.resize(threads);
    }
};

class FilterSeq
{
private:
    
    static const uint16_t BATCH_SIZE = 1024;
    
    // Args
    InputHandler    input;
    path            seq_file;
    path            output_prefix;

    double      threshold;
    bool        invert;
    bool        separate;
    bool        doStats;
    uint16_t    threads;
    uint16_t    merLen;
    bool        verbose;
    
    vector<uint32_t> keep;
    uint16_t recordsInBatch;
    size_t bucket_size, remaining; 
    uint32_t offset;
        
    seqan::StringSet<seqan::CharString> names;
    seqan::StringSet<seqan::CharString> seqs;
    string extension;
    
    ThreadedSeqStatsCounters stats;

    void init(const vector<path>& _input);

public:
    FilterSeq(const path& _seq_file, const path& _input);

    FilterSeq(const path& _seq_file, const vector<path>& _input);

    ~FilterSeq() {}

    bool isInvert() const {
        return invert;
    }

    void setInvert(bool invert) {
        this->invert = invert;
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
    
    double getThreshold() const {
        return threshold;
    }

    void setThreshold(double threshold) {
        this->threshold = threshold;
    }
    
    bool isDoStats() const {
        return doStats;
    }

    void setDoStats(bool doStats) {
        this->doStats = doStats;
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
    

protected:

    void processSeqFile();
        
    void analyseBatch();

    void analyseBatchSlice(int th_id);
    
    void processSeq(const size_t index, const uint16_t th_id);
        
    void save(vector<shared_ptr<SeqStats> >& stats);
    
    static string helpMessage() {            
        
        return string(  "Usage: kat filter seq [options] <seq_file_to_filter> <input>\n\n") +
                        "Filter sequences based on whether those sequences contain specific k-mers.\n\n" \
                        "The user loads a k-mer hash and then filters sequences (either in or out) depending on whether those\n" \
                        "sequences contain the k-mer or not.  The user can also apply a threshold requiring X% of k-mers to be\n" \
                        "in the sequence before filtering is applied.\n\n" \
                        "Options";

    }

public:

    static int main(int argc, char *argv[]);

};
}}

