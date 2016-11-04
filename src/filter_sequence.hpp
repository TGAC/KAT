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
#include <fstream>
using std::unique_ptr;
using std::stringstream;
using std::ofstream;

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <kat/input_handler.hpp>
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


class FilterSeq
{
private:
    
    // Args
    InputHandler    input;
    path            seq_file_1;
    path            seq_file_2;
    path            output_prefix;

    double      threshold;
    bool        invert;
    bool        separate;
    bool        doStats;
    uint16_t    threads;
    bool        verbose;
    
    uint64_t    keepers;
    uint64_t    total;
    
    seqan::CharString name;
    seqan::CharString seq;
    seqan::CharString qual;
    seqan::CharString name2;
    seqan::CharString seq2;
    seqan::CharString qual2;
    string extension;
    
    unique_ptr<seqan::SeqFileIn> reader = nullptr;
    unique_ptr<seqan::SeqFileIn> reader2 = nullptr;    
    
    unique_ptr<seqan::SeqFileOut> inWriter = nullptr;
    unique_ptr<seqan::SeqFileOut> outWriter = nullptr;
    unique_ptr<seqan::SeqFileOut> inWriter2 = nullptr;
    unique_ptr<seqan::SeqFileOut> outWriter2 = nullptr;
    
    unique_ptr<ofstream> stats_stream = nullptr;
    
    void init(const vector<path>& _input);

public:
    FilterSeq(const path& _seq_file_1, const path& _seq_file_2, const path& _input);

    FilterSeq(const path& _seq_file_1, const path& _seq_file_2, const vector<path>& _input);

    ~FilterSeq() {}
    
    uint16_t getThreads() const {
        return threads;
    }

    void setThreads(uint16_t threads) {
        this->threads = threads;
    }


    bool isInvert() const {
        return invert;
    }

    void setInvert(bool invert) {
        this->invert = invert;
    }

    uint16_t getMerLen() const {
        return input.merLen;
    }

    void setMerLen(uint16_t merLen) {
        this->input.merLen = merLen;
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
    
    bool isPaired() const {
        return !this->seq_file_2.empty(); 
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
    void processPairedSeqFile();
        
    void processSeq(uint64_t index);
    
    void getProfile(seqan::CharString& s, vector<bool>& hits);
        
    
    static string helpMessage() {            
        
        return string(  "Usage: kat filter seq [options] <seq_file_to_filter> [<seq_file_2>] <input>\n\n") +
                        "Filter sequences based on whether those sequences contain specific k-mers.\n\n" \
                        "The user loads a k-mer hash and then filters sequences (either in or out) depending on whether those\n" \
                        "sequences contain the k-mer or not.  The user can also apply a threshold requiring X% of k-mers to be\n" \
                        "in the sequence before filtering is applied.\n\n" \
                        "Should the user have paired-end data to filter the first two positional arguments represent the paired\n" \
                        "end read files to filter, and the remaining positional arguments are for loading the kmer hash.  If\n" \
                        "user wants filter paired end reads then the --paired option must be selected\n\n" \
                        "Options";

    }

public:

    static int main(int argc, char *argv[]);

};
}}

