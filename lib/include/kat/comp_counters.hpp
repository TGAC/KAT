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
#include <memory>
#include <vector>
using std::ostream;
using std::shared_ptr;
using std::vector;

#include <boost/filesystem/path.hpp>
using boost::filesystem::path;

namespace kat {
    
const uint32_t   DEFAULT_NB_BINS = 1001;

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

    vector<uint64_t> spectrum1;
    vector<uint64_t> spectrum2;
    vector<uint64_t> shared_spectrum1;
    vector<uint64_t> shared_spectrum2;

    path hash1_path;
    path hash2_path;
    path hash3_path;

    CompCounters();

    CompCounters(const path& _hash1_path, const path& _hash2_path, const path& _hash3_path, const size_t _dm_size);

    CompCounters(const CompCounters& o);

    void updateHash1Counters(const uint64_t hash1_count, const uint64_t hash2_count);

    void updateHash2Counters(const uint64_t hash1_count, const uint64_t hash2_count);

    void updateHash3Counters(const uint64_t hash3_count);

    void updateSharedCounters(const uint64_t hash1_count, const uint64_t hash2_count);

    static void updateSpectrum(vector<uint64_t>& spectrum, const uint64_t count);

    void printCounts(ostream &out);

    vector<uint64_t>& getSpectrum1() { return spectrum1; }

    vector<uint64_t>& getSpectrum2() { return spectrum2; }       

};

class ThreadedCompCounters {
private: 
    uint16_t threads;

    CompCounters final_matrix;
    vector<CompCounters> threaded_counters;

    static void merge_spectrum(vector<uint64_t>& spectrum, const vector<uint64_t>& threaded_spectrum);

public:

    ThreadedCompCounters();

    ThreadedCompCounters(const path& _hash1_path, const path& _hash2_path, const path& _hash3_path, const size_t _dm_size);

    void printCounts(ostream &out);

    void add(shared_ptr<CompCounters> cc);

    size_t size() {
        return threaded_counters.size();
    }

    CompCounters& getFinalMatrix() {
        return final_matrix;
    }

    CompCounters& getThreadedMatrixAt(uint16_t index) {
        return threaded_counters[index];
    }

    void merge();

};
}