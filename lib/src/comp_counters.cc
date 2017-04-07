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

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <cmath>
using std::ostream;
using std::shared_ptr;
using std::string;
using std::unique_ptr;
using std::vector;
using std::endl;

#include <boost/filesystem/path.hpp>
using boost::filesystem::path;

#include <kat/distance_metrics.hpp>
using kat::DistanceMetric;

#include <kat/comp_counters.hpp>

// ********** CompCounters ***********

kat::CompCounters::CompCounters() : CompCounters("", "", "", DEFAULT_NB_BINS) {}

kat::CompCounters::CompCounters(const size_t _dm_size) : CompCounters("", "", "", _dm_size) {}

kat::CompCounters::CompCounters(const path& _hash1_path, const path& _hash2_path, const path& _hash3_path, const size_t _dm_size) :
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
    
    spectrum1.resize(_dm_size, 0);
    spectrum2.resize(_dm_size, 0);
    
    shared_spectrum1.resize(_dm_size, 0);
    shared_spectrum2.resize(_dm_size, 0);
}
        
kat::CompCounters::CompCounters(const CompCounters& o) {
    hash1_path = path(o.hash1_path);
    hash2_path = path(o.hash2_path);
    hash3_path = path(o.hash3_path);
    hash1_total = o.hash1_total;
    hash2_total = o.hash2_total;
    hash3_total = o.hash3_total;
    hash1_distinct = o.hash1_distinct;
    hash2_distinct = o.hash2_distinct;
    hash3_distinct = o.hash3_distinct;
    hash1_only_total = o.hash1_only_total;
    hash2_only_total = o.hash2_only_total;
    hash1_only_distinct = o.hash1_only_distinct;
    hash2_only_distinct = o.hash2_only_distinct;
    shared_hash1_total = o.shared_hash1_total;
    shared_hash2_total = o.shared_hash2_total;
    shared_distinct = o.shared_distinct;
    spectrum1 = o.spectrum1;
    spectrum2 = o.spectrum2;
    shared_spectrum1 = o.shared_spectrum1;
    shared_spectrum2 = o.shared_spectrum2;
}

void kat::CompCounters::updateHash1Counters(const uint64_t hash1_count, const uint64_t hash2_count) {
    hash1_total += hash1_count;
    hash1_distinct++;
    updateSpectrum(spectrum1, hash1_count);

    if (!hash2_count) {
        hash1_only_total += hash1_count;
        hash1_only_distinct++;
    }
}

void kat::CompCounters::updateHash2Counters(const uint64_t hash1_count, const uint64_t hash2_count) {
    hash2_total += hash2_count;
    hash2_distinct++;
    updateSpectrum(spectrum2, hash2_count);

    if (!hash1_count) {
        hash2_only_total += hash2_count;
        hash2_only_distinct++;
    }
}

void kat::CompCounters::updateHash3Counters(const uint64_t hash3_count) {

    hash3_total += hash3_count;
    hash3_distinct++;
}

void kat::CompCounters::updateSharedCounters(const uint64_t hash1_count, const uint64_t hash2_count) {

    if (hash1_count && hash2_count) {
        shared_hash1_total += hash1_count;
        shared_hash2_total += hash2_count;
        shared_distinct++;
        updateSpectrum(shared_spectrum1, hash1_count);
        updateSpectrum(shared_spectrum2, hash2_count);
    }
}

void kat::CompCounters::updateSpectrum(vector<uint64_t>& spectrum, const uint64_t count) {
    
    size_t s_size = spectrum.size();
    
    if (count <= 0)
        ++(spectrum)[0];
    else if (count >= s_size)
        ++(spectrum)[s_size - 1];
    else
        ++(spectrum)[count];
}



void kat::CompCounters::printCounts(ostream &out) {

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
    
    vector<unique_ptr<DistanceMetric>> dms;
    dms.push_back(unique_ptr<DistanceMetric>(new ManhattanDistance()));
    dms.push_back(unique_ptr<DistanceMetric>(new EuclideanDistance()));
    dms.push_back(unique_ptr<DistanceMetric>(new CosineDistance()));
    dms.push_back(unique_ptr<DistanceMetric>(new CanberraDistance()));
    dms.push_back(unique_ptr<DistanceMetric>(new JaccardDistance()));
        
    out << "Distance between spectra 1 and 2 (all k-mers):" << endl;
    for(auto& dm : dms) {
        out << " - " << dm->getName() << " distance: " << dm->calcDistance(spectrum1, spectrum2) << endl;
    }
    out << endl;
    
        
    out << "Distance between spectra 1 and 2 (shared k-mers):" << endl;
    for(auto& dm : dms) {
        out << " - " << dm->getName() << " distance: " << dm->calcDistance(shared_spectrum1, shared_spectrum2) << endl;
    }
    out << endl;
    
}
     

// ******** ThreadedCompCounters *********

kat::ThreadedCompCounters::ThreadedCompCounters() : ThreadedCompCounters("", "", "", DEFAULT_NB_BINS) {}

kat::ThreadedCompCounters::ThreadedCompCounters(const size_t _dm_size) : ThreadedCompCounters("", "", "", _dm_size) {}

kat::ThreadedCompCounters::ThreadedCompCounters(const path& _hash1_path, const path& _hash2_path, const path& _hash3_path, const size_t _dm_size) {
    final_matrix = CompCounters(_hash1_path, _hash2_path, _hash3_path, _dm_size);            
}
                
void kat::ThreadedCompCounters::printCounts(ostream &out) {
    final_matrix.printCounts(out);
}
        
void kat::ThreadedCompCounters::add(shared_ptr<CompCounters> cc) {
    cc->hash1_path = final_matrix.hash1_path;
    cc->hash2_path = final_matrix.hash2_path;
    cc->hash3_path = final_matrix.hash3_path;
    threaded_counters.push_back(*cc);
}
        
void kat::ThreadedCompCounters::merge() {

    // Merge counters
    for (const auto& itp : threaded_counters) {

        final_matrix.hash1_total += itp.hash1_total;
        final_matrix.hash2_total += itp.hash2_total;
        final_matrix.hash3_total += itp.hash3_total;
        final_matrix.hash1_distinct += itp.hash1_distinct;
        final_matrix.hash2_distinct += itp.hash2_distinct;
        final_matrix.hash3_distinct += itp.hash3_distinct;
        final_matrix.hash1_only_total += itp.hash1_only_total;
        final_matrix.hash2_only_total += itp.hash2_only_total;
        final_matrix.hash1_only_distinct += itp.hash1_only_distinct;
        final_matrix.hash2_only_distinct += itp.hash2_only_distinct;
        final_matrix.shared_hash1_total += itp.shared_hash1_total;
        final_matrix.shared_hash2_total += itp.shared_hash2_total;
        final_matrix.shared_distinct += itp.shared_distinct;
        
        merge_spectrum(final_matrix.spectrum1, itp.spectrum1);
        merge_spectrum(final_matrix.spectrum2, itp.spectrum2);
        merge_spectrum(final_matrix.shared_spectrum1, itp.shared_spectrum1);
        merge_spectrum(final_matrix.shared_spectrum2, itp.shared_spectrum2);
    }
}

        
void kat::ThreadedCompCounters::merge_spectrum(vector<uint64_t>& spectrum, const vector<uint64_t>& threaded_spectrum) {
    
    for(size_t i = 0; i < spectrum.size(); i++) {
        spectrum[i] += threaded_spectrum[i];
    }
}
