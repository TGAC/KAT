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
#include <string.h>
#include <stdint.h>
#include <vector>
#include <math.h>
#include <memory>
#include <mutex>
#include <thread>
using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::ostream;
using std::ofstream;
using std::shared_ptr;
using std::unique_ptr;
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

#include "inc/matrix/matrix_metadata_extractor.hpp"
#include "inc/matrix/sparse_matrix.hpp"
#include "inc/matrix/threaded_sparse_matrix.hpp"
#include "inc/distance_metrics.hpp"

#include "input_handler.hpp"
#include "plot_spectra_cn.hpp"
#include "plot_density.hpp"
using kat::JellyfishHelper;
using kat::HashLoader;
using kat::InputHandler;
using kat::PlotSpectraCn;
using kat::PlotDensity;


#include "comp.hpp"

// ********** CompCounters ***********

kat::CompCounters::CompCounters() : CompCounters("", "", "", DEFAULT_NB_BINS) {}

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

    
// ********* Comp **********

kat::Comp::Comp() :
    kat::Comp::Comp(path(), path()) {}

kat::Comp::Comp(const path& _input1, const path& _input2) : 
    kat::Comp::Comp(_input1, _input2, path()) {}
        
kat::Comp::Comp(const path& _input1, const path& _input2, const path& _input3) {
    vector<path> vecInput1;
    vecInput1.push_back(_input1);
    vector<path> vecInput2;
    vecInput2.push_back(_input2);
    vector<path> vecInput3;
    if (!_input3.empty())
        vecInput3.push_back(_input3);    
    
    init(vecInput1, vecInput2, vecInput3);
}
    
    
kat::Comp::Comp(const vector<path>& _input1, const vector<path>& _input2) : 
    kat::Comp::Comp(_input1, _input2, vector<path>()) {}
        
kat::Comp::Comp(const vector<path>& _input1, const vector<path>& _input2, const vector<path>& _input3) {    
    init(_input1, _input2, _input3);
}

void kat::Comp::init(const vector<path>& _input1, const vector<path>& _input2, const vector<path>& _input3) {
    input = !_input3.empty() ? vector<InputHandler>(3) : vector<InputHandler>(2);
    
    input[0].setMultipleInputs(_input1);
    input[1].setMultipleInputs(_input2);
    input[0].index = 1;
    input[1].index = 2;
    if (!_input3.empty()) {
        input[2].setMultipleInputs(_input3);
        input[2].index = 3;
    }
    
    outputPrefix = "kat-comp";
    d1Scale = 1.0;
    d2Scale = 1.0;
    d1Bins = DEFAULT_NB_BINS;
    d2Bins = DEFAULT_NB_BINS;
    ioThreads = 1;
    analysisThreads = 1;
    merLen = DEFAULT_MER_LEN; 
    densityPlot = false;
    verbose = false;      
}

void kat::Comp::execute() {

    // Check input files exist and determine input mode
    for(uint16_t i = 0; i < input.size(); i++) {
        input[i].validateInput();
    }
        
    // Create output directory
    path parentDir = bfs::absolute(outputPrefix).parent_path();
    if (!bfs::exists(parentDir) || !bfs::is_directory(parentDir)) {
        if (!bfs::create_directories(parentDir)) {
            BOOST_THROW_EXCEPTION(CompException() << CompErrorInfo(string(
                    "Could not create output directory: ") + parentDir.string()));
        }
    }
    
    
    // Create the final K-mer counter matrices
    main_matrix = ThreadedSparseMatrix(d1Bins, d2Bins, analysisThreads);

    // Initialise extra matrices for hash3 (only allocates space if required)
    if (doThirdHash()) {

        ends_matrix = ThreadedSparseMatrix(d1Bins, d2Bins, analysisThreads);
        middle_matrix = ThreadedSparseMatrix(d1Bins, d2Bins, analysisThreads);
        mixed_matrix = ThreadedSparseMatrix(d1Bins, d2Bins, analysisThreads);
    }

    // Create the comp counters for each thread
    comp_counters = ThreadedCompCounters(
            input[0].getSingleInput(), 
            input[1].getSingleInput(), 
            input.size() == 3 ? input[2].getSingleInput() : path(),
            std::min(d1Bins, d2Bins));

    std::ostream* out_stream = verbose ? &cerr : (std::ostream*)0;

    string merLenStr = lexical_cast<string>(merLen);

    // Count kmers in sequence files if necessary (sets load and hashes and hashcounters as appropriate)
    for(int i = 0; i < input.size(); i++) {
        if (input[i].mode == InputHandler::InputHandler::InputMode::COUNT) {
            input[i].count(merLen, ioThreads);
        }
    }
    
    // Check to see if user specified any hashes to load
    bool anyLoad = false;
    bool allLoad = true;
    bool anyDump = false;
    for(uint16_t i = 0; i < input.size(); i++) {
        if (input[i].mode == InputHandler::InputMode::LOAD) {
            input[i].loadHeader();
            anyLoad = true;            
        }
        else if (input[i].mode == InputHandler::InputHandler::InputMode::COUNT) {
            anyDump = true;
            allLoad = false;
        }
    }
    
    // If all hashes are loaded directly there is no requirement that the user needs
    // to specify the merLen, so just set it to the merLen found in the header of the first input
    if (allLoad) merLen = input[0].header->key_len() / 2;
    
    for(uint16_t i = 0; i < input.size(); i++) {
        input[i].validateMerLen(merLen);
    }
    
    // Load any hashes if necessary
    if (anyLoad) loadHashes();
     
    
    // Run the threads
    compare();

    // Dump any hashes that were previously counted to disk if requested
    // NOTE: MUST BE DONE AFTER COMPARISON AS THIS CLEARS ENTRIES FROM HASH ARRAY!
    if (dumpHashes()) {
        for(uint16_t i = 0; i < input.size(); i++) {        
            path outputPath(outputPrefix.string() + "-hash" + lexical_cast<string>(input[i].index) + ".jf" + lexical_cast<string>(merLen));
            input[i].dump(outputPath, ioThreads, true);
        }    
    }    

    // Merge results
    merge();    
}

void kat::Comp::save() {
    
    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        

    cout << "Saving results to disk ...";
    cout.flush();
    
    // Send main matrix to output file
    ofstream main_mx_out_stream(string(outputPrefix.string() + "-main.mx").c_str());
    printMainMatrix(main_mx_out_stream);
    main_mx_out_stream.close();

    // Output ends matrices if required
    if (doThirdHash()) {
        // Ends matrix
        ofstream ends_mx_out_stream(string(outputPrefix.string() + "-ends.mx").c_str());
        printEndsMatrix(ends_mx_out_stream);
        ends_mx_out_stream.close();

        // Middle matrix
        ofstream middle_mx_out_stream(string(outputPrefix.string() + "-middle.mx").c_str());
        printMiddleMatrix(middle_mx_out_stream);
        middle_mx_out_stream.close();

        // Mixed matrix
        ofstream mixed_mx_out_stream(string(outputPrefix.string() + "-mixed.mx").c_str());
        printMixedMatrix(mixed_mx_out_stream);
        mixed_mx_out_stream.close();
    }

    // Send K-mer statistics to file
    ofstream stats_out_stream(string(outputPrefix.string() + ".stats").c_str());
    printCounters(stats_out_stream);
    stats_out_stream.close();

    cout << " done.";
    cout.flush();
}

void kat::Comp::merge() {
    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        

    cout << "Merging results ...";
    cout.flush();

    // Merge results from the threads
    main_matrix.mergeThreadedMatricies();
    if (doThirdHash()) {
        ends_matrix.mergeThreadedMatricies();
        middle_matrix.mergeThreadedMatricies();
        mixed_matrix.mergeThreadedMatricies();
    }

    comp_counters.merge();

    cout << " done.";
    cout.flush();
}

void kat::Comp::loadHashes() {
    
    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        
    
    cout << "Loading hashes into memory...";
    cout.flush();    

    // If using parallel IO load hashes in parallel, otherwise do one at a time
    if (ioThreads > 1) {
        
        thread threads[input.size()];
        
        void (kat::InputHandler::*memfunc)() = &kat::InputHandler::loadHash;
        
        for(int i = 0; i < input.size(); i++) {
            if (input[i].mode == InputHandler::InputMode::LOAD) {
                threads[i] = thread(memfunc, &input[i]);                
            }
        }
        
        for(int i = 0; i < input.size(); i++) {
            if (input[i].mode == InputHandler::InputMode::LOAD) {
                threads[i].join();                 
            }
        }        
    }
    else {
        for(int i = 0; i < input.size(); i++) {
            if (input[i].mode == InputHandler::InputMode::LOAD) {
               input[i].loadHash(false);                
            }
        }        
    }
    
    cout << " done.";
    cout.flush();
}


// Print K-mer comparison matrix

void kat::Comp::printMainMatrix(ostream &out) {

    const SM64& mx = main_matrix.getFinalMatrix();

    out << mme::KEY_TITLE << "K-mer comparison plot" << endl
            << mme::KEY_X_LABEL << "K-mer multiplicity for: " << input[0].getSingleInput().string() << endl
            << mme::KEY_Y_LABEL << "K-mer multiplicity for: " << input[1].getSingleInput().string() << endl
            << mme::KEY_Z_LABEL << "Distinct K-mers per bin" << endl
            << mme::KEY_NB_COLUMNS << mx.height() << endl
            << mme::KEY_NB_ROWS << mx.width() << endl
            << mme::KEY_MAX_VAL << mx.getMaxVal() << endl
            << mme::KEY_TRANSPOSE << "1" << endl
            << mme::MX_META_END << endl;

    mx.printMatrix(out);
}

// Print K-mer comparison matrix

void kat::Comp::printEndsMatrix(ostream &out) {

    out << "# Each row represents K-mer multiplicity for: " << input[0].getSingleInput().string() << endl;
    out << "# Each column represents K-mer multiplicity for sequence ends: " << input[2].getSingleInput().string() << endl;

    ends_matrix.getFinalMatrix().printMatrix(out);
}

// Print K-mer comparison matrix

void kat::Comp::printMiddleMatrix(ostream &out) {

    out << "# Each row represents K-mer multiplicity for: " << input[0].getSingleInput().string() << endl;
    out << "# Each column represents K-mer multiplicity for sequence middles: " << input[1].getSingleInput().string() << endl;

    middle_matrix.getFinalMatrix().printMatrix(out);
}

// Print K-mer comparison matrix

void kat::Comp::printMixedMatrix(ostream &out) {

    out << "# Each row represents K-mer multiplicity for hash file 1: " << input[0].getSingleInput().string() << endl;
    out << "# Each column represents K-mer multiplicity for mixed: " << input[1].getSingleInput().string() << " and " << input[2].getSingleInput().string() << endl;

    mixed_matrix.getFinalMatrix().printMatrix(out);
}

// Print K-mer statistics

void kat::Comp::printCounters(ostream &out) {

    comp_counters.printCounts(out);
}

        
void kat::Comp::compare() {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        

    cout << "Comparing hashes ...";
    cout.flush();
    
    thread t[analysisThreads];

    for(int i = 0; i < analysisThreads; i++) {
        t[i] = thread(&Comp::compareSlice, this, i);
    }

    for(int i = 0; i < analysisThreads; i++){
        t[i].join();
    }
    
    cout << " done.";
    cout.flush();
}

void kat::Comp::compareSlice(int th_id) {

    shared_ptr<CompCounters> cc = make_shared<CompCounters>();

    // Setup iterator for this thread's chunk of hash1
    LargeHashArray::eager_iterator hash1Iterator = input[0].hash->eager_slice(th_id, analysisThreads);

    // Go through this thread's slice for hash1
    while (hash1Iterator.next()) {
        
        // Get the current K-mer count for hash1
        uint64_t hash1_count = hash1Iterator.val();

        // Get the count for this K-mer in hash2 (assuming it exists... 0 if not)
        uint64_t hash2_count = JellyfishHelper::getCount(input[1].hash, hash1Iterator.key(), input[1].canonical);

        // Get the count for this K-mer in hash3 (assuming it exists... 0 if not)
        uint64_t hash3_count = doThirdHash() ? JellyfishHelper::getCount(input[2].hash, hash1Iterator.key(), input[2].canonical) : 0;

        // Increment hash1's unique counters
        cc->updateHash1Counters(hash1_count, hash2_count);

        // Increment shared counters
        cc->updateSharedCounters(hash1_count, hash2_count);

        // Scale counters to make the matrix look pretty
        uint64_t scaled_hash1_count = scaleCounter(hash1_count, d1Scale);
        uint64_t scaled_hash2_count = scaleCounter(hash2_count, d2Scale);
        uint64_t scaled_hash3_count = scaleCounter(hash3_count, d2Scale);

        // Modifies hash counts so that K-mer counts larger than MATRIX_SIZE are dumped in the last slot
        if (scaled_hash1_count >= d1Bins) scaled_hash1_count = d1Bins - 1;
        if (scaled_hash2_count >= d2Bins) scaled_hash2_count = d2Bins - 1;
        if (scaled_hash3_count >= d2Bins) scaled_hash3_count = d2Bins - 1;

        // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
        main_matrix.incTM(th_id, scaled_hash1_count, scaled_hash2_count, 1);

        // Update hash 3 related matricies if hash 3 was provided
        if (doThirdHash()) {
            if (scaled_hash2_count == scaled_hash3_count)
                ends_matrix.incTM(th_id, scaled_hash1_count, scaled_hash3_count, 1);
            else if (scaled_hash3_count > 0)
                mixed_matrix.incTM(th_id, scaled_hash1_count, scaled_hash3_count, 1);
            else
                middle_matrix.incTM(th_id, scaled_hash1_count, scaled_hash3_count, 1);
        }
    }

    // Setup iterator for this thread's chunk of hash2
    // We setup hash2 for random access, so hopefully performance isn't too bad here...
    // Hash2 should be smaller than hash1 in most cases so hopefully we can get away with this.
    LargeHashArray::eager_iterator hash2Iterator = input[1].hash->eager_slice(th_id, analysisThreads);

    // Iterate through this thread's slice of hash2
    while (hash2Iterator.next()) {
        // Get the current K-mer count for hash2
        uint64_t hash2_count = hash2Iterator.val();

        // Get the count for this K-mer in hash1 (assuming it exists... 0 if not)
        uint64_t hash1_count = JellyfishHelper::getCount(input[0].hash, hash2Iterator.key(), input[0].hash);

        // Increment hash2's unique counters (don't bother with shared counters... we've already done this)
        cc->updateHash2Counters(hash1_count, hash2_count);

        // Only bother updating thread matrix with K-mers not found in hash1 (we've already done the rest)
        if (hash1_count == 0) {
            // Scale counters to make the matrix look pretty
            uint64_t scaled_hash2_count = scaleCounter(hash2_count, d2Scale);

            // Modifies hash counts so that K-mer counts larger than MATRIX_SIZE are dumped in the last slot
            if (scaled_hash2_count >= d2Bins) scaled_hash2_count = d2Bins - 1;

            // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
            main_matrix.incTM(th_id, 0, scaled_hash2_count, 1);
        }
    }

    // Only update hash3 counters if hash3 was provided
    if (doThirdHash()) {
        // Setup iterator for this thread's chunk of hash3
        LargeHashArray::eager_iterator hash3Iterator = input[2].hash->eager_slice(th_id, analysisThreads);

        // Iterate through this thread's slice of hash2
        while (hash3Iterator.next()) {
            // Get the current K-mer count for hash2

            uint64_t hash3_count = hash3Iterator.val();

            // Increment hash3's unique counters (don't bother with shared counters... we've already done this)
            cc->updateHash3Counters(hash3_count);
        }
    }

    mu.lock();
    comp_counters.add(cc);
    mu.unlock();
}


void kat::Comp::plot(const string& output_type) {
    
    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        

    cout << "Creating plot ...";
    cout.flush();
    
    bool res = true;
    
    // Plot results
    if (densityPlot) {
        PlotDensity pd(getMxOutPath(), path(getMxOutPath().string() + ".density." + output_type));
        pd.setXLabel(string("# Distinct kmers for ") + input[0].pathString());
        pd.setYLabel(string("# Distinct kmers for ") + input[1].pathString());
        pd.setZLabel("Kmer multiplicity");
        pd.setTitle(string("Spectra Density Plot for: ") + input[0].pathString() + " vs " + input[1].pathString());
        pd.setOutputType(output_type);
        res = pd.plot();
    }
    else {
        PlotSpectraCn pscn(getMxOutPath(), path(getMxOutPath().string() + ".spectra-cn." + output_type));
        pscn.setTitle(string("Spectra CN Plot for: ") + input[0].pathString() + " vs " + input[1].pathString());
        pscn.setYLabel("# Distinct kmers");
        pscn.setXLabel("Kmer multiplicity");
        pscn.setOutputType(output_type);
        res = pscn.plot();
    }
    
    if (!res) {
        cout << "WARNING: gnuplot session not valid.  Probably gnuplot is not installed correctly on your system.  No plots produced.";
    }
    else {    
        cout << " done.";
    }
    
    cout.flush();
}



int kat::Comp::main(int argc, char *argv[]) {

    string input1;
    string input2;
    string input3;
    path output_prefix;
    double d1_scale;
    double d2_scale;
    uint16_t d1_bins;
    uint16_t d2_bins;
    uint16_t io_threads;
    uint16_t analysis_threads;
    uint16_t mer_len;
    bool canonical_1;
    bool canonical_2;
    bool canonical_3;
    uint64_t hash_size_1;
    uint64_t hash_size_2;
    uint64_t hash_size_3;
    bool dump_hashes;
    bool disable_hash_grow;
    bool density_plot;
    string plot_output_type;
    bool verbose;
    bool help;

    // Declare the supported options.
    po::options_description generic_options(Comp::helpMessage(), 100);
    generic_options.add_options()
            ("output_prefix,o", po::value<path>(&output_prefix)->default_value("kat-comp"), 
                "Path prefix for files generated by this program.")
            ("threads,t", po::value<uint16_t>(&io_threads)->default_value(1),
                "The number of threads to use.  Keep in mind initially threads will be used to IO purposes (counting kmers from file), so don't raise this above what your storage system can manage.")
            ("analysis_threads,T", po::value<uint16_t>(&analysis_threads)->default_value(1),
                "The number of threads to use for analysis stage.  This option only applies if you set it higher than option 't', otherwise it will have the same value as option 't'.  This option is intended to allow users to get extra performance on high core systems which might by bottlenecked by IO.")
            ("d1_scale,x", po::value<double>(&d1_scale)->default_value(1.0),
                "Scaling factor for the first dataset - float multiplier")
            ("d2_scale,y", po::value<double>(&d2_scale)->default_value(1.0),
                "Scaling factor for the second dataset - float multiplier")
            ("d1_bins,i", po::value<uint16_t>(&d1_bins)->default_value(1001),
                "Number of bins for the first dataset.  i.e. number of rows in the matrix")
            ("d2_bins,j", po::value<uint16_t>(&d2_bins)->default_value(1001),
                "Number of bins for the second dataset.  i.e. number of rows in the matrix")
            ("canonical1,C", po::bool_switch(&canonical_1)->default_value(false),
                "If counting fast(a/q) for input 1, this option specifies whether the jellyfish hash represents K-mers produced for both strands (canonical), or only the explicit kmer found.")
            ("canonical2,D", po::bool_switch(&canonical_2)->default_value(false),
                "If counting fast(a/q) for input 2, this option specifies whether the jellyfish hash represents K-mers produced for both strands (canonical), or only the explicit kmer found.")
            ("canonical3,E", po::bool_switch(&canonical_3)->default_value(false),
                "If counting fast(a/q) for input 3, this option specifies whether the jellyfish hash represents K-mers produced for both strands (canonical), or only the explicit kmer found.")
            ("mer_len,m", po::value<uint16_t>(&mer_len)->default_value(DEFAULT_MER_LEN),
                "The kmer length to use in the kmer hashes.  Larger values will provide more discriminating power between kmers but at the expense of additional memory and lower coverage.")
            ("hash_size_1,H", po::value<uint64_t>(&hash_size_1)->default_value(DEFAULT_HASH_SIZE),
                "If kmer counting is required for input 1, then use this value as the hash size.  If this hash size is not large enough for your dataset then the default behaviour is to double the size of the hash and recount, which will increase runtime and memory usage.")
            ("hash_size_2,I", po::value<uint64_t>(&hash_size_2)->default_value(DEFAULT_HASH_SIZE),
                "If kmer counting is required for input 2, then use this value as the hash size.  If this hash size is not large enough for your dataset then the default behaviour is to double the size of the hash and recount, which will increase runtime and memory usage.")
            ("hash_size_3,J", po::value<uint64_t>(&hash_size_3)->default_value(DEFAULT_HASH_SIZE),
                "If kmer counting is required for input 3, then use this value as the hash size.  If this hash size is not large enough for your dataset then the default behaviour is to double the size of the hash and recount, which will increase runtime and memory usage.")
            ("dump_hashes,d", po::bool_switch(&dump_hashes)->default_value(false), 
                "Dumps any jellyfish hashes to disk that were produced during this run.")
            ("disable_hash_grow,g", po::bool_switch(&disable_hash_grow)->default_value(false), 
                "By default jellyfish will double the size of the hash if it gets filled, and then attempt to recount.  Setting this option to true, disables automatic hash growing.  If the hash gets filled an error is thrown.  This option is useful if you are working with large genomes, or have strict memory limits on your system.")   
            ("density_plot,p", po::bool_switch(&density_plot)->default_value(false),
                "Makes a spectra_mx plot.  By default we create a spectra_cn plot.")
            ("output_type,p", po::value<string>(&plot_output_type)->default_value(DEFAULT_COMP_PLOT_OUTPUT_TYPE), 
                "The plot file type to create: png, ps, pdf.  Warning... if pdf is selected please ensure your gnuplot installation can export pdf files.")
            ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                "Print extra information.")
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
            ("input_1,i1", po::value<string>(&input1), "Path to the first input file.  Can be either FastA, FastQ or a jellyfish hash (non bloom filtered)")
            ("input_2,i2", po::value<string>(&input2), "Path to the second input file.  Can be either FastA, FastQ or a jellyfish hash (non bloom filtered)")
            ("input_3,i3", po::value<string>(&input3), "Path to the third input file.  Can be either FastA, FastQ or a jellyfish hash (non bloom filtered)")
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


    // Make sure we have at least the same number of cores allocated for analysis as we have for io.
    if (analysis_threads < io_threads) {
        analysis_threads = io_threads;
    }
    

    auto_cpu_timer timer(1, "KAT COMP completed.\nTotal runtime: %ws\n\n");        

    cout << "Running KAT in COMP mode" << endl
         << "------------------------" << endl << endl;

    // Glob input files
    vector<path> vecinput1 = InputHandler::globFiles(input1);
    vector<path> vecinput2 = InputHandler::globFiles(input2);

    vector<path> vecinput3;
    if ( !input3.empty() ){
        vecinput3 = InputHandler::globFiles(input3);
    }
    
    // Create the sequence coverage object
    Comp comp(vecinput1, vecinput2, vecinput3);
    comp.setOutputPrefix(output_prefix);
    comp.setD1Scale(d1_scale);
    comp.setD2Scale(d2_scale);
    comp.setD1Bins(d1_bins);
    comp.setD2Bins(d2_bins);
    comp.setIOThreads(io_threads);
    comp.setAnalysisThreads(analysis_threads);
    comp.setMerLen(mer_len);
    comp.setCanonical(0, canonical_1);
    comp.setCanonical(1, canonical_2);
    comp.setCanonical(2, canonical_3);
    comp.setHashSize(0, hash_size_1);
    comp.setHashSize(1, hash_size_2);
    comp.setHashSize(2, hash_size_3);
    comp.setDumpHashes(dump_hashes);
    comp.setDisableHashGrow(disable_hash_grow);
    comp.setDensityPlot(density_plot);
    comp.setVerbose(verbose);
    
    // Do the work
    comp.execute();

    // Save results to disk
    comp.save();
    
    // Plot results
    comp.plot(plot_output_type);
    
    // Send K-mer statistics to stdout as well
    comp.printCounters(cout);
    
    return 0;
}
