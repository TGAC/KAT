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

#include "jellyfish_helper.hpp"
#include "input_handler.hpp"
using kat::JellyfishHelper;
using kat::HashLoader;
using kat::InputHandler;


#include "comp.hpp"

// ********** CompCounters ***********

kat::CompCounters::CompCounters() : CompCounters("", "", "") {}

kat::CompCounters::CompCounters(const path& _hash1_path, const path& _hash2_path, const path& _hash3_path) :
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
}

void kat::CompCounters::updateHash1Counters(uint64_t hash1_count, uint64_t hash2_count) {
    hash1_total += hash1_count;
    hash1_distinct++;

    if (!hash2_count) {
        hash1_only_total += hash1_count;
        hash1_only_distinct++;
    }
}

void kat::CompCounters::updateHash2Counters(uint64_t hash1_count, uint64_t hash2_count) {
    hash2_total += hash2_count;
    hash2_distinct++;

    if (!hash1_count) {
        hash2_only_total += hash2_count;
        hash2_only_distinct++;
    }
}

void kat::CompCounters::updateHash3Counters(uint64_t hash3_count) {

    hash3_total += hash3_count;
    hash3_distinct++;
}

void kat::CompCounters::updateSharedCounters(uint64_t hash1_count, uint64_t hash2_count) {

    if (hash1_count && hash2_count) {
        shared_hash1_total += hash1_count;
        shared_hash2_total += hash2_count;
        shared_distinct++;
    }
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
}
     

// ******** ThreadedCompCounters *********


kat::ThreadedCompCounters::ThreadedCompCounters(const path& _hash1_path, const path& _hash2_path, const path& _hash3_path) {

    threaded_counters = make_shared<vector<shared_ptr<CompCounters>>>();
    final_matrix = make_shared<CompCounters>(_hash1_path, _hash2_path, _hash3_path);            
}
                
void kat::ThreadedCompCounters::printCounts(ostream &out) {
    final_matrix->printCounts(out);
}
        
void kat::ThreadedCompCounters::add(shared_ptr<CompCounters> cc) {
    cc->hash1_path = final_matrix->hash1_path;
    cc->hash2_path = final_matrix->hash2_path;
    cc->hash3_path = final_matrix->hash3_path;
    threaded_counters->push_back(cc);
}
        
void kat::ThreadedCompCounters::merge() {

    // Merge counters
    for (vector<shared_ptr<CompCounters>>::iterator itp = threaded_counters->begin();
            itp != threaded_counters->end(); ++itp) {

        shared_ptr<CompCounters> it = *itp;
        
        final_matrix->hash1_total += it->hash1_total;
        final_matrix->hash2_total += it->hash2_total;
        final_matrix->hash3_total += it->hash3_total;
        final_matrix->hash1_distinct += it->hash1_distinct;
        final_matrix->hash2_distinct += it->hash2_distinct;
        final_matrix->hash3_distinct += it->hash3_distinct;
        final_matrix->hash1_only_total += it->hash1_only_total;
        final_matrix->hash2_only_total += it->hash2_only_total;
        final_matrix->hash1_only_distinct += it->hash1_only_distinct;
        final_matrix->hash2_only_distinct += it->hash2_only_distinct;
        final_matrix->shared_hash1_total += it->shared_hash1_total;
        final_matrix->shared_hash2_total += it->shared_hash2_total;
        final_matrix->shared_distinct += it->shared_distinct;
    }
}

    
// ********* Comp **********

kat::Comp::Comp() :
    kat::Comp::Comp(path(), path()) {}

kat::Comp::Comp(const path& _input1, const path& _input2) : 
    kat::Comp::Comp(_input1, _input2, path()) {}
        
kat::Comp::Comp(const path& _input1, const path& _input2, const path& _input3) {
    input = !_input3.empty() ? vector<InputHandler>(3) : vector<InputHandler>(2);
    
    input[0].setSingleInput(_input1);
    input[1].setSingleInput(_input2);
    input[0].index = 1;
    input[1].index = 2;
    if (!_input3.empty()) {
        input[2].setSingleInput(_input3);
        input[2].index = 3;
    }
    
    outputPrefix = "kat-comp";
    d1Scale = 1.0;
    d2Scale = 1.0;
    d1Bins = 1001;
    d2Bins = 1001;
    threads = 1;
    merLen = DEFAULT_MER_LEN;    
    verbose = false;      
}

void kat::Comp::execute() {

    // Check input files exist and determine input mode
    for(uint16_t i = 0; i < input.size(); i++) {
        input[i].validateInput();
    }
    
    // Create the final K-mer counter matrices
    main_matrix = make_shared<ThreadedSparseMatrix>(d1Bins, d2Bins, threads);

    // Initialise extra matrices for hash3 (only allocates space if required)
    if (doThirdHash()) {

        ends_matrix = make_shared<ThreadedSparseMatrix>(d1Bins, d2Bins, threads);
        middle_matrix = make_shared<ThreadedSparseMatrix>(d1Bins, d2Bins, threads);
        mixed_matrix = make_shared<ThreadedSparseMatrix>(d1Bins, d2Bins, threads);
    }

    // Create the comp counters for each thread
    comp_counters = make_shared<ThreadedCompCounters>(
            input[0].getSingleInput(), 
            input[1].getSingleInput(), 
            input.size() == 3 ? input[2].getSingleInput() : path());

    std::ostream* out_stream = verbose ? &cerr : (std::ostream*)0;

    string merLenStr = lexical_cast<string>(merLen);

    // Count kmers in sequence files if necessary (sets load and hashes and hashcounters as appropriate)
    for(int i = 0; i < input.size(); i++) {
        if (input[i].mode == InputHandler::InputHandler::InputMode::COUNT) {
            input[i].count(merLen, threads);
        }
    }
    
    // Check to see if user specified any hashes to load
    bool anyLoad = false;
    bool allLoad = true;
    bool anyDump = false;
    for(InputHandler i : input) {
        if (i.mode == InputHandler::InputMode::LOAD) {
            anyLoad = true;            
        }
        else if (i.mode == InputHandler::InputHandler::InputMode::COUNT) {
            anyDump = true;
            allLoad = false;
        }
    }

    for(InputHandler i : input) {
        i.loadHeader();
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
    if (isDumpHashes()) {
        for(uint16_t i = 0; i < input.size(); i++) {        
            path outputPath(outputPrefix.string() + "-hash" + lexical_cast<string>(input[i].index) + ".jf" + lexical_cast<string>(merLen));
            input[i].dump(outputPath, threads, true);
        }    
    }    
    
    
    
    // Send main matrix to output file
    ofstream main_mx_out_stream(string(outputPrefix.string() + "_main.mx").c_str());
    printMainMatrix(main_mx_out_stream);
    main_mx_out_stream.close();

    // Output ends matricies if required
    if (doThirdHash()) {
        // Ends matrix
        ofstream ends_mx_out_stream(string(outputPrefix.string() + "_ends.mx").c_str());
        printEndsMatrix(ends_mx_out_stream);
        ends_mx_out_stream.close();

        // Middle matrix
        ofstream middle_mx_out_stream(string(outputPrefix.string() + "_middle.mx").c_str());
        printMiddleMatrix(middle_mx_out_stream);
        middle_mx_out_stream.close();

        // Mixed matrix
        ofstream mixed_mx_out_stream(string(outputPrefix.string() + "_mixed.mx").c_str());
        printMixedMatrix(mixed_mx_out_stream);
        mixed_mx_out_stream.close();
    }

    // Send K-mer statistics to file
    ofstream stats_out_stream(string(outputPrefix.string() + ".stats").c_str());
    printCounters(stats_out_stream);
    stats_out_stream.close();

    // Send K-mer statistics to stdout as well
    printCounters(cout);

}

void kat::Comp::merge() {
    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        

    cout << "Merging data returned from each thread ...";
    cout.flush();

    // Merge results from the threads
    main_matrix->mergeThreadedMatricies();
    if (doThirdHash()) {
        ends_matrix->mergeThreadedMatricies();
        middle_matrix->mergeThreadedMatricies();
        mixed_matrix->mergeThreadedMatricies();
    }

    comp_counters->merge();

    cout << " done.";
    cout.flush();
}

void kat::Comp::loadHashes() {
    
    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        
    
    cout << "Loading hashes into memory...";
    cout.flush();    

    // If using parallel IO load hashes in parallel, otherwise do one at a time
    if (threads > 1) {
        
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

    SM64 mx = main_matrix->getFinalMatrix();

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

    ends_matrix->getFinalMatrix().printMatrix(out);
}

// Print K-mer comparison matrix

void kat::Comp::printMiddleMatrix(ostream &out) {

    out << "# Each row represents K-mer multiplicity for: " << input[0].getSingleInput().string() << endl;
    out << "# Each column represents K-mer multiplicity for sequence middles: " << input[1].getSingleInput().string() << endl;

    middle_matrix->getFinalMatrix().printMatrix(out);
}

// Print K-mer comparison matrix

void kat::Comp::printMixedMatrix(ostream &out) {

    out << "# Each row represents K-mer multiplicity for hash file 1: " << input[0].getSingleInput().string() << endl;
    out << "# Each column represents K-mer multiplicity for mixed: " << input[1].getSingleInput().string() << " and " << input[2].getSingleInput().string() << endl;

    mixed_matrix->getFinalMatrix().printMatrix(out);
}

// Print K-mer statistics

void kat::Comp::printCounters(ostream &out) {

    comp_counters->printCounts(out);
}

        
void kat::Comp::compare() {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        

    cout << "Comparing hashes ...";
    cout.flush();
    
    thread t[threads];

    for(int i = 0; i < threads; i++) {
        t[i] = thread(&Comp::compareSlice, this, i);
    }

    for(int i = 0; i < threads; i++){
        t[i].join();
    }
    
    cout << " done.";
    cout.flush();
}

void kat::Comp::compareSlice(int th_id) {

    shared_ptr<CompCounters> cc = make_shared<CompCounters>();

    // Setup iterator for this thread's chunk of hash1
    LargeHashArray::region_iterator hash1Iterator = input[0].hash->region_slice(th_id, threads);

    // Go through this thread's slice for hash1
    while (hash1Iterator.next()) {
        // Get the current K-mer count for hash1
        uint64_t hash1_count = hash1Iterator.val();

        // Get the count for this K-mer in hash2 (assuming it exists... 0 if not)
        uint64_t hash2_count = JellyfishHelper::getCount(input[0].hash, hash1Iterator.key(), input[0].canonical);

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
        main_matrix->incTM(th_id, scaled_hash1_count, scaled_hash2_count, 1);

        // Update hash 3 related matricies if hash 3 was provided
        if (doThirdHash()) {
            if (scaled_hash2_count == scaled_hash3_count)
                ends_matrix->incTM(th_id, scaled_hash1_count, scaled_hash3_count, 1);
            else if (scaled_hash3_count > 0)
                mixed_matrix->incTM(th_id, scaled_hash1_count, scaled_hash3_count, 1);
            else
                middle_matrix->incTM(th_id, scaled_hash1_count, scaled_hash3_count, 1);
        }
    }

    // Setup iterator for this thread's chunk of hash2
    // We setup hash2 for random access, so hopefully performance isn't too bad here...
    // Hash2 should be smaller than hash1 in most cases so hopefully we can get away with this.
    LargeHashArray::region_iterator hash2Iterator = input[1].hash->region_slice(th_id, threads);

    // Iterate through this thread's slice of hash2
    while (hash2Iterator.next()) {
        // Get the current K-mer count for hash2
        uint64_t hash2_count = hash2Iterator.val();

        // Get the count for this K-mer in hash1 (assuming it exists... 0 if not)
        uint64_t hash1_count = JellyfishHelper::getCount(input[0].hash, hash2Iterator.key(), input[0].hash);

        // Increment hash2's unique counters (don't bother with shared counters... we've already done this)
        cc->updateHash2Counters(hash1_count, hash2_count);

        // Only bother updating thread matrix with K-mers not found in hash1 (we've already done the rest)
        if (!hash1_count) {
            // Scale counters to make the matrix look pretty
            uint64_t scaled_hash2_count = scaleCounter(hash2_count, d2Scale);

            // Modifies hash counts so that K-mer counts larger than MATRIX_SIZE are dumped in the last slot
            if (scaled_hash2_count >= d2Bins) scaled_hash2_count = d2Bins - 1;

            // Increment the position in the matrix determined by the scaled counts found in hash1 and hash2
            main_matrix->incTM(th_id, 0, scaled_hash2_count, 1);
        }
    }

    // Only update hash3 counters if hash3 was provided
    if (doThirdHash()) {
        // Setup iterator for this thread's chunk of hash3
        LargeHashArray::region_iterator hash3Iterator = input[2].hash->region_slice(th_id, threads);

        // Iterate through this thread's slice of hash2
        while (hash3Iterator.next()) {
            // Get the current K-mer count for hash2

            uint64_t hash3_count = hash3Iterator.val();

            // Increment hash3's unique counters (don't bother with shared counters... we've already done this)
            cc->updateHash3Counters(hash3_count);
        }
    }

    mu.lock();
    comp_counters->add(cc);
    mu.unlock();
}

int kat::Comp::main(int argc, char *argv[]) {

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
    bool parallel_io;
    bool dump_hashes;
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
            ("canonical1,C", po::bool_switch(&canonical_1)->default_value(false),
                "Whether the jellyfish hash for input 1 contains K-mers produced for both strands.  If this is not set to the same value as was produced during jellyfish counting then output from sect will be unpredicatable.")
            ("canonical2,D", po::bool_switch(&canonical_2)->default_value(false),
                "Whether the jellyfish hash for input 2 contains K-mers produced for both strands.  If this is not set to the same value as was produced during jellyfish counting then output from sect will be unpredicatable.")
            ("canonical3,E", po::bool_switch(&canonical_3)->default_value(false),
                "Whether the jellyfish hash for input 3 contains K-mers produced for both strands.  If this is not set to the same value as was produced during jellyfish counting then output from sect will be unpredicatable.  Only applicable if you are using a third input.")
            ("mer_len,m", po::value<uint16_t>(&mer_len)->default_value(DEFAULT_MER_LEN),
                "The kmer length to use in the kmer hashes.  Larger values will provide more discriminating power between kmers but at the expense of additional memory and lower coverage.")
            ("hash_size_1,H", po::value<uint64_t>(&hash_size_1)->default_value(DEFAULT_HASH_SIZE),
                "If kmer counting is required for input 1, then use this value as the hash size.  It is important this is larger than the number of distinct kmers in your set.  We do not try to merge kmer hashes in this version of KAT.")
            ("hash_size_2,J", po::value<uint64_t>(&hash_size_2)->default_value(DEFAULT_HASH_SIZE),
                "If kmer counting is required for input 2, then use this value as the hash size.  It is important this is larger than the number of distinct kmers in your set.  We do not try to merge kmer hashes in this version of KAT.")
            ("hash_size_3,K", po::value<uint64_t>(&hash_size_3)->default_value(DEFAULT_HASH_SIZE),
                "If kmer counting is required for input 3, then use this value as the hash size.  It is important this is larger than the number of distinct kmers in your set.  We do not try to merge kmer hashes in this version of KAT.")
            ("dump_hashes,d", po::bool_switch(&dump_hashes)->default_value(false), 
                "Dumps any jellyfish hashes to disk that were produced during this run.")        
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
    comp.setCanonical(0, canonical_1);
    comp.setCanonical(1, canonical_2);
    comp.setCanonical(2, canonical_3);
    comp.setHashSize(0, hash_size_1);
    comp.setHashSize(1, hash_size_2);
    comp.setHashSize(2, hash_size_3);
    comp.setDumpHashes(dump_hashes);
    comp.setVerbose(verbose);
    
    // Do the work
    comp.execute();

    return 0;
}
