
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

#include "inc/str_utils.hpp"

#include "input_handler.hpp"
#include "jellyfish_helper.hpp"
using kat::InputHandler;
using kat::JellyfishHelper;

#include "filter_kmer.hpp"


unique_ptr<kat::filter::Counter> kat::filter::ThreadedCounter::merge() {
    
    unique_ptr<kat::filter::Counter> merged( new Counter() );
    
    for(const auto& c : counter) {
        merged->distinct += c.distinct;
        merged->total += c.total;
    }
    
    return merged;
}

void kat::filter::ThreadedCounter::increment(const uint16_t th_id, const uint64_t total_inc) {
    
    counter[th_id].increment(total_inc);
}


kat::filter::FilterKmer::FilterKmer(const path& _input) {
    vector<path> vecInput;
    vecInput.push_back(_input);
    init(vecInput);
}
    

kat::filter::FilterKmer::FilterKmer(const vector<path>& _input) {
    
    init(_input);
}

void kat::filter::FilterKmer::init(const vector<path>& _input) {
    
    input.setMultipleInputs(_input);
    output_prefix = "kat.filter-kmer";
    
    threads = 1;
    input.canonical = false;
    merLen = DEFAULT_MER_LEN; 
    verbose = false; 
    
    low_count = DEFAULT_FILT_KMER_LOW_COUNT;
    high_count = DEFAULT_FILT_KMER_HIGH_COUNT;
    low_gc = DEFAULT_FILT_KMER_LOW_GC;
    high_gc = DEFAULT_FILT_KMER_HIGH_GC;
    invert = DEFAULT_FILT_KMER_INVERT;
    separate = DEFAULT_FILT_KMER_SEPARATE;
    
    all.resize(threads);
    in.resize(threads);
    out.resize(threads);
}

void kat::filter::FilterKmer::execute() {
    
    // Some validation first
    if(high_count < low_count) {
        BOOST_THROW_EXCEPTION(FilterKmerException() << FilterKmerErrorInfo(string(
                "High kmer count value must be >= to low kmer count value")));
    }
        
    if(high_gc < low_gc) {
        BOOST_THROW_EXCEPTION(FilterKmerException() << FilterKmerErrorInfo(string(
                "High GC count value must be >= to low GC count value")));
    }

    // Validate input
    input.validateInput();
    
    // Create output directory
    path parentDir = bfs::absolute(output_prefix).parent_path();
    if (!bfs::exists(parentDir) || !bfs::is_directory(parentDir)) {
        if (!bfs::create_directories(parentDir)) {
            BOOST_THROW_EXCEPTION(FilterKmerException() << FilterKmerErrorInfo(string(
                    "Could not create output directory: ") + parentDir.string()));
        }
    }
    // Either count or load input
    if (input.mode == InputHandler::InputHandler::InputMode::COUNT) {
        input.count(merLen, threads);
    }
    else {
        input.loadHeader();
        input.loadHash(true);                
    }
    
    size_t size = input.header->size();
    unsigned int key_len = input.header->key_len();
    unsigned int val_len = input.header->val_len();
    unsigned int max_reprobe_index = input.header->max_reprobe();
    
    
    if(verbose)
    {
        cerr << "Attempting to create output hash with the following settings: " << endl
             << " key length        = " << key_len << endl
             << " val length        = " << val_len << endl
             << " mer len           - " << key_len / 2 << endl
             << " hash size         = " << size << endl
             << " max reprobe index = " << input.header->max_reprobe() << endl
             << " nb mers           = " << input.header->nb_hashes() << endl << endl;
    }
    
    HashCounterPtr inCounter = make_shared<HashCounter>(size, key_len, val_len, threads);
    inCounter->do_size_doubling(false);   // We know the size of the hash
    
    HashCounterPtr outCounter = separate ? make_shared<HashCounter>(size, key_len, val_len, threads) : nullptr;
    if (separate) {
        outCounter->do_size_doubling(false);
    }
    
    // Do the work
    filter(*inCounter, *outCounter);

    // Output to disk
    path in_path(output_prefix.string() + "-in.jf" + lexical_cast<string>(merLen));
    path out_path(output_prefix.string() + "-out.jf" + lexical_cast<string>(merLen));
    
    dump(in_path, inCounter);    
    if (separate) {
        dump(out_path, outCounter);
    }    
}

void kat::filter::FilterKmer::dump(path& out_path, HashCounterPtr hash) {
     // Remove anything that exists at the target location
    if (bfs::is_symlink(out_path) || bfs::exists(out_path)) {
        bfs::remove(out_path.c_str());
    }

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n"); 
    cout << "Dumping hash to " << out_path.string() << " ...";
    cout.flush();

    JellyfishHelper::dumpHash(hash->ary(), *(input.header), threads, out_path);

    cout << " done.";
    cout.flush();
}

void kat::filter::FilterKmer::merge() {
    
    unique_ptr<Counter> all_counts = all.merge();
    unique_ptr<Counter> in_counts = in.merge();
    unique_ptr<Counter> out_counts = out.merge();    
    
    if (separate) {
       cout << "Distinct kmers in input      : " << all_counts->distinct << endl
            << "Distinct kmers in bounds     : " << in_counts->distinct << endl
            << "Distinct kmers out of bounds : " << out_counts->distinct << endl
            << "Total kmers in input         : " << all_counts->total << endl
            << "Total kmers in bounds        : " << in_counts->total << endl
            << "Total kmers out of bounds    : " << out_counts->total << endl; 
    }
}

void kat::filter::FilterKmer::filter(HashCounter& inCounter, HashCounter& outCounter) {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        

    cout << "Filtering kmers ...";
    cout.flush();

    thread t[threads];

    for(int i = 0; i < threads; i++) {
        t[i] = thread(&FilterKmer::filterSlice, this, i, std::ref(inCounter), std::ref(outCounter));
    }

    for(int i = 0; i < threads; i++){
        t[i].join();
    }

    cout << " done.";
    cout.flush();
}

void kat::filter::FilterKmer::filterSlice(int th_id, HashCounter& inCounter, HashCounter& outCounter) {
    
    LargeHashArray::region_iterator it = input.hash->region_slice(th_id, threads);
    while (it.next()) {        
        
        bool in_bounds = inBounds(it.key().to_str(), it.val());
        
        all.increment(th_id, it.val());
        
        // Logic to allocate kmer into correct counter
        if (!separate) {            
            if((in_bounds && !invert) || (!in_bounds && invert)) {
                inCounter.add(it.key(), it.val());
                in.increment(th_id, it.val());
            }
        }
        else {
            // We are just separating the kmers so update whichever hash is appropriate
            if (in_bounds) {
                inCounter.add(it.key(), it.val());
                in.increment(th_id, it.val());
            }
            else
            {
                outCounter.add(it.key(), it.val());
                out.increment(th_id, it.val());
            }
        }
    }

    inCounter.done();
    
    if (separate)
        outCounter.done();
}



bool kat::filter::FilterKmer::inBounds(const string& kmer_seq, const uint64_t& kmer_count) {
    
    // Calculate GC
    uint32_t gc_count = kat::gcCount(kmer_seq);

    // Are we within the limits
    bool in_gc_limits = low_gc <= gc_count && gc_count <= high_gc;

    // Are we within the limits
    bool in_cvg_limits = low_count <= kmer_count && kmer_count <= high_count;

    // Return if we want to keep this sequence or not based on the coverage (and, if we've got this far, GC)
    return in_gc_limits && in_cvg_limits;
}

int kat::filter::FilterKmer::main(int argc, char *argv[]) {

    vector<path>    inputs;
    path            output_prefix;
    uint16_t        threads;
    uint64_t        low_count;
    uint64_t        high_count;
    uint16_t        low_gc;
    uint16_t        high_gc;
    bool            canonical;
    uint16_t        mer_len;
    uint64_t        hash_size; 
    bool            verbose;
    bool            help;

    // Declare the supported options.
    po::options_description generic_options(FilterKmer::helpMessage(), 100);
    generic_options.add_options()
            ("output_prefix,o", po::value<path>(&output_prefix)->default_value("kat.hist"), 
                "Path prefix for files generated by this program.")
            ("threads,t", po::value<uint16_t>(&threads)->default_value(1),
                "The number of threads to use")
            ("low_count,c", po::value<uint64_t>(&low_count)->default_value(1),
                "Low count threshold")    
            ("high_count,d", po::value<uint64_t>(&high_count)->default_value(10000),
                "High count threshold")    
            ("low_gc,g", po::value<uint16_t>(&low_gc)->default_value(1),
                "Low GC count threshold") 
            ("high_gc,h", po::value<uint16_t>(&high_gc)->default_value(100),
                "Low GC count threshold") 
            ("canonical,C", po::bool_switch(&canonical)->default_value(false),
                "If counting fast(a/q) input, this option specifies whether the jellyfish hash represents K-mers produced for both strands (canonical), or only the explicit kmer found.")
            ("mer_len,m", po::value<uint16_t>(&mer_len)->default_value(DEFAULT_MER_LEN),
                "The kmer length to use in the kmer hashes.  Larger values will provide more discriminating power between kmers but at the expense of additional memory and lower coverage.")
            ("hash_size,H", po::value<uint64_t>(&hash_size)->default_value(DEFAULT_HASH_SIZE),
                "If kmer counting is required for the input, then use this value as the hash size.  If this hash size is not large enough for your dataset then the default behaviour is to double the size of the hash and recount, which will increase runtime and memory usage.")
            ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                "Print extra information.")
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
            ("inputs,i", po::value<std::vector<path>>(&inputs), "Path to the input file(s) to process.")
            ;

    // Positional option for the input bam file
    po::positional_options_description p;
    p.add("inputs", 100);

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



    auto_cpu_timer timer(1, "KAT filter kmer completed.\nTotal runtime: %ws\n\n");        

    cout << "Running KAT in filter kmer mode" << endl
         << "-------------------------------" << endl << endl;

    // Create the sequence coverage object
    FilterKmer filter(inputs);
    filter.setLow_count(low_count);
    filter.setHigh_count(high_count);
    filter.setLow_gc(low_gc);
    filter.setHigh_gc(high_gc);
    filter.setOutput_prefix(output_prefix);
    filter.setThreads(threads);
    filter.setCanonical(canonical);
    filter.setMerLen(mer_len);
    filter.setHashSize(hash_size);
    filter.setVerbose(verbose);

    // Do the work
    filter.execute();
    
    return 0;
}