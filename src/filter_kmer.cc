
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <sstream>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <math.h>
#include <memory>
#include <mutex>
#include <thread>
#include <sys/ioctl.h>
using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::ostream;
using std::ofstream;
using std::stringstream;
using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;
using std::thread;

#include <boost/filesystem/operations.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include <kat/str_utils.hpp>
#include <kat/input_handler.hpp>
#include <kat/jellyfish_helper.hpp>
#include <kat/kat_fs.hpp>
using kat::InputHandler;
using kat::JellyfishHelper;
using kat::KatFS;

#include "filter_kmer.hpp"


string kat::filter::Counter::toString() const {
    stringstream ss;
    ss << distinct << " distinct; " << total << " total.";
    return ss.str();
}

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
    verbose = false;

    low_count = DEFAULT_FILT_KMER_LOW_COUNT;
    high_count = DEFAULT_FILT_KMER_HIGH_COUNT;
    low_gc = DEFAULT_FILT_KMER_LOW_GC;
    high_gc = DEFAULT_FILT_KMER_HIGH_GC;
    invert = DEFAULT_FILT_KMER_INVERT;
    separate = DEFAULT_FILT_KMER_SEPARATE;
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
    KatFS::ensureDirectoryExists(parentDir);

    // Either count or load input
    if (input.mode == InputHandler::InputHandler::InputMode::COUNT) {
        input.count(threads);
    }
    else {
        input.loadHeader();
        input.loadHash();
    }

    size_t size = input.header->size();
    unsigned int key_len = input.header->key_len();
    unsigned int val_len = input.header->val_len();


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

    file_header out_header;
    out_header.fill_standard();
    out_header.canonical(input.header->canonical());
    out_header.counter_len(input.header->counter_len());
    out_header.format(input.header->format());
    out_header.fpr(input.header->fpr());
    out_header.key_len(input.header->key_len());
    out_header.max_reprobe(input.header->max_reprobe());
    out_header.nb_hashes(input.header->nb_hashes());
    out_header.size(input.header->size());
    out_header.val_len(input.header->val_len());

    HashCounter* inCounter = new HashCounter(size, key_len, val_len, threads);
    inCounter->do_size_doubling(false);   // We know the size of the hash

    HashCounter* outCounter = separate ? new HashCounter(size, key_len, val_len, threads) : nullptr;
    if (separate) {
        outCounter->do_size_doubling(false);
    }

    // Resize all the counters to the requested number of threads
    all.resize(threads);
    in.resize(threads);
    out.resize(threads);

    // Do the work
    filter(*inCounter, *outCounter);

    // Merge (and print) results
    merge();

    // Output to disk
    path in_path(output_prefix.string() + "-in.jf" + lexical_cast<string>(input.merLen));
    path out_path(output_prefix.string() + "-out.jf" + lexical_cast<string>(input.merLen));

    // Dumping automatically destroys hash counters
    dump(in_path, inCounter, out_header);
    if (separate) {
        dump(out_path, outCounter, out_header);
    }
}

void kat::filter::FilterKmer::dump(path& out_path, HashCounter* hash, file_header& header) {
     // Remove anything that exists at the target location
    if (bfs::is_symlink(out_path) || bfs::exists(out_path)) {
        bfs::remove(out_path.c_str());
    }

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");
    cout << "Dumping hash to " << out_path.string() << " ...";
    cout.flush();

    JellyfishHelper::dumpHash(hash->ary(), header, threads, out_path);

    cout << " done.";
    cout.flush();
}

void kat::filter::FilterKmer::merge() {

    unique_ptr<Counter> all_counts = all.merge();
    unique_ptr<Counter> in_counts = in.merge();

    cout << "K-mers in input   : " << all_counts->toString() << endl
         << "K-mers to keep    : " << in_counts->toString() << endl;

    if (separate) {
        unique_ptr<Counter> out_counts = out.merge();
        cout << "K-mers to discard : " << out_counts->toString() << endl;
    }

    cout << endl;
}

void kat::filter::FilterKmer::filter(HashCounter& inCounter, HashCounter& outCounter) {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");

    cout << "Filtering kmers ...";
    cout.flush();

    vector<thread> t(threads);

    for(uint16_t i = 0; i < threads; i++) {
        t[i] = thread(&FilterKmer::filterSlice, this, i, std::ref(inCounter), std::ref(outCounter));
    }

    for(uint16_t i = 0; i < threads; i++){
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
    bool            invert;
    bool            separate;
    bool            non_canonical;
    uint16_t        mer_len;
    uint64_t        hash_size;
    bool            verbose;
    bool            help;

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);


    // Declare the supported options.
    po::options_description generic_options(FilterKmer::helpMessage(), w.ws_col);
    generic_options.add_options()
            ("output_prefix,o", po::value<path>(&output_prefix)->default_value("kat.filter.kmer"),
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
            ("invert,i", po::bool_switch(&invert)->default_value(false),
                "Whether to take k-mers outside region as selected content, rather than those inside.")
            ("separate,s", po::bool_switch(&separate)->default_value(false),
                "Whether to partition the k-mers into two sets, those inside region and those outside.  Works in combination with \"invert\".")
            ("non_canonical,N", po::bool_switch(&non_canonical)->default_value(false),
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
            ("inputs", po::value<std::vector<path>>(&inputs), "Path to the input file(s) to process.")
            ;

    // Positional option for the input bam file
    po::positional_options_description p;
    p.add("inputs", -1);

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
    filter.setCanonical(!non_canonical);
    filter.setInvert(invert);
    filter.setSeparate(separate);
    filter.setMerLen(mer_len);
    filter.setHashSize(hash_size);
    filter.setVerbose(verbose);

    // Do the work
    filter.execute();

    return 0;
}
