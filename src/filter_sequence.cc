
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

#include <algorithm>
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

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <kat/str_utils.hpp>
#include <kat/input_handler.hpp>
#include <kat/jellyfish_helper.hpp>
#include <kat/kat_fs.hpp>
using kat::InputHandler;
using kat::JellyfishHelper;
using kat::KatFS;

#include "plot_density.hpp"
using kat::PlotDensity;

#include "filter_sequence.hpp"




unique_ptr<kat::filter::SeqFilterCounter> kat::filter::ThreadedSeqStatsCounters::merge() {
    
    unique_ptr<kat::filter::SeqFilterCounter> merged( new SeqFilterCounter() );
    
    for(const auto& c : counters) {        
        for(const auto& s : c.seq_stats) {
            merged->add(s);
        }
    }

    // Sort by index
    merged->sort();
        
    return merged;
}


kat::filter::FilterSeq::FilterSeq(const path& _seq_file, const path& _input) {
    this->seq_file = _seq_file;
    vector<path> vecInput;
    vecInput.push_back(_input);
    init(vecInput);
}
    

kat::filter::FilterSeq::FilterSeq(const path& _seq_file, const vector<path>& _input) {
    this->seq_file = _seq_file;
    init(_input);
}

void kat::filter::FilterSeq::init(const vector<path>& _input) {
    
    input.setMultipleInputs(_input);
    output_prefix = "kat.filter-kmer";
    
    threads = 1;
    input.canonical = false;
    verbose = false; 
    
    threshold = DEFAULT_FILT_SEQ_THRESHOLD;
    invert = DEFAULT_FILT_SEQ_INVERT;
    separate = DEFAULT_FILT_SEQ_SEPARATE;   
    doStats = false;
}

void kat::filter::FilterSeq::execute() {
    
    if (bfs::is_symlink(seq_file)) {
        if (bfs::symbolic_link_exists(seq_file)) {
            seq_file = bfs::canonical(seq_file);
        }
        else {
            BOOST_THROW_EXCEPTION(FilterSeqException() << FilterSeqErrorInfo(string(
                "Could not find input file at: ") + seq_file.string() + "; please check the path and try again."));
        }
    }

    if (!bfs::exists(seq_file)) {
        BOOST_THROW_EXCEPTION(FilterSeqException() << FilterSeqErrorInfo(string(
                "Could not find input file at: ") + seq_file.string() + "; please check the path and try again."));
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
    
    // Resize all the counters to the requested number of threads
    stats.resize(threads);
    
    // Do the core of the work here
    processSeqFile();
    
    // Merge (and print) results
    unique_ptr<SeqFilterCounter> merged_stats = stats.merge();
    
    cout << "Found " << merged_stats->calcKeepers(threshold) << " / " << merged_stats->size() << " to keep" << endl << endl;

    // Filter seqs
    save(merged_stats->seq_stats);
}


void kat::filter::FilterSeq::processSeqFile() {
    
    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");     
    
    cout << "Analysing sequences ...";
    cout.flush();
    
    bucket_size = BATCH_SIZE / threads;
    remaining = BATCH_SIZE % (bucket_size < 1 ? 1 : threads);
    
    // Setup space for storing output
    offset = 0;
    recordsInBatch = 0;

    names = seqan::StringSet<seqan::CharString>();
    seqs = seqan::StringSet<seqan::CharString>();
    
    // Open file, create RecordReader and check all is well
    seqan::SeqFileIn reader(seq_file.c_str());
    
    // Processes sequences in batches of records to reduce memory requirements
    while (!seqan::atEnd(reader)) {
        if (verbose)
            cout << "Loading Batch of sequences... ";

        seqan::clear(names);
        seqan::clear(seqs);

        seqan::readRecords(names, seqs, reader, BATCH_SIZE);

        recordsInBatch = seqan::length(names);

        if (verbose)
            cout << "Loaded " << recordsInBatch << " records.  Processing batch... ";

        // Process batch with worker threads
        // Process each sequence is processed in a different thread.
        // In each thread lookup each K-mer in the hash
        analyseBatch();

        // Increment batch management vars
        offset += recordsInBatch;

        if (verbose)
            cout << "done" << endl;
    }

    seqan::close(reader);
    
    cout << " done.";
    cout.flush();
}


void kat::filter::FilterSeq::analyseBatch() {

    vector<thread> t(threads);

    for(uint16_t i = 0; i < threads; i++) {
        t[i] = thread(&FilterSeq::analyseBatchSlice, this, i);
    }

    for(uint16_t i = 0; i < threads; i++){
        t[i].join();
    }            
}

void kat::filter::FilterSeq::analyseBatchSlice(int th_id) {
    // Check to see if we have useful work to do for this thread, return if not
    if (bucket_size < 1 && th_id >= recordsInBatch) {
        return;
    }
    
    for (size_t i = th_id; i < recordsInBatch; i += threads) {
        processSeq(i, th_id);
    }
}


void kat::filter::FilterSeq::processSeq(const size_t index, const uint16_t th_id) {

    // There's no substring functionality in SeqAn in this version (2.0.0).  So we'll just
    // use regular c++ string's for this bit.  This conversion of strings:
    // {CharString -> c++ string -> substring -> jellyfish mer_dna} is 
    // inefficient. Reducing the number of conversions necessary will 
    // make a big performance improvement here
    stringstream ssSeq;
    ssSeq << seqs[index];
    string seq = ssSeq.str();

    uint64_t seqLength = seq.length();
    uint64_t nbCounts = seqLength - input.merLen + 1;
    uint64_t nbInvalid = 0;
    
    vector<bool> kFound(nbCounts, 0);
    
    if (nbCounts <= 0) {

        // Can't analyse this sequence because it's too short
        //cerr << names[index] << ": " << seq << " is too short to compute coverage.  Sequence length is "
        //       << seqLength << " and K-mer length is " << merLen << ". Setting sequence coverage to 0." << endl;
                
    } else {

        for (uint64_t i = 0; i < nbCounts; i++) {

            string merstr = seq.substr(i, input.merLen);

            // Jellyfish compacted hash does not support Ns so if we find one set this kmer to false
            if (!validKmer(merstr)) {
                kFound[i] = false;
                nbInvalid++;
            } else {                
                mer_dna mer(merstr);
                uint64_t count = JellyfishHelper::getCount(input.hash, mer, input.canonical);
                kFound[i] = count > 0;
            }
        }
    }

    uint32_t nbFound = 0;
    for(const auto& b : kFound) {
        if (b) nbFound ++;
    }
    
    unique_ptr<SeqStats> k( new SeqStats(offset+index, nbFound, nbCounts) );        
    stats.add(th_id, std::move(k));
}

void kat::filter::FilterSeq::save(vector<shared_ptr<SeqStats> >& stats) {
    
    if (stats.empty()) {
        cout << "Nothing to filter" << endl << endl;
    }
    else {
    
        auto_cpu_timer timer(1, "  Time taken: %ws\n\n");     

        path ext = seq_file.extension();
        path output_path_in(output_prefix.string() + ".in" + ext.string());
        path output_path_out(output_prefix.string() + ".out" + ext.string());
        path stats_path_out(output_prefix.string() + ".stats");

        cout << "Saving kept sequences to " << output_path_in.string() << " ...";
        cout.flush();

        // Open file, create RecordReader and check all is well
        // Open file and create RecordReader.
        seqan::SeqFileIn reader(seq_file.c_str());
        seqan::SeqFileOut inWriter(output_path_in.c_str());
        seqan::SeqFileOut outWriter(output_path_out.c_str());
        
        ofstream stats_stream(stats_path_out.c_str());
        
        if (doStats) {
            stats_stream << "header\tnb_bases\tnb_kmers\tnb_hits\tratio" << endl;
        }

        // Processes sequences in batches of records to reduce memory requirements
        uint32_t i = 0;
        
        // Loop until end
        while (!atEnd(reader)) {

            // Read record
            seqan::CharString id;
            seqan::Dna5String seq;
            seqan::CharString qual;
            seqan::readRecord(id, seq, qual, reader);
            
            double ratio = stats[i]->calcRatio();

            if ((ratio >= threshold && !invert) || (invert && ratio < threshold)) {
                seqan::writeRecord(inWriter, id, seq, qual);
            }
            else if (separate) {
                seqan::writeRecord(outWriter, id, seq, qual);
            }

            if (doStats) {
                stats_stream << id << "\t" << seqan::length(seq) << "\t" << stats[i]->nb_kmers 
                        << "\t" << stats[i]->matches << "\t" << ratio << endl;
            }
            
            i++;
        }

        seqan::close(reader);
        seqan::close(inWriter);
        seqan::close(outWriter);
        stats_stream.close();  

        if (!separate) {
            boost::filesystem::remove(output_path_out);
        }
        
        if (!doStats) {
            boost::filesystem::remove(stats_path_out);
        }

        cout << " done.";
        cout.flush();
    }
}

int kat::filter::FilterSeq::main(int argc, char *argv[]) {

    vector<path>    inputs;
    path            seq_file;
    path            output_prefix;
    uint16_t        threads;
    double          threshold;
    bool            invert;
    bool            separate;
    bool            stats;
    bool            non_canonical;
    uint16_t        mer_len;
    uint64_t        hash_size;
    bool            verbose;
    bool            help;
    
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    

    // Declare the supported options.
    po::options_description generic_options(FilterSeq::helpMessage(), w.ws_col);
    generic_options.add_options()
            ("output_prefix,o", po::value<path>(&output_prefix)->default_value("kat.filter.kmer"), 
                "Path prefix for files generated by this program.")
            ("threads,t", po::value<uint16_t>(&threads)->default_value(1),
                "The number of threads to use")
            ("threshold,T", po::value<double>(&threshold)->default_value(DEFAULT_FILT_SEQ_THRESHOLD),
                "What percentage of the sequence needs to be covered with target k-mers to keep the sequence")
            ("invert,i", po::bool_switch(&invert)->default_value(false),
                "Whether to take k-mers outside region as selected content, rather than those inside.")
            ("separate,s", po::bool_switch(&separate)->default_value(false),
                "Whether to partition the k-mers into two sets, those inside region and those outside.  Works in combination with \"invert\".")
            ("stats", po::bool_switch(&stats)->default_value(false),
                "Whether to emit statistics about quantity of found k-mers in each sequence.")
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
            ("seq_file", po::value<path>(&seq_file), "Path to file containing sequences to filter.")
            ("inputs", po::value<std::vector<path>>(&inputs), "Path to the input file(s) to process.")
            ;

    // Positional option for the input bam file
    po::positional_options_description p;
    p.add("seq_file", 1);
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

    auto_cpu_timer timer(1, "KAT filter seq completed.\nTotal runtime: %ws\n\n");        

    cout << "Running KAT in filter sequence mode" << endl
         << "-----------------------------------" << endl << endl;

    // Create the sequence coverage object
    FilterSeq filter(seq_file, inputs);
    filter.setThreshold(threshold);
    filter.setOutput_prefix(output_prefix);
    filter.setThreads(threads);
    filter.setCanonical(!non_canonical);
    filter.setInvert(invert);
    filter.setSeparate(separate);
    filter.setDoStats(stats);
    filter.setMerLen(mer_len);
    filter.setHashSize(hash_size);
    filter.setVerbose(verbose);

    // Do the work
    filter.execute();
        
    return 0;
}
