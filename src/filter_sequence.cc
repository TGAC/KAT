
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

#include <boost/filesystem/operations.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
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
#include "comp.hpp"




kat::filter::FilterSeq::FilterSeq(const path& _seq_file_1, const path& _seq_file_2, const path& _input) {
    this->seq_file_1 = _seq_file_1;
    this->seq_file_2 = _seq_file_2;
    vector<path> vecInput;
    vecInput.push_back(_input);
    init(vecInput);
}


kat::filter::FilterSeq::FilterSeq(const path& _seq_file_1, const path& _seq_file_2, const vector<path>& _input) {
    this->seq_file_1 = _seq_file_1;
    this->seq_file_2 = _seq_file_2;
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

    keepers = 0;
    total = 0;
}

void kat::filter::FilterSeq::execute() {

    if (bfs::is_symlink(seq_file_1)) {
        if (bfs::symbolic_link_exists(seq_file_1)) {
            seq_file_1 = bfs::canonical(seq_file_1);
        }
        else {
            BOOST_THROW_EXCEPTION(FilterSeqException() << FilterSeqErrorInfo(string(
                "Could not find input file at: ") + seq_file_1.string() + "; please check the path and try again."));
        }
    }

    if (!bfs::exists(seq_file_1)) {
        BOOST_THROW_EXCEPTION(FilterSeqException() << FilterSeqErrorInfo(string(
                "Could not find input file at: ") + seq_file_1.string() + "; please check the path and try again."));
    }

    if (this->isPaired()) {
        if (bfs::is_symlink(seq_file_2)) {
            if (bfs::symbolic_link_exists(seq_file_2)) {
                seq_file_2 = bfs::canonical(seq_file_2);
            }
            else {
                BOOST_THROW_EXCEPTION(FilterSeqException() << FilterSeqErrorInfo(string(
                    "Could not find sequence file at: ") + seq_file_2.string() + "; please check the path and try again."));
            }
        }

        if (!bfs::exists(seq_file_2)) {
            BOOST_THROW_EXCEPTION(FilterSeqException() << FilterSeqErrorInfo(string(
                    "Could not find sequence file at: ") + seq_file_2.string() + "; please check the path and try again."));
        }
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


    // Do the work
    processSeqFile();

    // Summary
    cout << "Found " << keepers << " / " << total << " to keep" << endl << endl;

}



void kat::filter::FilterSeq::processSeqFile() {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");

    cout << "Filtering sequences ..." << endl;

    // Temporary storage for sequence data
    reader = unique_ptr<seqan::SeqFileIn>(new seqan::SeqFileIn(seq_file_1.c_str()));

    if (this->isPaired()) {
        reader2 = unique_ptr<seqan::SeqFileIn>(new seqan::SeqFileIn(seq_file_2.c_str()));
    }

    // Setup output file for statistics and output header if requested
    if (doStats) {
        path stats_path_out(output_prefix.string() + ".stats");
        stats_stream = unique_ptr<ofstream>(new ofstream(stats_path_out.c_str()));
        (*stats_stream) << "index\tnb_bases\tnb_kmers\tnb_hits\tratio" << endl;
    }

    // Setup file paths
    path ext = seq_file_1.extension();

    path output_path_in(output_prefix.string() + ".in" + (this->isPaired() ? ".R1" : "") + ext.string());
    inWriter = unique_ptr<seqan::SeqFileOut>(new seqan::SeqFileOut(output_path_in.c_str()));
    if (separate) {
        path output_path_out(output_prefix.string() + ".out" + (this->isPaired() ? ".R1" : "") + ext.string());
        outWriter = unique_ptr<seqan::SeqFileOut>(new seqan::SeqFileOut(output_path_out.c_str()));
    }


    if (this->isPaired()) {
        path output_path_in2(output_prefix.string() + ".in.R2" + ext.string());
        inWriter2 = unique_ptr<seqan::SeqFileOut>(new seqan::SeqFileOut(output_path_in2.c_str()));
        if (separate) {
            path output_path_out2(output_prefix.string() + ".out.R2" + ext.string());
            outWriter2 = unique_ptr<seqan::SeqFileOut>(new seqan::SeqFileOut(output_path_out2.c_str()));
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> urd;

    // Processes sequences in batches of records to reduce memory requirements
    uint64_t index = 0;
    while (!seqan::atEnd(*reader)) {

        seqan::readRecord(name, seq, qual, *reader);

        if (this->isPaired()) {
            seqan::readRecord(name2, seq2, qual2, *reader2);
        }

        // Generate a random value for this sequence between 0 and 1 (we may use
        // this for subsampling later, if requested by the user)
        double val = urd(gen);

        // Process batch with worker threads
        // Process each sequence is processed in a different thread.
        // In each thread lookup each K-mer in the hash

        processSeq(index++, val);

        if (index % 100000 == 0) {
            cout << "Processed " << index << (this->isPaired() ? " pairs" : " entries") << endl;
        }

    }

    if (this->isPaired() && !seqan::atEnd(*reader2)) {
        BOOST_THROW_EXCEPTION(FilterSeqException() << FilterSeqErrorInfo(string(
                    "Second sequence file appears to be longer than the first.")));
    }

    seqan::close(*reader);

    seqan::close(*inWriter);
    if (separate) {
        seqan::close(*outWriter);
    }

    if (this->isPaired()) {
        seqan::close(*reader2);
        seqan::close(*inWriter2);
        if (separate) {
            seqan::close(*outWriter2);
        }
    }

    if (this->doStats) {
        stats_stream->close();
    }

    cout << "Finished filtering.";
    cout.flush();
}


void kat::filter::FilterSeq::processSeq(uint64_t index, double random_val) {


    vector<bool> kFound;
    this->getProfile(seq, kFound);

    if (this->isPaired()) {
        vector<bool> kFound2;
        this->getProfile(seq2, kFound2);
        kFound.insert(kFound.end(), kFound2.begin(), kFound2.end());
    }

    uint32_t nbFound = 0;
    for(const auto& b : kFound) {
        if (b) nbFound ++;
    }

    SeqStats stats(index, nbFound, kFound.size());

    double ratio = stats.calcRatio();


    bool keep = true;

    // Check to see if seq stats are within limits, if so keep the sequence
    if ((ratio >= threshold && !invert) || (invert && ratio < threshold)) {

        // Also check to see if we have exceeded the threshold for subsampling
        if (this->frequency > 0.0 && this->frequency < random_val) {
            keep = false;
        }
        else {
            // Increase keeper count and output sequence
            keepers++;
            seqan::writeRecord(*inWriter, name, seq, qual);
            if (this->isPaired()) {
                seqan::writeRecord(*inWriter2, name2, seq2, qual2);
            }
        }
    }
    else {
        keep = false;
    }

    // If the user's requested to seperate the dataset and we are not keeping
    // this record then output it to the discard file(s)
    if (separate && !keep) {
        seqan::writeRecord(*outWriter, name, seq, qual);
        if (this->isPaired()) {
            seqan::writeRecord(*outWriter2, name2, seq2, qual2);
        }
    }

    if (doStats) {
        size_t len = this->isPaired() ? seqan::length(seq) + seqan::length(seq2) : seqan::length(seq);

        (*stats_stream) << index << "\t" << len << "\t" << stats.nb_kmers
                << "\t" << stats.matches << "\t" << ratio << endl;
    }

    total++;
}

void kat::filter::FilterSeq::getProfile(seqan::CharString& sequence, vector<bool>& hits) {
    // There's no substring functionality in SeqAn in this version (2.0.0).  So we'll just
    // use regular c++ string's for this bit.  This conversion of strings:
    // {CharString -> c++ string -> substring -> jellyfish mer_dna} is
    // inefficient. Reducing the number of conversions necessary will
    // make a big performance improvement here
    stringstream ssSeq;
    ssSeq << sequence;
    string s = ssSeq.str();

    uint64_t seqLength = s.length();
    uint64_t nbCounts = seqLength - input.merLen + 1;
    uint64_t nbInvalid = 0;


    if (nbCounts <= 0) {

        // Can't analyse this sequence because it's too short
        //cerr << names[index] << ": " << seq << " is too short to compute coverage.  Sequence length is "
        //       << seqLength << " and K-mer length is " << merLen << ". Setting sequence coverage to 0." << endl;

    } else {

        for (uint64_t i = 0; i < nbCounts; i++) {

            string merstr = s.substr(i, input.merLen);

            // Jellyfish compacted hash does not support Ns so if we find one set this kmer to false
            if (!validKmer(merstr)) {
                hits.push_back(false);
                nbInvalid++;
            } else {
                mer_dna mer(merstr);
                uint64_t count = JellyfishHelper::getCount(input.hash, mer, input.canonical);
                hits.push_back(count > 0);
            }
        }
    }
}


int kat::filter::FilterSeq::main(int argc, char *argv[]) {

    vector<path>    inputs;
    path            seq_file_1;
    path            seq_file_2;
    path            output_prefix;
    uint16_t        threads;
    double          threshold;
    double          frequency;
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
                "Whether to output sequences not found in the kmer hash, rather than those inside.")
            ("separate,s", po::bool_switch(&separate)->default_value(false),
                "Whether to partition the sequences into two sets, those with k-mers detected and those without.  Works in combination with \"invert\".")
            ("seq", po::value<path>(&seq_file_1),
                "The sequence file to filter")
            ("seq2", po::value<path>(&seq_file_2),
                "The second sequence file to filter (use this if you want to filter paired end reads)")
            ("frequency,f", po::value<double>(&frequency)->default_value(DEFAULT_FILT_SEQ_FREQUENCY),
                "If a value is set here then only keep the sequence if matching the kmer dataset and a random number is generated between 0 and 1 that exceeds this threshold.  The default is 0.0 which means keep every hit.")
            ("stats", po::bool_switch(&stats)->default_value(false),
                "Whether to emit statistics about quantity of found k-mers in each sequence.  If the user specifies seq2, then each entry will represent both sequences combined.")
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

    if (seq_file_1.empty()) {
        BOOST_THROW_EXCEPTION(FilterSeqException() << FilterSeqErrorInfo(string(
                    "You must specify at least one sequence file to filter")));
    }

    auto_cpu_timer timer(1, "KAT filter seq completed.\nTotal runtime: %ws\n\n");

    cout << "Running KAT in filter sequence mode" << endl
         << "-----------------------------------" << endl << endl;


    // Create the sequence coverage object
    FilterSeq filter(seq_file_1, seq_file_2, inputs);
    filter.setThreshold(threshold);
    filter.setOutput_prefix(output_prefix);
    filter.setThreads(threads);
    filter.setCanonical(!non_canonical);
    filter.setInvert(invert);
    filter.setSeparate(separate);
    filter.setFrequency(frequency);
    filter.setDoStats(stats);
    filter.setMerLen(mer_len);
    filter.setHashSize(hash_size);
    filter.setVerbose(verbose);

    // Do the work
    filter.execute();

    return 0;
}
