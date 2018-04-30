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
#include <string.h>
#include <stdint.h>
#include <vector>
#include <math.h>
#include <memory>
#include <thread>
#include <sys/ioctl.h>
using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::stringstream;
using std::shared_ptr;
using std::make_shared;
using std::thread;
using std::ofstream;

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;
using boost::lexical_cast;

#include <jellyfish/mer_dna.hpp>

#include <kat/jellyfish_helper.hpp>
#include <kat/matrix_metadata_extractor.hpp>
#include <kat/kat_fs.hpp>
using kat::KatFS;

#include "plot.hpp"
using kat::Plot;

#include "cold.hpp"

kat::Cold::Cold(const vector<path> _reads_files, const path _asm_file) {
    reads.setMultipleInputs(_reads_files);
    reads.index=1;
    assembly.setSingleInput(_asm_file);
    assembly.index=1;
    outputPrefix = "kat-cold";
    gcBins = 1001;
    cvgBins = 1001;
    threads = 1;
    verbose = false;
}


void kat::Cold::execute() {

    bucket_size = BATCH_SIZE / threads;
    remaining = BATCH_SIZE % (bucket_size < 1 ? 1 : threads);

    // Validate input
    reads.validateInput();
    assembly.validateInput();

    // Create output directory
    path parentDir = bfs::absolute(outputPrefix).parent_path();
    KatFS::ensureDirectoryExists(parentDir);

    // Either count or load reads
    if (reads.mode == InputHandler::InputHandler::InputMode::COUNT) {
        reads.count(threads);
    }
    else {
        reads.loadHeader();
        reads.loadHash();
    }

    // Either count or load assembly
    if (assembly.mode == InputHandler::InputHandler::InputMode::COUNT) {
        assembly.count(threads);
    }
    else {
        assembly.loadHeader();
        assembly.loadHash();
    }

    // Do the core of the work here
    processSeqFile();

    // Dump any hashes that were previously counted to disk if requested
    // NOTE: MUST BE DONE AFTER COMPARISON AS THIS CLEARS ENTRIES FROM HASH ARRAY!
    if (this->dumpHashes()) {
        path readsOutputPath(outputPrefix.string() + "-reads_hash.jf" + lexical_cast<string>(reads.merLen));
        reads.dump(readsOutputPath, threads);
        path asmOutputPath(outputPrefix.string() + "-asm_hash.jf" + lexical_cast<string>(assembly.merLen));
        assembly.dump(asmOutputPath, threads);
    }
}


void kat::Cold::processSeqFile() {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");

    cout << "Calculating kmer coverage across sequences ...";
    cout.flush();

    // Setup space for storing output
    offset = 0;
    recordsInBatch = 0;

    names = seqan::StringSet<seqan::CharString>();
    seqs = seqan::StringSet<seqan::CharString>();

    // Open file, create RecordReader and check all is well
    seqan::SeqFileIn reader(assembly.pathString().c_str());

    // Setup output stream for jellyfish initialisation
    std::ostream* out_stream = verbose ? &cerr : (std::ostream*)0;


    // Setup output streams for files
    if (verbose)
        *out_stream << endl;

    // Average sequence coverage and GC% scores output stream
    ofstream cvg_gc_stream(string(outputPrefix.string() + "-stats.tsv").c_str());
    cvg_gc_stream << "seq_name\tread_median_cvg\tread_mean_cvg\tasm_cn\tgc%\tseq_length\tkmers_in_seq\tinvalid_kmers\t%_invalid\tnon_zero_kmers\t%_non_zero\t%_non_zero_corrected" << endl;

    // Processes sequences in batches of records to reduce memory requirements
    while (!seqan::atEnd(reader)) {
        if (verbose)
            *out_stream << "Loading Batch of sequences... ";

        seqan::clear(names);
        seqan::clear(seqs);

        seqan::readRecords(names, seqs, reader, BATCH_SIZE);

        recordsInBatch = seqan::length(names);

        if (verbose)
            *out_stream << "Loaded " << recordsInBatch << " records.  Processing batch... ";

        // Allocate memory for output produced by this batch
        createBatchVars(recordsInBatch);

        // Process batch with worker threads
        // Process each sequence is processed in a different thread.
        // In each thread lookup each K-mer in the hash
        analyseBatch();

        // Output stats
        printStatTable(cvg_gc_stream);

        // Remove any batch specific variables from memory
        destroyBatchVars();

        // Increment batch management vars
        offset += recordsInBatch;

        if (verbose)
            *out_stream << "done" << endl;
    }

    seqan::close(reader);

    cvg_gc_stream.close();

    cout << " done.";
    cout.flush();
}

void kat::Cold::analyseBatch() {

    vector<thread> t(threads);

    for(uint16_t i = 0; i < threads; i++) {
        t[i] = thread(&Cold::analyseBatchSlice, this, i);
    }

    for(uint16_t i = 0; i < threads; i++){
        t[i].join();
    }
}

void kat::Cold::analyseBatchSlice(int th_id) {
    // Check to see if we have useful work to do for this thread, return if not
    if (bucket_size < 1 && th_id >= recordsInBatch) {
        return;
    }

    //processInBlocks(th_id);
    processInterlaced(th_id);
}

void kat::Cold::destroyBatchVars() {
    medians->clear();
    means->clear();
    asmCns->clear();
    gcs->clear();
    lengths->clear();
    invalid->clear();
    percentInvalid->clear();
    nonZero->clear();
    percentNonZero->clear();
    percentNonZeroCorrected->clear();

}

void kat::Cold::createBatchVars(uint16_t batchSize) {
    medians = make_shared<vector<uint32_t>>(batchSize);
    means = make_shared<vector<double>>(batchSize);
    asmCns = make_shared<vector<uint32_t>>(batchSize);
    gcs = make_shared<vector<double>>(batchSize);
    lengths = make_shared<vector<uint32_t>>(batchSize);
    invalid = make_shared<vector<uint32_t>>(batchSize);
    percentInvalid = make_shared<vector<double>>(batchSize);
    nonZero = make_shared<vector<uint32_t>>(batchSize);
    percentNonZero = make_shared<vector<double>>(batchSize);
    percentNonZeroCorrected = make_shared<vector<double>>(batchSize);
}

double kat::Cold::gcCountToPercentage(int16_t count) {
    return count == -1 ? -0.1 : (((double)count / (double)this->getMerLen()) * 100.0);
}

void kat::Cold::printStatTable(std::ostream &out) {

    out << std::fixed << std::setprecision(5);

    for (uint32_t i = 0; i < recordsInBatch; i++) {
        out << names[i] << "\t"
            << (*medians)[i] << "\t"
            << (*means)[i] << "\t"
            << (*asmCns)[i] << "\t"
            << (*gcs)[i] << "\t"
            << (*lengths)[i] << "\t"
            << (*lengths)[i] - this->assembly.merLen + 1 << "\t"
            << (*invalid)[i] << "\t"
            << (*percentInvalid)[i] << "\t"
            << (*nonZero)[i] << "\t"
            << (*percentNonZero)[i] << "\t"
            << (*percentNonZeroCorrected)[i] << endl;
    }
}


// This method won't be optimal in most cases... Fasta files are normally sorted by length (largest first)
// So first thread will be asked to do more work than the rest

void kat::Cold::processInBlocks(uint16_t th_id) {
    size_t start = bucket_size < 1 ? th_id : th_id * bucket_size;
    size_t end = bucket_size < 1 ? th_id : start + bucket_size - 1;
    for (size_t i = start; i <= end; i++) {
        processSeq(i, th_id);
    }

    // Process a remainder if required
    if (th_id < remaining) {
        size_t rem_idx = (threads * bucket_size) + th_id;
        processSeq(rem_idx, th_id);
    }
}

// This method is probably makes more efficient use of multiple cores on a length sorted fasta file

void kat::Cold::processInterlaced(uint16_t th_id) {
    size_t start = th_id;
    size_t end = recordsInBatch;
    for (size_t i = start; i < end; i += threads) {
        processSeq(i, th_id);
    }
}



void kat::Cold::processSeq(const size_t index, const uint16_t th_id) {

    // There's no substring functionality in SeqAn in this version (2.0.0).  So we'll just
    // use regular c++ string's for this bit.  This conversion of strings:
    // {CharString -> c++ string -> substring -> jellyfish mer_dna} is
    // inefficient. Reducing the number of conversions necessary will
    // make a big performance improvement here
    stringstream ssSeq;
    ssSeq << seqs[index];
    string seq = ssSeq.str();

    uint64_t seqLength = seq.length();
    int64_t nbCounts = seqLength - reads.merLen + 1;
    double average_cvg = 0.0;
    uint64_t nbNonZero = 0;
    uint64_t nbInvalid = 0;

    if (nbCounts <= 0) {

        // Can't analyse this sequence because it's too short
        // Hopefully this kind of thing doesn't happen too often in an assembly
        //cerr << names[index] << ": " << seq << " is too short to compute coverage.  Sequence length is "
        //       << seqLength << " and K-mer length is " << merLen << ". Setting sequence coverage to 0." << endl;

        (*medians)[index] = 0;
        (*means)[index] = 0.0;
        (*asmCns)[index] = 0;

    } else {

        shared_ptr<vector<uint64_t>> readsCounts = make_shared<vector<uint64_t>>(nbCounts, 0);
        shared_ptr<vector<uint64_t>> asmCounts = make_shared<vector<uint64_t>>(nbCounts, 0);

        uint64_t sum = 0;

        for (int64_t i = 0; i < nbCounts; i++) {

            string merstr = seq.substr(i, reads.merLen);

            // Jellyfish compacted hash does not support Ns so if we find one set this mer count to 0
            if (!validKmer(merstr)) {
                (*readsCounts)[i] = 0;
                (*asmCounts)[i] = 0;
                nbInvalid++;
            } else {
                mer_dna mer(merstr);
                uint64_t readcount = JellyfishHelper::getCount(reads.hash, mer, reads.canonical);
                sum += readcount;
                uint64_t asmcount = JellyfishHelper::getCount(assembly.hash, mer, assembly.canonical);
                (*readsCounts)[i] = readcount;
                (*asmCounts)[i] = asmcount;
                if (readcount != 0) nbNonZero++;
            }
        }

        // Create a copy of the counts, and sort it first, then take median value
        vector<uint64_t> sortedSeqCounts = *readsCounts;
        std::sort(sortedSeqCounts.begin(), sortedSeqCounts.end());
        (*medians)[index] = (double)(sortedSeqCounts[sortedSeqCounts.size() / 2]);

        // Calculate the mean
        (*means)[index] = (double)sum / (double)nbCounts;

        // Create a copy of the counts, and sort it first, then take median value
        vector<uint64_t> sortedAsmCounts = *asmCounts;
        std::sort(sortedAsmCounts.begin(), sortedAsmCounts.end());
        (*asmCns)[index] = (double)(sortedAsmCounts[sortedAsmCounts.size() / 2]);
    }

    // Add length
    (*lengths)[index] = seqLength;
    (*nonZero)[index] = nbNonZero;
    (*percentNonZero)[index] = nbNonZero == 0 || nbCounts <= 0 ?
        0.0 :
        ((double)nbNonZero / (double)nbCounts) * 100.0;
    (*invalid)[index] = nbInvalid;
    (*percentInvalid)[index] = nbInvalid == 0 || nbCounts <= 0 ?
        0.0 :
        ((double)nbInvalid / (double)nbCounts) * 100.0;

    uint64_t notInvalid = nbCounts - nbInvalid;
    (*percentNonZeroCorrected)[index] = nbNonZero == 0 || notInvalid <= 0 ?
        0.0 :
        ((double)nbNonZero / (double)notInvalid) * 100.0;


    // Calc GC%
    uint64_t gs = 0;
    uint64_t cs = 0;
    uint64_t ns = 0;

    for (uint64_t i = 0; i < seqLength; i++) {
        char c = seq[i];

        if (c == 'G' || c == 'g')
            gs++;
        else if (c == 'C' || c == 'c')
            cs++;
        else if (c == 'N' || c == 'n')
            ns++;
    }

    double gc_perc = ((double) (gs + cs)) / ((double) (seqLength - ns));
    (*gcs)[index] = gc_perc;
}


void kat::Cold::plot(const string& output_type) {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");

    cout << "Creating plot ...";
    cout.flush();

    string kstr = lexical_cast<string>(this->getMerLen());

    string outputFile = outputPrefix.string() + "." + output_type;

#ifdef HAVE_PYTHON
        vector<string> args;
        args.push_back("kat/plot/cold.py");
        args.push_back(string("--output=") + outputFile);
        if (verbose) {
            args.push_back("--verbose");
        }
        args.push_back(outputPrefix.string() + "-stats.tsv");
        Plot::executePythonPlot(Plot::PlotMode::COLD, args);
#endif

    cout << " done.";
    cout.flush();
}



int kat::Cold::main(int argc, char *argv[]) {

    vector<path>    reads_files;
    path            asm_file;
    path            output_prefix;
    uint16_t        gc_bins;
    uint16_t        cvg_bins;
    uint16_t        threads;
    string          trim5p;
    uint16_t        mer_len;
    uint64_t        hash_size;
    bool            dump_hash;
    bool            disable_hash_grow;
    string          plot_output_type;
    bool            verbose;
    bool            help;

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);


    // Declare the supported options.
    po::options_description generic_options(Cold::helpMessage(), w.ws_col);
    generic_options.add_options()
            ("output_prefix,o", po::value<path>(&output_prefix)->default_value("kat-cold"),
                "Path prefix for files generated by this program.")
            ("gc_bins,x", po::value<uint16_t>(&gc_bins)->default_value(1001),
                "Number of bins for the gc data when creating the contamination matrix.")
            ("cvg_bins,y", po::value<uint16_t>(&cvg_bins)->default_value(1001),
                "Number of bins for the cvg data when creating the contamination matrix.")
            ("threads,t", po::value<uint16_t>(&threads)->default_value(1),
                "The number of threads to use")
            ("5ptrim", po::value<string>(&trim5p)->default_value("0"),
                "Ignore the first X bases from reads.  If more that one file is provided you can specify different values for each file by seperating with commas.")
            ("mer_len,m", po::value<uint16_t>(&mer_len)->default_value(DEFAULT_MER_LEN),
                "The kmer length to use in the kmer hashes.  Larger values will provide more discriminating power between kmers but at the expense of additional memory and lower coverage.")
            ("hash_size,H", po::value<uint64_t>(&hash_size)->default_value(DEFAULT_HASH_SIZE),
                "If kmer counting is required, then use this value as the hash size for the reads.  We assume the assembly should use half this value.  If this hash size is not large enough for your dataset then the default behaviour is to double the size of the hash and recount, which will increase runtime and memory usage.")
            ("dump_hashes,d", po::bool_switch(&dump_hash)->default_value(false),
                "Dumps any jellyfish hashes to disk that were produced during this run. Normally, this is not recommended, and will likely consume a significant amount of disk space.")
            ("disable_hash_grow,g", po::bool_switch(&disable_hash_grow)->default_value(false),
                "By default jellyfish will double the size of the hash if it gets filled, and then attempt to recount.  Setting this option to true, disables automatic hash growing.  If the hash gets filled an error is thrown.  This option is useful if you are working with large genomes, or have strict memory limits on your system.")
            ("output_type,p", po::value<string>(&plot_output_type)->default_value(DEFAULT_COLD_PLOT_OUTPUT_TYPE),
                "The plot file type to create: png, ps, pdf.")
            ("verbose,v", po::bool_switch(&verbose)->default_value(false),
                "Print extra information.")
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
            ("asm_file", po::value<path>(&asm_file), "Path to the sequnce file to analyse for kmer coverage.")
            ("reads_files", po::value<std::vector<path>>(&reads_files), "Path(s) to the input files containing kmer counts.")
            ;

    // Positional option for the input bam file
    po::positional_options_description p;
    p.add("asm_file", 1);
    p.add("reads_files", -1);


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

    vector<string> d1_5ptrim_strs;
    vector<uint16_t> d1_5ptrim_vals;
    boost::split(d1_5ptrim_strs,trim5p,boost::is_any_of(","));
    for (auto& v : d1_5ptrim_strs) d1_5ptrim_vals.push_back(boost::lexical_cast<uint16_t>(v));

    auto_cpu_timer timer(1, "KAT CoLD completed.\nTotal runtime: %ws\n\n");

    cout << "Running KAT in Cold mode" << endl
         << "------------------------" << endl << endl;

    // Create the sequence coverage object
    Cold cold(reads_files, asm_file);
    cold.setOutputPrefix(output_prefix);
    cold.setGcBins(gc_bins);
    cold.setCvgBins(cvg_bins);
    cold.setThreads(threads);
    cold.setReadsTrim(d1_5ptrim_vals);
    cold.setMerLen(mer_len);
    cold.setHashSize(hash_size);
    cold.setDumpHashes(dump_hash);
    cold.setVerbose(verbose);

    // Do the work (outputs data to files as it goes)
    cold.execute();

#ifdef HAVE_PYTHON
    cold.plot(plot_output_type);
#endif

    return 0;
}
