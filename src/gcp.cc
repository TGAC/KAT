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

#include <stdint.h>
#include <iostream>
#include <math.h>
#include <memory>
#include <thread>
#include <vector>
using std::shared_ptr;
using std::make_shared;
using std::ostream;
using std::ofstream;
using std::thread;
using std::vector;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;
using boost::lexical_cast;

#include <jellyfish/mer_dna.hpp>
#include <jellyfish_helper.hpp>
using kat::HashLoader;

#include <matrix/matrix_metadata_extractor.hpp>
#include <matrix/threaded_sparse_matrix.hpp>

#include "input_handler.hpp"
using kat::InputHandler;

#include "gcp.hpp"

    
kat::Gcp::Gcp(vector<path>& _inputs) {
    
    input.input = _inputs;
    input.index = 1;
    outputPrefix = "kat-gcp";
    cvgScale = 1.0;
    cvgBins = 1000;
    merLen = DEFAULT_MER_LEN;
    threads = 1;
}
 
void kat::Gcp::execute() {

    // Validate input
    input.validateInput();
    
    // Setup output stream for jellyfish initialisation
    std::ostream* out_stream = verbose ? &cerr : (std::ostream*)0;

    // Either count or load input
    if (input.mode == InputHandler::InputHandler::InputMode::COUNT) {
        input.count(merLen, threads);
    }
    else {
        input.loadHeader();
        input.loadHash(true);                
    }
    
    // Create matrix of appropriate size (adds 1 to cvg bins to account for 0)
    gcp_mx = make_shared<ThreadedSparseMatrix>(input.header->key_len() / 2, cvgBins + 1, threads);

    // Process batch with worker threads
    // Process each sequence is processed in a different thread.
    // In each thread lookup each K-mer in the hash
    analyse();

    // Dump any hashes that were previously counted to disk if requested
    // NOTE: MUST BE DONE AFTER COMPARISON AS THIS CLEARS ENTRIES FROM HASH ARRAY!
    if (input.dumpHash) {
        path outputPath(outputPrefix.string() + "-hash.jf" + lexical_cast<string>(merLen));
        input.dump(outputPath, threads, true);     
    }
    
    // Merge results
    merge();
    
}
   
void kat::Gcp::save() {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        

    cout << "Saving results to disk ...";
    cout.flush();
    
    // Send main matrix to output file
    ofstream main_mx_out_stream(string(outputPrefix.string() + ".mx").c_str());
    printMainMatrix(main_mx_out_stream);
    main_mx_out_stream.close();
    
    cout << " done.";
    cout.flush();
}

void kat::Gcp::merge() {
    
    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        

    cout << "Merging matrices ...";
    cout.flush(); 
    gcp_mx->mergeThreadedMatricies();
    
    cout << "done.";
    cout.flush();
}

void kat::Gcp::printMainMatrix(ostream &out) {
    SM64 mx = gcp_mx->getFinalMatrix();

    out << mme::KEY_TITLE << "K-mer coverage vs GC count plot for: " << input.pathString() << endl;
    out << mme::KEY_X_LABEL << "K-mer multiplicity" << endl;
    out << mme::KEY_Y_LABEL << "GC count" << endl;
    out << mme::KEY_Z_LABEL << "Distinct K-mers per bin" << endl;
    out << mme::KEY_NB_COLUMNS << mx.height() << endl;
    out << mme::KEY_NB_ROWS << mx.width() << endl;
    out << mme::KEY_MAX_VAL << mx.getMaxVal() << endl;
    out << mme::KEY_TRANSPOSE << "0" << endl;
    out << mme::MX_META_END << endl;

    mx.printMatrix(out);
}

void kat::Gcp::analyse() {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        

    cout << "Analysing kmers in hash ...";
    cout.flush();
    
    thread t[threads];

    for(int i = 0; i < threads; i++) {
        t[i] = thread(&Gcp::analyseSlice, this, i);
    }

    for(int i = 0; i < threads; i++){
        t[i].join();
    }
    
    cout << "done.";
    cout.flush();
}

void kat::Gcp::analyseSlice(int th_id) {
   
    LargeHashArray::region_iterator it = input.hash->region_slice(th_id, threads);
    while (it.next()) {
        string kmer = it.key().to_str();
        uint64_t kmer_count = it.val();

        uint16_t g_or_c = 0;

        for (uint16_t i = 0; i < kmer.length(); i++) {
            char c = kmer[i];

            if (c == 'G' || c == 'g' || c == 'C' || c == 'c')
                g_or_c++;
        }

        // Apply scaling factor
        uint64_t cvg_pos = kmer_count == 0 ? 0 : ceil((double) kmer_count * cvgScale);

        if (cvg_pos > cvgBins)
            gcp_mx->incTM(th_id, g_or_c, cvgBins, 1);
        else
            gcp_mx->incTM(th_id, g_or_c, cvg_pos, 1);
    }
}

int kat::Gcp::main(int argc, char *argv[]) {

    vector<path>    inputs;
    path            output_prefix;
    uint16_t        threads;
    double          cvg_scale;
    uint16_t        cvg_bins;
    bool            canonical;
    uint16_t        mer_len;
    uint64_t        hash_size;
    bool            dump_hash;
    bool            verbose;
    bool            help;

    // Declare the supported options.
    po::options_description generic_options(Gcp::helpMessage(), 100);
    generic_options.add_options()
            ("output_prefix,o", po::value<path>(&output_prefix)->default_value(path("kat-gcp")), 
                "Path prefix for files generated by this program.")
            ("threads,t", po::value<uint16_t>(&threads)->default_value(1),
                "The number of threads to use")
            ("cvg_scale,x", po::value<double>(&cvg_scale)->default_value(1.0),
                "Number of bins for the gc data when creating the contamination matrix.")
            ("cvg_bins,y", po::value<uint16_t>(&cvg_bins)->default_value(1000),
                "Number of bins for the cvg data when creating the contamination matrix.")
            ("canonical,C", po::bool_switch(&canonical)->default_value(false),
                "IMPORTANT: Whether the jellyfish hashes contains K-mers produced for both strands.  If this is not set to the same value as was produced during jellyfish counting then output from sect will be unpredicatable.")
            ("mer_len,m", po::value<uint16_t>(&mer_len)->default_value(DEFAULT_MER_LEN),
                "The kmer length to use in the kmer hashes.  Larger values will provide more discriminating power between kmers but at the expense of additional memory and lower coverage.")
            ("hash_size,H", po::value<uint64_t>(&hash_size)->default_value(DEFAULT_HASH_SIZE),
                "If kmer counting is required for the input, then use this value as the hash size.  It is important this is larger than the number of distinct kmers in your set.  We do not try to merge kmer hashes in this version of KAT.")
            ("dump_hash,d", po::bool_switch(&dump_hash)->default_value(false), 
                        "Dumps any jellyfish hashes to disk that were produced during this run.") 
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



    auto_cpu_timer timer(1, "KAT GCP completed.\nTotal runtime: %ws\n\n");        

    cout << "Running KAT in GCP mode" << endl
         << "------------------------" << endl << endl;

    // Create the sequence coverage object
    Gcp gcp(inputs);
    gcp.setThreads(threads);
    gcp.setCanonical(canonical);
    gcp.setCvgBins(cvg_bins);
    gcp.setCvgScale(cvg_scale);
    gcp.setHashSize(hash_size);
    gcp.setMerLen(mer_len);
    gcp.setOutputPrefix(output_prefix);
    gcp.setDumpHash(dump_hash);
    gcp.setVerbose(verbose);

    // Do the work (outputs data to files as it goes)
    gcp.execute();

    // Save results
    gcp.save();
    
    return 0;
}
  