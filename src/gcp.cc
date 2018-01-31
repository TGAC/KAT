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

#include <stdint.h>
#include <iostream>
#include <math.h>
#include <memory>
#include <thread>
#include <vector>
#include <sys/ioctl.h>
using std::shared_ptr;
using std::make_shared;
using std::ostream;
using std::ofstream;
using std::thread;
using std::vector;

#include <boost/algorithm/string.hpp>
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
#include <kat/sparse_matrix.hpp>
#include <kat/input_handler.hpp>
#include <kat/pyhelper.hpp>
using kat::InputHandler;
using kat::HashLoader;
using kat::ThreadedSparseMatrix;
using kat::SparseMatrix;

#include "plot.hpp"
using kat::Plot;

#include "gcp.hpp"


kat::Gcp::Gcp(vector<path>& _inputs) {

    input.setMultipleInputs(_inputs);
    input.index = 1;
    outputPrefix = "kat-gcp";
    cvgScale = 1.0;
    cvgBins = 1000;
    threads = 1;
}

void kat::Gcp::execute() {

    // Validate input
    input.validateInput();

    // Create output directory
    path parentDir = bfs::absolute(outputPrefix).parent_path();
    KatFS::ensureDirectoryExists(parentDir);

    // Either count or load input
    if (input.mode == InputHandler::InputHandler::InputMode::COUNT) {
        input.count(threads);
    }
    else {
        input.loadHeader();
        input.loadHash();
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
        path outputPath(outputPrefix.string() + "-hash.jf" + lexical_cast<string>(input.merLen));
        input.dump(outputPath, threads);
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

    out << mme::KEY_TITLE << "K-mer coverage vs GC count plot for: " << input.fileName() << endl;
    out << mme::KEY_X_LABEL << input.merLen << "-mer frequency" << endl;
    out << mme::KEY_Y_LABEL << "GC count" << endl;
    out << mme::KEY_Z_LABEL << "# distinct " << input.merLen << "-mers" << endl;
    out << mme::KEY_NB_COLUMNS << mx.height() << endl;
    out << mme::KEY_NB_ROWS << mx.width() << endl;
    out << mme::KEY_MAX_VAL << mx.getMaxVal() << endl;
    out << mme::KEY_TRANSPOSE << "0" << endl;
    out << mme::KEY_KMER << input.merLen << endl;
    out << mme::KEY_INPUT_1 << input.pathString() << endl;
    out << mme::MX_META_END << endl;

    mx.printMatrix(out);
}

void kat::Gcp::analyse() {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");

    cout << "Analysing kmers in hash ...";
    cout.flush();

    vector<thread> t(threads);

    for(uint16_t i = 0; i < threads; i++) {
        t[i] = thread(&Gcp::analyseSlice, this, i);
    }

    for(uint16_t i = 0; i < threads; i++){
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

        // Count gs and cs
        uint16_t g_or_c = gcCount(kmer);

        // Apply scaling factor
        uint64_t cvg_pos = kmer_count == 0 ? 0 : ceil((double) kmer_count * cvgScale);

        if (cvg_pos > cvgBins)
            gcp_mx->incTM(th_id, g_or_c, cvgBins, 1);
        else
            gcp_mx->incTM(th_id, g_or_c, cvg_pos, 1);
    }
}

void kat::Gcp::plot(const string& output_type) {

    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");

    cout << "Creating plot ...";
    cout.flush();

    string kstr = lexical_cast<string>(this->getMerLen());

    string outputFile = outputPrefix.string() + ".mx." + output_type;

#ifdef HAVE_PYTHON
        vector<string> args;
        args.push_back("kat_plot_density.py");
        args.push_back(string("--output=") + outputFile);
        if (verbose) {
            args.push_back("--verbose");
        }
        args.push_back(outputPrefix.string() + ".mx");
        Plot::executePythonPlot(Plot::PlotMode::DENSITY, args);
#endif

    cout << " done.";
    cout.flush();
}

void kat::Gcp::analysePeaks() {
#ifdef HAVE_PYTHON
    cout << "Analysing peaks ... ";
    cout.flush();

    vector<string> args;
    args.push_back("kat_distanalysis.py");
    if (verbose) {
        args.push_back("--verbose");
    }
    args.push_back(outputPrefix.string() + ".mx");

    char* char_args[50];

    for(size_t i = 0; i < args.size(); i++) {
        char_args[i] = strdup(args[i].c_str());
    }

    PyHelper::getInstance().execute("kat_distanalysis.py", (int)args.size(), char_args);

    for(size_t i = 0; i < args.size(); i++) {
        free(char_args[i]);
    }

    cout << endl;
#endif
}

int kat::Gcp::main(int argc, char *argv[]) {

    vector<path>    inputs;
    path            output_prefix;
    uint16_t        threads;
    double          cvg_scale;
    uint16_t        cvg_bins;
    string          trim5p;
    bool            non_canonical;
    uint16_t        mer_len;
    uint64_t        hash_size;
    bool            dump_hash;
    string          plot_output_type;
    bool            verbose;
    bool            help;

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);


    // Declare the supported options.
    po::options_description generic_options(Gcp::helpMessage(), w.ws_col);
    generic_options.add_options()
            ("output_prefix,o", po::value<path>(&output_prefix)->default_value(path("kat-gcp")),
                "Path prefix for files generated by this program.")
            ("threads,t", po::value<uint16_t>(&threads)->default_value(1),
                "The number of threads to use")
            ("cvg_scale,x", po::value<double>(&cvg_scale)->default_value(1.0),
                "Number of bins for the gc data when creating the contamination matrix.")
            ("cvg_bins,y", po::value<uint16_t>(&cvg_bins)->default_value(1000),
                "Number of bins for the cvg data when creating the contamination matrix.")
            ("5ptrim", po::value<string>(&trim5p)->default_value("0"),
                "Ignore the first X bases from reads.  If more that one file is provided you can specify different values for each file by seperating with commas.")
            ("non_canonical,N", po::bool_switch(&non_canonical)->default_value(false),
                "If counting fast(a/q), store explicit kmer as found.  By default, we store 'canonical' k-mers, which means we count both strands.")
            ("mer_len,m", po::value<uint16_t>(&mer_len)->default_value(DEFAULT_MER_LEN),
                "The kmer length to use in the kmer hashes.  Larger values will provide more discriminating power between kmers but at the expense of additional memory and lower coverage.")
            ("hash_size,H", po::value<uint64_t>(&hash_size)->default_value(DEFAULT_HASH_SIZE),
                "If kmer counting is required for the input, then use this value as the hash size.  If this hash size is not large enough for your dataset then the default behaviour is to double the size of the hash and recount, which will increase runtime and memory usage.")
            ("dump_hash,d", po::bool_switch(&dump_hash)->default_value(false),
                        "Dumps any jellyfish hashes to disk that were produced during this run.")
            ("output_type,p", po::value<string>(&plot_output_type)->default_value(DEFAULT_GCP_PLOT_OUTPUT_TYPE),
                "The plot file type to create: png, ps, pdf.  Warning... if pdf is selected please ensure your gnuplot installation can export pdf files.")
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

    vector<string> d1_5ptrim_strs;
    vector<uint16_t> d1_5ptrim_vals;
    boost::split(d1_5ptrim_strs,trim5p,boost::is_any_of(","));
    for (auto& v : d1_5ptrim_strs) d1_5ptrim_vals.push_back(boost::lexical_cast<uint16_t>(v));

    auto_cpu_timer timer(1, "KAT GCP completed.\nTotal runtime: %ws\n\n");

    cout << "Running KAT in GCP mode" << endl
         << "------------------------" << endl << endl;

    // Create the sequence coverage object
    Gcp gcp(inputs);
    gcp.setThreads(threads);
    gcp.setCanonical(!non_canonical);
    gcp.setCvgBins(cvg_bins);
    gcp.setTrim(d1_5ptrim_vals);
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

#ifdef HAVE_PYTHON

    // Plot results
    gcp.plot(plot_output_type);

    // Run the distribution analysis script
    gcp.analysePeaks();

#endif

    return 0;
}
