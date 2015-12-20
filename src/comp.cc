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

#include <kat/matrix_metadata_extractor.hpp>
#include <kat/sparse_matrix.hpp>
#include <kat/threaded_sparse_matrix.hpp>
#include <kat/distance_metrics.hpp>
#include <kat/input_handler.hpp>
#include <kat/comp_counters.hpp>
using kat::JellyfishHelper;
using kat::InputHandler;
using kat::HashLoader;
using kat::CompCounters;
using kat::ThreadedCompCounters;

#include "plot.hpp"
#include "plot_spectra_cn.hpp"
#include "plot_density.hpp"
using kat::PlotSpectraCn;
using kat::PlotDensity;
using kat::Plot;


#include "comp.hpp"



    
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
    threads = 1;
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
    main_matrix = ThreadedSparseMatrix(d1Bins, d2Bins, threads);

    // Initialise extra matrices for hash3 (only allocates space if required)
    if (doThirdHash()) {

        ends_matrix = ThreadedSparseMatrix(d1Bins, d2Bins, threads);
        middle_matrix = ThreadedSparseMatrix(d1Bins, d2Bins, threads);
        mixed_matrix = ThreadedSparseMatrix(d1Bins, d2Bins, threads);
    }

    // Create the comp counters for each thread
    comp_counters = ThreadedCompCounters(
            input[0].getSingleInput(), 
            input[1].getSingleInput(), 
            input.size() == 3 ? input[2].getSingleInput() : path(),
            std::min(d1Bins, d2Bins));

    string merLenStr = lexical_cast<string>(this->getMerLen());

    // Count kmers in sequence files if necessary (sets load and hashes and hashcounters as appropriate)
    for(size_t i = 0; i < input.size(); i++) {
        if (input[i].mode == InputHandler::InputHandler::InputMode::COUNT) {
            input[i].count(threads);
        }
    }
    
    // Check to see if user specified any hashes to load
    bool anyLoad = false;
    bool allLoad = true;
    for(size_t i = 0; i < input.size(); i++) {
        if (input[i].mode == InputHandler::InputMode::LOAD) {
            input[i].loadHeader();
            anyLoad = true;            
        }
        else if (input[i].mode == InputHandler::InputHandler::InputMode::COUNT) {
            allLoad = false;
        }
    }
    
    // If all hashes are loaded directly there is no requirement that the user needs
    // to specify the merLen, so just set it to the merLen found in the header of the first input
    if (allLoad) this->setMerLen(input[0].header->key_len() / 2);
    
    for(uint16_t i = 0; i < input.size(); i++) {
        input[i].validateMerLen(this->getMerLen());
    }
    
    // Load any hashes if necessary
    if (anyLoad) loadHashes();
     
    
    // Run the threads
    compare();

    // Dump any hashes that were previously counted to disk if requested
    // NOTE: MUST BE DONE AFTER COMPARISON AS THIS CLEARS ENTRIES FROM HASH ARRAY!
    if (dumpHashes()) {
        for(uint16_t i = 0; i < input.size(); i++) {        
            path outputPath(outputPrefix.string() + "-hash" + lexical_cast<string>(input[i].index) + ".jf" + lexical_cast<string>(this->getMerLen()));
            input[i].dump(outputPath, threads);
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
    
    if (outputHists) {
        
        ofstream hist1_out_stream(string(outputPrefix.string() + ".1.hist").c_str());
        printHist(hist1_out_stream, input[0], comp_counters.getFinalMatrix().getSpectrum1());
        hist1_out_stream.close();
        
        ofstream hist2_out_stream(string(outputPrefix.string() + ".2.hist").c_str());
        printHist(hist2_out_stream, input[1], comp_counters.getFinalMatrix().getSpectrum2());
        hist2_out_stream.close();        
    }

    cout << " done.";
    cout.flush();
}

void kat::Comp::printHist(std::ostream &out, InputHandler& input, vector<uint64_t>& hist) {
    
    // Output header
    out << mme::KEY_TITLE << input.merLen << "-mer spectra for: " << input.pathString() << endl;
    out << mme::KEY_X_LABEL << input.merLen << "-mer frequency" << endl;
    out << mme::KEY_Y_LABEL << "# distinct " << input.merLen << "-mers" << endl;
    out << mme::MX_META_END << endl;

    for (uint64_t i = 0; i < hist.size(); i++) {
        out << i << " " << hist[i] << "\n";
    }
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
    if (threads > 1) {
        
        vector<thread> threads(input.size());
        
        void (kat::InputHandler::*memfunc)() = &kat::InputHandler::loadHash;
        
        for(size_t i = 0; i < input.size(); i++) {
            if (input[i].mode == InputHandler::InputMode::LOAD) {
                threads[i] = thread(memfunc, &input[i]);                
            }
        }
        
        for(size_t i = 0; i < input.size(); i++) {
            if (input[i].mode == InputHandler::InputMode::LOAD) {
                threads[i].join();                 
            }
        }        
    }
    else {
        for(size_t i = 0; i < input.size(); i++) {
            if (input[i].mode == InputHandler::InputMode::LOAD) {
               input[i].loadHash();                
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
    
    vector<thread> t(threads);

    for(uint16_t i = 0; i < threads; i++) {
        t[i] = thread(&Comp::compareSlice, this, i);
    }

    for(uint16_t i = 0; i < threads; i++){
        t[i].join();
    }
    
    cout << " done.";
    cout.flush();
}

void kat::Comp::compareSlice(int th_id) {

    shared_ptr<CompCounters> cc = make_shared<CompCounters>();

    // Setup iterator for this thread's chunk of hash1
    LargeHashArray::eager_iterator hash1Iterator = input[0].hash->eager_slice(th_id, threads);

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
    LargeHashArray::eager_iterator hash2Iterator = input[1].hash->eager_slice(th_id, threads);

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
        LargeHashArray::eager_iterator hash3Iterator = input[2].hash->eager_slice(th_id, threads);

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

    cout << "Creating plot(s) ...";
    cout.flush();
    
    bool res = true;
    string kstr = lexical_cast<string>(this->getMerLen());
    
    // Plot results
    if (densityPlot) {
        
        path outputFile = path(getMxOutPath().string() + ".density." + output_type);
        string xLabel = kstr + "-mer frequency for " + input[0].pathString();
        string yLabel = kstr + "-mer frequency for " + input[1].pathString();
        string zLabel = string("# distinct ") + kstr + "-mers";
        string title = string("Spectra Density Plot for: ") + input[0].pathString() + " vs " + input[1].pathString();
        
#if HAVE_PYTHON
        vector<string> args;
        args.push_back("kat_plot_density.py");
        args.push_back(string("--output=") + outputFile.string());
        args.push_back(string("--x_label=") + xLabel);
        args.push_back(string("--y_label=") + yLabel);
        args.push_back(string("--z_label=") + zLabel);
        args.push_back(string("--title=") + title);
        args.push_back(getMxOutPath().string());            
        Plot::executePythonPlot(Plot::PlotMode::DENSITY, args);
#elif HAVE_GNUPLOT
        PlotDensity pd(getMxOutPath(), outputFile);
        pd.setXLabel(xLabel);
        pd.setYLabel(yLabel);
        pd.setZLabel(zLabel);
        pd.setTitle(title);
        pd.setOutputType(output_type);
        res = pd.plot();
        
        if (!res) {
            cout << endl << "WARNING: gnuplot session not valid.  Probably gnuplot is not installed correctly on your system.  No plots produced.";
            return;
        }
        
#endif

    }
    else {
        
        path outputFile = path(getMxOutPath().string() + ".density." + output_type);
        string xLabel = kstr + "-mer frequency for " + input[0].pathString();
        string yLabel = string("# distinct ") + kstr + "-mers";
        string zLabel = kstr + "mer frequency for " + input[1].pathString();
        string title = string("Spectra CN Plot for: ") + input[0].pathString() + " vs " + input[1].pathString();
        
#if HAVE_PYTHON

        vector<string> args;
        args.push_back("kat_plot_spectra-cn.py");
        args.push_back(string("--output=") + outputFile.string());
        args.push_back(string("--x_label=") + xLabel);
        args.push_back(string("--y_label=") + yLabel);
        args.push_back(string("--z_label=") + zLabel);
        args.push_back(string("--title=") + title);
        args.push_back(getMxOutPath().string());
        Plot::executePythonPlot(Plot::PlotMode::SPECTRA_CN, args);
        
#elif HAVE_GNUPLOT
        PlotSpectraCn pscn(getMxOutPath(), path(getMxOutPath().string() + ".spectra-cn." + output_type));
        pscn.setTitle(title);
        pscn.setYLabel(yLabel);
        pscn.setXLabel(xLabel);
        pscn.setOutputType(output_type);
        res = pscn.plot();
        
        if (!res) {
            cout << endl << "WARNING: gnuplot session not valid.  Probably gnuplot is not installed correctly on your system.  No plots produced.";
            return;
        }       
#endif
    }
/*    
    if (outputHists) {
       
        path outputFile1 = path(outputPrefix.string() + ".1.hist." + output_type);
        path outputFile2 = path(outputPrefix.string() + ".2.hist." + output_type);
        string xLabel = kstr + "-mer frequency";
        string yLabel = "# distinct " + kstr + "-mers";
        string title1 = kstr + "-mer spectra for: " + input[0].pathString();
        string title2 = kstr + "-mer spectra for: " + input[1].pathString();
        
#if HAVE_PYTHON        
       
        vector<string> args1;
        args1.push_back("kat_plot_spectra-hist.py");
        args1.push_back(string("--output=") + outputFile1.string());
        args1.push_back(string("--x_label=") + xLabel);
        args1.push_back(string("--y_label=") + yLabel);
        args1.push_back(string("--title=") + title1);
        args1.push_back(outputPrefix.string() + ".1.hist");
        Plot::executePythonPlot(Plot::PlotMode::SPECTRA_HIST, args1);
        
        vector<string> args2;
        args2.push_back("kat_plot_spectra-hist.py");
        args2.push_back(string("--output=") + outputFile2.string());
        args2.push_back(string("--x_label=\"") + xLabel + "\"");
        args2.push_back(string("--y_label=\"") + yLabel + "\"");
        args2.push_back(string("--title=\"") + title2 + "\"");
        args2.push_back(outputPrefix.string() + ".2.hist");
        Plot::executePythonPlot(Plot::PlotMode::SPECTRA_HIST, args2);

#elif HAVE_GNUPLOT
        
        PlotSpectraHist psh(input[0].input, outputFile1);
        psh.setXLabel(xLabel);
        psh.setYLabel(yLabel);
        psh.setTitle(title1);
        psh.setOutputType(output_type);
        bool res1 = psh.plot(); 
        
        if (!res1) {
            cout << endl << "WARNING: gnuplot session not valid.  Probably gnuplot is not installed correctly on your system.  No plots produced.";
            return;
        }     
        
        PlotSpectraHist psh(input[1].input, outputFile2);
        psh.setXLabel(xLabel);
        psh.setYLabel(yLabel);
        psh.setTitle(title2);
        psh.setOutputType(output_type);
        bool res2 = psh.plot(); 
        
        if (!res2) {
            cout << endl << "WARNING: gnuplot session not valid.  Probably gnuplot is not installed correctly on your system.  No plots produced.";
            return;
        }     
        
#endif
    }
*/    
    cout << " done.";
    cout.flush();
}



int kat::Comp::main(int argc, char *argv[]) {

    string input1;
    string input2;
    string input3;
    string output_prefix;
    double d1_scale;
    double d2_scale;
    uint16_t d1_bins;
    uint16_t d2_bins;
    uint16_t threads;
    uint16_t mer_len;
    bool canonical_1;       // Deprecated... for removal in KAT 3.0
    bool canonical_2;       // Deprecated... for removal in KAT 3.0
    bool canonical_3;       // Deprecated... for removal in KAT 3.0
    bool non_canonical_1;
    bool non_canonical_2;
    bool non_canonical_3;
    uint64_t hash_size_1;
    uint64_t hash_size_2;
    uint64_t hash_size_3;
    bool dump_hashes;
    bool disable_hash_grow;
    bool density_plot;
    string plot_output_type;
    bool output_hists;
    bool verbose;
    bool help;

    // Declare the supported options.
    po::options_description generic_options(Comp::helpMessage(), 100);
    generic_options.add_options()
            ("output_prefix,o", po::value<string>(&output_prefix)->default_value("kat-comp"), 
                "Path prefix for files generated by this program.")
            ("threads,t", po::value<uint16_t>(&threads)->default_value(1),
                "The number of threads to use.")
            ("d1_scale,x", po::value<double>(&d1_scale)->default_value(1.0),
                "Scaling factor for the first dataset - float multiplier")
            ("d2_scale,y", po::value<double>(&d2_scale)->default_value(1.0),
                "Scaling factor for the second dataset - float multiplier")
            ("d1_bins,i", po::value<uint16_t>(&d1_bins)->default_value(1001),
                "Number of bins for the first dataset.  i.e. number of rows in the matrix")
            ("d2_bins,j", po::value<uint16_t>(&d2_bins)->default_value(1001),
                "Number of bins for the second dataset.  i.e. number of rows in the matrix")
            ("canonical1,C", po::bool_switch(&canonical_1)->default_value(false),
                "(DEPRECATED) If counting fast(a/q) for input 1, this option specifies whether the jellyfish hash represents K-mers produced for both strands (canonical), or only the explicit kmer found.")
            ("canonical2,D", po::bool_switch(&canonical_2)->default_value(false),
                "(DEPRECATED) If counting fast(a/q) for input 2, this option specifies whether the jellyfish hash represents K-mers produced for both strands (canonical), or only the explicit kmer found.")
            ("canonical3,E", po::bool_switch(&canonical_3)->default_value(false),
                "(DEPRECATED) If counting fast(a/q) for input 3, this option specifies whether the jellyfish hash represents K-mers produced for both strands (canonical), or only the explicit kmer found.")
            ("non_canonical_1,N", po::bool_switch(&non_canonical_1)->default_value(false),
                "If counting fast(a/q) for input 1, this option specifies whether the jellyfish hash represents K-mers produced for both strands (canonical), or only the explicit kmer found.")
            ("non_canonical_2,O", po::bool_switch(&non_canonical_2)->default_value(false),
                "If counting fast(a/q) for input 2, this option specifies whether the jellyfish hash represents K-mers produced for both strands (canonical), or only the explicit kmer found.")
            ("non_canonical_3,P", po::bool_switch(&non_canonical_3)->default_value(false),
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
            ("density_plot,n", po::bool_switch(&density_plot)->default_value(false),
                "Makes a spectra_mx plot.  By default we create a spectra_cn plot.")
            ("output_type,p", po::value<string>(&plot_output_type)->default_value(DEFAULT_COMP_PLOT_OUTPUT_TYPE), 
                "The plot file type to create: png, ps, pdf.  Warning... if pdf is selected please ensure your gnuplot installation can export pdf files.")
            ("output_hists,h", po::bool_switch(&output_hists)->default_value(false), 
                "Whether or not to output histogram data and plots for input 1 and input 2")
            ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                "Print extra information.")
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
            ("input_1,X", po::value<string>(&input1), "Path to the first input file.  Can be either FastA, FastQ or a jellyfish hash (non bloom filtered)")
            ("input_2,Y", po::value<string>(&input2), "Path to the second input file.  Can be either FastA, FastQ or a jellyfish hash (non bloom filtered)")
            ("input_3,Z", po::value<string>(&input3), "Path to the third input file.  Can be either FastA, FastQ or a jellyfish hash (non bloom filtered)")
            ;

    // Positional options for the input file groups
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

    // Glob input files
    vector<path> vecinput1, vecinput2;
    InputHandler::globFiles(input1, vecinput1);
    InputHandler::globFiles(input2, vecinput2);

    vector<path> vecinput3;
    if ( !input3.empty() ){
        InputHandler::globFiles(input3, vecinput3);
    }
    
    // Create the sequence coverage object
    Comp comp(vecinput1, vecinput2, vecinput3);
    comp.setOutputPrefix(path(output_prefix));
    comp.setD1Scale(d1_scale);
    comp.setD2Scale(d2_scale);
    comp.setD1Bins(d1_bins);
    comp.setD2Bins(d2_bins);
    comp.setThreads(threads);
    comp.setMerLen(mer_len);
    comp.setCanonical(0, non_canonical_1 ? non_canonical_1 : canonical_1 ? canonical_1 : true);
    comp.setCanonical(1, non_canonical_2 ? non_canonical_2 : canonical_2 ? canonical_2 : true);
    comp.setCanonical(2, non_canonical_3 ? non_canonical_3 : canonical_3 ? canonical_3 : true);
    comp.setHashSize(0, hash_size_1);
    comp.setHashSize(1, hash_size_2);
    comp.setHashSize(2, hash_size_3);
    comp.setDumpHashes(dump_hashes);
    comp.setDisableHashGrow(disable_hash_grow);
    comp.setDensityPlot(density_plot);
    comp.setOutputHists(output_hists);
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
