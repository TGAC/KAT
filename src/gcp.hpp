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

#pragma once

#include <stdint.h>
#include <iostream>
#include <math.h>
#include <memory>
#include <thread>
#include <vector>
using std::shared_ptr;
using std::make_shared;
using std::thread;
using std::vector;

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include <jellyfish/mer_dna.hpp>
#include <jellyfish_helper.hpp>

#include <matrix/matrix_metadata_extractor.hpp>
#include <matrix/threaded_sparse_matrix.hpp>

using std::ostream;

typedef boost::error_info<struct GcpError,string> GcpErrorInfo;
struct GcpException: virtual boost::exception, virtual std::exception { };


namespace kat {

    class Gcp {
    private:

        // Input args
        vector<path>    inputs;
        uint16_t        threads;
        double          cvg_scale;
        uint16_t        cvg_bins;
        bool            canonical;
        uint64_t        hash_size;            
        bool            verbose;
        
        // 
        const uint64_t nb_slices;
        uint64_t slice_id;

        // Variables that live for the lifetime of this object
        shared_ptr<JellyfishHelper> jfh;
        
        // Stores results
        shared_ptr<ThreadedSparseMatrix> gcp_mx; // Stores cumulative base count for each sequence where GC and CVG are binned

    public:

        Gcp(vector<path>& _inputs, uint16_t _threads) :
        inputs(_inputs), threads(_threads), nb_slices(threads * 100) {
            
        }

        ~Gcp() {            
        }

        void execute() {
            
            if (inputs.empty()) {
                BOOST_THROW_EXCEPTION(GcpException() << GcpErrorInfo(string(
                    "No input files provided"));
            }
            
            for (path p : inputs) {
                if (!JellyfishHelper::isSequenceFile(p.extension())) {
                    BOOST_THROW_EXCEPTION(GcpException() << GcpErrorInfo(string(
                        "You provided sequence files to generate a kmer hash from, however some of the input files do not have a recognised sequence file extension: \".fa,.fasta,.fq,.fastq,.fna\""));
                }

                // Check input files exist
                if (!bfs::exists(countsFile) && !bfs::symbolic_link_exists(countsFile)) {
                    BOOST_THROW_EXCEPTION(GcpException() << GcpErrorInfo(string(
                        "Could not find input file at: ") + p.string() + "; please check the path and try again."));                        
                }
            }

            cout << "Provided one or more sequence files.  Executing jellyfish to count kmers." << endl;

            hashFile = path(outputPrefix + string(".jf") + merLen;

            JellyfishHelper::jellyfishCount(countsFiles, canonical, hashFile, merLen, hashSize1, threads);
            
            // Setup output stream for jellyfish initialisation
            std::ostream* out_stream = verbose ? &cerr : (std::ostream*)0;

            // Setup handle to jellyfish hash
            jfh = make_shared<JellyfishHelper>(hash, AccessMethod::SEQUENTIAL);
            
            // Create matrix of appropriate size (adds 1 to cvg bins to account for 0)
            gcp_mx = make_shared<ThreadedSparseMatrix>(jfh->getKeyLen(), cvg_bins + 1, threads);

            // Process batch with worker threads
            // Process each sequence is processed in a different thread.
            // In each thread lookup each K-mer in the hash
            startAndJoinThreads();

            // Merge the contamination matrix
            gcp_mx->mergeThreadedMatricies();
        }
        
        void printVars(ostream &out) {
            out << "GCP parameters:" << endl;
            out << " - Hash File Path: " << db_arg << endl;
            out << " - Hash: " << (jfh != nullptr ? "loaded" : "not loaded") << endl;
            out << " - Threads: " << threads << endl;
        }

        // Print K-mer comparison matrix

        void printMainMatrix(ostream &out) {
            SM64 mx = gcp_mx->getFinalMatrix();

            out << mme::KEY_TITLE << "K-mer coverage vs GC count plot for: " << db_arg << endl;
            out << mme::KEY_X_LABEL << "K-mer multiplicity" << endl;
            out << mme::KEY_Y_LABEL << "GC count" << endl;
            out << mme::KEY_Z_LABEL << "Distinct K-mers per bin" << endl;
            out << mme::KEY_NB_COLUMNS << mx->height() << endl;
            out << mme::KEY_NB_ROWS << mx->width() << endl;
            out << mme::KEY_MAX_VAL << mx->getMaxVal() << endl;
            out << mme::KEY_TRANSPOSE << "0" << endl;
            out << mme::MX_META_END << endl;

            mx->printMatrix(out);
        }
        
    protected:
        
         void startAndJoinThreads() {
            
            thread t[threads];
            
            for(int i = 0; i < threads; i++) {
                t[i] = thread(&Gcp::start, this, i);
            }
            
            for(int i = 0; i < threads; i++){
                t[i].join();
            }
        }

        void start(int th_id) {
            SM64 mx = gcp_mx->getThreadMatrix(th_id);

            for (size_t i = slice_id++; i < nb_slices; i = slice_id++) {
                lha::region_iterator it = jfh->getSlice(i, nb_slices);
                while (it.next()) {
                    string kmer = ((jellyfish::mer_dna)it.key()).to_str();
                    uint64_t kmer_count = it.val();

                    uint16_t g_or_c = 0;

                    for (uint16_t i = 0; i < kmer.length(); i++) {
                        char c = kmer[i];

                        if (c == 'G' || c == 'g' || c == 'C' || c == 'c')
                            g_or_c++;
                    }

                    // Apply scaling factor
                    uint64_t cvg_pos = kmer_count == 0 ? 0 : ceil((double) kmer_count * cvg_scale);

                    if (cvg_pos > cvg_bins)
                        mx->inc(g_or_c, cvg_bins, 1);
                    else
                        mx->inc(g_or_c, cvg_pos, 1);
                }
            }
        }

        const static string helpMessage() const {
             return string("Usage: kat gcp <jellyfish_hash>\n\n") +
                            "Compares GC content and K-mer coverage within a single jellyfish hash." +
                            "This tool takes a single jellyfish hash as input and then counts the GC nucleotides for each distinct K-mer " \
                            "in the hash.  For each GC count and K-mer coverage level, the number of distinct K-mers are counted and " \
                            "stored in a matrix.  This matrix can be used to analyse biological content within the hash.  For example, " \
                            "it can be used to distinguish legitimate content from contamination, or unexpected content.";
        }
        
    public:
        
        static const string DEFAULT_OUTPUT_PREFIX = "kat-gcp";
        static const uint16_t DEFAULT_THREADS = 1;
        static const double DEFAULT_CVG_SCALE = 1.0;
        static const uint16_t DEFAULT_CVG_BINS = 1000;
        static const uint64_t DEFAULT_HASH_SIZE = 10000000000;
        
        static int main(int argc, char *argv[]) {
            
            vector<path>    inputs;
            string          output_prefix;
            uint16_t        threads;
            double          cvg_scale;
            uint16_t        cvg_bins;
            bool            canonical;
            uint64_t        hash_size;            
            bool            verbose;
            bool            help;
        
            // Declare the supported options.
            po::options_description generic_options(Gcp::helpMessage());
            generic_options.add_options()
                    ("output_prefix,o", po::value<path>(&output_prefix)->default_value(DEFAULT_OUTPUT_PREFIX), 
                        "Path prefix for files generated by this program.")
                    ("threads,t", po::value<uint16_t>(&threads)->default_value(DEFAULT_THREADS),
                        "The number of threads to use")
                    ("cvg_scale,x", po::value<double>(&cvg_scale)->default_value(DEFAULT_CVG_SCALE),
                        "Number of bins for the gc data when creating the contamination matrix.")
                    ("cvg_bins,y", po::value<uint16_t>(&cvg_bins)->default_value(DEFAULT_CVG_BINS),
                        "Number of bins for the cvg data when creating the contamination matrix.")
                    ("canonical,c", po::bool_switch(&canonical)->default_value(false),
                        "IMPORTANT: Whether the jellyfish hashes contains K-mers produced for both strands.  If this is not set to the same value as was produced during jellyfish counting then output from sect will be unpredicatable.")
                    ("hash_size,s", po::value<uint64_t>(&hash_size)->default_value(DEFAULT_HASH_SIZE),
                        "If kmer counting is required for the input, then use this value as the hash size.  It is important this is larger than the number of distinct kmers in your set.  We do not try to merge kmer hashes in this version of KAT.")
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
            Gcp gcp(inputs, threads);

            // Output seqcvg parameters to stderr if requested
            if (verbose)
                gcp.printVars(cerr);

            // Do the work (outputs data to files as it goes)
            gcp.execute();

            // Send main matrix to output file
            std::ostringstream main_mx_out_path;
            main_mx_out_path << output_prefix << ".mx";
            ofstream_default main_mx_out_stream(main_mx_out_path.str().c_str(), cout);
            gcp.printMainMatrix(main_mx_out_stream);
            main_mx_out_stream.close();

            return 0;
        }
    };
}
