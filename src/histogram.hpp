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

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <thread>
using std::shared_ptr;
using std::make_shared;
using std::thread;

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include <matrix/matrix_metadata_extractor.hpp>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish_helper.hpp>

typedef boost::error_info<struct HistogramError,string> HistogramErrorInfo;
struct HistogramException: virtual boost::exception, virtual std::exception { };

namespace kat {

    class Histogram : public jellyfish::thread_exec {
    private:
        
        // Arguments from user
        vector<path>    inputs;
        uint16_t        threads;
        uint64_t        low;
        uint64_t        high;
        bool            canonical;
        uint64_t        hash_size;            
        bool            verbose;

        // Jellyfish mapped file hash vars
        shared_ptr<JellyfishHelper> jfh;
        path hashFile;
        
        // Internal vars
        uint64_t base, ceil, inc, nb_buckets, nb_slices;
        uint64_t *data;
        uint64_t slice_id;

    public:

        Histogram(vector<path> _inputs) : inputs(_inputs) {
            
        }

        virtual ~Histogram() {
            if (data)
                delete [] data;
        }
        
        uint64_t getHash_size() const {
            return hash_size;
        }

        void setHash_size(uint64_t hash_size) {
            this->hash_size = hash_size;
        }

        uint64_t getHigh() const {
            return high;
        }

        void setHigh(uint64_t high) {
            this->high = high;
        }

        uint64_t getInc() const {
            return inc;
        }

        void setInc(uint64_t inc) {
            this->inc = inc;
        }

        vector<path> getInputs() const {
            return inputs;
        }

        void setInputs(vector<path> inputs) {
            this->inputs = inputs;
        }

        uint64_t getLow() const {
            return low;
        }

        void setLow(uint64_t low) {
            this->low = low;
        }

        uint16_t getThreads() const {
            return threads;
        }

        void setThreads(uint16_t threads) {
            this->threads = threads;
        }

        bool isVerbose() const {
            return verbose;
        }

        void setVerbose(bool verbose) {
            this->verbose = verbose;
        }


        void execute() {
            
            // Some validation first
            if (high < low) {
                BOOST_THROW_EXCEPTION(HistogramException() << HistogramErrorInfo(string(
                    "High count value must be >= to low count value.  High: ") + high + "; Low: " + low)); 
            }
            
            if (inputs.empty()) {
                BOOST_THROW_EXCEPTION(HistogramException() << HistogramErrorInfo(string(
                    "No input files provided"));
            }
            
            if (countsFiles.size() == 1 && !JellyfishHelper::isSequenceFile(countsFiles[0].extension())) {
                hashFile = countsFiles[0];                
            }
            else {
                
                for (path p : inputs) {
                    if (!JellyfishHelper::isSequenceFile(p.extension())) {
                        BOOST_THROW_EXCEPTION(HistogramException() << HistogramErrorInfo(string(
                            "You provided multiple sequence files to generate a kmer hash from, however some of the input files do not have a recognised sequence file extension: \".fa,.fasta,.fq,.fastq,.fna\""));
                    }
                    
                    // Check input files exist
                    if (!bfs::exists(p) && !bfs::symbolic_link_exists(p)) {
                        BOOST_THROW_EXCEPTION(HistogramException() << HistogramErrorInfo(string(
                            "Could not find input file at: ") + p.string() + "; please check the path and try again."));                        
                    }
                }
                
                cout << "Provided one or more sequence files.  Executing jellyfish to count kmers." << endl;
                
                hashFile = path(outputPrefix + string(".jf") + merLen;
                
                JellyfishHelper::jellyfishCount(countsFiles, canonical, hashFile, merLen, hashSize1, threads);
            }
            
            
            // Setup handles to load hashes
            jfh = make_shared<JellyfishHelper>(hashFile, AccessMethod::SEQUENTIAL);

            // Calculate other vars required for this run
            base = calcBase();
            ceil = calcCeil();
            nb_buckets = ceil + 1 - base;
            nb_slices = threads * 100;

            data = new uint64_t[threads * nb_buckets];
            memset(data, '\0', threads * nb_buckets * sizeof (uint64_t));

            std::ostream* out_stream = verbose ? &cerr : (std::ostream*)0;

            if (verbose)
                cerr << endl
                    << "Hash loaded successfully." << endl
                    << "Starting threads...";

            // Do the work
            startAndJoinThreads();

            if (verbose)
                cerr << "done." << endl;
        }
        
        void print(std::ostream &out) {
            // Output header
            out << mme::KEY_TITLE << "K-mer spectra for: " << hashFile << endl;
            out << mme::KEY_X_LABEL << "K" << jfh->getKeyLen() << " multiplicity: " << hashFile << endl;
            out << mme::KEY_Y_LABEL << "Number of distinct K" << jfh->getKeyLen() << " mers" << endl;
            out << mme::MX_META_END << endl;

            uint64_t col = base;
            for (uint64_t i = 0; i < nb_buckets; i++, col += inc) {
                uint64_t count = 0;

                for (uint16_t j = 0; j < threads; j++)
                    count += data[j * nb_buckets + i];

                out << col << " " << count << "\n";
            }
        }

    protected:
        
        uint64_t calcBase() {
            return low > DEFAULT_HIST_LOW ? (1 >= low ? DEFAULT_HIST_LOW : low - 1) : DEFAULT_HIST_LOW;
        }

        uint64_t calcCeil() {
            return high + 1;
        }
        
        void startAndJoinThreads() {
            
            thread t[threads];
            
            for(int i = 0; i < threads; i++) {
                t[i] = thread(&Histogram::start, this, i);
            }
            
            for(int i = 0; i < threads; i++){
                t[i].join();
            }
        }
         
        void start(int th_id) {
            uint64_t *hist = &data[th_id * nb_buckets];

            for (size_t i = slice_id++; i < nb_slices; i = slice_id++) {
                lha::region_iterator it = jfh->getSlice(i, nb_slices);
                while (it.next()) {
                    if (it.val() < base)
                        ++hist[0];
                    else if (it.val() > ceil)
                        ++hist[nb_buckets - 1];
                    else
                        ++hist[(it.val() - base) / inc];
                }
            }
        }
        
        static const string helpMessage() const {
            
            return string("Usage: kat hist [options] <jellyfish_hash>\n\n") +
                            "Create an histogram of k-mer occurrences in a sequence file.\n\n" +
                            "Create an histogram with the number of k-mers having a given count. In bucket 'i' are tallied the k-mers " \
                            "which have a count 'c' satisfying 'low+i*inc <= c < low+(i+1)'. Buckets in the output are labeled by the " \
                            "low end point (low+i). </br> " \
                            "The last bucket in the output behaves as a catchall: it tallies all k-mers with a count greater or equal to " \
                            "the low end point of this bucket. </br> " \
                            "This tool is very similar to the \"histo\" tool in jellyfish itself.  The primary difference being that the " \
                            "output contains metadata that make the histogram easier for the user to plot.";

        }
      
    public:
        
        static const uint64_t DEFAULT_LOW = 1;
        static const uint64_t DEFAULT_HIGH = 10000;
        static const uint64_t DEFAULT_INC = 1;
        static const uint16_t DEFAULT_THREADS = 1;
        static const char* DEFAULT_OUTPUT = "kat.hist";
        
        static int main(int argc, char *argv[]) {
            
            vector<path>    inputs;
            path            output;
            uint16_t        threads;
            uint64_t        low;
            uint64_t        high;
            bool            canonical;
            uint64_t        hash_size;            
            bool            verbose;
            bool            help;
        
            // Declare the supported options.
            po::options_description generic_options(Histogram::helpMessage());
            generic_options.add_options()
                    ("output,o", po::value<path>(&output)->default_value(DEFAULT_OUTPUT), 
                        "Path prefix for files generated by this program.")
                    ("threads,t", po::value<uint16_t>(&threads)->default_value(DEFAULT_THREADS),
                        "The number of threads to use")
                    ("low,l", po::value<uint64_t>(&low)->default_value(DEFAULT_LOW),
                        "Low count value of histogram")    
                    ("high,h", po::value<uint64_t>(&high)->default_value(DEFAULT_HIGH),
                        "High count value of histogram")    
                    ("canonical,c", po::bool_switch(&canonical)->default_value(false),
                        "Whether the jellyfish hashes contains K-mers produced for both strands.  If this is not set to the same value as was produced during jellyfish counting then output from sect will be unpredicatable.")
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
        
        

            auto_cpu_timer timer(1, "KAT HIST completed.\nTotal runtime: %ws\n\n");        

            cout << "Running KAT in HIST mode" << endl
                 << "------------------------" << endl << endl;
        
            // Create the sequence coverage object
            Histogram histo(inputs);
            histo.setThread(threads);
            histo.setLow(low);
            histo.setHigh(high);
            histo.setInc(inc);
            histo.setCanonical(canonical);
            histo.setHashSize(hash_size);
            histo.setVerbose(verbose);

            // Do the work
            histo.execute();

            // Output the results
            histo.print(out);

            // Close the output channel
            out.close();

            return 0;
        }
    };
}
