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

#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include <stdint.h>
#include <vector>

#include <common_args.hpp>
#include <str_utils.hpp>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

namespace kat {
    const double DEFAULT_D1_SCALE = 1.0;
    const double DEFAULT_D2_SCALE = 1.0;
    const uint16_t DEFAULT_THREADS = 1;
    const string DEFAULT_OUTPUT_PREFIX = "kat-comp";
    const uint16_t DEFAULT_D1_BINS = 1001;
    const uint16_t DEFAULT_D2_BINS = 1001;

    const uint16_t COMP_MIN_ARGS = 2;

    class CompArgs : public BaseArgs {
    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const string usage() const {
            return "Usage: kat comp [options] <jellyfish_hash_1> <jellyfish_hash_2> [<jellyfish_hash_3>]";
        }

        const string shortDescription() const {
            return "Compares jellyfish K-mer count hashes.";
        }

        const string longDescription() const {
            string long_desc = "The most common use case for this tool is to compare two jellyfish K-mer hashes.  The typical use case for " \
                               "this tool is to compare K-mers from two jellyfish hashes both representing K-mer counts for reads.  However, " \
                               "it is also common to compare K-mers generated from reads to those generated from an assembly. </br> " \
                               "If comparing K-mers from reads to K-mers from an assembly, the larger (most likely the read) K-mer hash " \
                               "should be provided first, then the assembly K-mer hash second. </br> " \
                               "The third optional jellyfish hash acts as a filter, restricting the analysis to the K-mers present on that " \
                               "set.  The manual contains more details on specific use cases.";

            return lineBreakString(long_desc, 78, "  ");
        }

        const string optionsDescription() const {
            ostringstream help_str;

            help_str << " -o, --output_prefix=string  Path prefix for files produced by this program (\"" << DEFAULT_OUTPUT_PREFIX << "\")" << endl
                    << " -x, --d1_scale=double       Scaling factor for the first dataset  - float multiplier (" << DEFAULT_D1_SCALE << ").  Max value: 1.0." << endl
                    << " -y, --d2_scale=double       Scaling factor for the second dataset - float multiplier (" << DEFAULT_D2_SCALE << ").  Max value: 1.0." << endl
                    << " -i, --d1_bins=uint16        Number of bins for the first dataset.  i.e. number of rows in the matrix (" << DEFAULT_D1_BINS << ")." << endl
                    << " -j, --d2_bins=uint16        Number of bins for the second dataset. i.e. number of columns in the matrix" << endl
                    << "                             (" << DEFAULT_D2_BINS << ")." << endl
                    << " -t, --threads=uint16        The number of threads to use (" << DEFAULT_THREADS << ")";

            return help_str.str();
        }

        vector<option>* longOptions() {
            static struct option long_options_array[] ={
                {"output_prefix", required_argument, 0, 'o'},
                {"d1_scale", required_argument, 0, 'x'},
                {"d2_scale", required_argument, 0, 'y'},
                {"d1_bins", required_argument, 0, 'i'},
                {"d2_bins", required_argument, 0, 'j'},
                {"threads", required_argument, 0, 't'},
            };

            vector<option>* long_options = new vector<option>();

            for (uint8_t i = 0; i < 6; i++) {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        string shortOptions() {
            return "o:x:y:i:j:t:";
        }

        void setOption(int c, string& option_arg) {

            switch (c) {
                case 'o':
                    output_prefix = string(option_arg);
                    break;
                case 't':
                    threads = strToInt16(option_arg);
                    break;
                case 'x':
                    d1_scale = strToDouble(option_arg);
                    if (d1_scale > 1.0)
                        d1_scale = 1.0;
                    break;
                case 'y':
                    d2_scale = strToDouble(option_arg);
                    if (d2_scale > 1.0)
                        d2_scale = 1.0;
                    break;
                case 'i':
                    d1_bins = strToInt16(option_arg);
                    break;
                case 'j':
                    d2_bins = strToInt16(option_arg);
                    break;
            }
        }

        void processRemainingArgs(const vector<string>& remaining_args) {
            db1_path = remaining_args[0];
            db2_path = remaining_args[1];

            if (remaining_args.size() >= 3)
                db3_path = remaining_args[2];
        }

        const string currentStatus() const {
            ostringstream status;

            status << "Threads requested: " << threads << endl
                    << "Dataset 1 Scaling Factor: " << d1_scale << endl
                    << "Dataset 2 Scaling Factor: " << d2_scale << endl
                    << "Number of Dataset 1 bins: " << d1_bins << endl
                    << "Number of Dataset 2 bins: " << d2_bins << endl
                    << "Jellyfish hash 1: " << db1_path << endl
                    << "Jellyfish hash 2: " << db2_path << endl
                    << "Jellyfish hash 3: " << db3_path << endl
                    << "Output file path prefix: " << output_prefix << endl;

            return status.str().c_str();
        }

    public:
        string db1_path;
        string db2_path;
        string db3_path;
        string output_prefix;
        double d1_scale;
        double d2_scale;
        uint16_t d1_bins;
        uint16_t d2_bins;
        uint16_t threads;

        // Default constructor

        CompArgs() : BaseArgs(COMP_MIN_ARGS),
        output_prefix(DEFAULT_OUTPUT_PREFIX), d1_scale(DEFAULT_D1_SCALE), d2_scale(DEFAULT_D2_SCALE),
        d1_bins(DEFAULT_D1_BINS), d2_bins(DEFAULT_D2_BINS),
        threads(DEFAULT_THREADS) {
        }

        // Constructor that parses command line options

        CompArgs(int argc, char* argv[]) : BaseArgs(COMP_MIN_ARGS),
        output_prefix(DEFAULT_OUTPUT_PREFIX), d1_scale(DEFAULT_D1_SCALE), d2_scale(DEFAULT_D2_SCALE),
        d1_bins(DEFAULT_D1_BINS), d2_bins(DEFAULT_D2_BINS),
        threads(DEFAULT_THREADS) {
            parse(argc, argv);
        }

        ~CompArgs() {
        }
    };
}

