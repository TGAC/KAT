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
using std::string;
using std::vector;
using std::ostringstream;

namespace kat {
    const uint64_t DEFAULT_HIST_LOW = 1;
    const uint64_t DEFAULT_HIST_HIGH = 10000;
    const uint64_t DEFAULT_HIST_INC = 1;
    const uint16_t DEFAULT_HIST_THREADS = 1;
    const char* DEFAULT_HIST_OUTPUT = "kat.hist";

    const uint16_t HISTO_MIN_ARGS = 1;

    class HistArgs : public BaseArgs {
    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const string usage() const {
            return "Usage: kat hist [options] <jellyfish_hash>";
        }

        const string shortDescription() const {
            return "Create an histogram of k-mer occurrences in a sequence file.";
        }

        const string longDescription() const {
            string long_desc = "Create an histogram with the number of k-mers having a given count. In bucket 'i' are tallied the k-mers " \
                               "which have a count 'c' satisfying 'low+i*inc <= c < low+(i+1)'. Buckets in the output are labeled by the " \
                               "low end point (low+i). </br> " \
                               "The last bucket in the output behaves as a catchall: it tallies all k-mers with a count greater or equal to " \
                               "the low end point of this bucket. </br> " \
                               "This tool is very similar to the \"histo\" tool in jellyfish itself.  The primary difference being that the " \
                               "output contains metadata that make the histogram easier for the user to plot.";

            return lineBreakString(long_desc, 78, "  ");
        }

        const string optionsDescription() const {
            ostringstream help_str;

            help_str << " -o, --output=path           Output file. (\"" << DEFAULT_HIST_OUTPUT << "\")" << endl
                    << " -l, --low=uint64            Low count value of histogram (" << DEFAULT_HIST_LOW << ")" << endl
                    << " -h, --high=uint64           High count value of histogram (" << DEFAULT_HIST_HIGH << ")" << endl
                    << " -t, --threads=uint16        Number of threads (" << DEFAULT_HIST_THREADS << ")";

            return help_str.str();
        }

        vector<option>* longOptions() {
            static struct option long_options_array[] ={
                {"low", required_argument, 0, 'l'},
                {"high", required_argument, 0, 'h'},
                {"threads", required_argument, 0, 't'},
                {"output", required_argument, 0, 'o'},
            };

            vector<option>* long_options = new vector<option>();

            for (uint8_t i = 0; i < 4; i++) {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        string shortOptions() {
            return "l:h:t:o:";
        }

        void setOption(int c, string& option_arg) {

            switch (c) {
                case 'l':
                    low = strToInt64(optarg);
                    break;
                case 'h':
                    high = strToInt64(optarg);
                    break;
                case 't':
                    threads = strToInt16(optarg);
                    break;
                case 'o':
                    output = string(optarg);
                    break;
            }
        }

        void processRemainingArgs(const vector<string>& remaining_args) {
            db_path = remaining_args[0];
        }

        const string currentStatus() const {
            ostringstream status;

            status << "low: " << low << endl
                    << "high: " << high << endl
                    << "threads: " << threads << endl
                    << "output: " << output << endl
                    << "db_path: " << db_path << endl;

            return status.str().c_str();
        }

    public:

        uint64_t low;
        uint64_t high;
        uint16_t threads;
        string output;
        string db_path;

        HistArgs() : BaseArgs(HISTO_MIN_ARGS),
        low(DEFAULT_HIST_LOW),
        high(DEFAULT_HIST_HIGH),
        threads(DEFAULT_HIST_THREADS),
        output(DEFAULT_HIST_OUTPUT) {
        }

        HistArgs(int argc, char* argv[]) : BaseArgs(HISTO_MIN_ARGS),
        low(DEFAULT_HIST_LOW),
        high(DEFAULT_HIST_HIGH),
        threads(DEFAULT_HIST_THREADS),
        output(DEFAULT_HIST_OUTPUT) {
            parse(argc, argv);
        }

        virtual ~HistArgs() {
        }

        uint64_t calcBase() {
            return low > DEFAULT_HIST_LOW ? (1 >= low ? DEFAULT_HIST_LOW : low - 1) : DEFAULT_HIST_LOW;
        }

        uint64_t calcCeil() {
            return high + 1;
        }
    };
}




