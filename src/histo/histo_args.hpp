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

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::ostringstream;

namespace kat
{
    const uint64_t  DEFAULT_HISTO_LOW          = 1;
    const uint64_t  DEFAULT_HISTO_HIGH         = 10000;
    const uint64_t  DEFAULT_HISTO_INC          = 1;
    const uint32_t  DEFAULT_HISTO_THREADS      = 1;
    const bool      DEFAULT_HISTO_FULL         = false;
    const bool      DEFAULT_HISTO_BOTH_STRANDS = false;
    const char*     DEFAULT_HISTO_OUTPUT       = "kat.histo";

    const uint16_t  HISTO_MIN_ARGS = 1;


    class HistoArgs : public BaseArgs
    {
    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const char* usage() const
        {
            return "Usage: kat histo [options] <jellyfish_hash_file>\n";
        }

        const char* shortDescription() const
        {
            return "Create an histogram of k-mer occurrences in a sequence file";
        }

        const char* longDescription() const
        {
            return  "  Create an histogram with the number of k-mers having a given count. In bucket 'i' are tallied the k-mers\n" \
                    "  which have a count 'c' satisfying 'low+i*inc <= c < low+(i+1)'. Buckets in the output are labeled by\n" \
                    "  the low end point (low+i).\n\n" \
                    "  The last bucket in the output behaves as a catchall: it tallies all k-mers with a count greater or equal\n" \
                    "  to the low end point of this bucket.\n\n" \
                    "  This tool is very similar to the \"histo\" tool in jellyfish itself.  The primary difference being that\n" \
                    "  the output contains metadata that make the histogram easier for the user to plot.";
        }

        const string optionsDescription() const
        {
            ostringstream help_str;

            help_str << " -o, --output=path           Output file. (\"" << DEFAULT_HISTO_OUTPUT << "\")" << endl
                     << " -l, --low=uint64            Low count value of histogram (" << DEFAULT_HISTO_LOW << ")" << endl
                     << " -h, --high=uint64           High count value of histogram (" << DEFAULT_HISTO_HIGH << ")" << endl
                     << " -t, --threads=uint32        Number of threads (" << DEFAULT_HISTO_THREADS << ")" << endl
                     << " -f, --full                  Full histo. Don't skip count 0. (" << DEFAULT_HISTO_FULL << ")" << endl
                     << " -C, --both_strands          IMPORTANT: Whether the jellyfish hash contains K-mers produced for both" << endl
                     << "                             strands. If this is not set to the same value as was produced during jellyfish" << endl
                     << "                             counting then output from histo will be unpredicatable (" << DEFAULT_HISTO_BOTH_STRANDS << ").";

            return help_str.str();
        }

        vector<option>* longOptions()
        {
            static struct option long_options_array[] =
            {
                {"low",           required_argument,  0, 'l'},
                {"high",          required_argument,  0, 'h'},
                {"threads",       required_argument,  0, 't'},
                {"full",          no_argument,        0, 'f'},
                {"output",        required_argument,  0, 'o'},
                {"both_strands",  no_argument,        0, 'C'}
            };

            vector<option>* long_options = new vector<option>();

            for(uint8_t i = 0; i < 7; i++)
            {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        string shortOptions()
        {
            return "l:h:t:fo:C";
        }

        void setOption(int c, char* option_arg) {

            switch(c)
            {
            case 'l':
                low = atoi(optarg);
                break;
            case 'h':
                high = atoi(optarg);
                break;
            case 't':
                threads = atoi(optarg);
                break;
            case 'f':
                full = true;
                break;
            case 'o':
                output = optarg;
                break;
            case 'C':
                both_strands = true;
                break;
            }
        }

        void processRemainingArgs(const vector<string>& remaining_args)
        {
            db_path = remaining_args[0];
        }

        const char* currentStatus() const
        {
            ostringstream status;

            status  << "low: " << low << endl
                    << "high: " << high << endl
                    << "threads: " << threads << endl
                    << "full: " << full << endl
                    << "output: " << output << endl
                    << "db_path: " << db_path << endl
                    << "both_strands: " << both_strands << endl;

            return status.str().c_str();
        }

    public:

        uint64_t        low;
        uint64_t        high;
        uint32_t        threads;
        bool            full;
        bool            both_strands;
        string          output;
        string          db_path;

        HistoArgs() : BaseArgs(HISTO_MIN_ARGS),
            low(DEFAULT_HISTO_LOW),
            high(DEFAULT_HISTO_HIGH),
            threads(DEFAULT_HISTO_THREADS),
            full(DEFAULT_HISTO_FULL),
            both_strands(DEFAULT_HISTO_BOTH_STRANDS),
            output(DEFAULT_HISTO_OUTPUT)
        { }

        HistoArgs(int argc, char* argv[]) : BaseArgs(HISTO_MIN_ARGS),
            low(DEFAULT_HISTO_LOW),
            high(DEFAULT_HISTO_HIGH),
            threads(DEFAULT_HISTO_THREADS),
            full(DEFAULT_HISTO_FULL),
            both_strands(DEFAULT_HISTO_BOTH_STRANDS),
            output(DEFAULT_HISTO_OUTPUT)
        { parse(argc, argv); }


        uint64_t calcBase()
        {
            return low > DEFAULT_HISTO_LOW ? (1 >= low ? DEFAULT_HISTO_LOW : low - 1) : DEFAULT_HISTO_LOW;
        }

        uint64_t calcCeil()
        {
            return high + 1;
        }
    };
}




