//  ********************************************************************
//  This file is part of KAT - the Kmer Analysis Toolkit.
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

using std::cout;
using std::cerr;
using std::endl;
using std::string;

namespace kat
{
    const uint64_t  DEFAULT_HISTO_LOW          = 1;
    const uint64_t  DEFAULT_HISTO_HIGH         = 10000;
    const uint64_t  DEFAULT_HISTO_INC          = 1;
    const uint32_t  DEFAULT_HISTO_THREADS      = 1;
    const bool      DEFAULT_HISTO_FULL         = false;
    const bool      DEFAULT_HISTO_BOTH_STRANDS = false;
    const char*     DEFAULT_HISTO_OUTPUT       = "";
    const bool      DEFAULT_HISTO_VERBOSE      = false;


    class HistoArgs
    {
    public:
        uint64_t        low;
        uint64_t        high;
        uint64_t        increment;
        uint32_t        threads;
        bool            full;
        bool            both_strands;
        const char*     output;
        bool            verbose;
        const char*     db_path;

        HistoArgs() :
            low(DEFAULT_HISTO_LOW),
            high(DEFAULT_HISTO_HIGH),
            increment(DEFAULT_HISTO_INC),
            threads(DEFAULT_HISTO_THREADS),
            full(DEFAULT_HISTO_FULL),
            both_strands(DEFAULT_HISTO_BOTH_STRANDS),
            output(DEFAULT_HISTO_OUTPUT),
            verbose(DEFAULT_HISTO_VERBOSE)
        { }

        HistoArgs(int argc, char* argv[]) :
            low(DEFAULT_HISTO_LOW),
            high(DEFAULT_HISTO_HIGH),
            increment(DEFAULT_HISTO_INC),
            threads(DEFAULT_HISTO_THREADS),
            full(DEFAULT_HISTO_FULL),
            both_strands(DEFAULT_HISTO_BOTH_STRANDS),
            output(DEFAULT_HISTO_OUTPUT),
            verbose(DEFAULT_HISTO_VERBOSE)
        { parse(argc, argv); }


#define histo_args_USAGE "Usage: kat histo [options] db:path"
        const char * usage() const { return histo_args_USAGE; }

        void error(const char *msg)
        {
            cerr << endl
                 << "Error: " << msg << endl << endl
                 << usage() << endl
                 << "Use --help for more information" << endl << endl;
            exit(1);
        }

#define histo_args_HELP "Create an histogram of k-mer occurrences\n\nCreate an histogram with the number of k-mers having a given\n" \
    "count. In bucket 'i' are tallied the k-mers which have a count 'c'\n" \
    "satisfying 'low+i*inc <= c < low+(i+1)*inc'. Buckets in the output are\n" \
    "labeled by the low end point (low+i*inc).\n" \
    "\n" \
    "The last bucket in the output behaves as a catchall: it tallies all\n" \
    "k-mers with a count greater or equal to the low end point of this\n" \
    "bucket.\n\n" \
    "Options (default value in (), *required):\n" \
    " -l, --low=uint64            Low count value of histogram (1)\n" \
    " -h, --high=uint64           High count value of histogram (10000)\n" \
    " -i, --increment=uint64      Increment value for buckets (1)\n" \
    " -t, --threads=uint32        Number of threads (1)\n" \
    " -f, --full                  Full histo. Don't skip count 0. (false)\n" \
    " -C, --both_strands          IMPORTANT: Whether the jellyfish hash contains kmers produced for both strands.\n" \
    "                             If this is not set to the same value as was produced during jellyfish counting then output from histo will be unpredicatable.\n" \
    " -o, --output=string         Output file\n" \
    " -v, --verbose               Outputs additional information to stderr (false)\n" \
    "     --usage                 Usage\n" \
    "     --help                  This message\n"

        const char * help() const { return histo_args_HELP; }

        void parse(int argc, char* argv[])
        {
            int c;
            int help_flag = 0;
            int usage_flag = 0;

            static struct option long_options[] =
            {
                {"low",           required_argument,  0, 'l'},
                {"high",          required_argument,  0, 'h'},
                {"increment",     required_argument,  0, 'i'},
                {"threads",       required_argument,  0, 't'},
                {"full",          no_argument,        0, 'f'},
                {"output",        required_argument,  0, 'o'},
                {"both_strands",  required_argument,  0, 'C'},
                {"verbose",       no_argument,        0, 'v'},
                {"help",          no_argument,        &help_flag, 1},
                {"usage",         no_argument,        &usage_flag, 1},
                {0, 0, 0, 0}
            };
            static const char *short_options = "l:h:i:t:fo:vC";

            while(true)
            {
                /* getopt_long stores the option index here. */
                int index = -1;

                c = getopt_long (argc, argv, short_options, long_options, &index);

                /* Detect the end of the options. */
                if (c == -1)
                    break;

                switch (c)
                {
                case ':':
                    cerr << "Missing required argument for "
                              << (index == -1 ? std::string(1, (char)optopt) : std::string(long_options[index].name))
                              << endl;
                    exit(1);
                case '?':
                    cerr << "Use --usage or --help for some help" << endl << endl;
                    exit(1);
                case 'l':
                    low = atoi(optarg);
                    break;
                case 'h':
                    high = atoi(optarg);
                    break;
                case 'i':
                    increment = atoi(optarg);
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
                case 'v':
                    verbose = true;
                    break;
                case 'C':
                    both_strands = true;
                    break;
                }
            }

            if (help_flag)
            {
                cout << usage() << endl
                     << help() << endl;
                exit(0);
            }

            if (usage_flag)
            {
                cout << usage() << endl
                     << "Use --help for more information." << endl << endl;
                exit(0);
            }


            // Parse arguments
            int remaining_args = argc - optind;

            if (verbose)
                cerr << "Found " << remaining_args << " remaining arguments on the command line." << endl;

            if(remaining_args != 1)
                error("Requires exactly 1 argument.");

            db_path = argv[optind++];
        }

        uint64_t calcBase()
        {
            return low > DEFAULT_HISTO_LOW ? (increment >= low ? DEFAULT_HISTO_LOW : low - increment) : DEFAULT_HISTO_LOW;
        }

        uint64_t calcCeil()
        {
            return high + increment;
        }


        void print()
        {
            cerr << "low: " << low << "\n";
            cerr << "high: " << high << "\n";
            cerr << "increment: " << increment << "\n";
            cerr << "threads: " << threads << "\n";
            cerr << "full: " << full << "\n";
            cerr << "output: " << output << "\n";
            cerr << "verbose: " << verbose << "\n";
            cerr << "db_path: " << db_path << "\n";
            cerr << "both_strands: " << both_strands << endl;
        }
    private:
    };
}




