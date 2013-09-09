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

const double    DEFAULT_D1_SCALE        = 1.0;
const double    DEFAULT_D2_SCALE        = 1.0;
const uint16_t  DEFAULT_THREADS         = 1;
#define DEFAULT_OUTPUT_PREFIX "./kat_comp_output"
const uint16_t  DEFAULT_D1_BINS         = 1001;
const uint16_t  DEFAULT_D2_BINS         = 1001;
const bool      DEFAULT_BOTH_STRANDS    = false;

class CompArgs
{
public:
    const char * 	db1_path;
    const char * 	db2_path;
    const char * 	db3_path;

    const char *  output_prefix;
    double d1_scale;
    double d2_scale;
    uint16_t d1_bins;
    uint16_t d2_bins;
    uint16_t threads;
    bool both_strands;
    bool verbose;

    // Default constructor
    CompArgs() :
        output_prefix(DEFAULT_OUTPUT_PREFIX), d1_scale(DEFAULT_D1_SCALE), d2_scale(DEFAULT_D2_SCALE),
        d1_bins(DEFAULT_D1_BINS), d2_bins(DEFAULT_D2_BINS),
        threads(DEFAULT_THREADS), both_strands(DEFAULT_BOTH_STRANDS), verbose(false)
    {}

    // Constructor that parses command line options
    CompArgs(int argc, char* argv[]) :
        output_prefix(DEFAULT_OUTPUT_PREFIX), d1_scale(DEFAULT_D1_SCALE), d2_scale(DEFAULT_D2_SCALE),
        d1_bins(DEFAULT_D1_BINS), d2_bins(DEFAULT_D2_BINS),
        threads(DEFAULT_THREADS), both_strands(DEFAULT_BOTH_STRANDS), verbose(false)
    {
        parse(argc, argv);
    }




#define seqcvg_args_USAGE "Usage: kat comp [options] db1_path db2_path [db3_path]\n"
    const char * usage() const
    {
        return seqcvg_args_USAGE;
    }

    void error(const char *msg)
    {
        cerr << endl
             << "Error: " << msg << endl << endl
             << usage() << endl
             << "Use --help for more information" << endl << endl;
        exit(1);
    }


#define sect_args_HELP "Compares sets of jellyfish kmer counts.  Normally, this will be two sets of kmers, although " \
  "when comparing kmers from assemblies, it is possible to add a 3rd set of kmers which describe the ends of each assembly " \
  "sequence to determine issues at end of sequences that may have prevented the assembler joining sequences together.\n\n" \
  "If comparing kmers from reads to kmers from an assembly, the larger (most likely the read) kmer hash should be provided " \
  "first, then the assembly kmer hash second.  Finally, if comparing kmers at the ends of sequences this should be supplied last.\n\n" \
  "Options (default value in (), *required):\n" \
  " -o, --output_prefix=string  Path prefix for files produced by this program (" DEFAULT_OUTPUT_PREFIX ")\n" \
  " -x, --d1_scale=double       Scaling factor for the first dataset  - float multiplier (1.0).  Max value: 1.0\n" \
  " -y, --d2_scale=double       Scaling factor for the second dataset - float multiplier (1.0).  Max value: 1.0.\n" \
  " -i, --d1_bins=uint16        Number of bins for the first dataset.  i.e. number of rows in the matrix (1001)\n" \
  " -j, --d2_bins-uint16        Number of bins for the second dataset. i.e. number of columns in the matrix (1001)\n" \
  " -t, --threads=uint16        The number of threads to use (1)\n" \
  " -C, --both_strands          IMPORTANT: Whether the jellyfish hashes contains kmers produced for both strands.\n" \
  "                             If this is not set to the same value as was produced during jellyfish counting then output from sect will be unpredicatable.\n" \
  "                             Note that all hashes must be built with either both strands or with single strands.  Using a mix will also produce unpredictable results.\n" \
  " -v, --verbose               Outputs additional information to stderr (false)\n" \
  "     --usage                 Usage\n" \
  "     --help                  This message\n" \

    const char * help() const
    {
        return sect_args_HELP;
    }

#define sect_args_HIDDEN "Hidden options:"
    const char * hidden() const
    {
        return sect_args_HIDDEN;
    }


    void parse(int argc, char *argv[])
    {
        int c;
        int help_flag = 0;
        int usage_flag = 0;

        static struct option long_options[] =
        {
            {"verbose",         no_argument,        0, 'v'},
            {"output_prefix",   required_argument,  0, 'o'},
            {"d1_scale" ,       required_argument,  0, 'x'},
            {"d2_scale" ,       required_argument,  0, 'y'},
            {"d1_bins",         required_argument,  0, 'i'},
            {"d2_bins",         required_argument,  0, 'j'},
            {"threads",         required_argument,  0, 't'},
            {"both_strands",    required_argument,  0, 'C'},
            {"help",            no_argument,        &help_flag, 1},
            {"usage",           no_argument,        &usage_flag, 1},
            {0, 0, 0, 0}
        };

        static const char *short_options = "f:e:o:x:y:i:j:t:Cvuh";

        if (argc <= 1)
        {
            cerr << endl
                 << usage() << endl
                 << help() << endl;
            exit(1);
        }


        while (true)
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
            case 'v':
                verbose = true;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 't':
                threads = atoi(optarg);
                break;
            case 'x':
                d1_scale = atof(optarg);
                if (d1_scale > 1.0)
                    d1_scale = 1.0;
                break;
            case 'y':
                d2_scale = atof(optarg);
                if (d2_scale > 1.0)
                    d2_scale = 1.0;
                break;
            case 'i':
                d1_bins = atoi(optarg);
                break;
            case 'j':
                d2_bins = atoi(optarg);
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

        if(remaining_args < 2 || remaining_args > 3)
            error("Requires 2 or 3 arguments describing paths to jellyfish kmer counts.");

        db1_path = argv[optind++];
        db2_path = argv[optind++];
        db3_path = remaining_args == 3 ? argv[optind++] : NULL;
    }

    void print()
    {
        if (verbose)
            cerr << "Verbose flag set" << endl;

        if (threads)
            cerr << "Threads requested: " << threads << endl;

        if (d1_scale)
            cerr << "Dataset 1 Scaling Factor: " << d1_scale << endl;

        if (d2_scale)
            cerr << "Dataset 2 Scaling Factor: " << d2_scale << endl;

        if (d1_bins)
            cerr << "Number of Dataset 1 bins: " << d1_bins << endl;

        if (d2_bins)
            cerr << "Number of Dataset 2 bins: " << d2_bins << endl;

        if (db1_path)
            cerr << "Jellyfish hash 1: " << db1_path << endl;

        if (db2_path)
            cerr << "Jellyfish hash 2: " << db2_path << endl;

        if (db3_path)
            cerr << "Jellyfish hash 3: " << db3_path << endl;

        if (output_prefix)
            cerr << "Output file path prefix: " << output_prefix << endl;

        if (both_strands)
            cerr << "Jellyfish hash to be treated as containing double_stranded information." << endl;

        cerr << endl;
    }

private:
};


