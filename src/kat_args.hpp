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
#include <string.h>
#include <stdlib.h>
#include <iostream>

using std::string;
using std::cerr;
using std::cout;
using std::endl;

class KatArgs
{
private:
    string  mode_arg;
    int     mode_argc;
    char**   mode_argv;

public:


    // Default constructor
    KatArgs()
    {}

    // Constructor that parses command line options
    KatArgs(int argc, char* argv[])
    {
        parse(argc, argv);
    }

    string getMode() {
        return mode_arg;
    }

    int getModeArgC() {
        return mode_argc;
    }

    char** getModeArgV() {
        return mode_argv;
    }




#define kat_args_USAGE "Usage: kat <mode>\n"
    const char * usage() const
    {
        return kat_args_USAGE;
    }

    void error(const char *msg)
    {
        cerr << endl
             << "Error: " << msg << endl << endl
             << usage() << endl
             << "Use --help for more information" << endl;
        exit(1);
    }


#define kat_args_HELP "The Kmer Analysis Toolkist (KAT) contains a number of tools that analyse jellyfish kmer hashes\n\n" \
  "First argument should be the tool/mode you wish to use:\n\n" \
  "   - sect: SEquence Coverage estimator Tool.  Estimates the coverage of each sequence in a fasta file using kmers from a jellyfish hash\n" \
  "   - comp: Kmer comparison tool.  Creates a matrix of shared kmers between two jellyfish hashes.\n" \
  "   - gcp:  Kmer GC Processor.  Creates a matrix of the number of kmers found given a GC count and a kmer count.\n" \
  "   - plot: Plotting tool.  Creates useful plots to visualise kmer distributions.  Requires gnuplot \n\n" \
  "Options:\n" \
  "     --usage                              Usage\n" \
  "     --help                               This message\n" \
  " -V, --version                            Version\n"

    const char * help() const
    {
        return kat_args_HELP;
    }

#define kat_args_HIDDEN "Hidden options:"
    const char * hidden() const
    {
        return kat_args_HIDDEN;
    }

    void print_version(std::ostream &os = std::cout) const
    {
#ifndef PACKAGE_NAME
#define PACKAGE_NAME "Kmer Analysis Toolkit (KAT)
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.2.0"
#endif
        os << PACKAGE_NAME << " V" << PACKAGE_VERSION << "\n";
    }

    bool validMode(char* mode_str)
    {
        return (strcmp(mode_str, "sect") == 0 ||
                strcmp(mode_str, "comp") == 0 ||
                strcmp(mode_str, "gcp") == 0 ||
                strcmp(mode_str, "plot") == 0) ?
                    true : false;
    }

    void parse(int argc, char *argv[])
    {
        int c;

        static struct option long_options[] =
        {
            {"help",  no_argument,         0, 'h'},
            {"usage", no_argument,         0, 'u'},
            {"version", no_argument,       0, 'V'},
            {0, 0, 0, 0}
        };

        static const char *short_options = "Vuh";


        if (argc <= 1) {
            cerr << endl
                 << usage() << endl
                 << help() << endl;
            exit(1);
        }
        else if (validMode(argv[1])) {

            mode_arg = argv[1];
            mode_argc = argc - 1;
            mode_argv = argv + 1;
        }
        else {

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
                case 'h':
                    cout << usage() << endl
                         << help() << endl;
                    exit(0);
                case 'u':
                    cout << usage() << endl
                         << "Use --help for more information." << endl << endl;
                    exit(0);
                case 'V':
                    print_version();
                    exit(0);
                case '?':
                    cerr << "Use --usage or --help for some help" << endl << endl;
                    exit(1);

                }
            }

            error("Invalid command line arguments passed to \"kat\"\n\n");
            exit(1);
        }
    }
};

