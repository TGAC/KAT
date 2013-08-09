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
#include <iostream>
#include <stdlib.h>

using std::string;
using std::cerr;
using std::cout;
using std::endl;

class PlotArgs
{
private:
    string  mode_arg;
    int     mode_argc;
    char**   mode_argv;
public:

    // Default constructor
    PlotArgs()
    {}

    // Constructor that parses command line options
    PlotArgs(int argc, char* argv[])
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

#define seqcvg_args_USAGE "Usage: kat plot <mode>\n"
    const char * usage() const
    {
        return seqcvg_args_USAGE;
    }

    void error(const char *msg)
    {
        cerr << endl
             << "Error: " << msg << endl << endl
             << usage() << endl
             << "Use --help for more information" << endl;
        exit(1);
    }


#define sect_args_HELP "Create Kmer Plots\n\n" \
  "First argument should be the plot mode you wish to use:\n\n" \
  "   - flame:           Creates a flame plot from a matrix created with the \"comp\" tool.  Typically this\n" \
  "                      is used to compare two kmer hashes produced by different NGS reads.\n" \
  "   - asm:             Creates a histogram using a matrix created with the \"comp\" tool.  Typically\n" \
  "                      this is used to compare a jellyfish hash produced from a read set to a jellyfish\n" \
  "                      hash produced from an assembly. The plot shows the amount of distinct kmers absent\n" \
  "                      in the assembly, those that are found once, and those found more (up to 10\n" \
  "                      duplications plotted)\n" \
  "   - sect:            Creates a kmer coverage plot for a single sequence.  Takes in fasta coverage output\n" \
  "                      coverage from the \"sect\" tool\n\n" \
  "Options (default value in (), *required):\n" \
  "     --usage          Usage\n" \
  "     --help           This message\n" \

    const char * help() const
    {
        return sect_args_HELP;
    }

#define sect_args_HIDDEN "Hidden options:"
    const char * hidden() const
    {
        return sect_args_HIDDEN;
    }

    bool validMode(char* mode_str)
    {
        return (strcmp(mode_str, "flame") == 0 ||
                strcmp(mode_str, "asm") == 0 ||
                strcmp(mode_str, "sect") == 0) ?
                    true : false;
    }


    void parse(int argc, char *argv[])
    {
        int c;

        static struct option long_options[] =
        {
            {"help",            no_argument,        0, 'h'},
            {"usage",           no_argument,        0, 'u'},
            {0, 0, 0, 0}
        };

        static const char *short_options = "uh";

        if (argc <= 1)
        {
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
                case '?':
                    cerr << "Use --usage or --help for some help" << endl << endl;
                    exit(1);
                }
            }

            error("Invalid command line arguments passed to \"kat plot\"\n\n");
            exit(1);
        }
    }
};
