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
#include <string.h>
#include <stdlib.h>
#include <iostream>

#include <common_args.hpp>

using std::string;
using std::cerr;
using std::cout;
using std::endl;

namespace kat
{
    const string KAT_SECT_ID  = "sect";
    const string KAT_COMP_ID  = "comp";
    const string KAT_GCP_ID   = "gcp";
    const string KAT_HIST_ID  = "hist";
    const string KAT_PLOT_ID  = "plot";

    const uint16_t MIN_ARGS = 0;

    class KatArgs : public BaseArgs
    {
    private:
        string  mode_arg;
        int     mode_argc;
        char**  mode_argv;

    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const char* usage() const               { return "Usage: kat <mode>\n"; }
        const char* shortDescription() const    { return "The K-mer Analysis Toolkist (KAT) contains a number of tools that analyse jellyfish K-mer hashes."; }
        const char* longDescription() const
        {
            return  "First argument should be the tool/mode you wish to use:\n\n" \
                    "   - sect:  SEquence Coverage estimator Tool.  Estimates the coverage of each sequence in a fasta file using\n" \
                    "            K-mers from a jellyfish hash.\n" \
                    "   - comp:  K-mer comparison tool.  Creates a matrix of shared K-mers between two jellyfish hashes.\n" \
                    "   - gcp:   K-mer GC Processor.  Creates a matrix of the number of K-mers found given a GC count and a K-mer\n" \
                    "            count.\n" \
                    "   - histo: Create an histogram of k-mer occurrences from a jellyfish hash.  Adds metadata in output for easy\n" \
                    "            plotting.\n" \
                    "   - plot:  Plotting tool.  Contains several plotting tools to visualise K-mer and compare distributions.\n" \
                    "            Requires gnuplot.";
        }

        const string optionsDescription() const    { return " -V, --version               Version"; }

        vector<option>* longOptions()
        {
            static struct option long_options_array[] =
            {
                {"version", no_argument,       0, 'V'}
            };


            vector<option>* long_options = new vector<option>();

            long_options->push_back(long_options_array[0]);

            return long_options;
        }

        string shortOptions()                   { return "V"; }

        void setOption(int c, char* option_arg)
        {
            switch (c)
            {
            case 'V':
                printVersion();
                exit(0);
            }
        }

        void processRemainingArgs(const vector<string>& remaining_args) {}
        const char* currentStatus() const       { return ""; }

    public:


        // Default constructor
        KatArgs() : BaseArgs(MIN_ARGS)
        {}

        // Constructor that parses command line options
        KatArgs(int argc, char* argv[]) : BaseArgs(MIN_ARGS)
        {
            customParse(argc, argv);
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


        void printVersion(std::ostream &os = std::cout) const
        {
    #ifndef PACKAGE_NAME
    #define PACKAGE_NAME "K-mer Analysis Toolkit (KAT)"
    #endif

    #ifndef PACKAGE_VERSION
    #define PACKAGE_VERSION "0.3.1"
    #endif
            os << PACKAGE_NAME << " V" << PACKAGE_VERSION << "\n";
        }

        bool validMode(string mode_str)
        {
            return (mode_str.compare(KAT_SECT_ID) == 0 ||
                    mode_str.compare(KAT_COMP_ID) == 0 ||
                    mode_str.compare(KAT_GCP_ID) == 0 ||
                    mode_str.compare(KAT_HIST_ID) == 0 ||
                    mode_str.compare(KAT_PLOT_ID) == 0) ?
                        true : false;
        }

        void customParse(int argc, char *argv[])
        {
            if (argc <= 1)
            {
                error("No mode specified");
            }
            else if (validMode(string(argv[1]))) {

                mode_arg = argv[1];
                mode_argc = argc - 1;
                mode_argv = argv + 1;
            }
            else {

                // Let BaseArgs have a go, but make sure we fail after
                parse(argc, argv);

                error("Invalid command line arguments passed to \"kat plot\"");
                exit(1);
            }
        }
    };
}
