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
#include <iostream>
#include <stdlib.h>
#include <vector>

#include <common_args.hpp>

#include "plot_main.hpp"

using std::string;
using std::cerr;
using std::cout;
using std::endl;

namespace kat
{
    const string KAT_PLOT_DENSITY_ID        = "density";
    const string KAT_PLOT_PROFILE_ID        = "profile";
    const string KAT_PLOT_SPECTRA_CN_ID     = "spectra-cn";
    const string KAT_PLOT_SPECTRA_HIST_ID   = "spectra-hist";
    const string KAT_PLOT_SPECTRA_MX_ID   = "spectra-mx";

    const uint16_t MIN_ARGS = 0;

    class PlotArgs : public BaseArgs
    {
    private:
        string  mode_arg;
        int     mode_argc;
        char**  mode_argv;

    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const string usage() const               { return "Usage: kat plot <mode>"; }
        const string shortDescription() const    { return "Create K-mer Plots"; }
        const string longDescription() const
        {
            return  "First argument should be the plot mode you wish to use:\n" \
                    "  - density:         Creates a density plot from a matrix created with the \"comp\" tool.  Typically this is\n" \
                    "                     used to compare two K-mer hashes produced by different NGS reads.\n" \
                    "  - profile:         Creates a K-mer coverage plot for a single sequence.  Takes in fasta coverage output\n" \
                    "                     coverage from the \"sect\" tool\n" \
                    "  - spectra-cn:      Creates a stacked histogram using a matrix created with the \"comp\" tool.  Typically\n" \
                    "                     this is used to compare a jellyfish hash produced from a read set to a jellyfish hash\n" \
                    "                     produced from an assembly. The plot shows the amount of distinct K-mers absent, as well\n" \
                    "                     as the copy number variation present within the assembly.\n" \
                    "  - spectra-hist:    Creates a K-mer spectra plot for a set of K-mer histograms produced either by jellyfish-\n" \
                    "                     histo or kat-histo.\n" \
                    "  - spectra-mx:      Creates a K-mer spectra plot for a set of K-mer histograms that are derived from\n" \
                    "                     selected rows or columns in a matrix produced by the \"comp\".";
        }

        const string optionsDescription() const    { return ""; }

        vector<option>* longOptions()
        {
            vector<option>* long_options = new vector<option>();

            return long_options;
        }

        string shortOptions()                   { return ""; }
        void setOption(int c, string& option_arg) {}
        void processRemainingArgs(const vector<string>& remaining_args) {}
        const string currentStatus() const       { return ""; }

    public:

        // Default constructor
        PlotArgs() : BaseArgs(MIN_ARGS)
        {}

        // Constructor that parses command line options
        PlotArgs(int argc, char* argv[]) : BaseArgs(MIN_ARGS)
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

        bool validMode(string mode_str)
        {
            return (mode_str.compare(KAT_PLOT_DENSITY_ID) == 0 ||
                    mode_str.compare(KAT_PLOT_PROFILE_ID) == 0 ||
                    mode_str.compare(KAT_PLOT_SPECTRA_CN_ID) == 0 ||
                    mode_str.compare(KAT_PLOT_SPECTRA_HIST_ID) == 0 ||
                    mode_str.compare(KAT_PLOT_SPECTRA_MX_ID) == 0) ?
                        true : false;
        }


        void customParse(int argc, char *argv[])
        {
            if (argc <= 1)
            {
                error("No plot mode specified");
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
