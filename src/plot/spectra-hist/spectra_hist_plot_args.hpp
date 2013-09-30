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
#include <stdint.h>
#include <vector>

#include "../common_plot_args.hpp"
#include <str_utils.hpp>

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

namespace kat
{
    const uint32_t DEFAULT_X_MIN = 0;
    const uint32_t DEFAULT_Y_MIN = 0;
    const uint32_t DEFAULT_X_MAX = 10000;
    const uint32_t DEFAULT_Y_MAX = 1000000;
    const bool DEFAULT_X_LOGSCALE = false;
    const bool DEFAULT_Y_LOGSCALE = false;

    const uint16_t MIN_ARGS = 1;

    class SpectraHistPlotArgs : public BasePlotArgs
    {
    private:

        void init()
        {
            title = defaultTitle();
            x_label = defaultXLabel();
            y_label = defaultYLabel();
            width = defaultWidth();
            height = defaultHeight();
        }

    protected:


        // ***********************************************
        // These methods override BaseArgs virtual methods

        const string usage() const
        {
            return "Usage: kat plot spectra-hist [options] <histo_file> [<histo_file> ...]*";
        }

        const string shortDescription() const
        {
            return "Creates K-mer Spectra Plot from one or more histograms.";
        }

        const string longDescription() const
        {
            string long_desc = "Produces K-mer spectras from \"kat hist\" or \"jellyfish histo\" output.  This tool is designed to plot line " \
                   "graphs of one or more histograms.  The idea is to be able to compare total K-mer counts between different " \
                   "datasets.";

            return lineBreakString(long_desc, 78, "  ");
        }

        const string optionsDescription() const
        {
            ostringstream help_str;

            help_str << BasePlotArgs::optionsDescription() << endl
                     << " -r  --x_min=uint32          Minimum value for the x-axis (" << DEFAULT_X_MIN << ")" << endl
                     << " -s  --y_min=uint32          Minimum value for the y-axis (" << DEFAULT_Y_MIN << ")" << endl
                     << " -x  --x_max=uint32          Maximum value for the x-axis (" << DEFAULT_X_MAX << ")" << endl
                     << " -y  --y_max=uint32          Maximum value for the y-axis (Auto calculate \'--y_max\' from data)" << endl
                     << " -l  --x_logscale            X-axis is logscale.  This overrides the x_min and x_max limits." << endl
                     << " -m  --y_logscale            Y-axis is logscale.  This overrides the y_min and y_max limits.";

            return help_str.str();
        }

        vector<option>* longOptions()
        {
            static struct option long_options_array[] =
            {
                {"x_min",           required_argument,  0, 'r'},
                {"y_min",           required_argument,  0, 's'},
                {"x_max",           required_argument,  0, 'x'},
                {"y_max",           required_argument,  0, 'y'},
                {"x_logscale",      no_argument,        0, 'l'},
                {"y_logscale",      no_argument,        0, 'm'}
            };

            vector<option>* long_options = BasePlotArgs::longOptions();

            for(uint8_t i = 0; i < 6; i++)
            {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        string shortOptions()
        {
            return BasePlotArgs::shortOptions() + "r:s:x:y:lm";
        }

        void setOption(int c, string& option_arg) {

            BasePlotArgs::setOption(c, option_arg);

            switch(c)
            {
            case 'r':
                x_min = strToInt32(option_arg);
                break;
            case 's':
                y_min = strToInt32(option_arg);
                break;
            case 'x':
                x_max = strToInt32(option_arg);
                break;
            case 'y':
                y_max = strToInt32(option_arg);
                break;
            case 'l':
                x_logscale = true;
                break;
            case 'm':
                y_logscale = true;
                break;
            }
        }


        void processRemainingArgs(const vector<string>& remaining_args)
        {
            for(uint16_t i = 0; i < remaining_args.size(); i++)
            {
                histo_paths.push_back(string(remaining_args[i]));
            }
        }



        const string currentStatus() const
        {
            ostringstream status;

            status  << BasePlotArgs::currentStatus()
                    << "X Min: " << x_min << endl
                    << "Y Min: " << y_min << endl
                    << "X Max: " << x_max << endl
                    << "Y Max: " << y_max << endl
                    << "X Logscale: " << x_logscale << endl
                    << "Y Logscale: " << y_logscale << endl;

            return status.str().c_str();
        }

    public:

        vector<string> histo_paths;
        uint32_t x_min;
        uint32_t y_min;
        uint32_t x_max;
        uint32_t y_max;
        bool x_logscale;
        bool y_logscale;

        // Default constructor
        SpectraHistPlotArgs() : BasePlotArgs(MIN_ARGS),
            histo_paths(vector<string>()),
            x_min(DEFAULT_X_MIN), y_min(DEFAULT_Y_MIN), x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX),
            x_logscale(DEFAULT_X_LOGSCALE), y_logscale(DEFAULT_Y_LOGSCALE)
        {
            init();
        }

        // Constructor that parses command line options
        SpectraHistPlotArgs(int argc, char* argv[]) : BasePlotArgs(MIN_ARGS),
            histo_paths(vector<string>()),
            x_min(DEFAULT_X_MIN), y_min(DEFAULT_Y_MIN), x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX),
            x_logscale(DEFAULT_X_LOGSCALE), y_logscale(DEFAULT_Y_LOGSCALE)
        {
            init();

            parse(argc, argv);
        }

        ~SpectraHistPlotArgs()
        {}

        // ***************************************************
        // These methods override BasePlotArgs virtual methods

        const string defaultOutputPrefix() const    { return "kat-plot-spectra-hist"; }
        const string defaultTitle() const           { return "K-mer Spectra"; }
        const string defaultXLabel() const          { return "K-mer Multiplicity"; }
        const string defaultYLabel() const          { return "Distinct K-mers"; }
        const uint16_t defaultWidth() const         { return 1024; }
        const uint16_t defaultHeight() const        { return 1024; }

    };
}
