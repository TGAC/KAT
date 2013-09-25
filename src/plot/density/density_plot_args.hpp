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
#include <stdint.h>
#include <iostream>

#include "../common_plot_args.hpp"
#include <str_utils.hpp>

using std::string;
using std::cerr;
using std::cout;
using std::endl;

namespace kat
{
    const string DEFAULT_Z_LABEL    = "Z";
    const int32_t DEFAULT_X_MAX     = 1000;
    const int32_t DEFAULT_Y_MAX     = 1000;
    const int64_t DEFAULT_Z_MAX     = 10000;

    const uint16_t MIN_ARGS = 1;

    class DensityPlotArgs : public BasePlotArgs
    {
    private:
        bool        z_label_mod;
        bool        x_max_mod;
        bool        y_max_mod;
        bool        z_max_mod;

        void init()
        {
            title = defaultTitle();
            x_label = defaultXLabel();
            y_label = defaultYLabel();
            width = defaultWidth();
            height = defaultHeight();

            z_label_mod = false;
            x_max_mod = false;
            y_max_mod = false;
            z_max_mod = false;
        }

    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const string usage() const
        {
            return "Usage: kat plot density [options] <matrix_file>";
        }

        const string shortDescription() const
        {
            return "Create K-mer Density Plots.";
        }

        const string longDescription() const
        {
            string long_desc = "Creates a scatter plot, where the density or \"heat\" at each point represents the number of distinct K-mers " \
                               "at that point.  Typically this is used to visualise a matrix produced by the \"kat comp\" tool to compare " \
                               "multiplicities from two K-mer hashes produced by different NGS reads, or to visualise the GC vs K-mer " \
                               "multiplicity matricies produced by the \"kat gcp\" tool.";

            return lineBreakString(long_desc, 78, "  ");
        }

        const string optionsDescription() const
        {
            ostringstream help_str;

            help_str << BasePlotArgs::optionsDescription() << endl
                     << " -k, --z_label=string        Label for the z-axis (\"Z\", or value from matrix metadata if present)" << endl
                     << " -x, --x_max=uint32          Maximum value for the x-axis (1000, or value from matrix metadata if present)" << endl
                     << " -y  --y_max=uint32          Maximum value for the y-axis (1000, or value from matrix metadata if present)" << endl
                     << " -z, --z_max=uint64          Cap for matrix values.  Values greater than this cap will be displayed at " << endl
                     << "                             maximum intensity, i.e. white. (10000, or value from matrix metadata if present)";

            return help_str.str();
        }

        vector<option>* longOptions()
        {
            static struct option long_options_array[] =
            {
                {"z_label",         required_argument,  0, 'k'},
                {"x_max",           required_argument,  0, 'x'},
                {"y_max",           required_argument,  0, 'y'},
                {"z_max",           required_argument,  0, 'z'}
            };

            vector<option>* long_options = BasePlotArgs::longOptions();

            for(uint8_t i = 0; i < 4; i++)
            {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        string shortOptions()
        {
            return BasePlotArgs::shortOptions() + "k:x:y:z:";
        }

        void setOption(int c, char* option_arg) {

            BasePlotArgs::setOption(c, option_arg);

            switch(c)
            {
            case 'k':
                z_label = string(optarg);
                break;
            case 'x':
                x_max = atoi(optarg);
                x_max_mod = true;
                break;
            case 'y':
                y_max = atoi(optarg);
                y_max_mod = true;
                break;
            case 'z':
                z_max = atoi(optarg);
                z_max_mod = true;
                break;
            }
        }


        void processRemainingArgs(const vector<string>& remaining_args)
        {
            mx_arg = remaining_args[0];
        }



        const string currentStatus() const
        {
            ostringstream status;

            status  << BasePlotArgs::currentStatus()
                    << "Z Label: " << z_label << endl
                    << "X Max: " << x_max << endl
                    << "Y Max: " << y_max << endl
                    << "Z Max: " << z_max << endl;

            return status.str().c_str();
        }



    public:
        string      mx_arg;
        string      z_label;
        int16_t     x_max;
        int16_t     y_max;
        int64_t     z_max;

        // Default constructor
        DensityPlotArgs() : BasePlotArgs(MIN_ARGS),
            mx_arg(""), z_label(DEFAULT_Z_LABEL),
            x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX), z_max(DEFAULT_Z_MAX)
        {
            init();
        }

        // Constructor that parses command line options
        DensityPlotArgs(int argc, char* argv[]) : BasePlotArgs(MIN_ARGS),
            mx_arg(""), z_label(DEFAULT_Z_LABEL),
            x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX), z_max(DEFAULT_Z_MAX)
        {
            init();

            parse(argc, argv);
        }

        ~DensityPlotArgs()
        {}

        // ***************************************************
        // These methods override BasePlotArgs virtual methods

        const string defaultOutputPrefix() const    { return "kat-plot-density"; }
        const string defaultTitle() const           { return "Flame Plot"; }
        const string defaultXLabel() const          { return "X"; }
        const string defaultYLabel() const          { return "Y"; }
        const uint16_t defaultWidth() const         { return 1024; }
        const uint16_t defaultHeight() const        { return 1024; }



        const bool zLabelModified()     {return z_label_mod;}
        const bool xMaxModified()       {return x_max_mod;}
        const bool yMaxModified()       {return y_max_mod;}
        const bool zMaxModified()       {return z_max_mod;}
    };
}
