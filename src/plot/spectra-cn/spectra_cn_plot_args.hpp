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

#include "../common_plot_args.hpp"
#include <str_utils.hpp>

using std::string;
using std::cerr;
using std::cout;

namespace kat
{

    const int32_t   DEFAULT_X_MAX           = 1000;
    const int32_t   DEFAULT_Y_MAX           = 1000000;
    const uint16_t  DEFAULT_DUPLICATION     = 5;
    const bool      DEFAULT_IGNORE_ABSENT   = false;
    const bool      DEFAULT_CUMULATIVE      = false;

    const uint16_t MIN_ARGS = 1;


    class SpectraCnPlotArgs : public BasePlotArgs
    {
    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const string usage() const
        {
            return "Usage: kat plot spectra-cn [options] <matrix_file>";
        }

        const string shortDescription() const
        {
            return "Creates a stacked histogram showing the level of duplication in an assembly.";
        }

        const string longDescription() const
        {
            string long_desc =  "Shows K-mer duplication levels, which correspond to copy number variation within an assembly by comparing " \
                                "K-mers found in sequenced reads, to K-mers found in an assembly of those reads. Uses matrix output from the " \
                                "\"kat comp\" tool.";

            return lineBreakString(long_desc, 78, "  ");
        }

        const string optionsDescription() const
        {
            ostringstream help_str;

            help_str << BasePlotArgs::optionsDescription() << endl
                     << " -x  --x_max=uint16          Maximum value for the x-axis (1000)" << endl
                     << " -y  --y_max=uint64          Maximum value for the y-axis (10000000)" << endl
                     << " -a, --ignore_absent         Ignore K-mers in reads but absent from the assembly" << endl
                     << " -m, --max_dup=uint16        Maximum duplication level to show in plots (5)" << endl
                     << " -c, --columns=string        Comma separated string listing columns to show in plot.  If used, this" << endl
                     << "                             overrides \"--ignore_absent\" and \"--columns\"" << endl
                     << " -u, --cumulative            Plot cumulative distribution of kmers";

            return help_str.str();
        }

        vector<option>* longOptions()
        {
            static struct option long_options_array[] =
            {
                {"x_max",           required_argument,  0, 'x'},
                {"y_max",           required_argument,  0, 'y'},
                {"ignore_absent",   no_argument,        0, 'a'},
                {"max_dup",         required_argument,  0, 'm'},
                {"columns",         required_argument,  0, 'c'},
                {"cumulative",      no_argument,        0, 'u'}
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
            return BasePlotArgs::shortOptions() + "x:y:am:c:u";
        }

        void setOption(int c, string& option_arg) {

            BasePlotArgs::setOption(c, option_arg);

            switch(c)
            {
            case 'x':
                x_max = strToInt16(optarg);
                break;
            case 'y':
                y_max = strToInt64(optarg);
                break;
            case 'a':
                ignore_absent = true;
                break;
            case 'm':
                max_duplication = strToInt16(optarg);
                break;
            case 'c':
                columns = string(optarg);
                break;
            case 'u':
                cumulative = true;
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
                    << "X Max: " << x_max << endl
                    << "Y Max: " << y_max << endl
                    << "Max duplication level to plot: " << max_duplication << endl
                    << "Columns to plot: " << columns << endl
                    << "Ignore absent K-mers: " << ignore_absent << endl
                    << "Cumulative mode: " << cumulative << endl;

            return status.str().c_str();
        }


    public:
        string      mx_arg;
        uint16_t    x_max;
        uint64_t    y_max;
        bool        ignore_absent;
        uint16_t    max_duplication;
        string      columns;
        bool        cumulative;

        // Default constructor
        SpectraCnPlotArgs() : BasePlotArgs(MIN_ARGS),
            mx_arg(""), x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX), 
            ignore_absent(DEFAULT_IGNORE_ABSENT), max_duplication(DEFAULT_DUPLICATION), 
            columns(""), cumulative(DEFAULT_CUMULATIVE)
        {
            title = defaultTitle();
            x_label = defaultXLabel();
            y_label = defaultYLabel();
            width = defaultWidth();
            height = defaultHeight();
        }

        // Constructor that parses command line options
        SpectraCnPlotArgs(int argc, char* argv[]) : BasePlotArgs(MIN_ARGS),
            mx_arg(""), x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX), 
            ignore_absent(DEFAULT_IGNORE_ABSENT), max_duplication(DEFAULT_DUPLICATION), 
            columns(""), cumulative(DEFAULT_CUMULATIVE)
        {
            title = defaultTitle();
            x_label = defaultXLabel();
            y_label = defaultYLabel();
            width = defaultWidth();
            height = defaultHeight();

            parse(argc, argv);
        }

        ~SpectraCnPlotArgs()
        {}


        // ***************************************************
        // These methods override BasePlotArgs virtual methods

        const string defaultOutputPrefix() const    { return "kat-plot-asm"; }
        const string defaultTitle() const           { return "Assembly duplication histogram"; }
        const string defaultXLabel() const          { return "K-mer Multiplicity"; }
        const string defaultYLabel() const          { return "Distinct K-mer Count"; }
        const uint16_t defaultWidth() const         { return 1024; }
        const uint16_t defaultHeight() const        { return 1024; }

    };
}
