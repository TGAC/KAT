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

using std::string;
using std::cerr;
using std::cout;

namespace kat
{

    const int32_t   DEFAULT_X_MAX           = 1000;
    const int32_t   DEFAULT_Y_MAX           = 1000000;
    const uint16_t  DEFAULT_DUPLICATION     = 5;
    const bool      DEFAULT_IGNORE_ABSENT   = false;

    const uint16_t MIN_ARGS = 1;


    class AsmPlotArgs : public BasePlotArgs
    {
    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const char* usage() const
        {
            return "Usage: kat plot asm [options] <matrix_file>\n";
        }

        const char* shortDescription() const
        {
            return "Creates a stacked histogram showing level of duplication in an assembly";
        }

        const char* longDescription() const
        {
            return  "  Shows K-mer duplication levels within an assembly by comparing K-mers fround in sequenced reads, to K-mers found in an assembly of those reads.\n" \
                    "  Uses matrix output from the \"kat comp\" tool";
        }

        const string optionsDescription() const
        {
            ostringstream help_str;

            help_str << BasePlotArgs::optionsDescription() << endl
                     << " -x  --x_max=uint16          Maximum value for the x-axis (1000)" << endl
                     << " -y  --y_max=uint64          Maximum value for the y-axis (10000000)" << endl
                     << " -a, --ignore_absent         Ignore K-mers in reads but absent from the assembly" << endl
                     << " -m, --max_dup=uint16        Maximum duplication level to show in plots (5)" << endl
                     << " -c, --columns=string        Comma separated string listing columns to show in plot.  If used," << endl
                     << "                             this overrides \"--ignore_absent\" and \"--columns\"";

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
            };

            vector<option>* long_options = BasePlotArgs::longOptions();

            for(uint8_t i = 0; i < 5; i++)
            {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        string shortOptions()
        {
            return BasePlotArgs::shortOptions() + "x:y:a:m:c";
        }

        void setOption(int c, char* option_arg) {

            BasePlotArgs::setOption(c, option_arg);

            switch(c)
            {
            case 'x':
                x_max = atoi(optarg);
                break;
            case 'y':
                y_max = atol(optarg);
                break;
            case 'a':
                ignore_absent = true;
                break;
            case 'm':
                max_duplication = atoi(optarg);
                break;
            case 'c':
                columns = string(optarg);
                break;
            }
        }


        void processRemainingArgs(const vector<string>& remaining_args)
        {
            mx_arg = remaining_args[0];
        }



        const char* currentStatus() const
        {
            ostringstream status;

            status  << BasePlotArgs::currentStatus()
                    << "X Max: " << x_max << endl
                    << "Y Max: " << y_max << endl
                    << "Max duplication level to plot: " << max_duplication << endl
                    << "Columns to plot: " << columns << endl
                    << "Ignore absent K-mers: " << ignore_absent << endl;

            return status.str().c_str();
        }
    public:
        string      mx_arg;
        uint16_t    x_max;
        uint64_t    y_max;
        bool        ignore_absent;
        uint16_t    max_duplication;
        string      columns;

        // Default constructor
        AsmPlotArgs() : BasePlotArgs(MIN_ARGS),
            mx_arg(""), x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX), ignore_absent(DEFAULT_IGNORE_ABSENT), max_duplication(DEFAULT_DUPLICATION), columns("")
        {
            title = defaultTitle();
            x_label = defaultXLabel();
            y_label = defaultYLabel();
            width = defaultWidth();
            height = defaultHeight();
        }

        // Constructor that parses command line options
        AsmPlotArgs(int argc, char* argv[]) : BasePlotArgs(MIN_ARGS),
            mx_arg(""), x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX), ignore_absent(DEFAULT_IGNORE_ABSENT), max_duplication(DEFAULT_DUPLICATION), columns("")
        {
            title = defaultTitle();
            x_label = defaultXLabel();
            y_label = defaultYLabel();
            width = defaultWidth();
            height = defaultHeight();

            parse(argc, argv);
        }

        ~AsmPlotArgs()
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
