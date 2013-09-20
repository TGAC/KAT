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
using std::ostringstream;

namespace kat
{
    const uint32_t DEFAULT_Y_MAX    = 1000;
    const uint32_t DEFAULT_FASTA_INDEX = 0;

    const uint16_t MIN_ARGS = 1;


    class SectPlotArgs : public BasePlotArgs
    {

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const char* usage() const
        {
            return "Usage: kat plot sect [options] <sect_profile_file>\n";
        }

        const char* shortDescription() const
        {
            return "Create Sequence Coverage Plot";
        }

        const char* longDescription() const
        {
            return "  Shows K-mer coverage level across an sequence.";
        }

        const string optionsDescription() const
        {
            ostringstream help_str;

            help_str << BasePlotArgs::optionsDescription() << endl
                     << " -y  --y_max=uint32          Maximum value for the y-axis (" << DEFAULT_Y_MAX << ", or value from matrix metadata if present)" << endl
                     << " -n, --index=uint32          Index of fasta entry to plot.  First entry is 1. (" << DEFAULT_FASTA_INDEX << ")" << endl
                     << " -d, --header=string         Fasta header of fasta entry to plot.  NOTE: \'--header\' has priority over \'--index\'.";

            return help_str.str();
        }

        vector<option>* longOptions()
        {
            static struct option long_options_array[] =
            {
                {"y_max",           required_argument,  0, 'y'},
                {"index",           required_argument,  0, 'n'},
                {"header",          required_argument,  0, 'd'}
            };

            vector<option>* long_options = BasePlotArgs::longOptions();

            for(uint8_t i = 0; i < 3; i++)
            {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        string shortOptions()
        {
            return BasePlotArgs::shortOptions() + "y:n:d:";
        }

        void setOption(int c, char* option_arg) {

            BasePlotArgs::setOption(c, option_arg);

            switch(c)
            {
            case 'y':
                y_max = atoi(option_arg);
                break;
            case 'n':
                fasta_index = atoi(option_arg);
                break;
            case 'd':
                fasta_header = option_arg;
                break;
            }
        }


        void processRemainingArgs(const vector<string>& remaining_args)
        {
            sect_file_arg = remaining_args[0];
        }



        const char* currentStatus() const
        {
            ostringstream status;

            status  << BasePlotArgs::currentStatus()
                    << "Y Max: " << y_max << endl
                    << "Fasta index to plot: " << fasta_index << endl
                    << "Fasta header to plot: " << fasta_header << endl;

            return status.str().c_str();
        }

    public:
        string      sect_file_arg;
        uint32_t    y_max;
        uint32_t    fasta_index;
        string      fasta_header;

        // Default constructor
        SectPlotArgs() :
            BasePlotArgs(MIN_ARGS), sect_file_arg(""), y_max(DEFAULT_Y_MAX), fasta_index(DEFAULT_FASTA_INDEX), fasta_header("")
        {
            title = defaultTitle();
            x_label = defaultXLabel();
            y_label = defaultYLabel();
            width = defaultWidth();
            height = defaultHeight();
        }

        // Constructor that parses command line options
        SectPlotArgs(int argc, char* argv[]) :
            BasePlotArgs(MIN_ARGS), sect_file_arg(""), y_max(DEFAULT_Y_MAX), fasta_index(DEFAULT_FASTA_INDEX), fasta_header("")
        {
            title = defaultTitle();
            x_label = defaultXLabel();
            y_label = defaultYLabel();
            width = defaultWidth();
            height = defaultHeight();

            parse(argc, argv);
        }


        string autoTitle(string& header)
        {
            std::ostringstream output_str;
            output_str << defaultTitle() << ": " << sect_file_arg << " - " << header;
            return title.compare(defaultTitle()) == 0 ? output_str.str() : title;
        }


        ~SectPlotArgs()
        {}



        // ***************************************************
        // These methods override BasePlotArgs virtual methods

        const string defaultOutputPrefix() const    { return "kat-plot-sect"; }
        const string defaultTitle() const           { return "Sequence Coverage Plot"; }
        const string defaultXLabel() const          { return "Position (nt)"; }
        const string defaultYLabel() const          { return "K-mer Coverage"; }
        const uint16_t defaultWidth() const         { return 1024; }
        const uint16_t defaultHeight() const        { return 1024; }



    };

}
