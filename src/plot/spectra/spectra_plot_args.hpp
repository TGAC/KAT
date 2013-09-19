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

#include <common_args.hpp>

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

namespace kat
{
    const string DEFAULT_OUTPUT_TYPE = "png";
    const string DEFAULT_OUTPUT_FILE_PREFIX = "kat-spectra";
    const string DEFAULT_TITLE              = "K-mer Spectra";
    const string DEFAULT_X_LABEL            = "K-mer Multiplicity";
    const string DEFAULT_Y_LABEL            = "Distinct K-mers";

    const uint32_t DEFAULT_X_MIN = 0;
    const uint32_t DEFAULT_Y_MIN = 0;
    const uint32_t DEFAULT_X_MAX = 10000;
    const uint32_t DEFAULT_Y_MAX = 1000000;
    const bool DEFAULT_X_LOGSCALE = false;
    const bool DEFAULT_Y_LOGSCALE = false;
    const uint16_t DEFAULT_WIDTH = 1024;
    const uint16_t DEFAULT_HEIGHT = 1024;

    const uint16_t MIN_ARGS = 1;

    class SpectraPlotArgs : public BaseArgs
    {
    protected:

        // These methods override base class virtual methods

        const char* usage() const
        {
            return "Usage: kat plot spectra [options] -o <output_file_path> histo_file [histo_file ...]*";
        }

        const char* shortDescription() const
        {
            return "Create K-mer Spectra Plot";
        }

        const char* longDescription() const
        {
            return "  Shows K-mer spectras from kat-histo or jellyfish-histo output.";
        }

        const string optionsList() const
        {
            ostringstream help_str;

            help_str << " -p, --output_type=string    The plot file type to create: png, ps, pdf.  Warning... if pdf is selected" << endl
                 << "                             please ensure your gnuplot installation can export pdf files. (\"" << DEFAULT_OUTPUT_TYPE << "\")" << endl
                 << " -o, --output=string         Output file (\"" << DEFAULT_OUTPUT_FILE_PREFIX << "." << DEFAULT_OUTPUT_TYPE << "\")" << endl
                 << " -t, --title=string          Title for plot (\"" << DEFAULT_TITLE << "\")" << endl
                 << " -i, --x_label=string        Label for the x-axis (\"" << DEFAULT_X_LABEL << "\")" << endl
                 << " -j, --y_label=string        Label for the y-axis (\"" << DEFAULT_Y_LABEL << "\")" << endl
                 << " -r  --x_min=uint32          Minimum value for the x-axis (" << DEFAULT_X_MIN << ")" << endl
                 << " -s  --y_min=uint32          Minimum value for the y-axis (" << DEFAULT_Y_MIN << ")" << endl
                 << " -x  --x_max=uint32          Maximum value for the x-axis (" << DEFAULT_X_MAX << ")" << endl
                 << " -y  --y_max=uint32          Maximum value for the y-axis (Auto calculate \'--y_max\' from data)" << endl
                 << " -l  --x_logscale            X-axis is logscale.  This overrides the x_min and x_max limits." << endl
                 << " -m  --y_logscale            Y-axis is logscale.  This overrides the y_min and y_max limits." << endl
                 << " -w, --width=uint16          Width of canvas (" << DEFAULT_WIDTH << ")" << endl
                 << " -h, --height=uint16         Height of canvas (" << DEFAULT_HEIGHT << ")";

            return help_str.str();
        }

        vector<option>* longOptions()
        {
            static struct option long_options_array[] =
            {
                {"output_type",     required_argument,  0, 'p'},
                {"output",          required_argument,  0, 'o'},
                {"title",           required_argument,  0, 't'},
                {"x_label",         required_argument,  0, 'i'},
                {"y_label",         required_argument,  0, 'j'},
                {"x_min",           required_argument,  0, 'r'},
                {"y_min",           required_argument,  0, 's'},
                {"x_max",           required_argument,  0, 'x'},
                {"y_max",           required_argument,  0, 'y'},
                {"x_logscale",      no_argument,        0, 'l'},
                {"y_logscale",      no_argument,        0, 'm'},
                {"width",           required_argument,  0, 'w'},
                {"height",          required_argument,  0, 'h'},
                {"index",           required_argument,  0, 'n'},
                {"header",          required_argument,  0, 'd'}
            };

            vector<option>* long_options = new vector<option>();

            for(uint8_t i = 0; i < 15; i++)
            {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        const char* shortOptions() const
        {
            return "o:p:t:i:j:r:s:x:y:lmw:h:n:d:";
        }

        void setOption(int c, char* option_arg) {

            switch(c)
            {
            case 'o':
                output_arg = string(option_arg);
                break;
            case 'p':
                output_type = string(option_arg);
                break;
            case 't':
                title = string(option_arg);
                break;
            case 'i':
                x_label = string(option_arg);
                break;
            case 'j':
                y_label = string(option_arg);
                break;
            case 'r':
                x_min = atoi(option_arg);
                break;
            case 's':
                y_min = atoi(option_arg);
                break;
            case 'x':
                x_max = atoi(option_arg);
                break;
            case 'y':
                y_max = atoi(option_arg);
                break;
            case 'l':
                x_logscale = true;
                break;
            case 'm':
                y_logscale = true;
                break;
            case 'w':
                width = atoi(option_arg);
                break;
            case 'h':
                height = atoi(option_arg);
                break;
            }
        }


        void processRemainingArgs(const vector<string>& remaining_args)
        {
            histo_paths = remaining_args;
        }



        const char* currentStatus() const
        {
            ostringstream status;

            status << "Output type: " << output_type.c_str() << endl;
            status << "Output file specified: " << output_arg.c_str() << endl;
            status << "Plot title: " << title << endl;
            status << "X Label: " << x_label << endl;
            status << "Y Label: " << y_label << endl;
            status << "X Min: " << x_min << endl;
            status << "Y Min: " << y_min << endl;
            status << "X Max: " << x_max << endl;
            status << "Y Max: " << y_max << endl;
            status << "X Logscale: " << x_logscale << endl;
            status << "Y Logscale: " << y_logscale << endl;
            status << "Width: " << width << endl;
            status << "Height: " << height << endl;

            return status.str().c_str();
        }

    public:

        vector<string> histo_paths;
        string  output_type;
        string  output_arg;
        string  title;
        string  x_label;
        string  y_label;
        uint32_t x_min;
        uint32_t y_min;
        uint32_t x_max;
        uint32_t y_max;
        bool x_logscale;
        bool y_logscale;
        uint16_t width;
        uint16_t height;

        // Default constructor
        SpectraPlotArgs() :
            BaseArgs(MIN_ARGS), histo_paths(vector<string>()), output_type(DEFAULT_OUTPUT_TYPE), output_arg(""), title(DEFAULT_TITLE),
            x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL),
            x_min(DEFAULT_X_MIN), y_min(DEFAULT_Y_MIN), x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX),
            x_logscale(DEFAULT_X_LOGSCALE), y_logscale(DEFAULT_Y_LOGSCALE),
            width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT)
        {
        }

        // Constructor that parses command line options
        SpectraPlotArgs(int argc, char* argv[]) :
            BaseArgs(MIN_ARGS), histo_paths(vector<string>()), output_type(DEFAULT_OUTPUT_TYPE), output_arg(""), title(DEFAULT_TITLE),
            x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL),
            x_min(DEFAULT_X_MIN), y_min(DEFAULT_Y_MIN), x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX),
            x_logscale(DEFAULT_X_LOGSCALE), y_logscale(DEFAULT_Y_LOGSCALE),
            width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT)
        {
            parse(argc, argv);
        }

        ~SpectraPlotArgs()
        {}


        // Work out the output path to use (either user specified or auto generated)
        string determineOutputPath()
        {
            std::ostringstream output_str;
            output_str << DEFAULT_OUTPUT_FILE_PREFIX << "." << output_type;
            return output_arg.empty() ? output_str.str() : output_arg;
        }

    };
}
