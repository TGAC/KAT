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
#include <str_utils.hpp>

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

namespace kat
{
    const string DEFAULT_OUTPUT_TYPE        = "png";


    class BasePlotArgs : public BaseArgs
    {
    private:
        bool    title_mod;
        bool    x_label_mod;
        bool    y_label_mod;


    protected:

        // These methods override BaseArgs virtual methods

        const string optionsDescription() const
        {
            ostringstream help_str;

            help_str << " -p, --output_type=string    The plot file type to create: png, ps, pdf.  Warning... if pdf is selected" << endl
                     << "                             please ensure your gnuplot installation can export pdf files. (\"" << DEFAULT_OUTPUT_TYPE << "\")" << endl
                     << " -o, --output=string         Output file (\"" << defaultOutputPrefix() << "." << DEFAULT_OUTPUT_TYPE << "\")" << endl
                     << " -t, --title=string          Title for plot (\"" << defaultTitle() << "\")" << endl
                     << " -i, --x_label=string        Label for the x-axis (\"" << defaultXLabel() << "\")" << endl
                     << " -j, --y_label=string        Label for the y-axis (\"" << defaultYLabel() << "\")" << endl
                     << " -w, --width=uint16          Width of canvas (" << defaultWidth() << ")" << endl
                     << " -h, --height=uint16         Height of canvas (" << defaultHeight() << ")";

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
                {"width",           required_argument,  0, 'w'},
                {"height",          required_argument,  0, 'h'}
            };

            vector<option>* long_options = new vector<option>();

            for(uint8_t i = 0; i < 7; i++)
            {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        string shortOptions()
        {
            return "o:p:t:i:j:w:h:";
        }

        void setOption(int c, string& option_arg) {

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
                x_label_mod = true;
                break;
            case 'j':
                y_label = string(option_arg);
                y_label_mod = true;
                break;
            case 'w':
                width = strToInt16(option_arg);
                break;
            case 'h':
                height = strToInt16(option_arg);
                break;
            }
        }

        const string currentStatus() const
        {
            ostringstream status;

            status << "Output type: " << output_type.c_str() << endl;
            status << "Output file specified: " << output_arg.c_str() << endl;
            status << "Plot title: " << title << endl;
            status << "X Label: " << x_label << endl;
            status << "Y Label: " << y_label << endl;
            status << "Width: " << width << endl;
            status << "Height: " << height << endl;

            return status.str();
        }

    public:

        string  output_type;
        string  output_arg;
        string  title;
        string  x_label;
        string  y_label;
        uint16_t width;
        uint16_t height;

        /**
         * @brief BasePlotArgs Constructor.  Requires the number of trailing arguments to be specified for this plotting tool.
         * @param min_args
         */
        BasePlotArgs(uint16_t min_args) : BaseArgs(min_args), output_type(DEFAULT_OUTPUT_TYPE), output_arg("")
        {
            title_mod = false;
            x_label_mod = false;
            y_label_mod = false;
        }

        /**
         * @brief ~BasePlotArgs Virtual destructor makes this class abstract
         */
        virtual ~BasePlotArgs() {}


        // **************************************************
        // Default values for properties which must be defined
        // by the child class

        virtual const string defaultOutputPrefix() const = 0;
        virtual const string defaultTitle() const = 0;
        virtual const string defaultXLabel() const = 0;
        virtual const string defaultYLabel() const = 0;
        virtual const uint16_t defaultWidth() const = 0;
        virtual const uint16_t defaultHeight() const = 0;


        /**
         * @brief determineOutputPath Work out the output path to use (either user specified or auto generated)
         * @return
         */
        string determineOutputPath()
        {
            std::ostringstream output_str;
            output_str << defaultOutputPrefix() << "." << output_type;
            return output_arg.empty() ? output_str.str() : output_arg;
        }


        const bool titleModified()    { return title_mod; }
        const bool xLabelModified()   { return x_label_mod; }
        const bool yLabelModified()   { return y_label_mod; }

    };
}
