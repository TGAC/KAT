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
using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

#include <common_args.hpp>
#include <str_utils.hpp>

namespace kat {
    const string DEFAULT_OUTPUT_TYPE = "png";

    class BasePlotArgs {
    private:
        bool title_mod;
        bool x_label_mod;
        bool y_label_mod;


        const string currentStatus() const {
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

        string output_type;
        string output_arg;
        string title;
        string x_label;
        string y_label;
        uint16_t width;
        uint16_t height;

        /**
         * @brief BasePlotArgs Constructor.  Requires the number of trailing arguments to be specified for this plotting tool.
         * @param min_args
         */
        BasePlotArgs(uint16_t min_args) : BaseArgs(min_args), output_type(DEFAULT_OUTPUT_TYPE), output_arg("") {
            title_mod = false;
            x_label_mod = false;
            y_label_mod = false;
        }

        /**
         * @brief ~BasePlotArgs Virtual destructor makes this class abstract
         */
        virtual ~BasePlotArgs() {
        }


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
        string determineOutputPath() {
            std::ostringstream output_str;
            output_str << defaultOutputPrefix() << "." << output_type;
            return output_arg.empty() ? output_str.str() : output_arg;
        }

        const bool titleModified() {
            return title_mod;
        }

        const bool xLabelModified() {
            return x_label_mod;
        }

        const bool yLabelModified() {
            return y_label_mod;
        }

    };
}
