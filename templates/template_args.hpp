//  ********************************************************************
//  This file is part of KAT - the Kmer Analysis Toolkit.
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
#include <stdlib.h>
#include <iostream>
#include <stdint.h>

using std::cout;
using std::cerr;
using std::endl;

namespace kat
{
    // TOOL CONSTANTS GO HERE

    // Must specify the minimum number of arguments to the too here (modify as required)
    const uint16_t MIN_ARGS = 1;

    // Change class name to whatever is appropriate
    class TemplateArgs : public BaseArgs
    {
    private:
        // TODO add any private content here

    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const char* usage() const
        {
            return "Usage: kat TEMPLATETOOL \n";
        }

        const char* shortDescription() const
        {
            return "Short description";
        }

        const char* longDescription() const
        {
            return  "  Template long description";
        }

        const string optionsDescription() const
        {
            ostringstream help_str;

            // TODO add option descriptions here.  Replace this option with your own.
            help_str << " -o, --option=string         Put any custom option descriptions here";

            return help_str.str();
        }

        vector<option>* longOptions()
        {
            static struct option long_options_array[] =
            {
                //TODO long options here
                {"option",       required_argument,  0, 'o'}
            };

            vector<option>* long_options = new vector<option>();

            //TODO Modify the limit on this loop as appropriate
            for(uint8_t i = 0; i < 1; i++)
            {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        string shortOptions()
        {
            //TODO short options here
            return "o:";
        }

        void setOption(int c, char* option_arg) {

            // TODO set class variable here
            switch(c)
            {
            case 'o':
                option = string(option_arg);
                break;
            }
        }

        void processRemainingArgs(const vector<string>& remaining_args)
        {
            //TODO This is ok if you have only a single required argument, although you may want to modify "prog_arg" to
            // something more appropriate.
            prog_arg = remaining_args[0];
        }

        const char* currentStatus() const
        {
            ostringstream status;

            // TODO modify as appropriate
            status << "Option: " << option << endl
                   << "Program Argument: " << prog_arg << endl;

            return status.str().c_str();
        }

    public:
        // Add class variables here (or in a private section with accessors here if you want to do it properly)
        string option;
        string prog_arg;

        // Default constructor (add default settings here)
        TemplateArgs()
        {}

        // Constructor that parses command line options
        TemplateArgs(int argc, char* argv[])
        {
            parse(argc, argv);
        }
    };
}

