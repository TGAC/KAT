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


// Add any constants here
#define DEFAULT_TEMPLATE_CONST 10

// Change class name to whatever is appropriate
class TemplateArgs
{
public:
    // Add class variables here (or in a private section with accessors here if you want to do it properly)    
    bool verbose;

    // Default constructor (add default settings here)
    TemplateArgs() :
	verbose(false)
    {}

    // Constructor that parses command line options
    TemplateArgs(int argc, char* argv[]) :
	verbose(false)
    {
        parse(argc, argv);
    }



// Modify usage string to your requirements
#define seqcvg_args_USAGE "Usage: kat <template> <mandatory args>\n"
    const char * usage() const
    {
        return seqcvg_args_USAGE;
    }

    void error(const char *msg)
    {
        cerr << endl
             << "Error: " << msg << endl << endl
             << usage() << endl
             << "Use --help for more information" << endl << endl;
        exit(1);
    }

// Modify help string to your requirements
#define sect_args_HELP "Template Help\n\n" \
  "Options (default value in (), *required):\n" \
  " <ADD ARGS HERE>" \
  " -v, --verbose        Outputs additional information to stderr\n" \
  "     --usage          Usage\n" \
  "     --help           This message\n" \

    const char * help() const
    {
        return sect_args_HELP;
    }

#define sect_args_HIDDEN "Hidden options:"
    const char * hidden() const
    {
        return sect_args_HIDDEN;
    }


    void parse(int argc, char *argv[])
    {
        int c;
        int help_flag = 0;
        int usage_flag = 0;

        static struct option long_options[] =
        {
            // ADD ADDITIONAL LONG OPTIONS HERE
            {"verbose",         no_argument,        0, 'v'},
            {"help",            no_argument,        &help_flag, 1},
            {"usage",           no_argument,        &usage_flag, 1},
            {0, 0, 0, 0}
        };

        // ADD ADDITIONAL SHORT OPTIONS HERE
        static const char *short_options = "vuh";

        if (argc <= 1)
        {
            cerr << endl
                 << usage() << endl
                 << help() << endl;
            exit(1);
        }


        while (true)
        {
            /* getopt_long stores the option index here. */
            int index = -1;

            c = getopt_long (argc, argv, short_options, long_options, &index);

            /* Detect the end of the options. */
            if (c == -1)
                break;

            switch (c)
            {
            case ':':
                cerr << "Missing required argument for "
                          << (index == -1 ? std::string(1, (char)optopt) : std::string(long_options[index].name))
                          << endl;
                exit(1);
            case '?':
                cerr << "Use --usage or --help for some help" << endl << endl;
                exit(1);
            case 'v':
                verbose = true;
                break;
            
            // ADD OPTION HANDLING CODE HERE (i.e. set class variables)
            }
        }

        if (help_flag)
        {
            cout << usage() << endl
                 << help() << endl;
            exit(0);
        }

        if (usage_flag)
        {
            cout << usage() << endl
                 << "Use --help for more information." << endl << endl;
            exit(0);
        }


        // Parse arguments
        int remaining_args = argc - optind;

        if (verbose)
            cerr << "Found " << remaining_args << " remaining arguments on the command line." << endl;

        // VALIDATION CHECK FOR REMAINING ARGS HERE
        //if(remaining_args < 2 || remaining_args > 3)
        //    error("Requires 2 or 3 arguments describing paths to jellyfish kmer counts.");


        // SET REMAINING ARGS HERE
    }

    void print()
    {
        if (verbose)
            cerr << "Verbose flag set" << endl;

        // Add any set options here for printing to stderr

        cerr << endl;
    }

private:
};

