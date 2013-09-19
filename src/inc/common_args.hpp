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


namespace kat
{
    using std::cout;
    using std::cerr;
    using std::endl;
    using std::ostringstream;
    using std::vector;

    const bool DEFAULT_VERBOSE = false;

    /**
     * @brief The BaseArgs class
     */
    class BaseArgs
    {
    private:
        uint16_t min_args;
        vector<string> remaining_args;

    protected:

        // These methods must be overriden by inheriting class


        /**
         * @brief usage A usage string that briefly describes how to use this program
         * @return
         */
        virtual const char * usage() const = 0;

        /**
         * @brief shortDescription A very brief description of the tool
         * @return
         */
        virtual const char * shortDescription() const = 0;

        /**
         * @brief longDescription A longer description of the tool
         * @return
         */
        virtual const char * longDescription() const = 0;

        /**
         * @brief longDescription Program specific help for the available options
         * @return
         */
        virtual const string optionsList() const = 0;

        /**
         * @brief longOptions Long options specific to this program
         * @return
         */
        virtual vector<option>* longOptions() = 0;

        /**
         * @brief shortOptions Short options specific to this program
         * @return
         */
        virtual const char* shortOptions() const = 0;

        /**
         * @brief setOption Sets this particular option
         * @param c short option identifier
         * @param option_arg The argument associated with this option
         */
        virtual void setOption(int c, char* option_arg) = 0;

        /**
         * @brief processRemainingArgs The inheriting class should implement this method so that
         * it can process the remaining arguments left on the command line.  This method will be
         * called at the end of the "parse" method.
         * @param remaining_args The remaining arguments left on the command line
         */
        virtual void processRemainingArgs(const vector<string>& remaining_args) = 0;


        /**
         * @brief printValues Print any program specific values
         */
        virtual const char* currentStatus() const = 0;


    public:

        /**
         * @brief Every program will have a verbose argument
         */
        bool verbose;


        /**
         * Default constructor
         */
        BaseArgs(uint16_t _min_args) : verbose(DEFAULT_VERBOSE)
        {
            min_args = _min_args;
        }

        /**
         * @brief ~BaseArgs Virtual destructor makes this class abstract
         */
        virtual ~BaseArgs() = 0;


        /**
         * @brief error The error message to show for this program if the user didn't enter the correct syntax
         * to run the program.  Also displays usage information.
         * @param msg The error message to display
         */
        void error(const char *msg)
        {
            cerr << endl
                 << "Error: " << msg << endl << endl
                 << usage() << endl
                 << "Use --help for more information" << endl;
            exit(1);
        }


        /**
         * @brief help A help message describing the syntax for all the programs options and arguments
         * @return
         */
        const string help() const
        {
            ostringstream help_str;

            help_str << shortDescription() << endl << endl
                 << longDescription() << endl << endl
                 << "Options (default value in (), *required):" << endl
                 << optionsList() << endl
                 << " -v, --verbose               Outputs additional information to stderr" << endl
                 << "     --usage                 Usage" << endl
                 << "     --help                  This message" << endl;

            return help_str.str();
        }


        void parse(int argc, char *argv[])
        {
            int c;
            int help_flag = 0;
            int usage_flag = 0;

            static struct option long_common_options[] =
            {
                {"verbose",         no_argument,        0,          'v'},
                {"help",            no_argument,        &help_flag, 1},
                {"usage",           no_argument,        &usage_flag, 1},
                {0, 0, 0, 0}
            };


            vector<option>* long_options = longOptions();

            for(uint8_t i = 0; i < 4; i++)
            {
                long_options->push_back(long_common_options[i]);
            }

            option* long_options_array = &(*long_options)[0];

            // Create short options
            ostringstream short_options_str;
            short_options_str << shortOptions() << "vuh";
            static const char *short_options = short_options_str.str().c_str();

            // Check that something is on the command line
            if (argc <= 1)
            {
                cerr << endl
                     << usage() << endl
                     << help() << endl;
                exit(1);
            }

            // Loop through the options
            while (true)
            {
                /* getopt_long stores the option index here. */
                int index = -1;

                c = getopt_long (argc, argv, short_options, long_options_array, &index);

                /* Detect the end of the options. */
                if (c == -1)
                    break;

                switch (c)
                {
                case ':':
                    cerr << "Missing required argument for "
                              << (index == -1 ? std::string(1, (char)optopt) : std::string(long_options_array[index].name))
                              << endl;
                    exit(1);
                case '?':
                    cerr << "Use --usage or --help for some help" << endl << endl;
                    exit(1);
                case 'v':
                    verbose = true;
                    break;
                }

                setOption(c, optarg);
            }

            // If the help flag was set print usage and help then exit
            if (help_flag)
            {
                cout << usage() << endl
                     << help() << endl;
                exit(0);
            }

            // If the usage flag was set print usage then exit
            if (usage_flag)
            {
                cout << usage() << endl
                     << "Use --help for more information." << endl << endl;
                exit(0);
            }

            // Find out how many arguments are remaining
            int nb_remaining_args = argc - optind;

            // Check we still have the required number left, if not error then exit
            if(nb_remaining_args < min_args)
            {
                const char* suffix = min_args > 1 ? "s" : "";

                ostringstream error_msg;
                error_msg << "Requires at least " << min_args << " argument" << suffix << ".";

                // Will die here
                error(error_msg.str().c_str());
            }

            // Add all remaining arguments to file list
            for(uint16_t i = 0; i < nb_remaining_args; i++)
            {
                remaining_args.push_back(string(argv[optind++]));
            }

            // Ask child class to process the remaining arguments
            processRemainingArgs(remaining_args);

            // Clean up
            delete long_options;
        }

        void print()
        {
            if (verbose)
                cerr << "Verbose flag set" << endl;

            cerr << currentStatus() << endl;
        }
    };

    BaseArgs::~BaseArgs() {}
}
