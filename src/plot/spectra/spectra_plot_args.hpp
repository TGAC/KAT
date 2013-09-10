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
#include <string.h>
#include <iostream>
#include <stdint.h>
#include <vector>

using std::string;
using std::cerr;
using std::cout;

const string DEFAULT_OUTPUT_TYPE = "png";

#define DEFAULT_OUTPUT_FILE_PREFIX "kat-kmer-spectra"
#define DEFAULT_TITLE "K-mer Spectra"
#define DEFAULT_X_LABEL "K-mer Multiplicity"
#define DEFAULT_Y_LABEL "Distinct K-mers"

const uint32_t DEFAULT_X_MIN = 0;
const uint32_t DEFAULT_Y_MIN = 0;
const uint32_t DEFAULT_X_MAX = 10000;
const uint32_t DEFAULT_Y_MAX = 1000000;
const bool DEFAULT_X_LOGSCALE = false;
const bool DEFAULT_Y_LOGSCALE = false;
const uint16_t DEFAULT_WIDTH = 1024;
const uint16_t DEFAULT_HEIGHT = 1024;

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;

namespace kat
{
    class SpectraPlotArgs
    {
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
        bool verbose;

        // Default constructor
        SpectraPlotArgs() :
            histo_paths(vector<string>()), output_type(DEFAULT_OUTPUT_TYPE), output_arg(""), title(DEFAULT_TITLE),
            x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL),
            x_min(DEFAULT_X_MIN), y_min(DEFAULT_Y_MIN), x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX),
            x_logscale(DEFAULT_X_LOGSCALE), y_logscale(DEFAULT_Y_LOGSCALE),
            width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT),
            verbose(false)
        {
        }

        // Constructor that parses command line options
        SpectraPlotArgs(int argc, char* argv[]) :
            histo_paths(vector<string>()), output_type(DEFAULT_OUTPUT_TYPE), output_arg(""), title(DEFAULT_TITLE),
            x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL),
            x_min(DEFAULT_X_MIN), y_min(DEFAULT_Y_MIN), x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX),
            x_logscale(DEFAULT_X_LOGSCALE), y_logscale(DEFAULT_Y_LOGSCALE),
            width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT),
            verbose(false)
        {
            parse(argc, argv);
        }


    #define spectra_plot_args_USAGE "Usage: kat plot spectra [options] -o <output_file_path> [list of space separated histo files]\n"
        const char * usage() const
        {
            return spectra_plot_args_USAGE;
        }

        void error(const char *msg)
        {
            cerr << endl
                 << "Error: " << msg << endl << endl
                 << usage() << endl
                 << "Use --help for more information" << endl;
            exit(1);
        }


    #define spectra_plot_args_HELP "Create K-mer Spectra Plot\n\n" \
      "  Shows K-mer spectras from kat-histo or jellyfish-histo output.\n\n" \
      "Options (default value in (), *required):\n" \
      " -p, --output_type    The plot file type to create: png, ps, pdf.  Warning... if pdf is selected\n" \
      "                      please ensure your gnuplot installation can export pdf files. (png)\n" \
      " -o, --output         Output file (" DEFAULT_OUTPUT_FILE_PREFIX ".png)\n" \
      " -t, --title          Title for plot (" DEFAULT_TITLE ")\n" \
      " -i, --x_label        Label for the x-axis (" DEFAULT_X_LABEL ")\n" \
      " -j, --y_label        Label for the y-axis (" DEFAULT_Y_LABEL ")\n" \
      " -r  --x_min          Minimum value for the x-axis (0)\n" \
      " -s  --y_min          Minimum value for the y-axis (0)\n" \
      " -x  --x_max          Maximum value for the x-axis (1000)\n" \
      " -y  --y_max          Maximum value for the y-axis (Auto calculate max y in data)\n" \
      " -l  --x_logscale     X-axis is logscale.  This overrides the x_min and x_max limits. (false)\n" \
      " -m  --y_logscale     Y-axis is logscale.  This overrides the y_min and y_max limits. (false)\n" \
      " -w, --width          Width of canvas (1024)\n" \
      " -h, --height         Height of canvas (1024)\n" \
      " -v, --verbose        Outputs additional information to stderr\n" \
      "     --usage          Usage\n" \
      "     --help           This message\n" \

        const char * help() const
        {
            return spectra_plot_args_HELP;
        }

    #define spectra_plot_args_HIDDEN "Hidden options:"
        const char * hidden() const
        {
            return spectra_plot_args_HIDDEN;
        }


        void parse(int argc, char *argv[])
        {
            int c;
            int help_flag = 0;
            int usage_flag = 0;

            static struct option long_options[] =
            {
                {"verbose",         no_argument,        0, 'v'},
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
                {"header",          required_argument,  0, 'd'},
                {"help",            no_argument,        &help_flag, 1},
                {"usage",           no_argument,        &usage_flag, 1},
                {0, 0, 0, 0}
            };

            static const char *short_options = "o:p:t:i:j:r:s:x:y:lmw:h:n:d:vuh";

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
                case 'o':
                    output_arg = string(optarg);
                    break;
                case 'p':
                    output_type = string(optarg);
                    break;
                case 't':
                    title = string(optarg);
                    break;
                case 'i':
                    x_label = string(optarg);
                    break;
                case 'j':
                    y_label = string(optarg);
                    break;
                case 'r':
                    x_min = atoi(optarg);
                    break;
                case 's':
                    y_min = atoi(optarg);
                    break;
                case 'x':
                    x_max = atoi(optarg);
                    break;
                case 'y':
                    y_max = atoi(optarg);
                    break;
                case 'l':
                    x_logscale = true;
                    break;
                case 'm':
                    y_logscale = true;
                    break;
                case 'w':
                    width = atoi(optarg);
                    break;
                case 'h':
                    height = atoi(optarg);
                    break;
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

            int remaining_args = argc - optind;

            // Parse arguments
            if(remaining_args < 1)
                error("Requires at least 1 argument.");


            // Add all remaining arguments to file list
            for(uint8_t i = 0; i < remaining_args; i++)
            {
                histo_paths.push_back(string(argv[optind++]));
            }
        }

        // Work out the output path to use (either user specified or auto generated)
        string determineOutputPath()
        {
            std::ostringstream output_str;
            output_str << DEFAULT_OUTPUT_FILE_PREFIX << "." << output_type;
            return output_arg.empty() ? output_str.str() : output_arg;
        }


        void print()
        {
            if (verbose)
                cerr << "Verbose flag set" << endl;

            cerr << "Output type: " << output_type.c_str() << endl;
            cerr << "Output file specified: " << output_arg.c_str() << endl;
            //cerr << "SECT input file specified: " << sect_file_arg.c_str() << endl;
            cerr << "Plot title: " << title << endl;
            cerr << "X Label: " << x_label << endl;
            cerr << "Y Label: " << y_label << endl;
            cerr << "X Min: " << x_min << endl;
            cerr << "Y Min: " << y_min << endl;
            cerr << "X Max: " << x_max << endl;
            cerr << "Y Max: " << y_max << endl;
            cerr << "X Logscale: " << x_logscale << endl;
            cerr << "Y Logscale: " << y_logscale << endl;
            cerr << "Width: " << width << endl;
            cerr << "Height: " << height << endl;

            cerr << endl;
        }
    };
}
