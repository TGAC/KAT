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
#include <stdint.h>
#include <iostream>

using std::string;
using std::cerr;
using std::cout;
using std::endl;

namespace kat
{
    const string DEFAULT_OUTPUT_TYPE = "png";

    const string DEFAULT_TITLE      = "Title";
    const string DEFAULT_X_LABEL    = "X";
    const string DEFAULT_Y_LABEL    = "Y";
    const string DEFAULT_Z_LABEL    = "Z";
    const uint16_t DEFAULT_WIDTH    = 1024;
    const uint16_t DEFAULT_HEIGHT   = 1024;
    const int32_t DEFAULT_X_MAX     = 1000;
    const int32_t DEFAULT_Y_MAX     = 1000;
    const int64_t DEFAULT_Z_MAX     = 10000;


    class FlamePlotArgs
    {
    public:
        string      mx_arg;
        string      output_type;
        string      output_path;
        string      title;
        string      x_label;
        string      y_label;
        string      z_label;
        int16_t     x_max;
        int16_t     y_max;
        int64_t     z_max;
        uint16_t    width;
        uint16_t    height;
        bool        verbose;

        // Default constructor
        FlamePlotArgs() :
            mx_arg(""), output_type(DEFAULT_OUTPUT_TYPE), output_path(""), title(DEFAULT_TITLE), x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL), z_label(DEFAULT_Z_LABEL),
            x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX), z_max(DEFAULT_Z_MAX), width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT), verbose(false)
        {}

        // Constructor that parses command line options
        FlamePlotArgs(int argc, char* argv[]) :
            mx_arg(""), output_type(DEFAULT_OUTPUT_TYPE), output_path(""), title(DEFAULT_TITLE), x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL), z_label(DEFAULT_Z_LABEL),
            x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX), z_max(DEFAULT_Z_MAX), width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT), verbose(false)
        {
            parse(argc, argv);
        }

        ~FlamePlotArgs()
        {}


    #define flame_plot_args_USAGE "Usage: kat plot flame [options] -o <output_file_path> matrix_path\n"
        const char * usage() const
        {
            return flame_plot_args_USAGE;
        }

        void error(const char *msg)
        {
            cerr << endl
                 << "Error: " << msg << endl << endl
                 << usage() << endl
                 << "Use --help for more information" << endl;
            exit(1);
        }


    #define flame_plot_args_HELP "Create K-mer Flame Plots\n\n" \
      "  Creates a scatter plot, where the \"heat\" in each point represents the number of\n" \
      "  distinct K-mers at that point.  Typically this is used to visualise a matrix produced\n" \
      "  by the \"kat comp\" tool to compare multiplicities from two K-mer hashes produced by\n" \
      "  different NGS reads, or to visualise the GC vs K-mer multiplicity matricies produced\n" \
      "  by the \"kat gcp\" tool.\n\n" \
      "Options (default value in (), *required):\n" \
      " -p, --output_type    The plot file type to create: png, ps, pdf.  Warning... if pdf is selected\n" \
      "                      please ensure your gnuplot installation can export pdf files. (png)\n" \
      " -o, --output         Output file (<matrix_path>.<output_type>)\n" \
      " -t, --title          Title for plot (\"Title\", or value from matrix metadata if present)\n" \
      " -i, --x_label        Label for the x-axis (\"X\", or value from matrix metadata if present)\n" \
      " -j, --y_label        Label for the y-axis (\"Y\", or value from matrix metadata if present)\n" \
      " -k, --z_label        Label for the z-axis (\"Z\", or value from matrix metadata if present)\n" \
      " -x, --x_max          Maximum value for the x-axis (1000, or value from matrix metadata if present)\n" \
      " -y  --y_max          Maximum value for the y-axis (1000, or value from matrix metadata if present)\n" \
      " -z, --z_max          Cap for matrix values.  Values greater than this cap will be displayed at maximum intensity, i.e. white. (10000, or value from matrix metadata if present)\n" \
      " -w, --width          Width of canvas (1024)\n" \
      " -h, --height         Height of canvas (1024)\n" \
      " -v, --verbose        Outputs additional information to stderr\n" \
      "     --usage          Usage\n" \
      "     --help           This message\n"

        const char * help() const
        {
            return flame_plot_args_HELP;
        }

    #define flame_plot_args_HIDDEN "Hidden options:"
        const char * hidden() const
        {
            return flame_plot_args_HIDDEN;
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
                {"z_label",         required_argument,  0, 'k'},
                {"x_max",           required_argument,  0, 'x'},
                {"y_max",           required_argument,  0, 'y'},
                {"z_max",           required_argument,  0, 'z'},
                {"width",           required_argument,  0, 'w'},
                {"height",          required_argument,  0, 'h'},
                {"help",            no_argument,        &help_flag, 1},
                {"usage",           no_argument,        &usage_flag, 1},
                {0, 0, 0, 0}
            };

            static const char *short_options = "o:p:t:i:j:k:x:y:w:h:z:vuh";

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
                              << endl << endl;
                    exit(1);
                case '?':
                    cerr << "Use --usage or --help for some help" << endl << endl;
                    exit(1);
                case 'v':
                    verbose = true;
                    break;
                case 'o':
                    output_path = string(optarg);
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
                case 'k':
                    z_label = string(optarg);
                    break;
                case 'x':
                    x_max = atoi(optarg);
                    break;
                case 'y':
                    y_max = atoi(optarg);
                    break;
                case 'z':
                    z_max = atoi(optarg);
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

            // Parse arguments
            if(argc - optind != 1)
                error("Requires exactly 1 argument.");

            mx_arg = string(argv[optind++]);
        }

        // Work out the output path to use (either user specified or auto generated)
        string determineOutputPath()
        {
            std::ostringstream output_str;
            output_str << mx_arg << "." << output_type;
            return output_path.empty() ? output_str.str() : output_path;
        }


        void print()
        {
            if (verbose)
                cerr << "Verbose flag set\n";

            cerr << "Output type: " << output_type.c_str() << endl;
            cerr << "Output file specified: " << output_path.c_str() << endl;
            cerr << "K-mer matrix input file specified: " << mx_arg.c_str() << endl;
            cerr << "Plot title: " << title << endl;
            cerr << "X Label: " << x_label << endl;
            cerr << "Y Label: " << y_label << endl;
            cerr << "Z Label: " << z_label << endl;
            cerr << "X Max: " << x_max << endl;
            cerr << "Y Max: " << y_max << endl;
            cerr << "Z Max: " << z_max << endl;
            cerr << "Width: " << width << endl;
            cerr << "Height: " << height << endl;

            cerr << endl;
        }
    };
}
