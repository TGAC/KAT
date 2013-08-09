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

using std::string;
using std::cerr;
using std::cout;

const string DEFAULT_OUTPUT_TYPE = "png";

#define DEFAULT_TITLE "Sequence Coverage Plot"
#define DEFAULT_X_LABEL "Position (nt)"
#define DEFAULT_Y_LABEL "Kmer Coverage"

const uint32_t DEFAULT_Y_MAX = -1;
const uint16_t DEFAULT_WIDTH = 1024;
const uint16_t DEFAULT_HEIGHT = 1024;


using std::cout;
using std::cerr;
using std::endl;

class SectPlotArgs
{
public:
    string  sect_file_arg;
    string  output_type;
    string  output_arg;
    string  title;
    string  x_label;
    string  y_label;
    uint32_t y_max;
    uint16_t width;
    uint16_t height;
    uint32_t fasta_index;
    string fasta_header;
    bool verbose;

    // Default constructor
    SectPlotArgs() :
        sect_file_arg(""), output_type(DEFAULT_OUTPUT_TYPE), output_arg(""), title(DEFAULT_TITLE),
        x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL), y_max(DEFAULT_Y_MAX),
        width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT),
        fasta_index(0), fasta_header(""),
        verbose(false)
    {
    }

    // Constructor that parses command line options
    SectPlotArgs(int argc, char* argv[]) :
        sect_file_arg(""), output_type(DEFAULT_OUTPUT_TYPE), output_arg(""), title(DEFAULT_TITLE),
        x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL), y_max(DEFAULT_Y_MAX),
        width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT),
        fasta_index(0), fasta_header(""),
        verbose(false)
    {
        parse(argc, argv);
    }


#define sect_plot_args_USAGE "Usage: kat plot sect [options] -o <output_file_path> sect_output_file\n"
    const char * usage() const
    {
        return sect_plot_args_USAGE;
    }

    void error(const char *msg)
    {
        cerr << endl
             << "Error: " << msg << endl << endl
             << usage() << endl
             << "Use --help for more information" << endl;
        exit(1);
    }


#define sect_plot_args_HELP "Create Sequence Coverage Plot\n\n" \
  "  Shows kmer coverage level across an sequence.\n\n" \
  "Options (default value in (), *required):\n" \
  " -p, --output_type    The plot file type to create: png, ps, pdf.  Warning... if pdf is selected\n" \
  "                      please ensure your gnuplot installation can export pdf files. (png)\n" \
  " -o, --output         Output file (<sect_file>.<output_type>)\n" \
  " -t, --title          Title for plot (" DEFAULT_TITLE ")\n" \
  " -i, --x_label        Label for the x-axis (" DEFAULT_X_LABEL ")\n" \
  " -j, --y_label        Label for the y-axis (" DEFAULT_Y_LABEL ")\n" \
  " -y  --y_max          Maximum value for the y-axis (Auto calculate max coverage in data)\n" \
  " -w, --width          Width of canvas (1024)\n" \
  " -h, --height         Height of canvas (1024)\n" \
  " -n, --index          Index of fasta entry to plot.  First entry is 1. (0)\n" \
  " -d, --header         Fasta header of fasta entry to plot.  Has priority over \'--index\'.\n" \
  " -v, --verbose        Outputs additional information to stderr\n" \
  "     --usage          Usage\n" \
  "     --help           This message\n" \

    const char * help() const
    {
        return sect_plot_args_HELP;
    }

#define sect_plot_args_HIDDEN "Hidden options:"
    const char * hidden() const
    {
        return sect_plot_args_HIDDEN;
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
            {"y_max",           required_argument,  0, 'y'},
            {"width",           required_argument,  0, 'w'},
            {"height",          required_argument,  0, 'h'},
            {"index",           required_argument,  0, 'n'},
            {"header",          required_argument,  0, 'd'},
            {"help",            no_argument,        &help_flag, 1},
            {"usage",           no_argument,        &usage_flag, 1},
            {0, 0, 0, 0}
        };

        static const char *short_options = "o:p:t:i:j:y:w:h:n:d:vuh";

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
            case 'y':
                y_max = atoi(optarg);
                break;
            case 'w':
                width = atoi(optarg);
                break;
            case 'h':
                height = atoi(optarg);
                break;
            case 'n':
                fasta_index = atoi(optarg);
                break;
            case 'd':
                fasta_header = optarg;
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

        sect_file_arg = string(argv[optind++]);
    }

    // Work out the output path to use (either user specified or auto generated)
    string determineOutputPath()
    {
        std::ostringstream output_str;
        output_str << sect_file_arg << "." << output_type;
        return output_arg.empty() ? output_str.str() : output_arg;
    }

    string autoTitle(string& header)
    {
        std::ostringstream output_str;
        output_str << DEFAULT_TITLE << ": " << sect_file_arg << " - " << header;
        return title.compare(DEFAULT_TITLE) == 0 ? output_str.str() : title;
    }



    void print()
    {
        if (verbose)
            cerr << "Verbose flag set" << endl;

        cerr << "Output type: " << output_type.c_str() << endl;
        cerr << "Output file specified: " << output_arg.c_str() << endl;
        cerr << "SECT input file specified: " << sect_file_arg.c_str() << endl;
        cerr << "Plot title: " << title << endl;
        cerr << "X Label: " << x_label << endl;
        cerr << "Y Label: " << y_label << endl;
        cerr << "Y Max: " << y_max << endl;
        cerr << "Width: " << width << endl;
        cerr << "Height: " << height << endl;
        cerr << "Fasta index to plot: " << fasta_index << endl;
        cerr << "Fasta header to plot: " << fasta_header << endl;

        cerr << endl;
    }
};
