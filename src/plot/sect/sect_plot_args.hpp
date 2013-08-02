#pragma once

#include <getopt.h>
#include <string.h>
#include <iostream>
#include <stdint.h>

using std::string;
using std::cerr;
using std::cout;


#define DEFAULT_TITLE "Sequence Coverage Plot"
#define DEFAULT_X_LABEL "Position (nt)"
#define DEFAULT_Y_LABEL "Kmer Coverage"

#define DEFAULT_WIDTH 1024
#define DEFAULT_HEIGHT 1024


using std::cout;
using std::cerr;
using std::endl;

class SectPlotArgs
{
public:
    string*  sect_file_arg;
    string*  output_type;
    string*  output_arg;
    const char*  title;
    const char*  x_label;
    const char*  y_label;
    uint16_t width;
    uint16_t height;
    uint32_t fasta_index;
    const char* fasta_header;
    bool verbose;

    // Default constructor
    SectPlotArgs() :
        sect_file_arg(NULL), output_type(NULL), output_arg(NULL), title(DEFAULT_TITLE),
        x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL),
        width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT),
        fasta_index(0), fasta_header(NULL),
        verbose(false)
    {
        output_type = new string("png");
    }

    // Constructor that parses command line options
    SectPlotArgs(int argc, char* argv[]) :
        sect_file_arg(NULL), output_type(NULL), output_arg(NULL), title(DEFAULT_TITLE),
        x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL),
        width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT),
        fasta_index(0), fasta_header(NULL),
        verbose(false)
    {
        output_type = new string("png");
        parse(argc, argv);
    }

    ~SectPlotArgs()
    {
        delete sect_file_arg;
        delete output_type;
        delete output_arg;
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
            {"index",           required_argument,  0, 'n'},
            {"header",          required_argument,  0, 'd'},
            {"help",            no_argument,        &help_flag, 1},
            {"usage",           no_argument,        &usage_flag, 1},
            {0, 0, 0, 0}
        };

        static const char *short_options = "o:p:t:i:j:w:h:n:d:vuh";

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
                output_arg = new string(optarg);
                break;
            case 'p':
                delete output_type;
                output_type = new string(optarg);
                break;
            case 't':
                title = optarg;
                break;
            case 'i':
                x_label = optarg;
                break;
            case 'j':
                y_label = optarg;
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

        sect_file_arg = new string(argv[optind++]);
    }

    // Work out the output path to use (either user specified or auto generated)
    string determineOutputPath()
    {
        std::ostringstream output_str;
        output_str << *sect_file_arg << "." << *output_type;
        return output_arg == NULL ? output_str.str() : *output_arg;
    }


    void print()
    {
        if (verbose)
            cerr << "Verbose flag set" << endl;

        if (output_type != NULL)
            cerr << "Output type: " << output_type->c_str() << endl;

        if (output_arg != NULL)
            cerr << "Output file specified: " << output_arg->c_str() << endl;

        if (sect_file_arg != NULL)
            cerr << "SECT input file specified: " << sect_file_arg->c_str() << endl;

        if (title)
            cerr << "Plot title: " << title << endl;

        if (x_label)
            cerr << "X Label: " << x_label << endl;

        if (y_label)
            cerr << "Y Label: " << y_label << endl;

        if (width)
            cerr << "Width: " << width << endl;

        if (height)
            cerr << "Height: " << height << endl;

        if (fasta_index)
            cerr << "Fasta index to plot: " << fasta_index << endl;

        if (fasta_header)
            cerr << "Fasta header to plot: " << fasta_header << endl;

        cerr << endl;
    }
};
