#ifndef __ASM_PLOT_ARGS_HPP__
#define __ASM_PLOT_ARGS_HPP__

#include <getopt.h>
#include <string.h>
#include <iostream>
#include <stdint.h>

using std::string;
using std::cerr;
using std::cout;


#define DEFAULT_TITLE "Assembly duplication histogram"
#define DEFAULT_X_LABEL "Kmer Multiplicity"
#define DEFAULT_Y_LABEL "Distinct Kmer Count"

#define DEFAULT_X_MAX 1000
#define DEFAULT_Y_MAX 1000000

#define DEFAULT_WIDTH 1024
#define DEFAULT_HEIGHT 1024

#define DEFAULT_DUPLICATION 5

using std::cout;
using std::cerr;
using std::endl;

class AsmPlotArgs
{
public:
    string*  mx_arg;
    string*  output_type;
    string*  output_arg;
    const char*  title;
    const char*  x_label;
    const char*  y_label;
    uint16_t x_max;
    uint64_t y_max;
    uint16_t width;
    uint16_t height;
    bool ignore_absent;
    uint16_t max_duplication;
    string* columns;
    bool verbose;

    // Default constructor
    AsmPlotArgs() :
        mx_arg(NULL), output_type(NULL), output_arg(NULL), title(DEFAULT_TITLE),
        x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL),
        x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX),
        width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT),
        ignore_absent(false), max_duplication(DEFAULT_DUPLICATION),
        verbose(false)
    {
        output_type = new string("png");
        columns = NULL;
    }

    // Constructor that parses command line options
    AsmPlotArgs(int argc, char* argv[]) :
        mx_arg(NULL), output_type(NULL), output_arg(NULL), title(DEFAULT_TITLE),
        x_label(DEFAULT_X_LABEL), y_label(DEFAULT_Y_LABEL),
        x_max(DEFAULT_X_MAX), y_max(DEFAULT_Y_MAX),
        width(DEFAULT_WIDTH), height(DEFAULT_HEIGHT),
        ignore_absent(false), max_duplication(DEFAULT_DUPLICATION),
        verbose(false)
    {
        output_type = new string("png");
        columns = NULL;
        parse(argc, argv);
    }

    ~AsmPlotArgs()
    {
        delete mx_arg;
        delete output_type;
        delete output_arg;
    }


#define asm_plot_args_USAGE "Usage: kat plot flame [options] -o <output_file_path> matrix_path\n"
    const char * usage() const
    {
        return asm_plot_args_USAGE;
    }

    void error(const char *msg)
    {
        cerr << endl
             << "Error: " << msg << endl << endl
             << usage() << endl
             << "Use --help for more information" << endl;
        exit(1);
    }


#define asm_plot_args_HELP "Create Assembly Kmer Histograms\n\n" \
  "  Shows kmer duplication levels within an assembly.\n\n" \
  "Options (default value in (), *required):\n" \
  " -p, --output_type    The plot file type to create: png, ps, pdf.  Warning... if pdf is selected\n" \
  "                      please ensure your gnuplot installation can export pdf files. (png)\n" \
  " -o, --output         *Output file\n" \
  " -t, --title          Title for plot (" DEFAULT_TITLE ")\n" \
  " -i, --x_label        Label for the x-axis (" DEFAULT_X_LABEL ")\n" \
  " -j, --y_label        Label for the y-axis (" DEFAULT_Y_LABEL ")\n" \
  " -x  --x_max          Maximum value for the x-axis (1000)\n" \
  " -y  --y_max          Maximum value for the y-axis (10000000)\n" \
  " -w, --width          Width of canvas (1024)\n" \
  " -h, --height         Height of canvas (1024)\n" \
  " -a, --ignore_absent  Ignore kmers in reads but absent from the assembly\n" \
  " -m, --max_dup        Maximum duplication level to show in plots (5)\n" \
  " -c, --columns        Comma separated string listing columns to show in plot.  If used, this overrides \"--ignore_absent\" and \"--columns\"\n" \
  " -v, --verbose        Outputs additional information to stderr\n" \
  "     --usage          Usage\n" \
  "     --help           This message\n" \

    const char * help() const
    {
        return asm_plot_args_HELP;
    }

#define asm_plot_args_HIDDEN "Hidden options:"
    const char * hidden() const
    {
        return asm_plot_args_HIDDEN;
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
            {"x_max",           required_argument,  0, 'x'},
            {"y_max",           required_argument,  0, 'y'},
            {"ignore_absent",   no_argument,        0, 'a'},
            {"max_dup",         required_argument,  0, 'm'},
            {"columns",         required_argument,  0, 'c'},
            {"help",            no_argument,        &help_flag, 1},
            {"usage",           no_argument,        &usage_flag, 1},
            {0, 0, 0, 0}
        };

        static const char *short_options = "o:p:t:i:j:x:y:w:h:m:c:avuh";

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
            case 'x':
                x_max = atoi(optarg);
                break;
            case 'y':
                y_max = atol(optarg);
                break;
            case 'w':
                width = atoi(optarg);
                break;
            case 'h':
                height = atoi(optarg);
                break;
            case 'a':
                ignore_absent = true;
                break;
            case 'm':
                max_duplication = atoi(optarg);
                break;
            case 'c':
                columns = new string(optarg);
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

        mx_arg = new string(argv[optind++]);
    }

    bool outputGiven()
    {
        return output_arg == NULL ? false : true;
    }


    void print()
    {
        if (verbose)
            cerr << "Verbose flag set" << endl;

        if (output_type != NULL)
            cerr << "Output type: " << output_type->c_str() << endl;

        if (output_arg != NULL)
            cerr << "Output file specified: " << output_arg->c_str() << endl;

        if (mx_arg != NULL)
            cerr << "Kmer Matrix input file specified: " << mx_arg->c_str() << endl;

        if (title)
            cerr << "Plot title: " << title << endl;

        if (x_label)
            cerr << "X Label: " << x_label << endl;

        if (y_label)
            cerr << "Y Label: " << y_label << endl;

        if (x_max)
            cerr << "X Max: " << x_max << endl;

        if (y_max)
            cerr << "Y Max: " << y_max << endl;

        if (width)
            cerr << "Width: " << width << endl;

        if (height)
            cerr << "Height: " << height << endl;

        if (max_duplication)
            cerr << "Max duplication level to plot: " << max_duplication << endl;

        if (columns)
            cerr << "Columns to plot: " << columns->c_str() << endl;

        if (ignore_absent)
            cerr << "Ignore absent kmers: " << ignore_absent << endl;


        cerr << endl;
    }
};

#endif // __ASM_PLOT_ARGS_HPP__

