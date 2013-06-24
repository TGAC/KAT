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
#define DEFAULT_Y_MAX 10000000

#define DEFAULT_X_LOGSCALE 0
#define DEFAULT_Y_LOGSCALE 0

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
    const char*  xlabel;
    const char*  ylabel;
    uint16_t xmax;
    uint64_t ymax;
    bool xlogscale;
    bool ylogscale;
    bool ignoreAbsent;
    bool verbose;

    // Default constructor
    AsmPlotArgs() :
        mx_arg(NULL), output_type(NULL), output_arg(NULL), title(DEFAULT_TITLE),
        xlabel(DEFAULT_X_LABEL), ylabel(DEFAULT_Y_LABEL),
        xmax(DEFAULT_X_MAX), ymax(DEFAULT_Y_MAX),
        xlogscale(DEFAULT_X_LOGSCALE), ylogscale(DEFAULT_Y_LOGSCALE),
        ignoreAbsent(false), verbose(false)
    {
        output_type = new string("png");
    }

    // Constructor that parses command line options
    AsmPlotArgs(int argc, char* argv[]) :
        mx_arg(NULL), output_type(NULL), output_arg(NULL), title(DEFAULT_TITLE),
        xlabel(DEFAULT_X_LABEL), ylabel(DEFAULT_Y_LABEL),
        xmax(DEFAULT_X_MAX), ymax(DEFAULT_Y_MAX),
        xlogscale(DEFAULT_X_LOGSCALE), ylogscale(DEFAULT_Y_LOGSCALE),
        ignoreAbsent(false), verbose(false)
    {
        output_type = new string("png");
        parse(argc, argv);
    }

    ~AsmPlotArgs()
    {
        delete mx_arg;
        delete output_type;
        delete output_arg;
    }


    string getOutputWithoutExt()
    {
        size_t result = output_arg->find(output_type->c_str(), output_arg->length() - output_type->length() - 2);

        return result == string::npos ? *output_arg : output_arg->substr(0, result);
    }


#define asm_plot_args_USAGE "Usage: kat plot flame [options] matrix_path\n"
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
  " -k  --x_logscale     Use logscale for x-axis (false)\n" \
  " -l  --y_logscale     Use logscale for y-axis (false)\n" \
  " -a  --ignore_absent  Ignore kmers in reads but absent from the assembly\n" \
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
            {"x_logscale",      no_argument,        0, 'k'},
            {"y_logscale",      required_argument,  0, 'l'},
            {"ignore_absent",   no_argument,        0, 'a'},
            {"help",            no_argument,        0, 'h'},
            {"usage",           no_argument,        0, 'u'},
            {0, 0, 0, 0}
        };

        static const char *short_options = "o:p:t:i:j:x:y:klavuh";

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
            case 'h':
                cout << usage() << endl
                     << help() << endl;
                exit(0);
            case 'u':
                cout << usage() << endl
                     << "Use --help for more information." << endl << endl;
                exit(0);
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
                xlabel = optarg;
                break;
            case 'j':
                ylabel = optarg;
                break;
            case 'k':
                xlogscale = true;
                break;
            case 'l':
                ylogscale = true;
                break;
            case 'x':
                xmax = atoi(optarg);
                break;
            case 'y':
                ymax = atol(optarg);
                break;
            case 'a':
                ignoreAbsent = true;
                break;

            }
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

        if (xlabel)
            cerr << "X Label: " << xlabel << endl;

        if (ylabel)
            cerr << "Y Label: " << ylabel << endl;

        if (xmax)
            cerr << "X Max: " << xmax << endl;

        if (ymax)
            cerr << "Y Max: " << ymax << endl;

        if (xlogscale)
            cerr << "X Logscale: " << xlogscale << endl;

        if (ylogscale)
            cerr << "Y Logscale: " << ylogscale << endl;

        if (ignoreAbsent)
            cerr << "Ingore absent kmers: " << ignoreAbsent << endl;


        cerr << endl;
    }
};

#endif // __ASM_PLOT_ARGS_HPP__

