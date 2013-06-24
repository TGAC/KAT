#ifndef __FLAME_PLOT_ARGS_HPP__
#define __FLAME_PLOT_ARGS_HPP__

#include <getopt.h>
#include <string.h>
#include <iostream>

using std::string;
using std::cerr;
using std::cout;
using std::endl;


#define DEFAULT_TITLE "Flame plot"
#define DEFAULT_X_LABEL "Hash1 Multiplicity"
#define DEFAULT_Y_LABEL "Hash2 Multiplicity"



class FlamePlotArgs
{
public:
    string*  mx_arg;
    string*  output_type;
    string*  output_arg;
    const char*  title;
    const char*  xlabel;
    const char*  ylabel;
    bool verbose;

    // Default constructor
    FlamePlotArgs() :
        mx_arg(NULL), output_type(NULL), output_arg(NULL), title(DEFAULT_TITLE), xlabel(DEFAULT_X_LABEL), ylabel(DEFAULT_Y_LABEL), verbose(false)
    {
        output_type = new string("png");
    }

    // Constructor that parses command line options
    FlamePlotArgs(int argc, char* argv[]) :
        mx_arg(NULL), output_type(NULL), output_arg(NULL), title(DEFAULT_TITLE), xlabel(DEFAULT_X_LABEL), ylabel(DEFAULT_Y_LABEL), verbose(false)
    {
        output_type = new string("png");
        parse(argc, argv);
    }

    ~FlamePlotArgs()
    {
        delete mx_arg;
        delete output_type;
        delete output_arg;
    }


    string getOutputWithoutExt()
    {
        size_t result = output_arg->find(".png", output_arg->length()-5);

        return result == string::npos ? *output_arg : output_arg->substr(0, result);
    }


#define flame_plot_args_USAGE "\nUsage: kat plot flame [options] matrix_path\n"
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


#define flame_plot_args_HELP "Create Kmer Flame Plots\n\n" \
  "  Creates a flame plot from a matrix created with the \"comp\" tool.  Typically this\n" \
  "  is used to compare two kmer hashes produced by different NGS reads.\n\n" \
  "Options (default value in (), *required):\n" \
  " -p, --output_type    The plot file type to create: png, ps, pdf.  Warning... if pdf is selected\n" \
  "                      please ensure your gnuplot installation can export pdf files. (png)\n" \
  " -o, --output         *Output file\n" \
  " -t, --title          Title for plot (" DEFAULT_TITLE ")\n" \
  " -x, --x_label        Label for the x-axis (" DEFAULT_X_LABEL ")\n" \
  " -y, --y_label        Label for the y-axis (" DEFAULT_Y_LABEL ")\n" \
  " -v, --verbose        Outputs additional information to stderr\n" \
  "     --usage          Usage\n" \
  "     --help           This message\n" \

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

        static struct option long_options[] =
        {
            {"verbose",         no_argument,        0, 'v'},
            {"output_type",     required_argument,  0, 'p'},
            {"output",          required_argument,  0, 'o'},
            {"title",           required_argument,  0, 't'},
            {"x_label",         required_argument,  0, 'x'},
            {"y_label",         required_argument,  0, 'y'},
            {"help",            no_argument,        0, 'h'},
            {"usage",           no_argument,        0, 'u'},
            {0, 0, 0, 0}
        };

        static const char *short_options = "o:p:t:x:y:vuh";

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
            case 'x':
                xlabel = optarg;
                break;
            case 'y':
                ylabel = optarg;
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
            cerr << "Verbose flag set\n";

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

        cerr << endl;
    }
};

#endif // __FLAME_PLOT_ARGS_HPP__

