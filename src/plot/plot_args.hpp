#ifndef __PLOT_ARGS_HPP__
#define __PLOT_ARGS_HPP__

#include <getopt.h>
#include <string.h>
#include <iostream>

using std::string;
using std::cerr;
using std::cout;

class PlotArgs
{
public:
    string  mode_arg;
    const char *  output_arg;
    bool verbose;

    // Default constructor
    PlotArgs() :
         output_arg(NULL), verbose(false)
    {}

    // Constructor that parses command line options
    PlotArgs(int argc, char* argv[]) :
        output_arg(NULL), verbose(false)
    {
        parse(argc, argv);
    }


#define seqcvg_args_USAGE "\nUsage: kat plot <mode> [options] matrix_path"
    const char * usage() const
    {
        return seqcvg_args_USAGE;
    }

    void error(const char *msg)
    {
        cerr << "Error: " << msg << "\n" << usage()
                  << "\nUse --help for more information"
                  << std::endl;
        exit(1);
    }


#define sect_args_HELP "Create Kmer Plots\n\n" \
  "First argument should be the  plot mode you wish to use:\n\n" \
  "   - flame:           Creates a flame plot from a matrix created with the \"comp\" tool.  Typically this\n" \
  "                      is used to compare two kmer hashes produced by different NGS reads.\n" \
  "   - asm:             Creates a histogram using a matrix created with the \"comp\" tool.  Typically\n" \
  "                      this is used to compare a jellyfish hash produced from a read set to a jellyfish\n" \
  "                      hash produced from an assembly. The plot shows the amount of distinct kmers absent\n" \
  "                      in the assembly, those that are found once, and those found more (up to 10\n" \
  "                      duplications plotted)\n\n" \
  "Options (default value in (), *required):\n" \
  " -o, --output         Output file in PNG format\n" \
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

    bool validMode(char* mode_str)
    {
        return (strcmp(mode_str, "flame") == 0 ||
                strcmp(mode_str, "asm") == 0) ?
                    true : false;
    }


    void parse(int argc, char *argv[])
    {
        int c;

        static struct option long_options[] =
        {
            {"verbose",         no_argument,        0, 'v'},
            {"output",          required_argument,  0, 'o'},
            {"help",            no_argument,        0, 'h'},
            {"usage",           no_argument,        0, 'u'},
            {0, 0, 0, 0}
        };

        static const char *short_options = "o:vuh";

        if (argc <= 0)
        {
            cout << usage() << "\n\n" << help() << std::endl;
            exit(1);
        }
        else if (validMode(argv[1])) {
            mode_arg = argv[1];
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
                          << std::endl;
                exit(1);
            case 'h':
                cout << usage() << "\n\n" << help() << std::endl;
                exit(0);
            case 'u':
                cout << usage() << "\nUse --help for more information." << std::endl;
                exit(0);
            case '?':
                cerr << "Use --usage or --help for some help\n";
                exit(1);
            case 'v':
                verbose = true;
                break;
            case 'o':
                output_arg = optarg;
                break;

            }
        }
    }

    bool outputGiven()
    {
        return output_arg == NULL ? false : true;
    }


    void print()
    {
        if (verbose)
            cerr << "Verbose flag set\n";

        if (output_arg)
            cerr << "Matrix Output file provided: " << output_arg << "\n";

        cerr << "\n";
    }

private:
};

#endif // __PLOT_ARGS_HPP__

