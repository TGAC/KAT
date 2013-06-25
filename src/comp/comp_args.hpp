#ifndef __COMP_ARGS_HPP__
#define __COMP_ARGS_HPP__

#include <getopt.h>
#include <stdlib.h>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

#define DEFAULT_X_SCALE 1.0
#define DEFAULT_Y_SCALE 1.0
#define DEFAULT_THREADS 1
#define DEFAULT_OUTPUT_PREFIX "./kat_comp_output"


class CompArgs
{
public:
    const char * 	db1_arg;
    const char * 	db2_arg;
    const char * 	db3_arg;

    const char *  output_prefix_arg;
    float xscale_arg;
    float yscale_arg;
    uint_t threads_arg;
    bool verbose;

    // Default constructor
    CompArgs() :
        output_prefix_arg(DEFAULT_OUTPUT_PREFIX), xscale_arg(DEFAULT_X_SCALE), yscale_arg(DEFAULT_Y_SCALE), threads_arg(DEFAULT_THREADS), verbose(false)
    {}

    // Constructor that parses command line options
    CompArgs(int argc, char* argv[]) :
        output_prefix_arg(DEFAULT_OUTPUT_PREFIX), xscale_arg(DEFAULT_X_SCALE), yscale_arg(DEFAULT_Y_SCALE), threads_arg(DEFAULT_THREADS), verbose(false)
    {
        parse(argc, argv);
    }




#define seqcvg_args_USAGE "Usage: kat comp [options] db1_path db2_path [db3_path]\n"
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


#define sect_args_HELP "Compares sets of jellyfish kmer counts.  Normally, this will be two sets of kmers, although " \
  "when comparing kmers from assemblies, it is possible to add a 3rd set of kmers which describe the ends of each assembly " \
  "sequence to determine issues at end of sequences that may have prevented the assembler joining sequences together.\n\n" \
  "If comparing kmers from reads to kmers from an assembly, the larger (most likely the read) kmer hash should be provided " \
  "first, then the assembly kmer hash second.  Finally, if comparing kmers at the ends of sequences this should be supplied last.\n\n" \
  "Options (default value in (), *required):\n" \
  " -o, --output_prefix  Path prefix for files produced by this program (./kat_comp_output)\n" \
  "specified, a second matrix file with the .ends suffix will be created also.\n" \
  " -x, --x_scale        Scaling factor for the first dataset - float multiplier (1.0).\n" \
  " -y, --y_scale        Scaling factor for the second dataset - float multiplier (1.0).\n" \
  " -t, --threads        The number of threads to use (1)\n" \
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

        static struct option long_options[] =
        {
            {"verbose",         no_argument,        0, 'v'},
            {"output_prefix",   required_argument,  0, 'o'},
            {"x_scale" ,        required_argument,  0, 'x'},
            {"y_scale" ,        required_argument,  0, 'y'},
            {"threads",         required_argument,  0, 't'},
            {"help",            no_argument,        0, 'h'},
            {"usage",           no_argument,        0, 'u'},
            {0, 0, 0, 0}
        };

        static const char *short_options = "f:e:o:x:y:t:vuh";

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
                output_prefix_arg = optarg;
                break;
            case 't':
                threads_arg = atoi(optarg);
                break;
            case 'x':
                xscale_arg = atof(optarg);
                break;
            case 'y':
                yscale_arg = atof(optarg);
                break;

            }
        }

        // Parse arguments
        int remaining_args = argc - optind;

        if (verbose)
            cerr << "Found " << remaining_args << " remaining arguments on the command line." << endl;

        if(remaining_args < 2 || remaining_args > 3)
            error("Requires 2 or 3 arguments describing paths to jellyfish kmer counts.");

        db1_arg = argv[optind++];
        db2_arg = argv[optind++];        
        db3_arg = remaining_args == 3 ? argv[optind++] : NULL;
    }

    void print()
    {
        if (verbose)
            cerr << "Verbose flag set" << endl;

        if (threads_arg)
            cerr << "Threads requested: " << threads_arg << endl;

        if (xscale_arg)
            cerr << "X Scale Arg: " << xscale_arg << endl;

        if (yscale_arg)
            cerr << "Y Scale Arg: " << yscale_arg << endl;

        if (db1_arg)
            cerr << "Jellyfish hash 1: " << db1_arg << endl;

        if (db2_arg)
            cerr << "Jellyfish hash 2: " << db2_arg << endl;

        if (db3_arg)
            cerr << "Jellyfish hash 3: " << db3_arg << endl;

        if (output_prefix_arg)
            cerr << "Output file path prefix: " << output_prefix_arg << endl;

        cerr << endl;
    }

private:
};

#endif // __COMP_ARGS_HPP__

