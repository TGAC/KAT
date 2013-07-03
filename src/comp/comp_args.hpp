#ifndef __COMP_ARGS_HPP__
#define __COMP_ARGS_HPP__

#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include <stdint.h>

using std::cout;
using std::cerr;
using std::endl;

#define DEFAULT_D1_SCALE 1.0
#define DEFAULT_D2_SCALE 1.0
#define DEFAULT_THREADS 1
#define DEFAULT_OUTPUT_PREFIX "./kat_comp_output"
#define DEFAULT_D1_BINS 1001
#define DEFAULT_D2_BINS 1001

class CompArgs
{
public:
    const char * 	db1_arg;
    const char * 	db2_arg;
    const char * 	db3_arg;

    const char *  output_prefix_arg;
    double d1_scale_arg;
    double d2_scale_arg;
    uint16_t d1_bins;
    uint16_t d2_bins;
    uint16_t threads_arg;
    bool verbose;

    // Default constructor
    CompArgs() :
        output_prefix_arg(DEFAULT_OUTPUT_PREFIX), d1_scale_arg(DEFAULT_D1_SCALE), d2_scale_arg(DEFAULT_D2_SCALE),
        d1_bins(DEFAULT_D1_BINS), d2_bins(DEFAULT_D2_BINS),
        threads_arg(DEFAULT_THREADS), verbose(false)
    {}

    // Constructor that parses command line options
    CompArgs(int argc, char* argv[]) :
        output_prefix_arg(DEFAULT_OUTPUT_PREFIX), d1_scale_arg(DEFAULT_D1_SCALE), d2_scale_arg(DEFAULT_D2_SCALE),
        d1_bins(DEFAULT_D1_BINS), d2_bins(DEFAULT_D2_BINS),
        threads_arg(DEFAULT_THREADS), verbose(false)
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
  " -o, --output_prefix  Path prefix for files produced by this program (" DEFAULT_OUTPUT_PREFIX ")\n" \
  "specified, a second matrix file with the .ends suffix will be created also.\n" \
  " -x, --d1_scale       Scaling factor for the first dataset - float multiplier (1.0).  Max value: 1.0\n" \
  " -y, --d2_scale       Scaling factor for the second dataset - float multiplier (1.0).  Max value: 1.0.\n" \
  " -i, --d1_bins        Number of bins for the first dataset.  i.e. number of rows in the matrix (1001)\n" \
  " -j, --d2_bins        Number of bins for the second dataset.  i.e. number of columns in the matrix (1001)\n" \
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
        int help_flag = 0;
        int usage_flag = 0;

        static struct option long_options[] =
        {
            {"verbose",         no_argument,        0, 'v'},
            {"output_prefix",   required_argument,  0, 'o'},
            {"d1_scale" ,       required_argument,  0, 'x'},
            {"d2_scale" ,       required_argument,  0, 'y'},
            {"d1_bins",         required_argument,  0, 'i'},
            {"d2_bins",         required_argument,  0, 'j'},
            {"threads",         required_argument,  0, 't'},
            {"help",            no_argument,        &help_flag, 1},
            {"usage",           no_argument,        &usage_flag, 1},
            {0, 0, 0, 0}
        };

        static const char *short_options = "f:e:o:x:y:i:j:t:vuh";

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
                output_prefix_arg = optarg;
                break;
            case 't':
                threads_arg = atoi(optarg);
                break;
            case 'x':
                d1_scale_arg = atof(optarg);
                if (d1_scale_arg > 1.0)
                    d1_scale_arg = 1.0;
                break;
            case 'y':
                d2_scale_arg = atof(optarg);
                if (d2_scale_arg > 1.0)
                    d2_scale_arg = 1.0;
                break;
            case 'i':
                d1_bins = atoi(optarg);
                break;
            case 'j':
                d2_bins = atoi(optarg);
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

        if (d1_scale_arg)
            cerr << "Dataset 1 Scaling Factor: " << d1_scale_arg << endl;

        if (d2_scale_arg)
            cerr << "Dataset 2 Scaling Factor: " << d2_scale_arg << endl;

        if (d1_bins)
            cerr << "Number of Dataset 1 bins: " << d1_bins << endl;

        if (d2_bins)
            cerr << "Number of Dataset 2 bins: " << d2_bins << endl;

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

