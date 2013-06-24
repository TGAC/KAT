#ifndef __COMP_ARGS_HPP__
#define __COMP_ARGS_HPP__

#include <getopt.h>
#include <stdlib.h>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

class CompArgs
{
public:
    const char * 	fasta_arg;
    const char * 	db1_arg;
    const char * 	db2_arg;
    const char *  input_type;
    const char *  output_arg;
    uint_t endlimit_arg;
    float xscale_arg;
    float yscale_arg;
    bool noindex;
    uint_t threads_arg;
    bool verbose;

    // Default constructor
    CompArgs() :
        output_arg(NULL), endlimit_arg(-1), xscale_arg(1.0), yscale_arg(1.0), noindex(false), threads_arg(1), verbose(false)
    {}

    // Constructor that parses command line options
    CompArgs(int argc, char* argv[]) :
        output_arg(NULL), endlimit_arg(-1), xscale_arg(1.0), yscale_arg(1.0), noindex(false), threads_arg(1), verbose(false)
    {
        parse(argc, argv);
    }




#define seqcvg_args_USAGE "\nUsage: kat comp [options] db1_path db2_path\n"
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


#define sect_args_HELP "Compares two sets of jellyfish kmer counts\n\n" \
  "If comparing kmers from reads to kmers from an assembly, the larger (most likely the read) kmer hash should be provided first, then the assembly kmer hash second.\n\n" \
  "Options (default value in (), *required):\n" \
  " -f, --fasta          Assembly fasta file. If provided, the ends up to endsize from this sequences will be included in a separate frequency comparison matrix.\n" \
  " -e, --endsize        Size of the end section from the fasta file sequences that will be used as filters when creating the ends matrix.\n" \
  " -o, --output         *File that should contain the frequency comparison matrix from this program. If fasta file and endlimit specified, a second matrix file with the .ends suffix will be created also.\n" \
  " -x, --x_scale        Scaling factor for the first dataset - float multiplier (1.0).\n" \
  " -y, --y_scale        Scaling factor for the second dataset - float multiplier (1.0).\n" \
  " -i, --no_index       Removes the first column which contains the kmer multiplicity value for the first dataset\n" \
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
            {"fasta",           required_argument,  0, 'f'},
            {"endlimit",        optional_argument,  0, 'e'},
            {"output",          required_argument,  0, 'o'},
            {"x_scale" ,        required_argument,  0, 'x'},
            {"y_scale" ,        required_argument,  0, 'y'},
            {"noindex" ,        no_argument,        0, 'i'},
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
            case 'f':
                fasta_arg = optarg;
                break;
            case 'o':
                output_arg = optarg;
                break;
            case 't':
                threads_arg = atoi(optarg);
                break;
            case 'e':
                endlimit_arg = atoi(optarg);
                break;
            case 'x':
                xscale_arg = atof(optarg);
                break;
            case 'y':
                yscale_arg = atof(optarg);
                break;
            case 'i':
                noindex = true;
                break;

            }
        }

        // Parse arguments
        if(argc - optind != 2)
            error("Requires exactly 2 arguments.");
        db1_arg = argv[optind++];
        db2_arg = argv[optind++];
    }

    bool outputGiven()
    {
        return output_arg == NULL ? false : true;
    }


    void print()
    {
        if (verbose)
            cerr << "Verbose flag set" << endl;

        if (fasta_arg)
            cerr << "Fasta file: " << fasta_arg << endl;

        if (threads_arg)
            cerr << "Threads requested: " << threads_arg << endl;

        if (endlimit_arg)
            cerr << "Endlimit: " << endlimit_arg << endl;

        if (xscale_arg)
            cerr << "X Scale Arg: " << xscale_arg << endl;

        if (yscale_arg)
            cerr << "Y Scale Arg: " << yscale_arg << endl;

        if (noindex)
            cerr << "No index column in output: " << noindex << endl;

        if (db1_arg)
            cerr << "Jellyfish hash 1: " << db1_arg << endl;

        if (db2_arg)
            cerr << "Jellyfish hash 2: " << db2_arg << endl;

        if (output_arg)
            cerr << "Matrix Output file provided: " << output_arg << endl;

        cerr << endl;
    }

private:
};

#endif // __COMP_ARGS_HPP__

