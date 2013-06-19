#ifndef __COMP_ARGS_HPP__
#define __COMP_ARGS_HPP__

#include <getopt.h>
#include <stdlib.h>
#include <iostream>

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
    uint_t threads_arg;
    uint_t kmer_arg;
    bool verbose;

    // Default constructor
    CompArgs() :
        output_arg(NULL), endlimit_arg(-1), xscale_arg(1.0), yscale_arg(1.0), threads_arg(1), kmer_arg(31), verbose(false)
    {}

    // Constructor that parses command line options
    CompArgs(int argc, char* argv[]) :
        output_arg(NULL), endlimit_arg(-1), xscale_arg(1.0), yscale_arg(1.0), threads_arg(1), kmer_arg(31), verbose(false)
    {
        parse(argc, argv);
    }




#define seqcvg_args_USAGE "Usage: kat comp [options] -f <fasta_file_path> db1_path db2_path"
    const char * usage() const
    {
        return seqcvg_args_USAGE;
    }

    void error(const char *msg)
    {
        std::cerr << "Error: " << msg << "\n" << usage()
                  << "\nUse --help for more information"
                  << std::endl;
        exit(1);
    }


#define sect_args_HELP "Compares two sets of jellyfish kmer counts\n\n" \
  "If comparing kmers from reads to kmers from an assembly, the larger (most likely the read) kmer hash should be provided first, then the assembly kmer hash second." \
  "Options (default value in (), *required):\n" \
  " -f, --fasta          Assembly fasta file. If provided, the ends up to endsize from this sequences will be included in a separate frequency comparisson matrix.\n" \
  " -e, --endsize        Size of the end section from the fasta file sequences that will be used as filters when creating the ends matrix.\n" \
  " -o, --output         File that should contain the frequency comparison matrix from this program. If fasta file and endlimit specified, a second matrix file with the .ends suffix will be created also.\n" \
  " -x, --x_scale	     Scaling factor for the first dataset (float multiplier).\n"\
  " -y, --y_scale	     Scaling factor for the second dataset (float multiplier).\n"\
  " -t, --threads        The number of threads to use.  Default: 1\n" \
  " -k, --kmer           The kmer size to use (must be the same as kmer sized used for jellyfish hash).  Default: 31.\n" \
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
            {"threads",         required_argument,  0, 't'},
            {"kmer",            required_argument,  0, 'k'},
            {"help",            no_argument,        0, 'h'},
            {"usage",           no_argument,        0, 'u'},
            {0, 0, 0, 0}
        };

        static const char *short_options = "f:e:o:x:y:t:k:vuh";

        if (argc <= 0)
        {
            std::cout << usage() << "\n\n" << help() << std::endl;
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
                std::cerr << "Missing required argument for "
                          << (index == -1 ? std::string(1, (char)optopt) : std::string(long_options[index].name))
                          << std::endl;
                exit(1);
            case 'h':
                std::cout << usage() << "\n\n" << help() << std::endl;
                exit(0);
            case 'u':
                std::cout << usage() << "\nUse --help for more information." << std::endl;
                exit(0);
            case '?':
                std::cerr << "Use --usage or --help for some help\n";
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
            case 'k':
                kmer_arg = atoi(optarg);
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
            std::cerr << "Verbose flag set\n";

        if (fasta_arg)
            std::cerr << "Fasta file: " << fasta_arg << "\n";

        if (threads_arg)
            std::cerr << "Threads requested: " << threads_arg << "\n";

        if (kmer_arg)
            std::cerr << "Kmer size: " << kmer_arg << "\n";

        if (endlimit_arg)
            std::cerr << "Endlimit: " << endlimit_arg << "\n";

        if (xscale_arg)
            std::cerr << "X Scale Arg: " << xscale_arg << "\n";

        if (yscale_arg)
            std::cerr << "Y Scale Arg: " << yscale_arg << "\n";

        if (db1_arg)
            std::cerr << "Jellyfish hash 1: " << db1_arg << "\n";

        if (db2_arg)
            std::cerr << "Jellyfish hash 2: " << db2_arg << "\n";

        if (output_arg)
            std::cerr << "Matrix Output file provided: " << output_arg << "\n";
    }

private:
};

#endif // __COMP_ARGS_HPP__

