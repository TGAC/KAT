#ifndef __SECT_ARGS_HPP__
#define __SECT_ARGS_HPP__

#include <getopt.h>
#include <stdlib.h>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

class SectArgs
{
public:
    const char * 	fasta_arg;
    const char * 	db_arg;
    const char *  input_type;
    const char *  output_arg;
    uint_t threads_arg;
    bool verbose;

    // Default constructor
    SectArgs() :
        output_arg(NULL), threads_arg(1), verbose(false)
    {}

    // Constructor that parses command line options
    SectArgs(int argc, char* argv[]) :
        output_arg(NULL), threads_arg(1), verbose(false)
    {
        parse(argc, argv);
    }




#define seqcvg_args_USAGE "Usage: kat sect [options] -f <fasta_file_path> db_path\n"
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


#define sect_args_HELP "Estimates coverage for sequences in a fasta file using jellyfish kmer counts. Kmers containing any Ns derived from sequences in the fasta files will have 0 coverage.\n\n" \
  "Options (default value in (), *required):\n" \
  " -f, --fasta          *Fasta file contains sequences that should have coverage estimated.  \n" \
  " -o, --count_output   File that should contain the sequence count profiles produced from this program. If not specified count data is not output\n" \
  " -t, --threads        The number of threads to use (1).\n" \
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
            {"verbose", no_argument,       0, 'v'},
            {"fasta",   required_argument, 0, 'f'},
            {"count_output",  required_argument, 0, 'o'},
            {"threads", required_argument, 0, 't'},
            {"help",  no_argument,         0, 'h'},
            {"usage", no_argument,         0, 'u'},
            {0, 0, 0, 0}
        };

        static const char *short_options = "f:o:t:vuh";

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
            case 'f':
                fasta_arg = optarg;
                break;
            case 'o':
                output_arg = optarg;
                break;
            case 't':
                threads_arg = atoi(optarg);
                break;

            }
        }

        // Parse arguments
        if(argc - optind != 1)
            error("Requires exactly 1 argument.");
        db_arg = argv[optind++];
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

        if (db_arg)
            cerr << "Jellyfish hash: " << db_arg << endl;

        if (outputGiven())
        {
            cerr << "Count Output file provided: " << output_arg << endl;
        }
        else
        {
            cerr << "No output argument provided.  Not outputing count information" << endl;
        }

        cerr << endl;
    }

private:
};

#endif // __SECT_ARGS_HPP__

