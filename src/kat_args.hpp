#ifndef __SECT_ARGS_HPP__
#define __SECT_ARGS_HPP__

#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

using std::string;
using std::cerr;
using std::cout;

class KatArgs
{
private:
    string  mode_arg;
    int     mode_argc;
    char**   mode_argv;

public:


    // Default constructor
    KatArgs()
    {}

    // Constructor that parses command line options
    KatArgs(int argc, char* argv[])
    {
        parse(argc, argv);
    }

    string getMode() {
        return mode_arg;
    }

    int getModeArgC() {
        return mode_argc;
    }

    char** getModeArgV() {
        return mode_argv;
    }




#define kat_args_USAGE "\nUsage: kat <mode>"
    const char * usage() const
    {
        return kat_args_USAGE;
    }

    void error(const char *msg)
    {
        cerr << "Error: " << msg << "\n" << usage()
                  << "\nUse --help for more information"
                  << std::endl;
        exit(1);
    }


#define kat_args_HELP "The Kmer Analysis Toolkist (KAT) contains a number of tools that analyse jellyfish kmer hashes\n\n" \
  "First argument should be the tool/mode you wish to use:\n\n" \
  "   - sect: SEquence Coverage estimator Tool.  Estimates the coverage of each sequence in a fasta file using kmers from a jellyfish hash\n" \
  "   - comp: Kmer comparison tool.  Creates a matrix of shared kmers between two jellyfish hashes.\n" \
  "   - plot: PLOT assisting tool.  Creates useful plots to visualise kmer distributions\n\n" \
  "Options:\n" \
  "     --usage                              Usage\n" \
  "     --help                               This message\n" \
  " -V, --version                            Version"

    const char * help() const
    {
        return kat_args_HELP;
    }

#define kat_args_HIDDEN "Hidden options:"
    const char * hidden() const
    {
        return kat_args_HIDDEN;
    }

    void print_version(std::ostream &os = std::cout) const
    {
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1.0"
#endif
        os << PACKAGE_VERSION << "\n";
    }

    bool validMode(char* mode_str)
    {
        return (strcmp(mode_str, "sect") == 0 ||
                strcmp(mode_str, "comp") == 0 ||
                strcmp(mode_str, "plot")) ?
                    true : false;
    }

    void parse(int argc, char *argv[])
    {
        int c;

        static struct option long_options[] =
        {
            {"help",  no_argument,         0, 'h'},
            {"usage", no_argument,         0, 'u'},
            {"version", no_argument,       0, 'V'},
            {0, 0, 0, 0}
        };

        static const char *short_options = "Vuh";


        if (argc <= 1) {
            cout << usage() << "\n\n" << help() << std::endl;
            exit(0);
        }
        else if (validMode(argv[1])) {

            mode_arg = argv[1];
            mode_argc = argc - 1;
            mode_argv = argv + 1;
        }
        else {

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
                case 'V':
                    print_version();
                    exit(0);
                case '?':
                    cerr << "Use --usage or --help for some help\n";
                    exit(1);

                }
            }

            error("Invalid command line arguments passed to KAT\n\n");
            exit(1);
        }
    }


private:
};

#endif // __SECT_ARGS_HPP__

