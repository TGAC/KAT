#ifndef __SECT_ARGS_HPP__
#define __SECT_ARGS_HPP__

#include <getopt.h>
#include <string.h>

class KatArgs
{
private:
    string  mode_arg;
    int     mode_argc;
    char*   mode_argv[];

public:


    // Default constructor
    SectArgs() :
        mode_arg(NULL), mode_argc(0), mode_argv(NULL)
    {}

    // Constructor that parses command line options
    SectArgs(int argc, char* argv[]) :
        mode_arg(NULL), mode_argc(0), mode_argv(NULL)
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




#define kat_args_USAGE "Usage: kat <mode>"
    const char * usage() const
    {
        return kat_args_USAGE;
    }

    void error(const char *msg)
    {
        std::cerr << "Error: " << msg << "\n" << usage()
                  << "\nUse --help for more information"
                  << std::endl;
        exit(1);
    }


#define kat_args_HELP "The Kmer Analysis Toolkist (KAT) contains a number of tools that analyse jellyfish kmer hashes\n\n" \
  "First argument should be the tool you wish to use: " \
  "Options (default value in (), *required):\n" \
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


        if (argc == 0) {
            std::cout << usage() << "\n\n" << help() << std::endl;
            exit(0);
        }
        else if (validMode(argv[0])) {
            mode_arg = argv[0];
            mode_argc = argc -1;
            mode_argv = &argv[1];
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
                case 'V':
                    print_version();
                    exit(0);
                case '?':
                    std::cerr << "Use --usage or --help for some help\n";
                    exit(1);

                }
            }
        }
    }


    void print()
    {
    }

private:
};

#endif // __SECT_ARGS_HPP__

