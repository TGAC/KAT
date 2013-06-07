#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "kat_args.hpp";
#include "sect/sect.hpp";


// Start point
int main(int argc, char *argv[])
{
    // Parse args
    KatArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
    {
        args.print();
    }


    // Pass remaining args to relevant child tool
    if (strncmp(args.getMode(), "sect"))
    {
        sectStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (strncmp(args.getMode(), "comp"))
    {
        compStart(args.getModeArgC(), args.getModeArgV());
    }


    return 0;
}
