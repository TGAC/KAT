#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "kat_args.hpp"
#include "sect/sect_main.hpp"
#include "comp/comp_main.hpp"
#include "gcp/gcp_main.hpp"
#include "plot/plot_main.hpp"

/**
 * Start point for KAT.  Processes the start of the command line and then
 * delegates the rest of the command line to the child tool.
 */
int main(int argc, char *argv[])
{
    // Parse args
    KatArgs args(argc, argv);


    // Pass remaining args to relevant child tool
    if (args.getMode().compare("sect") == 0)
    {
        sectStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("comp") == 0)
    {
        compStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("gcp") == 0)
    {
        gcpStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("plot") == 0)
    {
        plotStart(args.getModeArgC(), args.getModeArgV());
    }

    return 0;
}
