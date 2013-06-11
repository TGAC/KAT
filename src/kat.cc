#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "kat_args.hpp"
#include "sect/sect_main.hpp"
#include "comp/comp_main.hpp"

// Start point
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

    return 0;
}
