#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include "flame/flame_plot_main.hpp"
#include "asm/asm_plot_main.hpp"
#include "sect/sect_plot_main.hpp"

#include "plot_args.hpp"
#include "plot_main.hpp"

using std::string;

// Start point
int plotStart(int argc, char *argv[])
{
    // Parse args
    PlotArgs args(argc, argv);

    // Pass remaining args to relevant child tool
    if (args.getMode().compare("flame") == 0)
    {
        flamePlotStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("asm") == 0)
    {
        asmPlotStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("sect") == 0)
    {
        sectPlotStart(args.getModeArgC(), args.getModeArgV());
    }

    return 0;
}
