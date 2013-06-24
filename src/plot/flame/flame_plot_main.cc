#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <string>
#include <iostream>

#include <gnuplot/gnuplot_i.hpp>

#include "flame_plot_args.hpp"
#include "flame_plot_main.hpp"


void setFlamePlotOutputType(Gnuplot* plot, string* type, const char* output_path)
{
    if (type->compare("png") == 0)
    {
        plot->savetopng(output_path);
    }
    else if (type->compare("ps") == 0)
    {
        plot->savetops(output_path);
    }
    else if (type->compare("pdf") == 0)
    {
        plot->savetopdf(output_path);
    }
    else
    {
        std::cerr << "Unknown file type, assuming PNG\n";
        plot->savetopng(output_path);
    }
}


// Start point
int flamePlotStart(int argc, char *argv[])
{
    // Parse args
    FlamePlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    Gnuplot* flame = new Gnuplot("lines");

    setFlamePlotOutputType(flame, args.output_type, args.getOutputWithoutExt().c_str());

    flame->set_title(args.title);
    flame->set_xlabel(args.xlabel);
    flame->set_ylabel(args.ylabel);

    flame->set_xrange(0, 1000);
    flame->set_yrange(0, 1000);

    std::ostringstream plotstr;
    plotstr << "plot '" << args.mx_arg->c_str() << "' matrix with image";

    flame->cmd(plotstr.str());

    delete flame;

    return 0;
}
