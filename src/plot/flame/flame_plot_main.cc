#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <string>
#include <iostream>

#include <gnuplot/gnuplot_i.hpp>

#include "flame_plot_args.hpp"
#include "flame_plot_main.hpp"


void configurePlot(Gnuplot* plot, string* type, const char* output_path)
{
    if (type->compare("png") == 0)
    {
        plot->cmd("set terminal png font \"Helvetica\" size 1024,1024");
    }
    else if (type->compare("ps") == 0)
    {
        plot->cmd("set terminal postscript color font \"Helvetica\" size 1024,1024");
    }
    else if (type->compare("pdf") == 0)
    {
        plot->cmd("set terminal pdf color font \"Helvetica\" size 1024,1024");
    }
    else
    {
        std::cerr << "Unknown file type, assuming PNG\n";
        plot->cmd("set terminal png font \"Helvetica\" size 1024,1024");
    }

    std::ostringstream cmdstr;
    cmdstr << "set output \"" << output_path << "\"";
    plot->cmd(cmdstr.str());
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

    configurePlot(flame, args.output_type, args.output_arg);

    flame->set_title(args.title);
    flame->set_xlabel(args.xlabel);
    flame->set_ylabel(args.ylabel);

    flame->set_xrange(-1, 1000);
    flame->set_yrange(-1, 1000);

    //flame->set_xlogscale();
    //flame->set_ylogscale();
    //flame->set_zlogscale();

    flame->cmd("set palette rgb 21,22,23");
    flame->cmd("set size ratio 1");

    std::ostringstream rangestr;
    rangestr << "set cbrange [0:" << args.zrange << "]";
    flame->cmd(rangestr.str());

    std::ostringstream plotstr;
    plotstr << "plot '" << args.mx_arg->c_str() << "' matrix with image";

    flame->cmd(plotstr.str());

    delete flame;

    return 0;
}
