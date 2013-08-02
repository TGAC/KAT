#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <string>
#include <iostream>

#include <gnuplot/gnuplot_i.hpp>

#include "contamination_plot_args.hpp"
#include "contamination_plot_main.hpp"


// Start point
int contaminationPlotStart(int argc, char *argv[])
{
    // Parse args
    ContaminationPlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    Gnuplot* contamination = new Gnuplot("lines");

    // Work out the output path to use (either user specified or auto generated)
    string output_path = args.determineOutputPath();

    contamination->configurePlot(*(args.output_type), output_path, args.width, args.height);

    contamination->set_title(args.title);
    contamination->set_xlabel(args.x_label);
    contamination->set_ylabel(args.y_label);

    contamination->set_xrange(-1, 1000);
    contamination->set_yrange(-1, 1000);

    //flame->set_xlogscale();
    //flame->set_ylogscale();
    //flame->set_zlogscale();

    contamination->cmd("set palette rgb 21,22,23");
    contamination->cmd("set size ratio 1");

    std::ostringstream rangestr;
    rangestr << "set cbrange [0:" << args.z_cap << "]";
    contamination->cmd(rangestr.str());

    std::ostringstream plotstr;
    plotstr << "plot '" << args.mx_arg->c_str() << "' matrix with image";

    contamination->cmd(plotstr.str());

    delete contamination;

    return 0;
}
