#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <string>
#include <iostream>

#include <gnuplot/gnuplot_i.hpp>

#include "flame_plot_args.hpp"
#include "flame_plot_main.hpp"


void configureFlamePlot(Gnuplot* plot, string* type, const char* output_path,
	uint canvas_width, uint canvas_height)
{
    
    std::ostringstream term_str;

    if (type->compare("png") == 0)
    {
        term_str << "set terminal png";
    }
    else if (type->compare("ps") == 0)
    {
        term_str << "set terminal postscript color";
    }
    else if (type->compare("pdf") == 0)
    {
        term_str << "set terminal pdf color";
    }
    else
    {
        std::cerr << "Unknown file type, assuming PNG\n";
        term_str << "set terminal png";
    }

    term_str << " large";
    term_str << " size " << canvas_width << "," << canvas_height;

    plot->cmd(term_str.str());

    std::ostringstream output_str;
    output_str << "set output \"" << output_path << "\"";
    plot->cmd(output_str.str());
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

    configureFlamePlot(flame, args.output_type, args.output_arg->c_str(), args.width, args.height);

    flame->set_title(args.title);
    flame->set_xlabel(args.x_label);
    flame->set_ylabel(args.y_label);

    flame->set_xrange(-1, 1000);
    flame->set_yrange(-1, 1000);

    //flame->set_xlogscale();
    //flame->set_ylogscale();
    //flame->set_zlogscale();

    flame->cmd("set palette rgb 21,22,23");
    flame->cmd("set size ratio 1");

    std::ostringstream rangestr;
    rangestr << "set cbrange [0:" << args.z_cap << "]";
    flame->cmd(rangestr.str());

    std::ostringstream plotstr;
    plotstr << "plot '" << args.mx_arg->c_str() << "' matrix with image";

    flame->cmd(plotstr.str());

    delete flame;

    return 0;
}
