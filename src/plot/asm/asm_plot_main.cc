#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <string>
#include <iostream>

#include <gnuplot/gnuplot_i.hpp>

#include "asm_plot_args.hpp"
#include "asm_plot_main.hpp"


void configurePlot(Gnuplot* plot, string* type, const char* output_path)
{
    if (type->compare("png") == 0)
    {
        plot->cmd("set terminal png font \"Helvetica\"");
    }
    else if (type->compare("ps") == 0)
    {
        plot->cmd("set terminal postscript color font \"Helvetica\"");
    }
    else if (type->compare("pdf") == 0)
    {
        plot->cmd("set terminal pdf color font \"Helvetica\"");
    }
    else
    {
        std::cerr << "Unknown file type, assuming PNG\n";
        plot->cmd("set terminal png font \"Helvetica\"");
    }

    std::ostringstream cmdstr;
    cmdstr << "set output \"" << output_path << "\"";
    plot->cmd(cmdstr.str());
}

string createLineStyleStr(int i, const char* colour)
{
    std::ostringstream lineStyleStr;

    lineStyleStr << "set style line " << i << " lc rgb \"#" << colour << "\"";

    return lineStyleStr.str();
}


// Start point
int asmPlotStart(int argc, char *argv[])
{
    // Parse args
    AsmPlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    Gnuplot* asmPlot = new Gnuplot("lines");

    configurePlot(asmPlot, args.output_type, args.output_arg);

    asmPlot->set_title(args.title);
    asmPlot->set_xlabel(args.xlabel);
    asmPlot->set_ylabel(args.ylabel);


    int i = 1;

    if (!args.ignoreAbsent)
        asmPlot->cmd(createLineStyleStr(i++, "FF8080"));

    asmPlot->cmd(createLineStyleStr(i++, "80FF80"));
    asmPlot->cmd(createLineStyleStr(i++, "80FFFF"));
    asmPlot->cmd(createLineStyleStr(i++, "8080FF"));
    asmPlot->cmd(createLineStyleStr(i++, "FFFF40"));


    asmPlot->cmd("set style increment user");

    asmPlot->cmd("set style fill solid 1 noborder");
    asmPlot->cmd("set style histogram rowstacked");
    asmPlot->cmd("set style data histograms");

    asmPlot->set_xrange(0, args.xmax);
    asmPlot->set_yrange(0, args.ymax);

    if (args.xlogscale)
        asmPlot->set_xlogscale();

    if (args.ylogscale)
        asmPlot->set_ylogscale();

    std::ostringstream absentstr;

    if (!args.ignoreAbsent)
        absentstr << "'" << args.mx_arg->c_str() << "' u 1 t \"Absent\", ";

    std::ostringstream plotstr;
    plotstr << "plot " << absentstr.str() << \
                    "'" << args.mx_arg->c_str() << "' u 2 t \"Distinct\", " \
                    "'" << args.mx_arg->c_str() << "' u 3 t \"Duplicate\", " \
                    "'" << args.mx_arg->c_str() << "' u 4 t \"Triplicate\", " \
                    "'" << args.mx_arg->c_str() << "' u 5 t \"Quadruplate\"";

    asmPlot->cmd(plotstr.str());

    delete asmPlot;

    return 0;
}
