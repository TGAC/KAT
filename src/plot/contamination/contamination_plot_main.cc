//  ********************************************************************
//  This file is part of KAT - the Kmer Analysis Toolkit.
//
//  KAT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  KAT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with KAT.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <string>
#include <iostream>

#include <gnuplot/gnuplot_i.hpp>

#include "contamination_plot_args.hpp"
#include "contamination_plot_main.hpp"

using kat::ContaminationPlotArgs;

// Start point
int kat::contaminationPlotStart(int argc, char *argv[])
{
    // Parse args
    ContaminationPlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    Gnuplot contamination("lines");

    // Work out the output path to use (either user specified or auto generated)
    string output_path = args.determineOutputPath();

    contamination.configurePlot(*(args.output_type), output_path, args.width, args.height);

    contamination.set_title(args.title);
    contamination.set_xlabel(args.x_label);
    contamination.set_ylabel(args.y_label);

    contamination.set_xrange(-1, 1000);
    contamination.set_yrange(-1, 1000);

    //flame->set_xlogscale();
    //flame->set_ylogscale();
    //flame->set_zlogscale();

    contamination.cmd("set palette rgb 21,22,23");
    contamination.cmd("set size ratio 1");

    std::ostringstream rangestr;
    rangestr << "set cbrange [0:" << args.z_cap << "]";
    contamination.cmd(rangestr.str());

    std::ostringstream plotstr;
    plotstr << "plot '" << args.mx_arg->c_str() << "' matrix with image";

    contamination.cmd(plotstr.str());

    return 0;
}
