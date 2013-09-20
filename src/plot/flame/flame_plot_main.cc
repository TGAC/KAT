//  ********************************************************************
//  This file is part of KAT - the K-mer Analysis Toolkit.
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
#include <fstream>
#include <iostream>

#include <gnuplot/gnuplot_i.hpp>

#include <matrix/matrix_metadata_extractor.hpp>

#include "flame_plot_args.hpp"
#include "flame_plot_main.hpp"

using std::string;
using std::ifstream;
using std::istringstream;

using kat::FlamePlotArgs;

// Start point
int kat::flamePlotStart(int argc, char *argv[])
{
    // Parse args
    FlamePlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();


    // Get plotting properties, either from file, or user.  User args have precedence.
    uint16_t x_range = args.xMaxModified() ? args.x_max : mme::getNumeric(args.mx_arg, mme::KEY_NB_COLUMNS);
    uint16_t y_range = args.yMaxModified() ? args.y_max : mme::getNumeric(args.mx_arg, mme::KEY_NB_ROWS);
    uint32_t z_range = args.zMaxModified() ? args.z_max : mme::getNumeric(args.mx_arg, mme::KEY_MAX_VAL) / 1000;        // Scale down from the max value to saturate the hot spots

    string x_label = args.xLabelModified() ? args.x_label : mme::getString(args.mx_arg, mme::KEY_X_LABEL);
    string y_label = args.yLabelModified() ? args.y_label : mme::getString(args.mx_arg, mme::KEY_Y_LABEL);
    string z_label = args.zLabelModified() ? args.z_label : mme::getString(args.mx_arg, mme::KEY_Z_LABEL);

    string title = args.titleModified() ? args.title : mme::getString(args.mx_arg, mme::KEY_TITLE);


    // If neither the user or the data file contain any ideas of what values to use then use defaults
    x_range = x_range < 0 ? DEFAULT_X_MAX : x_range;
    y_range = y_range < 0 ? DEFAULT_Y_MAX : y_range;
    z_range = z_range < 0 ? DEFAULT_Z_MAX : z_range;    // Saturate the hot spots a bit to show more detail around the edges

    x_label = x_label.empty() ? args.defaultXLabel() : x_label;
    y_label = y_label.empty() ? args.defaultYLabel() : y_label;
    z_label = z_label.empty() ? DEFAULT_Z_LABEL : z_label;

    title = title.empty() ? args.defaultTitle() : title;

    // Work out the output path to use (either user specified or auto generated)
    string output_path = args.determineOutputPath();


    if (args.verbose)
    {
        cerr << "Actual variables used to create plot:" << endl;
        cerr << "Output Path: " << output_path << endl;
        cerr << "X Range: " << x_range << endl;
        cerr << "Y Range: " << y_range << endl;
        cerr << "Z Range: " << z_range << endl;
        cerr << "X Label: " << x_label << endl;
        cerr << "Y Label: " << y_label << endl;
        cerr << "Z Label: " << z_label << endl;
        cerr << "Title: " << title << endl;
    }


    // Start defining the plot
    Gnuplot flame("lines");

    flame.configurePlot(args.output_type, output_path, args.width, args.height);

    flame.set_title(title);
    flame.set_xlabel(x_label);
    flame.set_ylabel(y_label);

    std::ostringstream cblabelstr;
    cblabelstr << "set cblabel \"" << z_label << "\"";
    flame.cmd(cblabelstr.str());

    flame.set_xrange(0, x_range);
    flame.set_yrange(0, y_range);

    //flame.set_xlogscale();
    //flame.set_ylogscale();
    //flame.set_zlogscale();

    flame.cmd("set palette rgb 21,22,23");
    flame.cmd("set size ratio 1");

    std::ostringstream rangestr;
    rangestr << "set cbrange [0:" << z_range << "]";
    flame.cmd(rangestr.str());

    std::ostringstream plotstr;
    plotstr << "plot '" << args.mx_arg.c_str() << "' matrix with image";

    flame.cmd(plotstr.str());

    return 0;
}
