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
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>

#include <gnuplot/gnuplot_i.hpp>
#include <str_utils.hpp>
#include <matrix/sparse_matrix.hpp>
#include <matrix/matrix_metadata_extractor.hpp>

#include "spectra_mx_plot_args.hpp"
#include "spectra_mx_plot_main.hpp"

using std::string;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::ifstream;
using std::vector;

using kat::SpectraMxPlotArgs;

string getDataFromList(string mx_path, string list)
{
    ostringstream data_str;

    // Split string by commas
    vector<string> parts = kat::splitString(list, ',');

    // Need to change colors and labels
    for(uint16_t i = 0; i < parts.size(); i++)
    {
        data_str << "'-' using 1:2 with linespoints ps 0.25 linetype 1 linecolor " << i+1 << " title '" << parts[i] << "'";
        if (i != parts.size() - 1)
        {
            data_str << ", ";
        }
    }

    data_str << "\n";

    SparseMatrix<uint64_t> mx = SparseMatrix<uint64_t>(mx_path);

    for(uint16_t i = 0; i < parts.size(); i++)
    {
        char c_or_r = parts[i][0];
        uint16_t index = atoi(parts[i].substr(1).c_str());

        if (c_or_r == 'c')
        {
            for(uint16_t j = 0; j < mx.width(); j++)
            {
                data_str << j << " " << mx.get(index,j) << endl;
            }
        }
        else if (c_or_r == 'r')
        {
            for(uint16_t j = 0; j < mx.height(); j++)
            {
                data_str << j << " " << mx.get(j,index) << endl;
            }
        }
        else
        {
            throw;
        }

        data_str << "e\n";
    }

    return data_str.str();
}


string getVennData(string mx_path)
{
    // Ready buffer for data
    ostringstream data_str;

    // Load matrix
    SparseMatrix<uint64_t> mx = SparseMatrix<uint64_t>(mx_path);

    cerr << "Matrix loaded:- Width: " << mx.width() << "; Height: " << mx.height() << ";" << endl;

    // Define datasets
    data_str << "'-' using 1:2 with linespoints ps 0.25 linetype 1 linecolor " << "1" << " title '" << "dataset 1 exclusive content" << "',";
    data_str << "'-' using 1:2 with linespoints ps 0.25 linetype 1 linecolor " << "2" << " title '" << "dataset 1 shared content" << "',";
    data_str << "'-' using 1:2 with linespoints ps 0.25 linetype 1 linecolor " << "3" << " title '" << "dataset 2 exclusive content" << "',";
    data_str << "'-' using 1:2 with linespoints ps 0.25 linetype 1 linecolor " << "4" << " title '" << "dataset 2 shared content" << "'" << endl;

    cerr << "Datasets defined" << endl;


    // Acquire data

    // Dataset 1 Exclusive content
    for(uint16_t j = 0; j < mx.width(); j++)
    {
        data_str << j << " " << mx.get(0,j) << endl;
    }
    data_str << "e\n";
    cerr << "Dataset 1 Exclusive content collected" << endl;


    // Dataset 1 Shared content
    for(uint16_t j = 1; j < mx.width(); j++)
    {
        uint64_t sum = mx.sumColumn(j, 1, mx.height() - 1);
        data_str << j << " " << sum << endl;
    }
    data_str << "e\n";
    cerr << "Dataset 1 Shared content calculated" << endl;

    // Dataset 2 Exclusive content
    for(uint16_t j = 0; j < mx.height(); j++)
    {
        data_str << j << " " << mx.get(j,0) << endl;
    }
    data_str << "e\n";
    cerr << "Dataset 2 Exclusive content collected" << endl;


    // Dataset 2 Shared content
    for(uint16_t j = 1; j < mx.height(); j++)
    {
        uint64_t sum = mx.sumRow(j, 1, mx.width() - 1);
        data_str << j << " " << sum << endl;
    }
    data_str << "e\n";
    cerr << "Dataset 2 Shared content calculated" << endl;

    return data_str.str();
}


// Start point
int kat::spectraMxPlotStart(int argc, char *argv[])
{
    // Parse args
    SpectraMxPlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
    {
        cerr << "Starting Spectra Matrix Plotter" << endl;
        cerr << "Arguments passed from command line are: " << endl;
        args.print();
    }

    // Modify variables as appropriate
    string title = args.titleModified() ? args.title : mme::getString(args.mx_path, mme::KEY_TITLE);
    title = title.empty() ? args.defaultTitle() : title;


    // Initialise gnuplot
    Gnuplot spectra_mx_plot = Gnuplot("lines");

    // Work out the output path to use (either user specified or auto generated)
    string output_path = args.determineOutputPath();

    spectra_mx_plot.configurePlot(args.output_type, output_path, args.width, args.height);

    spectra_mx_plot.set_title(title);
    spectra_mx_plot.set_xlabel(args.x_label);
    spectra_mx_plot.set_ylabel(args.y_label);
    spectra_mx_plot.set_xrange(args.x_logscale ? 1 : args.x_min, args.x_logscale ? DEFAULT_X_MAX : args.x_max);
    spectra_mx_plot.set_yrange(args.y_logscale ? 1 : args.y_min, args.y_max);

    if (args.x_logscale)
        spectra_mx_plot.set_xlogscale();
    if (args.y_logscale)
        spectra_mx_plot.set_ylogscale();

    spectra_mx_plot.cmd("set size ratio 1");
    spectra_mx_plot.cmd("set key font \",8\"");
    spectra_mx_plot.cmd("set tics font \", 8\"");
    spectra_mx_plot.cmd("set palette rgb 33,13,10");
    spectra_mx_plot.cmd("unset colorbox");

    spectra_mx_plot.cmd("set style data linespoints");

    ostringstream plot_str;

    string data;

    if (!args.list.empty())
    {
        if (args.verbose)
            cerr << "Extracting requested data from matrix..." << endl;

        data = getDataFromList(args.mx_path, args.list);

        if (args.verbose)
            cerr << "done." << endl;
    }
    else if (args.venn)
    {
        if (args.verbose)
            cerr << "Extracting Venn data from matrix... ";

        data = getVennData(args.mx_path);

        if (args.verbose)
            cerr << "done." << endl;
    }
    else
    {
        cerr << "Error: not data to process.  You did not select a list of content from the matrix (\"--list\"), or alternatively select venn mode (\"--venn\")";
        return 1;
    }


    plot_str << "plot " << data ;

    spectra_mx_plot.cmd(plot_str.str());

    if (args.verbose)
        //cerr << "Data Plotted." << endl;
        cerr << "Plotted data: " << plot_str.str() << endl;

    return 0;
}
