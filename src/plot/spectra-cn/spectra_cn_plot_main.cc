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
#include <iomanip>
#include <vector>

#include <gnuplot/gnuplot_i.hpp>
#include <file_utils.hpp>

#include "spectra_cn_plot_args.hpp"
#include "spectra_cn_plot_main.hpp"

using std::string;
using std::ostringstream;
using std::vector;

using kat::SpectraCnPlotArgs;

string createLineStyleStr(uint16_t i, const char* colour)
{
    ostringstream line_style_str;
    line_style_str << "set style line " << i << " lc rgb \"" << colour << "\"";
    return line_style_str.str();
}

string createSinglePlotString(const string data_file, uint16_t idx, uint16_t level_count)
{
    ostringstream col_str;
    col_str << idx+1;

    int limit = level_count+1;

    ostringstream sum_str;
    sum_str << "(sum[i=" << limit << ":900] column(i))";

    string col = (idx == limit) ? sum_str.str() : col_str.str();

    string plus = (idx == limit) ? " +" : "";

    ostringstream plot_str;
    plot_str << "'" << data_file << "' u " << col << " t \"" << idx << "x" << plus << "\"";
    return plot_str.str();
}


vector<uint16_t>* getStandardCols(SpectraCnPlotArgs* args)
{
    vector<uint16_t>* col_defs = new vector<uint16_t>();

    // To cover absent content if required
    if (!args->ignore_absent)
    {
        col_defs->push_back(0);
    }

    // To cover 1 - max_duplication
    for(uint16_t i = 1; i <= args->max_duplication; i++)
    {
        col_defs->push_back(i);
    }

    // Placeholder to cover everything else
    col_defs->push_back(args->max_duplication + 1);

    return col_defs;
}


vector<uint16_t>* getUserDefinedCols(SpectraCnPlotArgs* args)
{
    vector<uint16_t>* col_defs = new vector<uint16_t>();

    string delimiter = ",";
    string s(args->columns);

    size_t pos = 0;
    string token;
    while ((pos = s.find(delimiter)) != string::npos) {
        token = s.substr(0, pos);

        col_defs->push_back(atoi(s.c_str()));

        s.erase(0, pos + delimiter.length());
    }

    col_defs->push_back(atoi(s.c_str()));

    return col_defs;
}




// Start point
int kat::spectraCnPlotStart(int argc, char *argv[])
{
    // Parse args
    SpectraCnPlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    // Check input file exists
    if (!fileExists(args.mx_arg))
    {
        cerr << endl << "Could not find matrix file at: " << args.mx_arg << "; please check the path and try again." << endl << endl;
        return 1;
    }

    vector<uint16_t>* plot_cols = args.columns.empty() ? getStandardCols(&args) : getUserDefinedCols(&args);

    if (!plot_cols->empty())
    {
        // Determine configuration
        bool request_absent = (*plot_cols)[0] == 0 ? true : false;
        uint16_t level_count = request_absent ? plot_cols->size() - 2 : plot_cols->size() - 1;

        if (args.verbose)
        {
            cerr << "Request plot for absent K-mers" << endl
                 << level_count << " levels of present K-mers requested for plotting" << endl << endl;
        }

        // Initialise gnuplot
        Gnuplot spectra_cn_plot("lines");


        // Work out the output path to use (either user specified or auto generated)
        string output_path = args.determineOutputPath();


        spectra_cn_plot.configurePlot(args.output_type, output_path, args.width, args.height);

        spectra_cn_plot.set_title(args.title);
        spectra_cn_plot.set_xlabel(args.x_label);
        spectra_cn_plot.set_ylabel(args.y_label);


        // Get plot strings
        ostringstream plot_str;



        bool first = true;

        if (request_absent)
        {
            plot_str << createSinglePlotString(args.mx_arg, 0, level_count) << " lt rgb \"black\"";
            first = false;
        }

        for(uint16_t i = 1; i <= level_count; i++)
        {
            if (first)
                first = false;
            else
                plot_str << ", ";


            double col_frac = 1.0 - ((double)(i-1) / (double)(level_count-1));

            plot_str << createSinglePlotString(args.mx_arg, i, level_count) << " lt palette frac " << std::fixed << col_frac;
        }

        // Do the rest
        plot_str << ", " << createSinglePlotString(args.mx_arg, level_count + 1, level_count) << " lt rgb \"gray\"";

        spectra_cn_plot.cmd("set palette rgb 33,13,10");
        spectra_cn_plot.cmd("unset colorbox");

        spectra_cn_plot.cmd("set style fill solid 1 noborder");
        spectra_cn_plot.cmd("set style histogram rowstacked");
        spectra_cn_plot.cmd("set style data histograms");

        spectra_cn_plot.set_xrange(0, args.x_max);
        spectra_cn_plot.set_yrange(0, args.y_max);

        ostringstream plot_cmd;
        plot_cmd << "plot " << plot_str.str();
        spectra_cn_plot.cmd(plot_cmd.str());

        if (args.verbose)
            cerr << "Gnuplot command: " << plot_cmd.str() << endl;
    }


    delete plot_cols;

    return 0;
}
