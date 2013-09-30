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
#include <file_utils.hpp>

#include "spectra_hist_plot_args.hpp"
#include "spectra_hist_plot_main.hpp"

using std::string;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::ifstream;
using std::vector;

using kat::SpectraHistPlotArgs;



// Start point
int kat::spectraHistPlotStart(int argc, char *argv[])
{
    // Parse args
    SpectraHistPlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    // Check input files exists
    for(uint16_t i = 0; i < args.histo_paths.size(); i++)
    {
        if (!fileExists(args.histo_paths[i]))
        {
            cerr << endl << "Could not find the histogram file at index " << i << ": " << args.histo_paths[i] << "; please check the path and try again." << endl << endl;
            return 1;
        }
    }

    // Initialise gnuplot
    Gnuplot spectra_hist_plot = Gnuplot("lines");

    // Work out the output path to use (either user specified or auto generated)
    string output_path = args.determineOutputPath();

    spectra_hist_plot.configurePlot(args.output_type, output_path, args.width, args.height);

    spectra_hist_plot.set_title(args.title);
    spectra_hist_plot.set_xlabel(args.x_label);
    spectra_hist_plot.set_ylabel(args.y_label);
    spectra_hist_plot.set_xrange(args.x_logscale ? 1 : args.x_min, args.x_logscale ? DEFAULT_X_MAX : args.x_max);
    spectra_hist_plot.set_yrange(args.y_logscale ? 1 : args.y_min, args.y_max);

    if (args.x_logscale)
        spectra_hist_plot.set_xlogscale();
    if (args.y_logscale)
        spectra_hist_plot.set_ylogscale();

    spectra_hist_plot.cmd("set size ratio 1");
    spectra_hist_plot.cmd("set key font \",8\"");
    spectra_hist_plot.cmd("set tics font \", 8\"");
    spectra_hist_plot.cmd("set palette rgb 33,13,10");
    spectra_hist_plot.cmd("unset colorbox");

    spectra_hist_plot.cmd("set style data linespoints");

    ostringstream plot_str;
    ostringstream data_str;


    // Need to change colors and labels
    for(uint8_t i = 0; i < args.histo_paths.size(); i++)
    {
        data_str << "'-' using 1:2 with linespoints ps 0.25 linetype 1 linecolor " << i+1 << " title '" << args.histo_paths[i] << "'";
        if (i != args.histo_paths.size() - 1)
        {
            data_str << ", ";
        }
    }

    data_str << "\n";

    for(uint8_t i = 0; i < args.histo_paths.size(); i++)
    {
        ifstream infile;
        infile.open(args.histo_paths[i].c_str());

        string line("");
        while(!infile.eof())
        {
            getline(infile, line);

            // Only do something if the start of the line isn't a #
            if (line[0] != '#')
            {
                vector<uint32_t> parts = kat::splitUInt32(line, ' ');

                uint32_t kmer_multiplicity = parts[0];
                uint32_t distinct_kmers = parts[1];

                data_str << kmer_multiplicity << " " << distinct_kmers << "\n";
            }
        }

        infile.close();

        data_str << "e\n";
    }

    plot_str << "plot " << data_str.str() ;

    spectra_hist_plot.cmd(plot_str.str());

    if (args.verbose)
        cerr << "Plotted data: " << plot_str.str() << endl;

    return 0;
}
