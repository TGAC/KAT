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

#pragma once

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
using std::string;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::ifstream;
using std::vector;

#include <gnuplot/gnuplot_i.hpp>
#include <str_utils.hpp>

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

namespace kat {
    
    class PlotSpectraHist {
    
    protected:
        
        static string helpMessage() {
            return string("Usage: kat plot spectra-hist [options] <histo_file> [<histo_file> ...]*\n\n") +
                    "Creates K-mer Spectra Plot from one or more histograms.\n\n" +
                    "Produces K-mer spectras from \"kat hist\" or \"jellyfish histo\" output.  This tool is designed to plot line " \
                    "graphs of one or more histograms.  The idea is to be able to compare total K-mer counts between different " \
                    "datasets.";
        }
        
    public:
  
        // Start point
        static int main(int argc, char *argv[]) {
            
            const string DEFAULT_TITLE       = "Kmer histograms";
            const string DEFAULT_X_LABEL     = "X";
            const string DEFAULT_Y_LABEL     = "Y";
            const string DEFAULT_Z_LABEL     = "Z";
            
            vector<path> histo_paths;        
            string      output_type;
            path        output;
            string      title;
            string      x_label;
            string      y_label;
            uint16_t    width;
            uint16_t    height;
            uint32_t    x_min;
            uint32_t    y_min;
            uint32_t    x_max;
            uint32_t    y_max;
            bool        x_logscale;
            bool        y_logscale;
            bool        verbose;
            bool        help;
        
            // Declare the supported options.
            po::options_description generic_options(PlotSpectraHist::helpMessage());
            generic_options.add_options()
                    ("output_type,p", po::value<string>(&output_type)->default_value("png"), 
                        "The plot file type to create: png, ps, pdf.  Warning... if pdf is selected please ensure your gnuplot installation can export pdf files.")
                    ("output,o", po::value<path>(&output)->required(),
                        "The path to the output file")
                    ("title,t", po::value<string>(&title)->default_value(DEFAULT_TITLE),
                        "Title for plot")
                    ("x_label,a", po::value<string>(&x_label)->default_value(DEFAULT_X_LABEL),
                        "Label for the x-axis (value taken from matrix metadata if present)")
                    ("y_label,b", po::value<string>(&y_label)->default_value(DEFAULT_Y_LABEL),
                        "Label for the y-axis (value taken from matrix metadata if present)")
                    ("x_min,r", po::value<uint32_t>(&x_min)->default_value(0),
                        "Minimum value for the x-axis")
                    ("y_min,s", po::value<uint32_t>(&y_min)->default_value(0),
                        "Minimum value for the y-axis")
                    ("x_max,x", po::value<uint32_t>(&x_max)->default_value(10000),
                        "Maximum value for the x-axis (value taken from matrix metadata if present)")
                    ("y_max,y", po::value<uint32_t>(&y_max)->default_value(1000000),
                        "Maximum value for the y-axis (value taken from matrix metadata if present)")
                    ("width,w", po::value<uint16_t>(&width)->default_value(1024),
                        "Width of canvas")
                    ("height,h", po::value<uint16_t>(&height)->default_value(1024),
                        "Height of canvas")
                    ("x_logscale,l", po::bool_switch(&x_logscale)->default_value(false),
                        "X-axis is logscale.  This overrides the x_min and x_max limits.")
                    ("y_logscale,m", po::bool_switch(&y_logscale)->default_value(false),
                        "Y-axis is logscale.  This overrides the y_min and y_max limits.")
                    ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                        "Print extra information.")
                    ("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
                    ;

            // Hidden options, will be allowed both on command line and
            // in config file, but will not be shown to the user.
            po::options_description hidden_options("Hidden options");
            hidden_options.add_options()
                    ("histo_paths,i", po::value<vector<path>>(&histo_paths), "List of histogram files to plot")                    
                    ;

            // Positional option for the input bam file
            po::positional_options_description p;
            p.add("histo_paths", 100);
            

            // Combine non-positional options
            po::options_description cmdline_options;
            cmdline_options.add(generic_options).add(hidden_options);

            // Parse command line
            po::variables_map vm;
            po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
            po::notify(vm);

            // Output help information the exit if requested
            if (help || argc <= 1) {
                cout << generic_options << endl;
                return 1;
            }
            
            // Check input files exists
            for(uint16_t i = 0; i < histo_paths.size(); i++)
            {
                if (!bfs::exists(histo_paths[i]) && !bfs::symbolic_link_exists(histo_paths[i])) {
                    cerr << endl << "Could not find the histogram file at index " << i << ": " << histo_paths[i] << "; please check the path and try again." << endl << endl;
                    return 1;
                }
            }

            if (verbose)
                cerr << "Input validated." << endl << "Setting up plot...";

            // Initialise gnuplot
            Gnuplot spectra_hist_plot = Gnuplot("lines");

            spectra_hist_plot.configurePlot(output_type, output.string(), width, height);

            spectra_hist_plot.set_title(title);
            spectra_hist_plot.set_xlabel(x_label);
            spectra_hist_plot.set_ylabel(y_label);
            spectra_hist_plot.set_xrange(x_logscale ? 1 : x_min, x_max);
            spectra_hist_plot.set_yrange(y_logscale ? 1 : y_min, y_max);

            if (x_logscale)
                spectra_hist_plot.set_xlogscale();
            if (y_logscale)
                spectra_hist_plot.set_ylogscale();

            spectra_hist_plot.cmd("set size ratio 1");
            spectra_hist_plot.cmd("set key font \",8\"");
            spectra_hist_plot.cmd("set tics font \", 8\"");
            spectra_hist_plot.cmd("set palette rgb 33,13,10");
            spectra_hist_plot.cmd("unset colorbox");

            spectra_hist_plot.cmd("set style data linespoints");

            if (verbose)
                cerr << "done." << endl << "Setting up " << histo_paths.size() << " datasets...";

            ostringstream plot_str;
            ostringstream data_str;


            // Need to change colors and labels
            for(size_t i = 0; i < histo_paths.size(); i++)
            {
                data_str << "'-' using 1:2 with linespoints ps 0.25 linetype 1 linecolor " << i+1 << " title '" << histo_paths[i] << "'";
                if (i != histo_paths.size() - 1)
                {
                    data_str << ", ";
                }
            }

            data_str << "\n";

            if (verbose)
                cerr << "done." << endl << "Acquiring data...";

            for(size_t i = 0; i < histo_paths.size(); i++)
            {
                ifstream infile;
                infile.open(histo_paths[i].c_str());

                string line;
                while(infile.good() && !infile.eof())
                {
                    getline(infile, line);

                    // Only do something if the start of the line looks like a number
                    if (line[0] >= '0' && line[0] <= '9')
                    {
                        vector<uint64_t> parts = kat::splitUInt64(line, ' ');

                        uint64_t kmer_multiplicity = parts[0];
                        uint64_t distinct_kmers = parts[1];

                        data_str << kmer_multiplicity << " " << distinct_kmers << "\n";
                    }            
                }

                infile.close();

                data_str << "e\n";

                if (verbose)
                    cerr << i << " ";
            }

            plot_str << "plot " << data_str.str() ;

            if (verbose)
                cerr << "done." << endl << "Plotting...";

            spectra_hist_plot.cmd(plot_str.str());

            if (verbose)
                cerr << "done." << endl;

            //if (verbose)
            //    cerr << plot_str.str() << endl;

            return 0;
        }

    };
}
