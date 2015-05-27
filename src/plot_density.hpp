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
#include <fstream>
#include <iostream>

#include <gnuplot/gnuplot_i.hpp>

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include <matrix/sparse_matrix.hpp>
#include <matrix/matrix_metadata_extractor.hpp>

using std::string;
using std::ifstream;
using std::istringstream;
using std::ostringstream;


typedef boost::error_info<struct PlotDensityError,string> PlotDensityErrorInfo;
struct PlotDensityException: virtual boost::exception, virtual std::exception { };

namespace kat {
    
    class PlotDensity {
    private:
        
        static string helpMessage() {
            return string("Usage: kat plot density [options] <matrix_file>\n\n") +
                    "Create K-mer Density Plots.\n\n" +
                    "Creates a scatter plot, where the density or \"heat\" at each point represents the number of distinct K-mers " \
                    "at that point.  Typically this is used to visualise a matrix produced by the \"kat comp\" tool to compare " \
                    "multiplicities from two K-mer hashes produced by different NGS reads, or to visualise the GC vs K-mer " \
                    "multiplicity matricies produced by the \"kat gcp\" tool.\n\n" \
                    "Options";
        }
        

    public:
       
        static int main(int argc, char *argv[]) {
            
            const string DEFAULT_TITLE      = "Density plot";
            const string DEFAULT_X_LABEL    = "X";
            const string DEFAULT_Y_LABEL    = "Y";
            const string DEFAULT_Z_LABEL    = "Z";
            const uint32_t DEFAULT_X_MAX    = 1000;
            const uint32_t DEFAULT_Y_MAX    = 1000;
            const uint64_t DEFAULT_Z_MAX    = 10000;
            
            path        mx_file;           
            string      output_type;
            path        output;
            string      title;
            string      x_label;
            string      y_label;
            string      z_label;
            uint16_t    width;
            uint16_t    height;
            uint32_t    x_max;
            uint32_t    y_max;
            uint64_t    z_max;
            bool        verbose;
            bool        help;
        
            // Declare the supported options.
            po::options_description generic_options(PlotDensity::helpMessage(), 100);
            generic_options.add_options()
                    ("output_type,p", po::value<string>(&output_type)->default_value("png"), 
                        "The plot file type to create: png, ps, pdf.  Warning... if pdf is selected please ensure your gnuplot installation can export pdf files.")
                    ("output,o", po::value<path>(&output),
                        "The path to the output file")
                    ("title,t", po::value<string>(&title)->default_value(DEFAULT_TITLE),
                        "Title for plot")
                    ("x_label,a", po::value<string>(&x_label)->default_value(DEFAULT_X_LABEL),
                        "Label for the x-axis (value taken from matrix metadata if present)")
                    ("y_label,b", po::value<string>(&y_label)->default_value(DEFAULT_Y_LABEL),
                        "Label for the y-axis (value taken from matrix metadata if present)")
                    ("z_label,c", po::value<string>(&z_label)->default_value(DEFAULT_Z_LABEL),
                        "Label for the z-axis (value taken from matrix metadata if present)")
                    ("x_max,x", po::value<uint32_t>(&x_max)->default_value(DEFAULT_X_MAX),
                        "Maximum value for the x-axis (value taken from matrix metadata if present)")
                    ("y_max,y", po::value<uint32_t>(&y_max)->default_value(DEFAULT_Y_MAX),
                        "Maximum value for the y-axis (value taken from matrix metadata if present)")
                    ("z_max,z", po::value<uint64_t>(&z_max)->default_value(DEFAULT_Z_MAX),
                        "Maximum value for the z-axis (value taken from matrix metadata if present)")
                    ("width,w", po::value<uint16_t>(&width)->default_value(1024),
                        "Width of canvas")
                    ("height,h", po::value<uint16_t>(&height)->default_value(1024),
                        "Height of canvas")
                    ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                        "Print extra information.")
                    ("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
                    ;

            // Hidden options, will be allowed both on command line and
            // in config file, but will not be shown to the user.
            po::options_description hidden_options("Hidden options");
            hidden_options.add_options()
                    ("mx_file,s", po::value<path>(&mx_file), "Path to the matrix file to plot.")                    
                    ;

            // Positional option for the input bam file
            po::positional_options_description p;
            p.add("mx_file", 1);
            

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
        
            // Check input file exists
            if (!bfs::exists(mx_file) && !bfs::symbolic_link_exists(mx_file))
            {
                cerr << endl << "Could not find matrix file at: " << mx_file << "; please check the path and try again." << endl << endl;
                return 1;
            }

            // Get plotting properties, either from file, or user.  User args have precedence.
            uint16_t x_range = x_max != DEFAULT_X_MAX ? x_max : mme::getNumeric(mx_file, mme::KEY_NB_COLUMNS);
            uint16_t y_range = y_max != DEFAULT_Y_MAX ? y_max : mme::getNumeric(mx_file, mme::KEY_NB_ROWS);
            uint32_t z_range = z_max != DEFAULT_Z_MAX ? z_max : mme::getNumeric(mx_file, mme::KEY_MAX_VAL) / 1000;        // Scale down from the max value to saturate the hot spots

            string xl = !boost::equals(x_label, DEFAULT_X_LABEL) ? x_label : mme::getString(mx_file, mme::KEY_X_LABEL);
            string yl = !boost::equals(y_label, DEFAULT_Y_LABEL) ? y_label : mme::getString(mx_file, mme::KEY_Y_LABEL);
            string zl = !boost::equals(z_label, DEFAULT_Z_LABEL) ? z_label : mme::getString(mx_file, mme::KEY_Z_LABEL);

            string t = !boost::equals(title, DEFAULT_TITLE) ? title : mme::getString(mx_file, mme::KEY_TITLE);

            bool transpose = mme::getNumeric(mx_file, mme::KEY_TRANSPOSE) == 0 ? false : true;


            // If neither the user or the data file contain any ideas of what values to use then use defaults
            x_range = x_range < 0 ? DEFAULT_X_MAX : x_range;
            y_range = y_range < 0 ? DEFAULT_Y_MAX : y_range;
            z_range = z_range < 0 ? DEFAULT_Z_MAX : z_range;    // Saturate the hot spots a bit to show more detail around the edges

            xl = xl.empty() ? DEFAULT_X_LABEL : xl;
            yl = yl.empty() ? DEFAULT_Y_LABEL : yl;
            zl = zl.empty() ? DEFAULT_Z_LABEL : zl;

            t = t.empty() ? DEFAULT_TITLE : t;

            if (verbose)
            {
                cerr << "Actual variables used to create plot:" << endl;
                cerr << "Output Path: " << output << endl;
                cerr << "X Range: " << x_range << endl;
                cerr << "Y Range: " << y_range << endl;
                cerr << "Z Range: " << z_range << endl;
                cerr << "X Label: " << x_label << endl;
                cerr << "Y Label: " << y_label << endl;
                cerr << "Z Label: " << z_label << endl;
                cerr << "Title: " << title << endl;
            }


            // Start defining the plot
            Gnuplot density("lines");

            density.configurePlot(output_type, output.string(), width, height);

            density.set_title(title);
            density.set_xlabel(x_label);
            density.set_ylabel(y_label);

            std::ostringstream cblabelstr;
            cblabelstr << "set cblabel \"" << z_label << "\"";
            density.cmd(cblabelstr.str());

            density.set_xrange(0, x_range);
            density.set_yrange(0, y_range);

            //flame.set_xlogscale();
            //flame.set_ylogscale();
            //flame.set_zlogscale();

            density.cmd("set palette rgb 21,22,23");
            density.cmd("set size ratio 1");

            std::ostringstream rangestr;
            rangestr << "set cbrange [0:" << z_range << "]";
            density.cmd(rangestr.str());

            // Transpose the matrix and store in ostream
            ostringstream data;
            SparseMatrix<uint64_t> mx(mx_file);
            mx.printMatrix(data, transpose);

            // Plot the transposed matrix as image
            std::ostringstream plotstr;
            plotstr << "plot '-' matrix with image" << endl
                    << data.str()
                    << "e" << endl
                    << "e" << endl;

            density.cmd(plotstr.str());

            return 0;
        }
    };
}
