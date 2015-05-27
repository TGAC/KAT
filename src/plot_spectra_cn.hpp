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
#include <iomanip>
#include <memory>
#include <vector>
using std::string;
using std::ostringstream;
using std::shared_ptr;
using std::make_shared;
using std::vector;

#include <gnuplot/gnuplot_i.hpp>

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;


namespace kat {
    

    class PlotSpectraCn {
        
    protected:
    
        static string createLineStyleStr(uint16_t i, const char* colour)
        {
            ostringstream line_style_str;
            line_style_str << "set style line " << i << " lc rgb \"" << colour << "\"";
            return line_style_str.str();
        }

        static string createSinglePlotString(const path& data_file, uint16_t idx, uint16_t level_count, bool cumulative)
        {
            int16_t mxCol = idx+1;

            ostringstream col_str;
            col_str << mxCol;

            int limit = level_count+1;

            ostringstream sum_str;
            sum_str << "(sum[i=" << (limit+1) << ":900] column(i))";

            string col = "";
            if (cumulative) {
                ostringstream cumul_str;
                cumul_str << " (cum_sum($" << mxCol << "))";
                col = cumul_str.str();
            }
            else {
                col = (idx == limit) ? sum_str.str() : col_str.str();
            }

            string plus = (idx == limit) ? " +" : "";

            ostringstream plot_str;

            plot_str << "'" << data_file.string() << "' u " << col << " t \"" << idx << "x" << plus << "\"";
            return plot_str.str();
        }


        static shared_ptr<vector<uint16_t>> getStandardCols(bool ignoreAbsent, uint16_t maxDuplication)
        {
            shared_ptr<vector<uint16_t>> col_defs = make_shared<vector<uint16_t>>();

            // To cover absent content if required
            if (!ignoreAbsent)
            {
                col_defs->push_back(0);
            }

            // To cover 1 - max_duplication
            for(uint16_t i = 1; i <= maxDuplication; i++)
            {
                col_defs->push_back(i);
            }

            // Placeholder to cover everything else
            col_defs->push_back(maxDuplication + 1);

            return col_defs;
        }


        static shared_ptr<vector<uint16_t>> getUserDefinedCols(const string& columns)
        {
            shared_ptr<vector<uint16_t>> col_defs = make_shared<vector<uint16_t>>();

            string delimiter = ",";
            string s(columns);

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

        
        static string helpMessage() {
            return string("Usage: kat plot spectra-cn [options] <matrix_file>\n\n") +
                    "Creates a stacked histogram showing the level of duplication in an assembly.\n\n" \
                    "Shows K-mer duplication levels, which correspond to copy number variation within an assembly by comparing " \
                    "K-mers found in sequenced reads, to K-mers found in an assembly of those reads. Uses matrix output from the " \
                    "\"kat comp\" tool.\n\n" \
                    "Options";
        }

    public:

        static int main(int argc, char *argv[])
        {
            const string     DEFAULT_TITLE           = "Spectra Copy Number Plot";
            const string     DEFAULT_X_LABEL         = "X";
            const string     DEFAULT_Y_LABEL         = "Y";

            path        mx_file;           
            string      output_type;
            path        output;
            string      title;
            string      x_label;
            string      y_label;
            uint16_t    width;
            uint16_t    height;
            uint32_t    x_max;
            uint32_t    y_max;
            bool        ignore_absent;
            uint16_t    max_duplication;
            string      columns;
            bool        cumulative;
            bool        verbose;
            bool        help;
        
            // Declare the supported options.
            po::options_description generic_options(PlotSpectraCn::helpMessage(), 100);
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
                    ("x_max,x", po::value<uint32_t>(&x_max)->default_value(1000),
                        "Maximum value for the x-axis (value taken from matrix metadata if present)")
                    ("y_max,y", po::value<uint32_t>(&y_max)->default_value(1000000),
                        "Maximum value for the y-axis (value taken from matrix metadata if present)")
                    ("width,w", po::value<uint16_t>(&width)->default_value(1024),
                        "Width of canvas")
                    ("height,h", po::value<uint16_t>(&height)->default_value(1024),
                        "Height of canvas")
                    ("ignore_absent,a", po::bool_switch(&ignore_absent)->default_value(false),
                        "Ignore K-mers in reads but absent from the assembly")
                    ("max_dup,m", po::value<uint16_t>(&max_duplication)->default_value(5),
                        "Maximum duplication level to show in plots")
                    ("columns,c", po::value<string>(&columns),
                        "Comma separated string listing columns to show in plot.  If used, this overrides \"--ignore_absent\" and \"--columns\"")
                    ("cumulative,u", po::bool_switch(&cumulative)->default_value(false),
                        "Plot cumulative distribution of kmers")
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
                cerr << endl << "Could not find matrix file at: " << mx_file.string() << "; please check the path and try again." << endl << endl;
                return 1;
            }

            shared_ptr<vector<uint16_t>> plot_cols = columns.empty() ? getStandardCols(ignore_absent, max_duplication) : getUserDefinedCols(columns);

            if (!plot_cols->empty())
            {
                // Determine configuration
                bool request_absent = (*plot_cols)[0] == 0 ? true : false;
                uint16_t level_count = request_absent ? plot_cols->size() - 2 : plot_cols->size() - 1;

                if (verbose)
                {
                    cerr << "Request plot for absent K-mers" << endl
                         << level_count << " levels of present K-mers requested for plotting" << endl << endl;
                }

                // Initialise gnuplot
                Gnuplot spectra_cn_plot("lines");

                spectra_cn_plot.configurePlot(output_type, output.string(), width, height);

                spectra_cn_plot.set_title(title);
                spectra_cn_plot.set_xlabel(x_label);
                spectra_cn_plot.set_ylabel(cumulative ? "Cumulative " + y_label : y_label);


                // Get plot strings
                ostringstream plot_str;


                if (cumulative) {
                    plot_str << "a = 0" << endl 
                             << "cum_sum(x)=(a=a+x,a)" << endl;
                }

                plot_str << "plot ";

                bool first = true;

                if (request_absent && !cumulative)
                {
                    plot_str << createSinglePlotString(mx_file, 0, level_count, false) << " lt rgb \"black\"";
                    first = false;
                }

                for(uint16_t i = 1; i <= level_count; i++)
                {
                    if (first)
                        first = false;
                    else
                        plot_str << ", ";

                    double col_frac = 1.0 - ((double)(i-1) / (double)(level_count-1));

                    plot_str << "a=0,";
                    plot_str << createSinglePlotString(mx_file, i, level_count, cumulative) << " lt palette frac " << std::fixed << col_frac;
                }

                // Do the rest
                plot_str << ", " << createSinglePlotString(mx_file, level_count + 1, level_count, cumulative) << " lt rgb \"gray\"";


                spectra_cn_plot.cmd("set palette rgb 33,13,10");
                spectra_cn_plot.cmd("unset colorbox");

                spectra_cn_plot.cmd("set style fill solid 1 noborder");
                spectra_cn_plot.cmd("set style histogram rowstacked");
                spectra_cn_plot.cmd("set style data histograms");

                spectra_cn_plot.set_xrange(0, x_max);
                spectra_cn_plot.set_yrange(0, y_max);

                spectra_cn_plot.cmd(plot_str.str());

                if (verbose)
                    cerr << "Gnuplot command: " << plot_str.str() << endl;
            }

            return 0;
        }
    };
    
}