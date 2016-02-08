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
#include <sys/ioctl.h>
using std::string;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::vector;
using std::endl;

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include <kat/gnuplot_i.hpp>
#include <kat/str_utils.hpp>

namespace kat {
    
    class PlotProfile {
    
    protected:
        static int readRecord(std::ifstream& stream, string& id, string& counts)
        {
            std::string line;
            if (std::getline(stream, line)) {
                if ('>' == line[0]) {

                    id.assign(line.begin() + 1, line.end());

                    //std::string sequence;
                    /*while (stream.good() && '>' != stream.peek()) {
                        std::getline(stream, line);
                        sequence += line;
                    }*/

                    // We assume the counts are on a single line
                    std::getline(stream, line);
                    counts = line;

                    stream.clear();

                    return 0;
                }
                else {
                    stream.clear(std::ios_base::failbit);
                }
            }

            return -1;
        }


        /**
         *  Finds a particular fasta header in a fasta file and returns the associated sequence
         **/
        static int getEntryFromFasta(const path& fasta_path, const string& header, string& sequence)
        {
            // Setup stream to sequence file and check all is well
            std::ifstream inputStream (fasta_path.c_str());

            if (!inputStream.is_open())
            {
                std::cerr << "ERROR: Could not open the sequence file: " << fasta_path << endl;
                return 0;
            }

            string id;
            string seq;

            // Read a record
            while(inputStream.good() && readRecord(inputStream, id, seq) == 0)
            {
                stringstream ssHeader;
                ssHeader << id;

                if (header.compare(ssHeader.str()) == 0)
                {
                    inputStream.close();
                    sequence = seq;

                    return 0;
                }
                seq.clear();
            }

            inputStream.close();
            return -1;
        }

        /**
         *  Finds the nth entry from the fasta file and returns the associated sequence
         **/
        static int getEntryFromFasta(const path& fasta_path, uint32_t n, string& header, string& sequence)
        {
            // Setup stream to sequence file and check all is well
            std::ifstream inputStream (fasta_path.c_str());

            if (!inputStream.is_open())
            {
                std::cerr << "ERROR: Could not open the sequence file: " << fasta_path << endl;
                return 0;
            }


            std::string id;
            std::string seq;
            uint32_t i = 1;

            // Read a record
            while(inputStream.good() && readRecord(inputStream, id, seq) == 0)
            {
                if (i == n)
                {
                    inputStream.close();

                    header.swap(id);
                    sequence = seq;
                    return 0;
                }

                seq.clear();
                i++;
            }

            inputStream.close();
            return -1;
        }
        
        
        static string helpMessage() {
            return string("Usage: kat plot profile [options] <sect_profile_file>\n\n") + 
                    "Create Sequence Coverage Plot.\n\n" \
                    "Shows K-mer coverage level across an sequence.\n\n" \
                    "Options";
        }

        static string autoTitle(string& title, string& header)
        {
            std::ostringstream output_str;
            output_str << "Sequence Coverage Plot" << ": " << header;
            return title.compare("Sequence Coverage Plot") == 0 ? output_str.str() : title;
        }

    public:
        
        
        
        static int main(int argc, char *argv[]) {

            const string DEFAULT_TITLE      = "Sequence Coverage Plot";
            const string DEFAULT_X_LABEL    = "X";
            const string DEFAULT_Y_LABEL    = "Y";
            const int32_t DEFAULT_X_MAX     = 1000;
            const int32_t DEFAULT_Y_MAX     = 1000;
            const uint32_t DEFAULT_FASTA_INDEX = 0;
        
            path        sect_file;           
            string      output_type;
            path        output;
            string      title;
            string      x_label;
            string      y_label;
            uint16_t    width;
            uint16_t    height;
            uint32_t    x_max;
            uint32_t    y_max;
            bool        y_logscale;
            uint32_t    fasta_index;
            string      fasta_header;
            bool        verbose;
            bool        help;
            
            struct winsize w;
            ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    

            // Declare the supported options.
            po::options_description generic_options(PlotProfile::helpMessage(), w.ws_col);
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
                    ("x_max,x", po::value<uint32_t>(&x_max)->default_value(DEFAULT_X_MAX),
                        "Maximum value for the x-axis (value taken from matrix metadata if present)")
                    ("y_max,y", po::value<uint32_t>(&y_max)->default_value(DEFAULT_Y_MAX),
                        "Maximum value for the y-axis (value taken from matrix metadata if present)")
                    ("width,w", po::value<uint16_t>(&width)->default_value(1024),
                        "Width of canvas")
                    ("height,h", po::value<uint16_t>(&height)->default_value(1024),
                        "Height of canvas")
                    ("index,n", po::value<uint32_t>(&fasta_index)->default_value(DEFAULT_FASTA_INDEX),
                        "Index of fasta entry to plot.  First entry is 1.")
                    ("header,d", po::value<string>(&fasta_header),
                        "Fasta header of fasta entry to plot.  NOTE: \'--header\' has priority over index")
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
                    ("sect_file", po::value<path>(&sect_file), "Path to the sect profile file to plot.")                    
                    ;

            // Positional option for the input bam file
            po::positional_options_description p;
            p.add("sect_file", 1);


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
            if (!bfs::exists(sect_file) && !bfs::symbolic_link_exists(sect_file)) {
                cerr << endl << "Could not find matrix file at: " << sect_file << "; please check the path and try again." << endl << endl;
                return 1;
            }

            string header;
            string coverages;

            if (!fasta_header.empty()) {
                header.assign(fasta_header);
                getEntryFromFasta(sect_file, header, coverages);
            }
            else if (fasta_index > 0) {
                getEntryFromFasta(sect_file, fasta_index, header, coverages);
            }

            if (coverages.length() == 0) {
                cerr << "Could not find requested fasta header in sect coverages fasta file" << endl;
            }
            else {
                if (verbose)
                    cerr << "Found requested sequence : " << header << endl << coverages << endl << endl;

                // Split coverages
                vector<uint32_t> cvs = kat::splitUInt32(coverages, ' ');

                uint32_t maxCvgVal = y_max != DEFAULT_Y_MAX ? y_max : (*(std::max_element(cvs.begin(), cvs.end())) + 1);

                string t = autoTitle(title, header);

                if (verbose)
                    cerr << "Acquired K-mer counts" << endl;

                // Initialise gnuplot
                Gnuplot profile_plot = Gnuplot("lines");

                profile_plot.configurePlot(output_type, output.string(), width, height);

                profile_plot.set_title(t);
                profile_plot.set_xlabel(x_label);
                profile_plot.set_ylabel(y_label);
                profile_plot.set_xrange(0, cvs.size());
                profile_plot.set_yrange(0, maxCvgVal);
                
                profile_plot.cmd("set style data linespoints");

                std::ostringstream data_str;

                for(uint32_t i = 0; i < cvs.size(); i++)
                {
                    uint32_t index = i+1;
                    double val = y_logscale ? (double)std::log(cvs[i]) : (double)cvs[i];
                    data_str << index << " " << val << "\n";
                }

                std::ostringstream plot_str;

                plot_str << "plot '-'\n" << data_str.str() << "e\n";

                profile_plot.cmd(plot_str.str());

                if (verbose)
                    cerr << "Plotted data: " << plot_str.str() << endl;
            }

            return 0;
        }
        
    };
}
