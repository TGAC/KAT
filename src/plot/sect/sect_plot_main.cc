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
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>

#include <gnuplot/gnuplot_i.hpp>
#include <string/str_utils.hpp>

#include "sect_plot_args.hpp"
#include "sect_plot_main.hpp"

using std::string;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::vector;
using std::endl;

using kat::SectPlotArgs;

namespace kat
{
    int readRecord(std::ifstream& stream, string& id, string& counts)
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



    // Finds a particular fasta header in a fasta file and returns the associated sequence
    int getEntryFromFasta(const string& fasta_path, const string& header, string& sequence)
    {
        // Setup stream to sequence file and check all is well
        std::ifstream inputStream (fasta_path.c_str());

        if (!inputStream.is_open())
        {
            std::cerr << "ERROR: Could not open the sequence file: " << fasta_path << endl;
            return NULL;
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

    // Finds the nth entry from the fasta file and returns the associated sequence
    int getEntryFromFasta(const string& fasta_path, uint32_t n, string& header, string& sequence)
    {
        // Setup stream to sequence file and check all is well
        std::ifstream inputStream (fasta_path.c_str());

        if (!inputStream.is_open())
        {
            std::cerr << "ERROR: Could not open the sequence file: " << fasta_path << endl;
            return NULL;
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

}


// Start point
int kat::sectPlotStart(int argc, char *argv[])
{
    // Parse args
    SectPlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    string header;
    string coverages;

    if (!args.fasta_header.empty())
    {
        header.assign(args.fasta_header);
        getEntryFromFasta(args.sect_file_arg, header, coverages);
    }
    else if (args.fasta_index > 0)
    {
        getEntryFromFasta(args.sect_file_arg, args.fasta_index, header, coverages);
    }

    if (coverages.length() == 0)
    {
        cerr << "Could not find requested fasta header in sect coverages fasta file" << endl;
    }
    else
    {
        if (args.verbose)
            cerr << "Found requested sequence : " << header << endl << coverages << endl << endl;

        // Split coverages
        vector<uint32_t> cvs;
        kat::split(coverages, cvs, ' ');

        uint32_t maxCvgVal = args.y_max != kat::DEFAULT_Y_MAX ? args.y_max : (*(std::max_element(cvs.begin(), cvs.end())) + 1);

        string title = args.autoTitle(header);

        if (args.verbose)
            cerr << "Acquired kmer counts" << endl;

        // Initialise gnuplot
        Gnuplot sect_plot = Gnuplot("lines");

        // Work out the output path to use (either user specified or auto generated)
        string output_path = args.determineOutputPath();

        sect_plot.configurePlot(args.output_type, output_path, args.width, args.height);

        sect_plot.set_title(title);
        sect_plot.set_xlabel(args.x_label);
        sect_plot.set_ylabel(args.y_label);
        sect_plot.set_xrange(0, cvs.size());
        sect_plot.set_yrange(0, maxCvgVal);

        sect_plot.cmd("set style data linespoints");

        std::ostringstream data_str;

        for(uint32_t i = 0; i < cvs.size(); i++)
        {
            uint32_t index = i+1;
            data_str << index << " " << cvs[i] << "\n";
        }

        std::ostringstream plot_str;

        plot_str << "plot '-'\n" << data_str.str() << "e\n";

        sect_plot.cmd(plot_str.str());

        if (args.verbose)
            cerr << "Plotted data: " << plot_str.str() << endl;
    }

    return 0;
}
