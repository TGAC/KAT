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

#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::path;

#include <kat/gnuplot_i.hpp>
#include <kat/str_utils.hpp>

namespace kat {

    class PlotProfile {

    protected:
        static int readRecord(std::ifstream& stream, string& id, string& counts);


        /**
         *  Finds a particular fasta header in a fasta file and returns the associated sequence
         **/
        static int getEntryFromFasta(const path& fasta_path, const string& header, string& sequence);

        /**
         *  Finds the nth entry from the fasta file and returns the associated sequence
         **/
        static int getEntryFromFasta(const path& fasta_path, uint32_t n, string& header, string& sequence);

        static string autoTitle(string& title, string& header);


    protected:
        static string helpMessage() {
            return string("Usage: kat plot profile [options] <sect_profile>\n\n") +
                    "Create Sequence Profile Plots.\n\n" +
                    "Shows k-mer coverage across one or more sequences.\n\n" \
                    "Options";
        }


    public:

       static int main(int argc, char *argv[]);
    };

}
