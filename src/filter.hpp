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

#include <string>
#include <vector>
using std::string;
using std::vector;

#include <boost/exception/exception.hpp>
#include <boost/exception/info.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::path;

typedef boost::error_info<struct KatFilterError,string> KatFilterErrorInfo;
struct KatFilterException: virtual boost::exception, virtual std::exception { };


namespace kat {


    class Filter {

    public:

        enum FilterMode {
            KMER,
            SEQ
        };

        static bool validateFilterOutputType();

        static int main(int argc, char *argv[]);

    protected:

        static FilterMode parseMode(const string& mode);

        static string helpMessage() {
            return string("Usage: kat filter <mode>\n\n") +
                    "Filtering tools\n\n" +
                    "First argument should be the filter mode you wish to use:\n" \
                    "  * kmer:            Filters a jellyfish k-mer hash using user defined properties.\n" \
                    "  * seq:             Filters sequences in a file based on k-mers in a given hash\n\n" \
                    "Options";
        }


    };

}
