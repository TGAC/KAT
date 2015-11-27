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

#include <string>
#include <vector>
using std::string;
using std::vector;

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include "inc/kat_fs.hpp"
using kat::KatFS;

namespace kat {
    
    typedef boost::error_info<struct KatPlotError,string> KatPlotErrorInfo;
    struct KatPlotException: virtual boost::exception, virtual std::exception { };
    
    class Plot {
    
    public:
        
        enum PlotMode {
            DENSITY,
            PROFILE,
            SPECTRA_CN,
            SPECTRA_HIST,
            SPECTRA_MX
        };
        
        static bool validatePlotOutputType();
        
        static int main(int argc, char *argv[], const KatFS& fs);
        
    private:
        static wchar_t* convertCharToWideChar(const char* c);
        
    protected:

        static PlotMode parseMode(const string& mode);
                
        static void executePythonPlot(const PlotMode mode, int argc, char *argv[], const KatFS& fs);
        
        static void executeGnuplotPlot(const PlotMode mode, int argc, char *argv[]);        
        
        static path getPythonScript(const PlotMode mode);
        
        static string helpMessage() {
            return string("Usage: kat plot <mode>\n\n") +
                    "Create K-mer Plots\n\n" +
                    "First argument should be the plot mode you wish to use:\n" \
                    "  * density:         Creates a density plot from a matrix created with the \"comp\" tool or the \"GCP\"\n" \
                    "                     tool.  Typically this is used to compare two K-mer hashes produced by different NGS\n" \
                    "                     reads, or to represent the kmer coverage vs GC count plots.\n" \
                    "  * profile:         Creates a K-mer coverage plot for a single sequence.  Takes in fasta coverage output\n" \
                    "                     coverage from the \"sect\" tool\n" \
                    "  * spectra-cn:      Creates a stacked histogram using a matrix created with the \"comp\" tool.  Typically\n" \
                    "                     this is used to compare a jellyfish hash produced from a read set to a jellyfish hash\n" \
                    "                     produced from an assembly. The plot shows the amount of distinct K-mers absent, as well\n" \
                    "                     as the copy number variation present within the assembly.\n" \
                    "  * spectra-hist:    Creates a K-mer spectra plot for a set of K-mer histograms produced either by jellyfish-\n" \
                    "                     histo or kat-histo.\n" \
                    "  * spectra-mx:      Creates a K-mer spectra plot for a set of K-mer histograms that are derived from\n" \
                    "                     selected rows or columns in a matrix produced by the \"comp\".\n\n" \
                    "Options";
        }
        
    
    };
    
}

