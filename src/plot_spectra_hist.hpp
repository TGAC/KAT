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

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include <kat/gnuplot_i.hpp>
#include <kat/str_utils.hpp>
#include <kat/spectra_helper.hpp>
using kat::SpectraHelper;

typedef boost::error_info<struct PlotSpectraHistError,string> PlotSpectraHistErrorInfo;
struct PlotSpectraHistException: virtual boost::exception, virtual std::exception { };

namespace kat {
    
    const string DEFAULT_SH_TITLE       = "Kmer histograms";
    const string DEFAULT_SH_X_LABEL     = "X";
    const string DEFAULT_SH_Y_LABEL     = "Y";
    const string DEFAULT_SH_OUTPUT_TYPE   = "png"; 
    const uint32_t DEFAULT_SH_X_MAX       = 1000;
    const uint32_t DEFAULT_SH_Y_MAX       = 1000;
    const uint16_t DEFAULT_SH_WIDTH       = 1024;
    const uint16_t DEFAULT_SH_HEIGHT      = 1024;

    
    class PlotSpectraHist {
    
    private:
        vector<path> histoPaths;        
        string      outputType;
        path        output;
        string      title;
        string      xLabel;
        string      yLabel;
        uint16_t    width;
        uint16_t    height;
        uint32_t    xMin;
        uint32_t    yMin;
        uint32_t    xMax;
        uint32_t    yMax;
        bool        xLogscale;
        bool        yLogscale;
        bool        verbose;
           
    public:
        PlotSpectraHist(const vector<path>& input, const path& _output) {
            histoPaths = input;
            output = _output;
            outputType = DEFAULT_SH_OUTPUT_TYPE;
            title = DEFAULT_SH_TITLE;
            xLabel = DEFAULT_SH_X_LABEL;
            yLabel = DEFAULT_SH_Y_LABEL;
            width = DEFAULT_SH_WIDTH;
            height = DEFAULT_SH_HEIGHT;
            xMax = 0;
            yMax = 0;
            xMin = 0;
            yMin = 0;
            xLogscale = false;
            yLogscale = false;
            verbose = false;
        }
        
        uint16_t getHeight() const {
            return height;
        }

        void setHeight(uint16_t height) {
            this->height = height;
        }

        vector<path> getHistoPaths() const {
            return histoPaths;
        }

        void setHistoPaths(vector<path> histoPaths) {
            this->histoPaths = histoPaths;
        }

        path getOutput() const {
            return output;
        }

        void setOutput(path output) {
            this->output = output;
        }

        string getOutputType() const {
            return outputType;
        }

        void setOutputType(string outputType) {
            this->outputType = outputType;
        }

        string getTitle() const {
            return title;
        }

        void setTitle(string title) {
            this->title = title;
        }

        bool isVerbose() const {
            return verbose;
        }

        void setVerbose(bool verbose) {
            this->verbose = verbose;
        }

        uint16_t getWidth() const {
            return width;
        }

        void setWidth(uint16_t width) {
            this->width = width;
        }

        string getXLabel() const {
            return xLabel;
        }

        void setXLabel(string xLabel) {
            this->xLabel = xLabel;
        }

        bool isXLogscale() const {
            return xLogscale;
        }

        void setXLogscale(bool xLogscale) {
            this->xLogscale = xLogscale;
        }

        uint32_t getXMax() const {
            return xMax;
        }

        void setXMax(uint32_t xMax) {
            this->xMax = xMax;
        }

        uint32_t getXMin() const {
            return xMin;
        }

        void setXMin(uint32_t xMin) {
            this->xMin = xMin;
        }

        string getYLabel() const {
            return yLabel;
        }

        void setYLabel(string yLabel) {
            this->yLabel = yLabel;
        }

        bool isYLogscale() const {
            return yLogscale;
        }

        void setYLogscale(bool yLogscale) {
            this->yLogscale = yLogscale;
        }

        uint32_t getYMax() const {
            return yMax;
        }

        void setYMax(uint32_t yMax) {
            this->yMax = yMax;
        }

        uint32_t getYMin() const {
            return yMin;
        }

        void setYMin(uint32_t yMin) {
            this->yMin = yMin;
        }

        bool plot();
        
    protected:
        
        static string helpMessage() {
            return string("Usage: kat plot spectra-hist [options] <histo_file> [<histo_file> ...]*\n\n") +
                    "Creates K-mer Spectra Plot from one or more histograms.\n\n" +
                    "Produces K-mer spectras from \"kat hist\" or \"jellyfish histo\" output.  This tool is designed to plot line " \
                    "graphs of one or more histograms.  The idea is to be able to compare total K-mer counts between different " \
                    "datasets.\n\n" \
                    "Options";
        }
        
    public:
  
        // Start point
        static int main(int argc, char *argv[]);

    };
}
