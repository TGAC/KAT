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

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include <gnuplot/gnuplot_i.hpp>
#include "inc/spectra_helper.hpp"
#include "inc/matrix/sparse_matrix.hpp"
using kat::SpectraHelper;

typedef boost::error_info<struct PlotSpectraCnError,string> PlotSpectraCnErrorInfo;
struct PlotSpectraCnException: virtual boost::exception, virtual std::exception { };

namespace kat {
    
    const string     DEFAULT_PSCN_TITLE           = "Spectra Copy Number Plot";
    const string     DEFAULT_PSCN_OUTPUT_TYPE     = "png"; 
    const string     DEFAULT_PSCN_X_LABEL         = "X";
    const string     DEFAULT_PSCN_Y_LABEL         = "Y";
    const uint16_t   DEFAULT_PSCN_WIDTH           = 1024;
    const uint16_t   DEFAULT_PSCN_HEIGHT          = 1024;
    const uint16_t   DEFAULT_MAX_DUPLICATION      = 6;

    class PlotSpectraCn {
    private:
        
        path        mxFile;           
        string      outputType;
        path        output;
        string      title;
        string      xLabel;
        string      yLabel;
        uint16_t    width;
        uint16_t    height;
        uint32_t    xMax;
        uint32_t    yMax;
        bool        ignoreAbsent;
        uint16_t    maxDuplication;
        string      columns;
        bool        cumulative;
        bool        verbose;
    
            
    public:
            
        PlotSpectraCn(const path& _mxFile, const path& _outFile) {
            mxFile = _mxFile;
            output = _outFile;
            outputType = DEFAULT_PSCN_OUTPUT_TYPE;
            title = DEFAULT_PSCN_TITLE;
            xLabel = DEFAULT_PSCN_X_LABEL;
            yLabel = DEFAULT_PSCN_Y_LABEL;
            width = DEFAULT_PSCN_WIDTH;
            height = DEFAULT_PSCN_HEIGHT;
            xMax = 0;
            yMax = 0;
            ignoreAbsent = false;
            maxDuplication = DEFAULT_MAX_DUPLICATION;
            columns = "";
            cumulative = false;
            verbose = false;
        }
        
        string getColumns() const {
            return columns;
        }

        void setColumns(string columns) {
            this->columns = columns;
        }

        bool isCumulative() const {
            return cumulative;
        }

        void setCumulative(bool cumulative) {
            this->cumulative = cumulative;
        }

        uint16_t getHeight() const {
            return height;
        }

        void setHeight(uint16_t height) {
            this->height = height;
        }

        bool isIgnoreAbsent() const {
            return ignoreAbsent;
        }

        void setIgnoreAbsent(bool ignoreAbsent) {
            this->ignoreAbsent = ignoreAbsent;
        }

        uint16_t getMaxDuplication() const {
            return maxDuplication;
        }

        void setMaxDuplication(uint16_t maxDuplication) {
            this->maxDuplication = maxDuplication;
        }

        path getMxFile() const {
            return mxFile;
        }

        void setMxFile(path mxFile) {
            this->mxFile = mxFile;
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

        uint32_t getXMax() const {
            return xMax;
        }

        void setXMax(uint32_t xMax) {
            this->xMax = xMax;
        }

        string getYLabel() const {
            return yLabel;
        }

        void setYLabel(string yLabel) {
            this->yLabel = yLabel;
        }

        uint32_t getYMax() const {
            return yMax;
        }

        void setYMax(uint32_t yMax) {
            this->yMax = yMax;
        }

        bool plot();
        
    protected:

        string createSinglePlotString(const path& data_file, uint16_t idx, uint16_t level_count, bool cumulative);
        
        string createLineStyleStr(uint16_t i, const char* colour);

        shared_ptr<vector<uint16_t>> getStandardCols(bool ignoreAbsent, uint16_t maxDuplication);

        shared_ptr<vector<uint16_t>> getUserDefinedCols(const string& columns);
        
        static string helpMessage() {
            return string("Usage: kat plot spectra-cn [options] <matrix_file>\n\n") +
                    "Creates a stacked histogram showing the level of duplication in an assembly.\n\n" \
                    "Shows K-mer duplication levels, which correspond to copy number variation within an assembly by comparing " \
                    "K-mers found in sequenced reads, to K-mers found in an assembly of those reads. Uses matrix output from the " \
                    "\"kat comp\" tool.\n\n" \
                    "Options";
        }

    public:

        static int main(int argc, char *argv[]);
    };
    
}