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
using std::string;
using std::ifstream;
using std::istringstream;
using std::ostringstream;

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

#include "inc/spectra_helper.hpp"
using kat::SpectraHelper;

typedef boost::error_info<struct PlotDensityError,string> PlotDensityErrorInfo;
struct PlotDensityException: virtual boost::exception, virtual std::exception { };

namespace kat {
    
    const string DEFAULT_PD_TITLE      = "Density plot";
    const string DEFAULT_PD_X_LABEL    = "X";
    const string DEFAULT_PD_Y_LABEL    = "Y";
    const string DEFAULT_PD_Z_LABEL    = "Z";
    const uint32_t DEFAULT_PD_X_MAX    = 1000;
    const uint32_t DEFAULT_PD_Y_MAX    = 1000;
    const uint64_t DEFAULT_PD_Z_MAX    = 10000;
    const string DEFAULT_PD_OUTPUT_TYPE   = "png"; 
    const uint16_t DEFAULT_PD_WIDTH       = 1024;
    const uint16_t DEFAULT_PD_HEIGHT      = 1024;

    
    class PlotDensity {
    private:
        
        path        mxFile;           
        string      outputType;
        path        output;
        string      title;
        string      xLabel;
        string      yLabel;
        string      zLabel;
        uint16_t    width;
        uint16_t    height;
        uint32_t    xMax;
        uint32_t    yMax;
        uint64_t    zMax;
        bool        verbose;
        
    public:
        
        PlotDensity(const path& _mxFile, const path& _outFile) {
            mxFile = _mxFile;
            output = _outFile;
            outputType = DEFAULT_PD_OUTPUT_TYPE;
            title = DEFAULT_PD_TITLE;
            xLabel = DEFAULT_PD_X_LABEL;
            yLabel = DEFAULT_PD_Y_LABEL;
            zLabel = DEFAULT_PD_Z_LABEL;
            width = DEFAULT_PD_WIDTH;
            height = DEFAULT_PD_HEIGHT;
            xMax = 0;
            yMax = 0;
            zMax = 0;            
            verbose = false;
        }
        
        uint16_t getHeight() const {
            return height;
        }

        void setHeight(uint16_t height) {
            this->height = height;
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

        string getZLabel() const {
            return zLabel;
        }

        void setZLabel(string zLabel) {
            this->zLabel = zLabel;
        }

        uint64_t getZMax() const {
            return zMax;
        }

        void setZMax(uint64_t zMax) {
            this->zMax = zMax;
        }
        
        void plot();

            
    protected:
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
       
        static int main(int argc, char *argv[]);
    };
}
