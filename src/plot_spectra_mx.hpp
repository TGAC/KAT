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

#include <boost/exception/exception.hpp>
#include <boost/exception/info.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
namespace bfs = boost::filesystem;
using bfs::path;

#include <kat/gnuplot_i.hpp>
#include <kat/str_utils.hpp>
#include <kat/sparse_matrix.hpp>
#include <kat/matrix_metadata_extractor.hpp>

typedef boost::error_info<struct PlotSpectraMxError,string> PlotSpectraMxErrorInfo;
struct PlotSpectraMxException: virtual boost::exception, virtual std::exception { };

namespace kat {

    const string DEFAULT_PSMX_TITLE         = "Spectra MX plot";
    const string DEFAULT_PSMX_OUTPUT_TYPE   = "png";
    const string DEFAULT_PSMX_X_LABEL       = "X";
    const string DEFAULT_PSMX_Y_LABEL       = "Y";
    const uint32_t DEFAULT_PSMX_X_MAX       = 1000;
    const uint32_t DEFAULT_PSMX_Y_MAX       = 1000;
    const uint16_t DEFAULT_PSMX_WIDTH       = 1024;
    const uint16_t DEFAULT_PSMX_HEIGHT      = 1024;


    class PlotSpectraMx {

    private:

        path        mxFile;
        string      outputType;
        path        output;
        string      title;
        string      xLabel;
        string      yLabel;
        uint16_t    width;
        uint16_t    height;
        bool        intersection;
        string      list;
        uint16_t    excCutoffD1;
        uint16_t    excCutoffD2;
        uint32_t    xMin;
        uint32_t    yMin;
        uint32_t    xMax;
        uint64_t    yMax;
        bool        xLogscale;
        bool        yLogscale;
        bool        verbose;

    public:

        PlotSpectraMx(const path& _mxFile, const path& _outFile) {
            mxFile = _mxFile;
            output = _outFile;
            outputType = DEFAULT_PSMX_OUTPUT_TYPE;
            title = DEFAULT_PSMX_TITLE;
            xLabel = DEFAULT_PSMX_X_LABEL;
            yLabel = DEFAULT_PSMX_Y_LABEL;
            width = DEFAULT_PSMX_WIDTH;
            height = DEFAULT_PSMX_HEIGHT;
            xMax = 0;
            yMax = 0;
            intersection = false;
            list = "";
            excCutoffD1 = 1;
            excCutoffD2 = 1;
            xMin = 0;
            yMin = 0;
            xLogscale = false;
            yLogscale = false;
            verbose = false;
        }

        uint16_t getExcCutoffD1() const {
            return excCutoffD1;
        }

        void setExcCutoffD1(uint16_t excCutoffD1) {
            this->excCutoffD1 = excCutoffD1;
        }

        uint16_t getExcCutoffD2() const {
            return excCutoffD2;
        }

        void setExcCutoffD2(uint16_t excCutoffD2) {
            this->excCutoffD2 = excCutoffD2;
        }

        uint16_t getHeight() const {
            return height;
        }

        void setHeight(uint16_t height) {
            this->height = height;
        }

        bool isIntersection() const {
            return intersection;
        }

        void setIntersection(bool intersection) {
            this->intersection = intersection;
        }

        string getList() const {
            return list;
        }

        void setList(string list) {
            this->list = list;
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

        uint64_t getYMax() const {
            return yMax;
        }

        void setYMax(uint64_t yMax) {
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

        string getDataFromList(const path& mx_file, string list);


        string getIntersectionData(const path& mx_file, uint16_t exc_cutoff_d1, uint16_t exc_cutoff_d2);

        static string helpMessage() {
            return string("Usage: kat plot spectra-mx [options] {--list <comma_seperated_list> | --intersection} <mx_file>\n\n") +
                    "Creates K-mer Spectra Plot from selected rows and/or columns in a \"comp\" matrix.\n\n" +
                    "Produces K-mer spectras from rows or columns in a matrix generated by \"kat comp\".  This tool is " \
                    "designed to plot line graphs for one or more histograms, each histogram being represented by a single row or column " \
                    "in the matrix.\n" \
                    "This tool also has a special mode for showing shared and exclusive content between two different samples. This mode " \
                    "takes the first row and column of the matrix representing content which is found exclusively in " \
                    "each sample.  Two more lines are plotting, one which has each following row summed, and the other that has " \
                    "each following column summed.  These two plots represent the shared content for each sample.  This mode can " \
                    "be activated using the \"--intersection\" flag.\n" \
                    "Alternatively, you can select specific rows and columns from the matrix using a comma separated list " \
                    "identified with the \"--list\" option.  Each element in the list should start with either a 'c' or a 'r' " \
                    "indicating whether or not the column or row is requested.  Then the element should contain a number " \
                    "indicating which column or row to select.  For example: \"--list c0,r1\" will select column 0 and row 1. " \
                    "Note: spaces are not tolerated in this list.\n\n" \
                    "Options";
        }


    public:

        static int main(int argc, char *argv[]);
    };
}
