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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <sys/ioctl.h>
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::istringstream;
using std::ostringstream;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include <kat/gnuplot_i.hpp>
#include <kat/sparse_matrix.hpp>
#include <kat/matrix_metadata_extractor.hpp>
#include <kat/spectra_helper.hpp>
using kat::SpectraHelper;

#include "plot_density.hpp"

bool kat::PlotDensity::plot() {
    // Check input file exists
    if (!bfs::exists(mxFile) && !bfs::symbolic_link_exists(mxFile))
    {
        BOOST_THROW_EXCEPTION(PlotDensityException() << PlotDensityErrorInfo(string(
            "Could not find matrix file at: ") + mxFile.string() + "; please check the path and try again.")); 
    }
    
    // Determine auto ranges
    SparseMatrix<uint64_t> tmx(mxFile);
    vector<Pos> cumulativeSpectraX(tmx.height());
    vector<Pos> cumulativeSpectraY(tmx.width());

    for(size_t i = 1; i < tmx.height() - 1; i++) {
        cumulativeSpectraX[i].first = i;
        cumulativeSpectraX[i].second += tmx.sumRow(i);
    }

    Pos posX = SpectraHelper::findPeak(cumulativeSpectraX);

    for(size_t i = 1; i < tmx.width() - 1; i++) {
        cumulativeSpectraY[i].first = i;
        cumulativeSpectraY[i].second += tmx.sumColumn(i);
    }

    Pos posY = SpectraHelper::findPeak(cumulativeSpectraY, false);

    // We choose the min rather than max because on the Y axis we will include 
    // the error kmer peak in the data, which we want to avoid.  Also normally
    // there is no harm in oversaturating the image.
    uint32_t maxZ = std::min(posX.second, posY.second);

    // If possible estimate a reasonable x and y range directly from the data
    uint16_t autoXMax = posX.first > 0 ? posX.first * 3 : 1000;
    uint16_t autoYMax = posY.first > 0 ? posY.first * 3 : 1000;
    uint32_t autoZMax = posX.first > 0 && posY.first > 0 ? maxZ / 7 : 10000;    // 7 seems to work well

    // Don't go over any limits in the data for the X and Y axis
    autoXMax = std::min((uint16_t)mme::getNumeric(mxFile, mme::KEY_NB_COLUMNS), autoXMax);
    autoYMax = std::min((uint16_t)mme::getNumeric(mxFile, mme::KEY_NB_ROWS), autoYMax);            

    // Get plotting properties, either from file, or user.  User args have precedence.
    uint16_t x_range = xMax != 0 && xMax != DEFAULT_PD_X_MAX ? xMax : autoXMax;
    uint16_t y_range = yMax != 0 && yMax != DEFAULT_PD_Y_MAX ? yMax : autoYMax;
    uint32_t z_range = zMax != 0 && zMax != DEFAULT_PD_Z_MAX ? zMax : autoZMax;

    string xl = !boost::equals(xLabel, DEFAULT_PD_X_LABEL) ? xLabel : mme::getString(mxFile, mme::KEY_X_LABEL);
    string yl = !boost::equals(yLabel, DEFAULT_PD_Y_LABEL) ? yLabel : mme::getString(mxFile, mme::KEY_Y_LABEL);
    string zl = !boost::equals(zLabel, DEFAULT_PD_Z_LABEL) ? zLabel : mme::getString(mxFile, mme::KEY_Z_LABEL);

    string t = !boost::equals(title, DEFAULT_PD_TITLE) ? title : mme::getString(mxFile, mme::KEY_TITLE);

    bool transpose = mme::getNumeric(mxFile, mme::KEY_TRANSPOSE) == 0 ? false : true;

    xl = xl.empty() ? DEFAULT_PD_X_LABEL : xl;
    yl = yl.empty() ? DEFAULT_PD_Y_LABEL : yl;
    zl = zl.empty() ? DEFAULT_PD_Z_LABEL : zl;

    t = t.empty() ? DEFAULT_PD_TITLE : t;

    if (verbose)
    {
        cerr << "Actual variables used to create plot:" << endl;
        cerr << "Output Path: " << output << endl;
        cerr << "X Range: " << x_range << endl;
        cerr << "Y Range: " << y_range << endl;
        cerr << "Z Range: " << z_range << endl;
        cerr << "X Label: " << xLabel << endl;
        cerr << "Y Label: " << yLabel << endl;
        cerr << "Z Label: " << zLabel << endl;
        cerr << "Title: " << title << endl;
    }


    // Start defining the plot
    Gnuplot density("lines");

    density.configurePlot(outputType, output.string(), width, height);

    density.set_title(title);
    density.set_xlabel(xLabel);
    density.set_ylabel(yLabel);

    std::ostringstream cblabelstr;
    cblabelstr << "set cblabel \"" << zLabel << "\"";
    density.cmd(cblabelstr.str());

    density.set_xrange(0, x_range);
    density.set_yrange(0, y_range);

    //flame.set_xlogscale();
    //flame.set_ylogscale();
    //flame.set_zlogscale();

    density.cmd("set palette rgb 21,22,23");
    density.cmd("set size ratio 1");

    std::ostringstream rangestr;
    rangestr << "set cbrange [0:" << z_range << "]";
    density.cmd(rangestr.str());

    // Transpose the matrix and store in ostream
    ostringstream data;
    SparseMatrix<uint64_t> mx(mxFile);
    mx.printMatrix(data, transpose);

    // Plot the transposed matrix as image
    std::ostringstream plotstr;
    plotstr << "plot '-' matrix with image" << endl
            << data.str()
            << "e" << endl
            << "e" << endl;

    if (!density.is_valid()) {
        return false;
    }
    
    density.cmd(plotstr.str());
    
    return true;
}

int kat::PlotDensity::main(int argc, char *argv[]) {

    path        mx_file;           
    string      output_type;
    path        output;
    string      title;
    string      x_label;
    string      y_label;
    string      z_label;
    uint16_t    width;
    uint16_t    height;
    uint32_t    x_max;
    uint32_t    y_max;
    uint64_t    z_max;
    bool        verbose;
    bool        help;
    
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    

    // Declare the supported options.
    po::options_description generic_options(PlotDensity::helpMessage(), w.ws_col);
    generic_options.add_options()
            ("output_type,p", po::value<string>(&output_type)->default_value("png"), 
                "The plot file type to create: png, ps, pdf.  Warning... if pdf is selected please ensure your gnuplot installation can export pdf files.")
            ("output,o", po::value<path>(&output),
                "The path to the output file")
            ("title,t", po::value<string>(&title)->default_value(DEFAULT_PD_TITLE),
                "Title for plot")
            ("x_label,a", po::value<string>(&x_label)->default_value(DEFAULT_PD_X_LABEL),
                "Label for the x-axis (value taken from matrix metadata if present)")
            ("y_label,b", po::value<string>(&y_label)->default_value(DEFAULT_PD_Y_LABEL),
                "Label for the y-axis (value taken from matrix metadata if present)")
            ("z_label,c", po::value<string>(&z_label)->default_value(DEFAULT_PD_Z_LABEL),
                "Label for the z-axis (value taken from matrix metadata if present)")
            ("x_max,x", po::value<uint32_t>(&x_max)->default_value(DEFAULT_PD_X_MAX),
                "Maximum value for the x-axis (value taken from matrix metadata if present)")
            ("y_max,y", po::value<uint32_t>(&y_max)->default_value(DEFAULT_PD_Y_MAX),
                "Maximum value for the y-axis (value taken from matrix metadata if present)")
            ("z_max,z", po::value<uint64_t>(&z_max)->default_value(DEFAULT_PD_Z_MAX),
                "Maximum value for the z-axis (value taken from matrix metadata if present)")
            ("width,w", po::value<uint16_t>(&width)->default_value(1024),
                "Width of canvas")
            ("height,h", po::value<uint16_t>(&height)->default_value(1024),
                "Height of canvas")
            ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                "Print extra information.")
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
            ("mx_file", po::value<path>(&mx_file), "Path to the matrix file to plot.")                    
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
    
    if (output.empty()) {
        BOOST_THROW_EXCEPTION(PlotDensityException() << PlotDensityErrorInfo(string(
            "Output file not specified.  Please use the '-o' option."))); 
    }
    
        
    PlotDensity pd(mx_file, output);
    pd.setHeight(height);
    pd.setOutputType(output_type);
    pd.setTitle(title);
    pd.setVerbose(verbose);
    pd.setWidth(width);
    pd.setXLabel(x_label);
    pd.setXMax(x_max);
    pd.setYLabel(y_label);
    pd.setYMax(y_max);
    pd.setZLabel(z_label);
    pd.setZMax(z_max);            
    pd.plot();

    return 0;
}
