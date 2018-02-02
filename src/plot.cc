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

#include <iostream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <unistd.h>
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include <kat/kat_fs.hpp>
#include <kat/pyhelper.hpp>

#include "plot.hpp"
using kat::Plot;


Plot::PlotMode kat::Plot::parseMode(const string& mode) {

    string upperMode = boost::to_upper_copy(mode);

    if (upperMode == string("DENSITY")) {
        return DENSITY;
    }
    else if (upperMode == string("PROFILE")) {
        return PROFILE;
    }
    else if (upperMode == string("SPECTRA-CN")) {
        return SPECTRA_CN;
    }
    else if (upperMode == string("SPECTRA-HIST")) {
        return SPECTRA_HIST;
    }
    else if (upperMode == string("SPECTRA-MX")) {
        return SPECTRA_MX;
    }
    else if (upperMode == string("BLOB")) {
        return BLOB;
    }
    else {
        BOOST_THROW_EXCEPTION(KatPlotException() << KatPlotErrorInfo(string(
                    "Could not recognise mode string: ") + mode));
    }
}

path kat::Plot::getPythonScript(const PlotMode mode) {


    switch (mode) {
        case DENSITY:
            return "kat_plot_density.py";
        case PROFILE:
            return "kat_plot_profile.py";
        case SPECTRA_CN:
            return "kat_plot_spectra_cn.py";
        case SPECTRA_HIST:
            return "kat_plot_spectra_hist.py";
        case SPECTRA_MX:
            return "kat_plot_spectra_mx.py";
        case BLOB:
            return "kat_plot_blob.py";
        default:
            BOOST_THROW_EXCEPTION(KatPlotException() << KatPlotErrorInfo(string(
                    "Unrecognised KAT PLOT mode")));
    }
}


void kat::Plot::executePythonPlot(const PlotMode mode, vector<string>& args) {

    char* char_args[50];

    for(size_t i = 0; i < args.size(); i++) {
        char_args[i] = strdup(args[i].c_str());
    }

    PyHelper::getInstance().execute(getPythonScript(mode).string(), (int)args.size(), char_args);

    for(size_t i = 0; i < args.size(); i++) {
        free(char_args[i]);
    }
}

void kat::Plot::executePythonPlot(const PlotMode mode, int argc, char *argv[]) {
    PyHelper::getInstance().execute(getPythonScript(mode).string(), argc, argv);
}


int kat::Plot::main(int argc, char *argv[]) {

    string modeStr;
    vector<string> others;
    bool verbose;
    bool help;

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);


    // Declare the supported options.
    po::options_description generic_options(Plot::helpMessage(), w.ws_col);
    generic_options.add_options()
            ("verbose,v", po::bool_switch(&verbose)->default_value(false), "Print extra information")
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
            ("mode", po::value<string>(&modeStr), "KAT PLOT mode.")
            ("others", po::value< std::vector<string> >(&others), "Other options.")
            ;

    // Positional option for the input bam file
    po::positional_options_description p;
    p.add("mode", 1);
    p.add("others", 100);

    // Combine non-positional options
    po::options_description cmdline_options;
    cmdline_options.add(generic_options).add(hidden_options);

    // Parse command line
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).allow_unregistered().run(), vm);
    po::notify(vm);

    // Output help information the exit if requested
    // Output help information the exit if requested
    if (argc == 1 || (argc == 2 && verbose) || (argc == 2 && help) || (argc == 3 && verbose && help)) {
            cout << generic_options << endl;
        return 1;
    }

    // Get the plot mode
    Plot::PlotMode mode = Plot::parseMode(modeStr);

    // Remove the plot mode from the command line arguments
    const int modeArgC = argc-1;
    char** modeArgV = argv+1;

    // Execute via appropriate method (or error)
#ifdef HAVE_PYTHON
    executePythonPlot(mode, modeArgC, modeArgV);
#else
    BOOST_THROW_EXCEPTION(KatPlotException() << KatPlotErrorInfo(string(
                "No suitable plotting environment detected.  We recommend you install anaconda3 to get a python plotting environment setup.  Otherwise install gnuplot.")));
#endif

    return 0;
}
