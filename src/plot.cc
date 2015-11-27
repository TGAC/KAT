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

#include <iostream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

#include <Python.h>

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/predicate.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include "plot_density.hpp"
#include "plot_profile.hpp"
#include "plot_spectra_cn.hpp"
#include "plot_spectra_hist.hpp"
#include "plot_spectra_mx.hpp"
#include "plot.hpp"
#include "inc/kat_fs.hpp"
using kat::PlotDensity;
using kat::PlotProfile;
using kat::PlotSpectraCn;
using kat::PlotSpectraHist;
using kat::PlotSpectraMx;
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
    else {
        BOOST_THROW_EXCEPTION(KatPlotException() << KatPlotErrorInfo(string(
                    "Could not recognise mode string: ") + mode));
    }
}

path kat::Plot::getPythonScript(const PlotMode mode) {
    
    
    switch (mode) {
        case DENSITY:
            return "density.py";            
        case PROFILE:
            return "profile.py";
        case SPECTRA_CN:
            return "spectra-cn.py";
        case SPECTRA_HIST:
            return "spectra-hist.py";
        case SPECTRA_MX:
            return "spectra-mx.py";
        default:
            BOOST_THROW_EXCEPTION(KatPlotException() << KatPlotErrorInfo(string(
                    "Unrecognised KAT PLOT mode")));
    }
}

wchar_t* kat::Plot::convertCharToWideChar(const char* c) {
    const size_t cSize = strlen(c)+1;
    wchar_t* wc = new wchar_t[cSize];
    mbstowcs (wc, c, cSize);
    
    return wc;
}

void kat::Plot::executePythonPlot(const PlotMode mode, int argc, char *argv[], const KatFS& fs) {
    
    const path script_name = getPythonScript(mode);
    const path scripts_dir = fs.GetScriptsDir();
    const path full_script_path = path(scripts_dir.string() + "/" + script_name.string());
    
    stringstream ss;
    
    // Create wide char alternatives
    wchar_t* wsn = convertCharToWideChar(script_name.c_str());
    wchar_t* wsp = convertCharToWideChar(full_script_path.c_str());    
    wchar_t* wargv[50]; // Can't use variable length arrays!
    wargv[0] = wsp;
    ss << full_script_path.c_str();
    for(int i = 1; i < argc; i++) {
        wargv[i] = convertCharToWideChar(argv[i]);
        ss << " " << argv[i];
    }
    
    cout << "Effective command line: " << ss.str() << endl << endl;

    std::ifstream script_in(full_script_path.c_str());
    std::string contents((std::istreambuf_iterator<char>(script_in)), std::istreambuf_iterator<char>());
        
    // Run python script
    Py_Initialize();
    Py_SetProgramName(wsp);
    PySys_SetArgv(argc, wargv);
    PyRun_SimpleString(contents.c_str());
    Py_Finalize();
    
    // Cleanup
    delete wsn;
    // No need to free up "wsp" as it is element 0 in the array
    for(int i = 0; i < argc; i++) {
        delete wargv[i];
    }
}

void kat::Plot::executeGnuplotPlot(const PlotMode mode, int argc, char *argv[]) {
    switch (mode) {
    case DENSITY:
        PlotDensity::main(argc, argv);
        break;
    case PROFILE:
        PlotProfile::main(argc, argv);
        break;
    case SPECTRA_CN:
        PlotSpectraCn::main(argc, argv);
        break;
    case SPECTRA_HIST:
        PlotSpectraHist::main(argc, argv);
        break;
    case SPECTRA_MX:
        PlotSpectraMx::main(argc, argv);            
        break;
    default:
        BOOST_THROW_EXCEPTION(KatPlotException() << KatPlotErrorInfo(string(
                "Unrecognised KAT PLOT mode")));
    }
}

       
int kat::Plot::main(int argc, char *argv[], const KatFS& fs) {

    string modeStr;
    vector<string> others;
    bool help;

    // Declare the supported options.
    po::options_description generic_options(Plot::helpMessage(), 100);
    generic_options.add_options()
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
    if (argc == 1 || (argc == 2 && help)) {
        cout << generic_options << endl;
        return 1;
    }

    // Get the plot mode
    Plot::PlotMode mode = Plot::parseMode(modeStr);

    // Remove the plot mode from the command line arguments
    const int modeArgC = argc-1;
    char** modeArgV = argv+1;

    // Execute via appropriate method (or error)
#if HAVE_PYTHON
    executePythonPlot(mode, modeArgC, modeArgV, fs);
#elif HAVE_GNUPLOT
    executeGnuplotPlot(mode, modeArgC, modeArgV);
#else
    BOOST_THROW_EXCEPTION(KatPlotException() << KatPlotErrorInfo(string(
                "No suitable plotting environment detected.  We recommend you install anaconda3 to get a python plotting environment setup.  Otherwise install gnuplot.")));
#endif
    
    return 0;
}

