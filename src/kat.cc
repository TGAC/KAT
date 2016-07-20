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

#include <sys/ioctl.h>
#include <string.h>
#include <exception>
#include <iostream>
using std::exception;
using std::string;
using std::wstring;

#include <boost/exception/all.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <kat/kat_fs.hpp>
#include <kat/str_utils.hpp>
using kat::KatFS;

#include "comp.hpp"
#include "filter.hpp"
#include "gcp.hpp"
#include "histogram.hpp"
#include "plot.hpp"
#include "sect.hpp"
using kat::Comp;
using kat::Filter;
using kat::Gcp;
using kat::Histogram;
using kat::Plot;
using kat::Sect;


typedef boost::error_info<struct KatError,string> KatErrorInfo;
struct KatException: virtual boost::exception, virtual std::exception { };

namespace kat {
// Global variable! :(
KatFS katFileSystem;
}

enum Mode {
    COMP,
    FILTER,
    GCP,
    HIST,
    PLOT,
    SECT
};

Mode parseMode(string mode) {
    
    string upperMode = boost::to_upper_copy(mode);
    
    if (upperMode == string("COMP")) {
        return COMP;                
    }
    else if (upperMode == string("FILTER")) {
        return FILTER;
    }
    else if (upperMode == string("GCP")) {
        return GCP;
    }
    else if (upperMode == string("HIST")) {
        return HIST;
    }
    else if (upperMode == string("PLOT")) {
        return PLOT;
    }    
    else if (upperMode == string("SECT")) {
        return SECT;
    }
    else {
        BOOST_THROW_EXCEPTION(KatException() << KatErrorInfo(string(
                    "Could not recognise mode string: ") + mode));
    }
}

const string helpMessage() {
    return "The K-mer Analysis Toolkit (KAT) contains a number of tools that analyse jellyfish K-mer hashes. \n\n"
                   "The First argument should be the tool/mode you wish to use:\n\n" \
                   "   * sect:   SEquence Coverage estimator Tool.  Estimates the coverage of each sequence in\n" \
                   "             a file using K-mers from another sequence file.\n" \
                   "   * comp:   K-mer comparison tool.  Creates a matrix of shared K-mers between two (or three)\n" \
                   "             sequence files.\n" \
                   "   * gcp:    K-mer GC Processor.  Creates a matrix of the number of K-mers found given a GC\n" \
                   "             count and a K-mer count.\n" \
                   "   * hist:   Create an histogram of k-mer occurrences from a sequence file.  Similar to\n" \
                   "             jellyfish histogram sub command but adds metadata in output for easy plotting,\n" \
                   "             also actually runs multi-threaded.\n" \
                   "   * filter: Filtering tools.  Contains tools for filtering k-mers and sequences based on\n" \
                   "             user-defined GC and coverage limits.\n" \
                   "   * plot:   Plotting tools.  Contains several plotting tools to visualise K-mer and compare\n" \
                   "             distributions.\n\n" \
                   "Options";
}



/**
 * Start point for KAT.  Processes the start of the command line and then
 * delegates the rest of the command line to the child tool.
 */
int main(int argc, char *argv[])
{
    try {
        // KAT args
        string modeStr;
        std::vector<string> others;
        bool verbose;
        bool version;
        bool help;
        
        struct winsize w;
        ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    
        // Declare the supported options.
        po::options_description generic_options(helpMessage(), w.ws_col);
        generic_options.add_options()
                ("verbose,v", po::bool_switch(&verbose)->default_value(false), "Print extra information")
                ("version", po::bool_switch(&version)->default_value(false), "Print version string")
                ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
                ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden_options("Hidden options");
        hidden_options.add_options()
                ("mode", po::value<string>(&modeStr), "KAT mode.")
                ("others", po::value< std::vector<string> >(&others), "Other options.")
                ;

        // Positional options
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

        kat::katFileSystem = KatFS(argv[0]);
        
        if (verbose) {
            cout << kat::katFileSystem << endl << endl;
        }
        // Output help information the exit if requested
        if (argc == 1 || (argc == 2 && verbose) || (argc == 2 && help) || (argc == 3 && verbose && help)) {
            cout << generic_options << endl;
            return 1;
        }

        // Always output version information but exit if version info was all user requested
        
#ifndef PACKAGE_NAME
#define PACKAGE_NAME "KAT"
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "2.X.X"
#endif
        
        if (version) {
            // Output in GNU standard format
            cout << PACKAGE_NAME << " " << PACKAGE_VERSION << endl;
            return 0;
        }
        else {
            // Output for the first line in a normal KAT run
            cout << "Kmer Analysis Toolkit (KAT) V" << PACKAGE_VERSION << endl << endl;
        }
        
        
        Mode mode = parseMode(modeStr);
        
        const int modeArgC = argc-1;
        char** modeArgV = argv+1;
        
        switch(mode) {
            case COMP:
                Comp::main(modeArgC, modeArgV);
                break;
            case FILTER:
                Filter::main(modeArgC, modeArgV);
                break;
            case GCP:
                Gcp::main(modeArgC, modeArgV);
                break;
            case HIST:
                Histogram::main(modeArgC, modeArgV);
                break;
            case PLOT:
                Plot::main(modeArgC, modeArgV);            
                break;
            case SECT:
                Sect::main(modeArgC, modeArgV);            
                break;
            default:
                BOOST_THROW_EXCEPTION(KatException() << KatErrorInfo(string(
                    "Unrecognised KAT mode: ") + modeStr));
        }
                
    } catch(po::error& e) { 
        cerr << "Error: Parsing Command Line: " << e.what() << endl; 
        return 1; 
    } 
    catch (boost::exception &e) { 
        cerr << boost::diagnostic_information(e); 
        return 4;
    } catch (exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 5;
    } catch (const char* msg) {
        cerr << "Error: " << msg << endl;
        return 6;
    } catch (...) {
        cerr << "Error: Exception of unknown type!" << endl;
        return 7;
    }

    return 0;    
}
