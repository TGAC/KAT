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
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

#include <boost/algorithm/string.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/predicate.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

#include "filter_kmer.hpp"
#include "filter.hpp"
using kat::filter::FilterKmer;
using kat::Filter;


Filter::FilterMode kat::Filter::parseMode(const string& mode) {

    string upperMode = boost::to_upper_copy(mode);

    if (upperMode == string("KMER")) {
        return KMER;                
    }
    else if (upperMode == string("SEQUENCE")) {
        return SEQUENCE;
    }
    else {
        BOOST_THROW_EXCEPTION(KatFilterException() << KatFilterErrorInfo(string(
                    "Could not recognise mode string: ") + mode));
    }
}
       
int kat::Filter::main(int argc, char *argv[]) {

    string modeStr;
    vector<string> others;
    bool help;

    // Declare the supported options.
    po::options_description generic_options(Filter::helpMessage(), 100);
    generic_options.add_options()
            ("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
            ("mode", po::value<string>(&modeStr), "KAT FILTER mode.")
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
    if (argc == 1 || argc == 2 && help) {
        cout << generic_options << endl;
        return 1;
    }

    Filter::FilterMode mode = Filter::parseMode(modeStr);

    const int modeArgC = argc-1;
    char** modeArgV = argv+1;

    switch (mode) {
    case KMER:
        FilterKmer::main(modeArgC, modeArgV);
        break;
    case SEQUENCE:
        //FilterSequence::main(modeArgC, modeArgV);
        break;
    default:
        BOOST_THROW_EXCEPTION(KatFilterException() << KatFilterErrorInfo(string(
                "Unrecognised KAT FILTER mode: ") + modeStr));
    }

    return 0;
}

