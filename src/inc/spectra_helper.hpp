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

#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
using std::ifstream;
using std::istringstream;
using std::getline;
using std::pair;
using std::string;
using std::vector;

#include <boost/lexical_cast.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::path;
using boost::lexical_cast;

typedef pair<uint32_t, uint64_t> Pos;

namespace kat {
    
    typedef boost::error_info<struct SpectraHelperError,string> SpectraHelperErrorInfo;
    struct SpectraHelperException: virtual boost::exception, virtual std::exception { };
    
    class SpectraHelper {
    public:
        
        static uint32_t findFirstMin(const vector<Pos>& histo) {
            
            uint64_t previous = std::numeric_limits<uint64_t>::max();
            
            for(uint32_t i = 0; i < histo.size(); i++) {
                if (histo[i].second <= previous) {
                    previous = histo[i].second;                    
                }
                else {
                    return i;
                }
            }            
            
            return 0;
        }
        
        
        static Pos findPeak(const vector<Pos>& histo) {
            
            uint64_t previous = std::numeric_limits<uint64_t>::max();
            
            Pos bestMax;
            Pos lastMax;
            
            for(uint32_t i = findFirstMin(histo); i < histo.size(); i++) {
                if (histo[i].second > previous) {
                    lastMax = histo[i];
                    bestMax = lastMax.second > bestMax.second ? lastMax : bestMax;
                }
                else {
                    previous = histo[i].second;                    
                }
            }
            
            return bestMax;
        }
        
        
        static void loadHist(const path& histFile, vector<Pos>& histo) {
            
            ifstream in(histFile.c_str());
            string line;
            uint64_t linenb = 0;
            while (getline(in, line)) {
                istringstream iss(line);
                uint32_t bin;
                uint64_t val;
                linenb++;
                
                if (line[0] != '#') {
                    if (!(iss >> bin >> val)) { 
                        BOOST_THROW_EXCEPTION(SpectraHelperException() << SpectraHelperErrorInfo(string(
                            "Encountered unexpected syntax on line ") + lexical_cast<string>(linenb)));
                    }
                    
                    histo.push_back(Pos(bin,val));
                }
            }
        }
    };
}

