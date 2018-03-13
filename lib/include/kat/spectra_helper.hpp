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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::getline;
using std::pair;
using std::string;
using std::vector;

#include <boost/lexical_cast.hpp>
#include <boost/exception/exception.hpp>
#include <boost/exception/info.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::path;
using boost::lexical_cast;

typedef pair<uint32_t, uint64_t> Pos;
typedef pair<uint32_t, uint32_t> Coord;

namespace kat {

    typedef boost::error_info<struct SpectraHelperError,string> SpectraHelperErrorInfo;
    struct SpectraHelperException: virtual boost::exception, virtual std::exception { };

    class SpectraHelper {
    public:

        static uint32_t findFirstMin(const vector<Pos>& histo) {

            return findFirstMin(histo, false);
        }

        static uint32_t findFirstMin(const vector<Pos>& histo, bool skipFirst) {

            uint64_t previous = std::numeric_limits<uint64_t>::max();

            for(uint32_t i = skipFirst ? 1 : 0; i < histo.size(); i++) {
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
            return findPeak(histo, true);
        }

        static Pos findPeak(const vector<Pos>& histo, bool findMin) {

            uint64_t previous = std::numeric_limits<uint64_t>::max();

            Pos bestMax(0,0);
            Pos lastMax(0,0);

            for(uint32_t i = findMin ? findFirstMin(histo) : 1; i < histo.size(); i++) {
                if (histo[i].second > previous) {
                    lastMax = histo[i];
                    bestMax = lastMax.second > bestMax.second ? lastMax : bestMax;
                }

                previous = histo[i].second;
            }

            return bestMax;
        }

        static Pos lim97(const vector<Pos>& histo) {

            uint32_t xStart = findFirstMin(histo, true);

            //cout << "minima: " << xStart << endl;

            // No initial minima, so no way of easily working out what the xlim should be
            if (xStart == 0) {
                return Pos();
            }

            uint64_t total = 0;

            for (uint32_t i = xStart; i < histo.size(); i++) {
                total += histo[i].second;
            }

            //cout << "total: " << total << endl;


            uint64_t cumulative = 0;
            for (uint32_t i = xStart; i < histo.size(); i++) {
                cumulative += histo[i].second;
                double proportion = (double)cumulative / (double)total;

                //cout << "proportion: " << proportion << endl;
                if (proportion > 0.97) {
                    return Pos(histo[i].first, cumulative);
                }
            }

            return Pos();
        }

        /*static Coord findPeak(const SparseMatrix<uint64_t>& mx) {

            uint64_t previous = std::numeric_limits<uint64_t>::max();

            Coord bestMax;
            Coord lastMax;

            for(uint32_t i = findFirstMin(mx); i < mx.size(); i++) {
                if (mx[i].second > previous) {
                    lastMax = mx[i];
                    bestMax = lastMax.second > bestMax.second ? lastMax : bestMax;
                }
                else {
                    previous = mx[i].second;
                }
            }

            return bestMax;
        }*/


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
