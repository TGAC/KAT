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

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE KAT_JELLYFISH
#define BOOST_TEST_LOG_LEVEL all

#include <iostream>
using std::cout;
using std::endl;

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>
#include <boost/filesystem.hpp>
using boost::filesystem::remove;

#include <../src/inc/spectra_helper.hpp>
using kat::SpectraHelper;

BOOST_AUTO_TEST_SUITE(KAT_SPECTRA_HELPER)

BOOST_AUTO_TEST_CASE(TEST_LOAD_HIST) {
    
    vector<Pos> hist;
    SpectraHelper::loadHist("data/kat.hist", hist);
    
    Pos p1(1, 54015667);
    Pos p10(10, 18649);
    Pos p10000(10000, 0);
    Pos last(10001, 358);
    
    BOOST_CHECK_EQUAL( hist.size(), 10001 );
    BOOST_CHECK_EQUAL( p1.first, hist[0].first );
    BOOST_CHECK_EQUAL( p1.second, hist[0].second );
    BOOST_CHECK_EQUAL( p10.second, hist[9].second );
    BOOST_CHECK_EQUAL( last.first, hist[10000].first );
    BOOST_CHECK_EQUAL( last.second, hist[10000].second );
}

BOOST_AUTO_TEST_CASE(TEST_PEAK) {
    
    vector<Pos> hist;
    SpectraHelper::loadHist("data/kat.hist", hist);
    Pos p = SpectraHelper::findPeak(hist);
    
    BOOST_CHECK_EQUAL( 229, p.first );
    BOOST_CHECK_EQUAL( 9762, p.second );
}

BOOST_AUTO_TEST_SUITE_END()