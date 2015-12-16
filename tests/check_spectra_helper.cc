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

#include <gtest/gtest.h>

#include <iostream>
using std::cout;
using std::endl;

#include <boost/filesystem.hpp>
using boost::filesystem::remove;

#include <../src/inc/spectra_helper.hpp>
using kat::SpectraHelper;

namespace kat {

TEST(KAT_SPECTRA_HELPER, TEST_LOAD_HIST) {
    
    vector<Pos> hist;
    SpectraHelper::loadHist("data/kat.hist", hist);
    
    Pos p1(1, 54015667);
    Pos p10(10, 18649);
    Pos p10000(10000, 0);
    Pos last(10001, 358);
    
    EXPECT_EQ( hist.size(), 10001 );
    EXPECT_EQ( p1.first, hist[0].first );
    EXPECT_EQ( p1.second, hist[0].second );
    EXPECT_EQ( p10.second, hist[9].second );
    EXPECT_EQ( last.first, hist[10000].first );
    EXPECT_EQ( last.second, hist[10000].second );
}

TEST(KAT_SPECTRA_HELPER, TEST_PEAK) {
    
    vector<Pos> hist;
    SpectraHelper::loadHist("data/kat.hist", hist);
    Pos p = SpectraHelper::findPeak(hist);
    
    EXPECT_EQ( 229, p.first );
    EXPECT_EQ( 9762, p.second );
}

}