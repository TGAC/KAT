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

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#define BOOST_TEST_MODULE COMP
#endif
#include <boost/test/unit_test.hpp>

#include <../src/comp.hpp>

using kat::Comp;

BOOST_AUTO_TEST_SUITE(COMP)

BOOST_AUTO_TEST_CASE( COMP1 )
{
    Comp comp("data/comp_test_1.jf31_0", "data/comp_test_2.jf31_0");
    comp.setOutputPrefix("temp/comp_test");
    
    comp.execute();

    SM64 results = comp.getMainMatrix();

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()