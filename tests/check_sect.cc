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
#define BOOST_TEST_MODULE KAT_SECT
#define BOOST_TEST_LOG_LEVEL all
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>



#include "../src/sect.hpp"
using kat::Sect;

BOOST_AUTO_TEST_SUITE( KAT_SECT )

BOOST_AUTO_TEST_CASE( quick )
{
    Sect sect("data/sect_test.fa.jf31_0", "data/sect_test.fa");
    sect.setOutputPrefix("temp/sect_quick");
    sect.setCanonical(true);
    
    sect.execute();

    BOOST_CHECK( true );

    // Remove any generated files
    remove("temp/sect_quick_counts.cvg");
    remove("temp/sect_quick_contamination.mx");
    remove("temp/sect_quick_stats.csv");
}

BOOST_AUTO_TEST_CASE( length_check )
{
    Sect sect("data/sect_test.fa.jf31_0", "data/sect_length_test.fa");
    sect.setOutputPrefix("temp/sect_length");
    sect.setCanonical(true);

    sect.execute();

    BOOST_CHECK( true );

    // Remove any generated files
    remove("temp/sect_length_counts.cvg");
    remove("temp/sect_length_contamination.mx");
    remove("temp/sect_length_stats.csv");
}

BOOST_AUTO_TEST_SUITE_END()