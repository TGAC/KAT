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


#include <../src/sect/sect_args.hpp>
#include <../src/sect/sect.hpp>
#include <stdio.h>
#define BOOST_TEST_MODULE SECT
#include <boost/test/included/unit_test.hpp>

using kat::SectArgs;
using kat::Sect;

BOOST_AUTO_TEST_CASE( quick )
{
    SectArgs args;

    args.seq_file = "data/sect_test.fa";
    args.jellyfish_hash = "data/sect_test.fa.jf31_0";
    args.output_prefix = "temp/sect_quick";
    args.both_strands = true;
    args.threads_arg = 1;

    Sect<hash_query_t> sect(&args);

    sect.do_it();

    BOOST_CHECK( true );

    // Remove any generated files
    remove("temp/sect_quick_counts.cvg");
    remove("temp/sect_quick_contamination.mx");
    remove("temp/sect_quick_stats.csv");
}

BOOST_AUTO_TEST_CASE( length_check )
{
    SectArgs args;

    args.seq_file = "data/sect_length_test.fa";
    args.jellyfish_hash = "data/sect_test.fa.jf31_0";
    args.output_prefix = "temp/sect_length";
    args.both_strands = true;
    args.threads_arg = 1;

    Sect<hash_query_t> sect(&args);

    sect.do_it();

    BOOST_CHECK( true );

    // Remove any generated files
    remove("temp/sect_length_counts.cvg");
    remove("temp/sect_length_contamination.mx");
    remove("temp/sect_length_stats.csv");
}
