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

#include <../src/comp/comp_args.hpp>
#include <../src/comp/comp.hpp>
#define BOOST_TEST_MODULE COMP
#include <boost/test/included/unit_test.hpp>

using kat::CompArgs;
using kat::Comp;

BOOST_AUTO_TEST_CASE( COMP )
{
    CompArgs args;

    args.db1_path = "data/comp_test_1.jf31_0";
    args.db2_path = "data/comp_test_2.jf31_0";
    args.output_prefix = "temp/comp_test";
    args.threads = 1;

    Comp<hash_query_t> comp(&args);

    comp.do_it();

    SparseMatrix<uint64_t>* results = comp.getMainMatrix();

    BOOST_CHECK( true );
}
