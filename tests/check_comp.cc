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


#include <comp/comp_args.hpp>
#include <comp/comp.hpp>
#define BOOST_TEST_MODULE Comp
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE( comp )
{
    CompArgs args;

    args.db1_path = "resources/comp_test_1.jf17_0";
    args.db2_path = "resources/comp_test_2.jf17_0";
    args.both_strands = true;

    Comp comp = Comp(&args);

    comp.do_it();

    SparseMatrix<uint64_t>* results = comp.getMainMatrix();

    BOOST_CHECK( true );
}
