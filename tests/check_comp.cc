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
#define BOOST_TEST_MODULE KAT_COMP
#define BOOST_TEST_LOG_LEVEL all
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>

#include <thread>
using std::thread;

#include <../src/comp.hpp>
using kat::Comp;
using kat::ThreadedCompCounters;
using kat::CompCounters;

BOOST_AUTO_TEST_SUITE( KAT_COMP )

/*BOOST_AUTO_TEST_CASE( COMP1 )
{
    Comp comp("data/ecoli.header.jf27", "data/ecoli.header.2.jf27");
    comp.setOutputPrefix("temp/comp_test");
    
    comp.execute();

    SM64 results = comp.getMainMatrix();

    BOOST_CHECK( true );
}*/

void addTcc(ThreadedCompCounters& tcc) {
    shared_ptr<CompCounters> cc = make_shared<CompCounters>();
    
    cc->updateHash1Counters(10, 2);
    cc->updateHash1Counters(20, 4);
    cc->updateHash2Counters(0, 3);
    
    tcc.add(cc);
}

BOOST_AUTO_TEST_CASE( THREADED_COUNTERS )
{
    const uint16_t threads = 2;
    
    ThreadedCompCounters tcc("path1", "path2", "path3", 1001);
    
    shared_ptr<CompCounters> cc1 = make_shared<CompCounters>();
    
    cc1->updateHash1Counters(10, 2);
    cc1->updateHash1Counters(20, 4);
    cc1->updateHash2Counters(0, 3);
    
    tcc.add(cc1);
    
    shared_ptr<CompCounters> cc2 = make_shared<CompCounters>();
    
    cc2->updateHash1Counters(10, 2);
    cc2->updateHash1Counters(20, 4);
    cc2->updateHash2Counters(0, 3);
    
    tcc.add(cc2);
    
    tcc.merge();
    
    BOOST_CHECK( tcc.size() == 2 );
    BOOST_CHECK( tcc.getFinalMatrix().hash1_path == path("path1"));
    BOOST_CHECK( tcc.getThreadedMatrixAt(0).hash1_path == path("path1"));
    BOOST_CHECK( tcc.getFinalMatrix().hash1_distinct == 4);
    BOOST_CHECK( tcc.getThreadedMatrixAt(0).hash1_distinct == 2);
    BOOST_CHECK( tcc.getThreadedMatrixAt(1).hash1_distinct == 2);
    BOOST_CHECK( tcc.getFinalMatrix().hash1_total == 60);
    
}

BOOST_AUTO_TEST_SUITE_END()