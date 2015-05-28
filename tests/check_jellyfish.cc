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
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>

#include <../src/jellyfish_helper.hpp>
using kat::JellyfishHelper;
using kat::AccessMethod;

BOOST_AUTO_TEST_SUITE( KAT_JELLYFISH )

BOOST_AUTO_TEST_CASE( TEST_JFCMD ) {
    
    string cmd = JellyfishHelper::jellyfishCountCmd("input.fa", "output.jf27", 27, 100000, 4, true);
    
    const string expected = "jellyfish count -C -m 27 -s 100000 -t 4 -o output.jf27 input.fa ";
    
    //cout << cmd << endl;
    //cout << expected << endl;
    
    BOOST_CHECK_EQUAL( cmd, expected );
}


BOOST_AUTO_TEST_CASE( TEST_READER ) {
    
    JellyfishHelper jfh("data/ecoli.header.jf27", AccessMethod::SEQUENTIAL);
    
    shared_ptr<binary_reader> reader = jfh.getReader();
    
    uint32_t count = 0;
    while (reader->next()) {
        count++;
        //cout << reader->val() << endl;
    }
    
    BOOST_CHECK_EQUAL( count, 1889 );
}

BOOST_AUTO_TEST_CASE( TEST_HEADER ) {
    
    JellyfishHelper jfh("data/ecoli.header.jf27", AccessMethod::SEQUENTIAL);
    
    file_header header = jfh.getHeader();
    unsigned int klen = header.key_len();
    unsigned int vlen = header.val_len();    
    unsigned int clen = header.counter_len();
    unsigned int mr = header.max_reprobe();
    unsigned int offset = header.offset();
    unsigned int hashes = header.nb_hashes();
    size_t s = header.size();
    
    cout << "Key length: " << klen << endl;
    cout << "Value length: " << vlen << endl;
    cout << "Counter length: " << clen << endl;
    cout << "Max reprobe: " << mr << endl;
    cout << "Offset: " << offset << endl;
    cout << "Hashes: " << hashes << endl;
    cout << "Size: " << s << endl;
    
    BOOST_CHECK_EQUAL( klen, 54 );
    BOOST_CHECK_EQUAL( vlen, 7 );
    BOOST_CHECK_EQUAL( clen, 4 );
    BOOST_CHECK_EQUAL( mr, 126 );
    BOOST_CHECK_EQUAL( offset, 1368 );
    BOOST_CHECK_EQUAL( hashes, 0 );
    BOOST_CHECK_EQUAL( s, 131072 );
    
    
}

BOOST_AUTO_TEST_CASE( TEST_QUERY ) {
    
    JellyfishHelper jfh("data/ecoli.header.jf27", AccessMethod::SEQUENTIAL);
    unsigned int keylen = jfh.getKeyLen();
    
    mer_dna kStart("AGCTTTTCATTCTGACTGCAACGGGCA");
    mer_dna kEarly("GCATAGCGCACAGACAGATAAAAATTA");
    mer_dna kMiddle("AATGAAAAAGGCGAACTGGTGGTGCTT");
    mer_dna kEnd("CTCACCAATGTACATGGCCTTAATCTG");
    
    uint64_t countStart = jfh.getCount(kStart);
    uint64_t countEarly = jfh.getCount(kEarly);
    uint64_t countMiddle = jfh.getCount(kMiddle);
    uint64_t countEnd = jfh.getCount(kEnd);
    
    BOOST_CHECK_EQUAL( countStart, 3 );
    BOOST_CHECK_EQUAL( countEarly, 1 );
    BOOST_CHECK_EQUAL( countMiddle, 1 );
    BOOST_CHECK_EQUAL( countEnd, 1 );
    BOOST_CHECK_EQUAL( keylen, 54 );
}

BOOST_AUTO_TEST_CASE( TEST_SLICE ) {
    
    JellyfishHelper jfh("data/ecoli.header.jf27", AccessMethod::SEQUENTIAL);
    lha::region_iterator r1 = jfh.getSlice(1,2);
    lha::region_iterator r2 = jfh.getSlice(2,2);
    
    uint32_t r1Count = 0;
    while (r1.next()) {
        r1Count++;        
    }
    
    uint32_t r2Count = 0;
    /*while (r2.next()) {
        r2Count++;        
    }*/
   
    BOOST_CHECK_EQUAL( r1Count, 200 );
    BOOST_CHECK_EQUAL( r2Count, 200 );
}

BOOST_AUTO_TEST_SUITE_END()