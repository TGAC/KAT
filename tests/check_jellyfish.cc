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
#include <boost/filesystem.hpp>
using boost::filesystem::remove;

#include <chrono>
using std::chrono::system_clock;
using std::chrono::duration;
using std::chrono::duration_cast;
template<typename DtnType>
inline double as_seconds(DtnType dtn) { return duration_cast<duration<double>>(dtn).count(); }

#include <../src/jellyfish_helper.hpp>
using kat::JellyfishHelper;
using kat::HashLoader;

BOOST_AUTO_TEST_SUITE(KAT_JELLYFISH)


BOOST_AUTO_TEST_CASE(TEST_HEADER) {
    
    file_header header = *(JellyfishHelper::loadHashHeader("data/ecoli.header.jf27"));
    unsigned int klen = header.key_len();
    unsigned int vlen = header.val_len();    
    unsigned int clen = header.counter_len();
    unsigned int mr = header.max_reprobe();
    unsigned int offset = header.offset();
    unsigned int hashes = header.nb_hashes();
    size_t s = header.size();
    string format = header.format();
       
    BOOST_CHECK_EQUAL( klen, 54 );
    BOOST_CHECK_EQUAL( vlen, 7 );
    BOOST_CHECK_EQUAL( clen, 4 );
    BOOST_CHECK_EQUAL( mr, 126 );
    BOOST_CHECK_EQUAL( offset, 1368 );
    BOOST_CHECK_EQUAL( hashes, 0 );
    BOOST_CHECK_EQUAL( s, 131072 );
    BOOST_CHECK_EQUAL( format, "binary/sorted" );
    
    
}

BOOST_AUTO_TEST_CASE(TEST_QUERY) {
    
    HashLoader hl;
    LargeHashArrayPtr hash = hl.loadHash("data/ecoli.header.jf27", false);
    
    mer_dna kStart("AGCTTTTCATTCTGACTGCAACGGGCA");
    mer_dna kEarly("GCATAGCGCACAGACAGATAAAAATTA");
    mer_dna kMiddle("AATGAAAAAGGCGAACTGGTGGTGCTT");
    mer_dna kEnd("CTCACCAATGTACATGGCCTTAATCTG");
    
    uint64_t countStart = JellyfishHelper::getCount(hash, kStart, false);
    uint64_t countEarly = JellyfishHelper::getCount(hash, kEarly, false);
    uint64_t countMiddle = JellyfishHelper::getCount(hash, kMiddle, false);
    uint64_t countEnd = JellyfishHelper::getCount(hash, kEnd, false);
    
    uint64_t countStartCan = JellyfishHelper::getCount(hash, kStart, true);
    uint64_t countEarlyCan = JellyfishHelper::getCount(hash, kEarly, true);
    uint64_t countMiddleCan = JellyfishHelper::getCount(hash, kMiddle, true);
    uint64_t countEndCan = JellyfishHelper::getCount(hash, kEnd, true);
   
    BOOST_CHECK_EQUAL( countStart, 3 );
    BOOST_CHECK_EQUAL( countEarly, 1 );
    BOOST_CHECK_EQUAL( countMiddle, 1 );
    BOOST_CHECK_EQUAL( countEnd, 1 ); 
    
    BOOST_CHECK_EQUAL( countStartCan, 3 );
    BOOST_CHECK_EQUAL( countEarlyCan, 1 );
    BOOST_CHECK_EQUAL( countMiddleCan, 0 );
    BOOST_CHECK_EQUAL( countEndCan, 0 );  
}

BOOST_AUTO_TEST_CASE(TEST_SLICE) {
    
    HashLoader hl;
    LargeHashArrayPtr hash = hl.loadHash("data/ecoli.header.jf27", false);
    
    LargeHashArray::region_iterator r1 = hash->region_slice(0,2);
    LargeHashArray::region_iterator r2 = hash->region_slice(1,2);
    
    uint32_t r1Count = 0;
    while (r1.next()) {
        r1Count++;        
        //cout << "i1: pos - " << r1.pos() << "; id - " << r1.id() << "; key - " << r1.key() << "; val - " << r1.val() << endl;    
    }
    
    uint32_t r2Count = 0;
    while (r2.next()) {
        r2Count++;        
        //cout << "i2: pos - " << r2.pos() << "; id - " << r2.id() << "; key - " << r2.key() << "; val - " << r2.val() << endl;    
    }
   
    size_t nb_records = r1Count + r2Count;
    
    BOOST_CHECK_EQUAL( nb_records, 1889 );
}

BOOST_AUTO_TEST_CASE(TEST_COUNT) {
    
    cout << "Start" << endl;
    HashCounter hc(10000000, 27 * 2, 7, 1);
    
    cout << "HC created" << endl;
    LargeHashArrayPtr hash = JellyfishHelper::countSeqFile("data/EcoliK12.fasta", hc, true, 1);
    
    cout << "Counted" << endl;
    
    mer_dna kStart("AGCTTTTCATTCTGACTGCAACGGGCA");
    
    uint64_t count = JellyfishHelper::getCount(hash, kStart, false);
    
    cout << "Kmer found" << endl;
    
    BOOST_CHECK_EQUAL( count, 1 );
}

/*BOOST_AUTO_TEST_CASE(TEST_TIME) {
    
    auto before_hash_count_lib = system_clock::now();
    
    HashCounter libCounter(10000000, 27 * 2, 7, 1);
    LargeHashArrayPtr libHash = JellyfishHelper::countSeqFile("data/EcoliK12.fasta", libCounter, true, 1);
    
    auto after_hash_count_lib = system_clock::now();
    
    JellyfishHelper::executeJellyfishCount("data/EcoliK12.fasta", "temp.jf", 27, 10000000, 1, true, false);
    
    auto after_hash_count_sys = system_clock::now();
        
    cout << "Lib time: " << as_seconds(after_hash_count_lib - before_hash_count_lib) << endl
         << "Sys time: " << as_seconds(after_hash_count_sys - after_hash_count_lib) << endl << endl;
    
    BOOST_CHECK( true );
    
    remove("temp.jf");
}*/

BOOST_AUTO_TEST_CASE(TEST_DUMP) {
    
    HashLoader hlBefore;
    LargeHashArrayPtr hashBefore = hlBefore.loadHash("data/ecoli.header.jf27", false);
    file_header header = hlBefore.getHeader();
    
    mer_dna kStart("AGCTTTTCATTCTGACTGCAACGGGCA");
    uint64_t countBefore = JellyfishHelper::getCount(hashBefore, kStart, false);
    
    BOOST_CHECK_EQUAL( countBefore, 3 );
    
    JellyfishHelper::dumpHash(hashBefore, header, 2, "temp_dump.jf");
    
    BOOST_CHECK( boost::filesystem::exists("temp_dump.jf") );
    
    HashLoader hlAfter;
    LargeHashArrayPtr hashAfter = hlAfter.loadHash("temp_dump.jf", false);
    uint64_t countAfter = JellyfishHelper::getCount(hashAfter, kStart, false);
    
    BOOST_CHECK_EQUAL( countAfter, 3 );
    
    remove("temp_dump.jf");
}

BOOST_AUTO_TEST_SUITE_END()