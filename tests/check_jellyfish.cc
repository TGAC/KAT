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

namespace kat {

TEST(KAT_JELLYFISH, TEST_HEADER) {
    
    file_header header = *(JellyfishHelper::loadHashHeader("data/ecoli.header.jf27"));
    unsigned int klen = header.key_len();
    unsigned int vlen = header.val_len();    
    unsigned int clen = header.counter_len();
    unsigned int mr = header.max_reprobe();
    unsigned int offset = header.offset();
    unsigned int hashes = header.nb_hashes();
    size_t s = header.size();
    string format = header.format();
       
    EXPECT_EQ( klen, 54 );
    EXPECT_EQ( vlen, 7 );
    EXPECT_EQ( clen, 4 );
    EXPECT_EQ( mr, 126 );
    EXPECT_EQ( offset, 1368 );
    EXPECT_EQ( hashes, 0 );
    EXPECT_EQ( s, 131072 );
    EXPECT_EQ( format, "binary/sorted" );
    
    
}

TEST(KAT_JELLYFISH, TEST_QUERY) {
    
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
   
    EXPECT_EQ( countStart, 3 );
    EXPECT_EQ( countEarly, 1 );
    EXPECT_EQ( countMiddle, 1 );
    EXPECT_EQ( countEnd, 1 ); 
    
    EXPECT_EQ( countStartCan, 3 );
    EXPECT_EQ( countEarlyCan, 1 );
    EXPECT_EQ( countMiddleCan, 0 );
    EXPECT_EQ( countEndCan, 0 );  
}

TEST(KAT_JELLYFISH, TEST_SLICE) {
    
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
    
    EXPECT_EQ( nb_records, 1889 );
}

TEST(KAT_JELLYFISH, TEST_COUNT) {
    
    cout << "Start" << endl;
    HashCounter hc(10000000, 27 * 2, 7, 1);
    
    cout << "HC created" << endl;
    LargeHashArrayPtr hash = JellyfishHelper::countSeqFile("data/EcoliK12.fasta", hc, true, 1);
    
    cout << "Counted" << endl;
    
    mer_dna kStart("AGCTTTTCATTCTGACTGCAACGGGCA");
    
    uint64_t count = JellyfishHelper::getCount(hash, kStart, false);
    
    cout << "Kmer found" << endl;
    
    EXPECT_EQ( count, 1 );
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

TEST(KAT_JELLYFISH, TEST_DUMP) {
    
    HashLoader hlBefore;
    LargeHashArrayPtr hashBefore = hlBefore.loadHash("data/ecoli.header.jf27", false);
    file_header header = hlBefore.getHeader();
    
    mer_dna kStart("AGCTTTTCATTCTGACTGCAACGGGCA");
    uint64_t countBefore = JellyfishHelper::getCount(hashBefore, kStart, false);
    
    EXPECT_EQ( countBefore, 3 );
    
    JellyfishHelper::dumpHash(hashBefore, header, 2, "temp_dump.jf");
    
    EXPECT_EQ( boost::filesystem::exists("temp_dump.jf"), true );
    
    HashLoader hlAfter;
    LargeHashArrayPtr hashAfter = hlAfter.loadHash("temp_dump.jf", false);
    uint64_t countAfter = JellyfishHelper::getCount(hashAfter, kStart, false);
    
    EXPECT_EQ( countAfter, 3 );
    
    remove("temp_dump.jf");
}


}