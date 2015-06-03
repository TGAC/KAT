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

#pragma once

#include <string>
#include <iostream>
#include <memory>
#include <vector>
using std::ifstream;
using std::ostream;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::shared_ptr;
using std::make_shared;
using std::vector;

#include <boost/exception/all.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>
using boost::filesystem::path;
using boost::lexical_cast;
using boost::timer::auto_cpu_timer;

#include <jellyfish/err.hpp>
#include <jellyfish/file_header.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/large_hash_array.hpp>
#include <jellyfish/large_hash_iterator.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/storage.hpp>
#include <jellyfish/stream_manager.hpp>
using jellyfish::mer_dna;
using jellyfish::file_header;
using jellyfish::mapped_file;

typedef boost::error_info<struct JellyfishError,string> JellyfishErrorInfo;
struct JellyfishException: virtual boost::exception, virtual std::exception { };

typedef shared_ptr<file_header> HashHeaderPtr;
typedef shared_ptr<binary_reader> HashReaderPtr;
typedef jellyfish::stream_manager<vector<const char*>::const_iterator> StreamManager;
typedef jellyfish::mer_overlap_sequence_parser<StreamManager> SequenceParser;
typedef jellyfish::mer_iterator<SequenceParser, mer_dna> MerIterator;
typedef jellyfish::cooperative::hash_counter<mer_dna> HashCounter;
typedef shared_ptr<HashCounter> HashCounterPtr;
typedef HashCounter::array LargeHashArray;
typedef LargeHashArray* LargeHashArrayPtr;

namespace kat {

    const uint64_t DEFAULT_HASH_SIZE = 10000000000;
    const uint16_t DEFAULT_MER_LEN = 27;
    
    class HashLoader {
        
    private:
        
        LargeHashArrayPtr hash;
        bool canonical;
        uint16_t merLen;
        file_header header;
        
    public:
        
        HashLoader() {
            hash = nullptr;
            canonical = false;
            merLen = 0;
        }
        
        virtual ~HashLoader() {
            
            if (hash != nullptr)
                delete hash;
        }
        
        /**
         * Loads an entire jellyfish hash into memory.  Results stored at the "hash" pointer variable,
         * which is also returned from this function.
         * @param jfHashPath Path to the jellyfish hash file
         * @param verbose Output additional information to cout
         * @return The hash array
         */
        LargeHashArrayPtr loadHash(const path& jfHashPath, bool verbose);
        
        LargeHashArrayPtr getHash() { return hash; }
        
        bool getCanonical() { return header.canonical(); }
        
        uint16_t getMerLen() { return merLen; }
        
        file_header getHeader() { return header; }
    };
    
    
    class JellyfishHelper {

        
    public:

        

        /**
         * This can be set be the client as a record of the path to the jellyfish
         * executable
         */
        static path jellyfishExe;
        
        static uint64_t getCount(LargeHashArrayPtr hash, const mer_dna& kmer, bool canonical);
        
        /**
        * Simple count routine
        * @param ary Hash array which contains the counted kmers
        * @param parser The parser that handles the input stream and chunking
        * @param canonical whether or not the kmers should be treated as canonical or not
        */
        static void countSlice(HashCounter& ary, SequenceParser& parser, bool canonical);

        /**
         * Counts kmers in the given sequence file (Fasta or Fastq) returning
         * a hash array of those kmers
         * @param seqFile Sequence file to count
         * @return The hash array counter
         */
        static LargeHashArrayPtr countSeqFile(const path& p, HashCounter& hashCounter, bool canonical, uint16_t threads) {
            vector<path> paths;
            paths.push_back(p);
            return countSeqFile(paths, hashCounter, canonical, threads);
        }

        /**
         * Counts kmers in the given sequence file (Fasta or Fastq) returning
         * a hash array of those kmers
         * @param seqFile Sequence file to count
         * @return The hash array counter
         */
        static LargeHashArrayPtr countSeqFile(const vector<path>& seqFiles, HashCounter& hashCounter, bool canonical, uint16_t threads);

        
        
        static void dumpHash(LargeHashArrayPtr ary, file_header& header,  uint16_t threads, path outputFile);
        
        /**
        * Extracts the jellyfish hash file header
        * @param jfHashPath Path to the jellyfish hash file
        * @return The hash header
        */
        static file_header loadHashHeader(const path& jfHashPath);
        
        /**
        * Output header in human-readable format
        * @param header Jellyfish hash header
        * @param out Output stream
        */
        static void printHeader(const file_header& header, ostream& out);
        
 
        /**
         * Returns whether or not the specified file path looks like it belongs to
         * a sequence file (either FastA or FastQ).  Gzipped sequence files are
         * also supported.
         * @param filename Path to file
         * @return Whether or not the file is a seqeunce file
         */
        static bool isSequenceFile(const path& filename);
        
        static string createJellyfishCountCmd(const path& input, const path& output, uint16_t merLen, uint64_t hashSize, uint16_t threads, bool canonical) {
            
            vector<path> paths;
            paths.push_back(input);
            return createJellyfishCountCmd(paths, output, merLen, hashSize, threads, canonical);
        }
        
        static string createJellyfishCountCmd(const vector<path>& input, const path& output, uint16_t merLen, uint64_t hashSize, uint16_t threads, bool canonical);
        
        static path executeJellyfishCount(const path& input, const path& output, uint16_t merLen, uint64_t hashSize, uint16_t threads, bool canonical, bool verbose) {
            
            vector<path> paths;
            paths.push_back(input);
            
            return executeJellyfishCount(paths, output, merLen, hashSize, threads, canonical, verbose);
        }
        
        static path executeJellyfishCount(const vector<path>& inputs, const path& output, uint16_t merLen, uint64_t hashSize, uint16_t threads, bool canonical, bool verbose);
        
    protected:

        
        
        static void executeJellyfishCount(const string& cmd, bool verbose);
    };
}
