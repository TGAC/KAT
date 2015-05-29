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

#include <string.h>
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
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/large_hash_array.hpp>
#include <jellyfish/large_hash_iterator.hpp>
#include <jellyfish/storage.hpp>

using jellyfish::mer_dna;
using jellyfish::file_header;
using jellyfish::mapped_file;
typedef jellyfish::large_hash::array<mer_dna> lha;

#include <fstream_default.hpp>

typedef boost::error_info<struct JellyfishError,string> JellyfishErrorInfo;
struct JellyfishException: virtual boost::exception, virtual std::exception { };

namespace kat {

    const uint64_t DEFAULT_HASH_SIZE = 10000000000;
    const uint16_t DEFAULT_MER_LEN = 27;
    
    enum AccessMethod {
        SEQUENTIAL,
        RANDOM
    };
    
    struct HashProperties {
        size_t fileSizeBytes;
        size_t arraySizeBytes;
        size_t blockSize;
        size_t nbBlocks;
    };
    
    class JellyfishHelper {
    private:
        path jfHashPath;
        AccessMethod accessMethod;
        shared_ptr<ifstream> in;
        file_header header;
        shared_ptr<binary_reader> reader;
        shared_ptr<mapped_file> map;
        shared_ptr<binary_query> query;
        shared_ptr<lha> hash; 
        
    public:

        static path jellyfishExe;
        
        JellyfishHelper(const path& _jfHashPath, AccessMethod _accessMethod);
        
        JellyfishHelper(const path& _jfHashPath, AccessMethod _accessMethod, bool verbose);

        virtual ~JellyfishHelper();

        void load();
        
        unsigned int getKeyLen() {
            return header.key_len();
        }

        uint64_t getCount(const mer_dna& kmer);
        
        shared_ptr<binary_reader> getReader() {
            return reader;
        }
        
        file_header getHeader() {
            return header;
        }
        
        lha::eager_iterator getEagerSlice(int index, uint16_t nbSlices) {
            return hash->eager_slice(index, nbSlices);            
        }
        
        lha::region_iterator getRegionSlice(int index, uint16_t nbSlices) {
            return hash->region_slice(index, nbSlices);            
        }


        static bool isSequenceFile(const string& p) {
            return  boost::iequals(p, ".fastq") || 
                    boost::iequals(p, ".fq") || 
                    boost::iequals(p, ".fasta") ||
                    boost::iequals(p, ".fa") || 
                    boost::iequals(p, ".fna");
        }
        
        static string jellyfishCountCmd(const path& input, const path& output, uint16_t merLen, uint64_t hashSize, uint16_t threads, bool canonical) {
            
            vector<path> paths;
            paths.push_back(input);
            return jellyfishCountCmd(paths, output, merLen, hashSize, threads, canonical);
        }
        
        static string jellyfishCountCmd(const vector<path>& input, const path& output, uint16_t merLen, uint64_t hashSize, uint16_t threads, bool canonical);
        
        static path jellyfishCount(const path& input, const path& output, uint16_t merLen, uint64_t hashSize, uint16_t threads, bool canonical, bool verbose) {
            
            vector<path> paths;
            paths.push_back(input);
            
            return jellyfishCount(paths, output, merLen, hashSize, threads, canonical, verbose);
        }
        
        static path jellyfishCount(const vector<path>& inputs, const path& output, uint16_t merLen, uint64_t hashSize, uint16_t threads, bool canonical, bool verbose);
        
    protected:
        
        static void executeJellyfishCount(const string& cmd, bool verbose);
    };
}
