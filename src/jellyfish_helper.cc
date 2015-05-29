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

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>
namespace bfs = boost::filesystem;
using bfs::path;
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
using jellyfish::large_hash::reprobe_limit_t;
using jellyfish::Offsets;
using jellyfish::quadratic_reprobes;
typedef jellyfish::large_hash::array_raw<mer_dna> lha;

#include <fstream_default.hpp>

#include "jellyfish_helper.hpp"
using kat::JellyfishHelper;

path kat::JellyfishHelper::jellyfishExe = "";

kat::JellyfishHelper::JellyfishHelper(const path& _jfHashPath, AccessMethod _accessMethod) :
    JellyfishHelper(_jfHashPath, _accessMethod, false) {}

kat::JellyfishHelper::JellyfishHelper(const path& _jfHashPath, AccessMethod _accessMethod, bool verbose) :
    jfHashPath(_jfHashPath), accessMethod(_accessMethod) {

    in = make_shared<ifstream>(jfHashPath.c_str(), std::ios::in | std::ios::binary);
    header = file_header(*in);

    if (!in->good()) {
        BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
            "Failed to parse header of file: ") + jfHashPath.string()));
    }

    mer_dna::k(header.key_len() / 2);

    // Output jellyfish hash details if requested
    if (verbose) {
        header.write(cerr);
        cerr << endl;
    }

    if (header.format() == "bloomcounter") {
        BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
            "KAT does not currently support bloom counted kmer hashes.  Please create a binary hash with jellyfish and use that instead.")));
    } else if (header.format() == binary_dumper::format) {

        // Create a binary reader for the input file, configured using the header properties
        reader = make_shared<binary_reader>(*in, &header);

        // Create a binary map for the input file
        map = make_shared<mapped_file>(jfHashPath.c_str());

        if (accessMethod == SEQUENTIAL) {
            map->sequential();
        }
        else {
            map->random();
        }

        map->load();

        size_t dataSize = map->length() - header.offset();
        cout << "Datasize: " << dataSize << endl;
        
        query = make_shared<binary_query>(
                map->base() + header.offset(), 
                header.key_len(), 
                header.counter_len(), 
                header.matrix(),
                header.size() - 1, 
                dataSize);
        
        size_t lsize = jellyfish::ceilLog2(dataSize);
        size_t size_ = (size_t)1 << lsize;
        size_t size_mask_ = size_ - 1;
        
        reprobe_limit_t reprobe_limit_(header.max_reprobe(), quadratic_reprobes, size_);
        unsigned int raw_key_len_ = (header.key_len() > lsize ? header.key_len() - lsize : 0);
        Offsets<uint64_t> offsets_(raw_key_len_ + jellyfish::bitsize(reprobe_limit_.val() + 1), header.val_len(), reprobe_limit_.val() + 1);
        size_t size_bytes_(jellyfish::div_ceil(size_, (size_t)offsets_.block_len()) * offsets_.block_word_len() * sizeof(uint64_t));
        
        
        hash = make_shared<lha>(
                map->base() + header.offset(),
                size_bytes_,
                dataSize,
                header.key_len(),
                header.val_len(),
                header.max_reprobe(),
                header.matrix());


    } else if (header.format() == text_dumper::format) {
        BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
            "Processing a text format hash will be painfully slow, so we don't support it.  Please create a binary hash with jellyfish and use that instead.")));
    } else {
        BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
            "Unknown format '") + header.format() + "'"));
    }

}

kat::JellyfishHelper::~JellyfishHelper() {

    if (in)
        in->close();
}

uint64_t kat::JellyfishHelper::getCount(const mer_dna& kmer) {
    mer_dna k = kmer;
    if(header.canonical())
        k.canonicalize();            
    return (*query)[k];
}
                
string kat::JellyfishHelper::jellyfishCountCmd(const vector<path>& input, const path& output, uint16_t merLen, uint64_t hashSize, uint16_t threads, bool canonical) {

    string i;
    for (path p : input) {
        i += p.string();                
        i += " ";
    }

    string jellyfishCmd =   kat::JellyfishHelper::jellyfishExe.empty() ? 
                                "jellyfish" : 
                                kat::JellyfishHelper::jellyfishExe.string();
    
    return jellyfishCmd + " count " +
            (canonical ? "-C " : "") +
            "-m " + lexical_cast<string>(merLen) + 
            " -s " + lexical_cast<string>(hashSize) + 
            " -t " + lexical_cast<string>(threads) + 
            " -o " + output.string() + 
            " " + i;
}
    

path kat::JellyfishHelper::jellyfishCount(const vector<path>& inputs, const path& output, uint16_t merLen, uint64_t hashSize, uint16_t threads, bool canonical, bool verbose) {

    if (inputs.empty()) {
        BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
            "No input files provided")));
    }

    if (inputs.size() == 1 && !JellyfishHelper::isSequenceFile(inputs[0].extension().string())) {
        // No need to jellyfish count, input is already a jellyfish hash
        return inputs[0];                
    }
    else {

        for (path p : inputs) {
            if (!JellyfishHelper::isSequenceFile(p.extension().string())) {
                BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
                    "You provided multiple sequence files to generate a kmer hash from, however some of the input files do not have a recognised sequence file extension: \".fa,.fasta,.fq,.fastq,.fna\"")));
            }

            // Check input files exist
            if (!bfs::exists(p) && !bfs::symbolic_link_exists(p)) {
                BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
                    "Could not find input file at: ") + p.string() + "; please check the path and try again."));                        
            }
        }

        if (verbose)
            cout << "Provided one or more sequence files.  Executing jellyfish to count kmers." << endl;

        const string cmd = jellyfishCountCmd(inputs, output, merLen, hashSize, threads, canonical);
        
        if (verbose) {
            cout << "Command issued: " << cmd << endl;
        }
        
        executeJellyfishCount(cmd, verbose);            

        return output;
    }
}
     
void kat::JellyfishHelper::executeJellyfishCount(const string& cmd, bool verbose) {

    auto_cpu_timer timer(1, "Kmer counting total runtime: %ws\n\n");

    if (verbose) {                
        cout << "Counting kmers...";
        cout.flush();
    }

    int res = system(cmd.c_str());

    if (res != 0) {
        BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
                "Problem executing jellyfish count.  Non-0 return code.  Return code: ") + lexical_cast<string>(res)));
    }

    if (verbose)
        cout << " done" << endl;
}
    
