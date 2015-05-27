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
using std::ifstream;
using std::ostream;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::shared_ptr;
using std::make_shared;

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
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
typedef jellyfish::large_hash::array_raw<mer_dna> lha;

#include <fstream_default.hpp>


namespace kat {

    const uint64_t DEFAULT_HASH_SIZE = 10000000000;
    const uint16_t DEFAULT_MER_LEN = 27;
        
    
    typedef boost::error_info<struct JellyfishError,string> JellyfishErrorInfo;
    struct JellyfishException: virtual boost::exception, virtual std::exception { };
    
    enum AccessMethod {
        SEQUENTIAL,
        RANDOM
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
        ostream* out;

    public:

        
        JellyfishHelper(path _jfHashPath, AccessMethod _accessMethod) :
            jfHashPath(_jfHashPath), accessMethod(_accessMethod) {

            in = make_shared<ifstream>(jfHashPath.c_str(), std::ios::in | std::ios::binary);
            header = file_header(*in);

            if (!in->good())
                cerr << "Failed to parse header of file '" << jfHashPath << "'";
                throw;

            mer_dna::k(header.key_len() / 2);

            // Output jellyfish hash details if requested
            if (out) {
                header.write(*out);
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
                
                query = make_shared<binary_query>(
                        map->base() + header.offset(), 
                        header.key_len(), 
                        header.counter_len(), 
                        header.matrix(1),
                        header.size() - 1, 
                        map->length() - header.offset());

                hash = make_shared<lha>(
                        map->base() + header.offset(),
                        map->end() - map->base(),
                        header.size(),
                        header.key_len(),
                        header.counter_len(),
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

        virtual ~JellyfishHelper() {

            if (in)
                in->close();
        }

        unsigned int getKeyLen() {
            return header.key_len();
        }

        uint64_t getCount(const mer_dna& kmer) {
            return (*query)[kmer.get_canonical()];
        }
        
        void setOut(ostream* out) {
            this->out = out;
        }

        shared_ptr<binary_reader> getReader() {
            return reader;
        }

        lha::region_iterator getSlice(int index, uint16_t nbSlices) {
            return hash->region_slice(index, nbSlices);            
        }

        static bool isSequenceFile(const path& p) {
            return p == "fastq" || p == "fq" || p == "fasta" || p == "fa" || p == "fna";
        }
        
        static string jellyfishCountCmd(const path& input, const path& output, uint8_t merLen, uint64_t hashSize, uint16_t threads, bool canonical) {
            
            vector<path> paths;
            paths.push_back(input);
            return jellyfishCountCmd(paths, output, merLen, hashSize, threads, canonical);
        }
        
        static string jellyfishCountCmd(const vector<path>& input, const path& output, uint8_t merLen, uint64_t hashSize, uint16_t threads, bool canonical) {
            
            string i;
            for (path p : input) {
                i += p.string();                
                i += " ";
            }
            
            return string("jellyfish count ") 
                    + (canonical ? "-C " : "") 
                    + "-m " + lexical_cast<string>(merLen) + 
                    " -s " + lexical_cast<string>(hashSize) + 
                    " -t " + lexical_cast<string>(threads) + 
                    " -o " + output.string() + 
                    " " + i;
        }
        
        static path jellyfishCount(const path& input, const path& output, uint8_t merLen, uint64_t hashSize, uint16_t threads, bool canonical, bool verbose) {
            
            vector<path> paths;
            paths.push_back(input);
            
            return jellyfishCount(paths, output, merLen, hashSize, threads, canonical, verbose);
        }
        
        static path jellyfishCount(const vector<path>& inputs, const path& output, uint8_t merLen, uint64_t hashSize, uint16_t threads, bool canonical, bool verbose) {
            
            if (inputs.empty()) {
                BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
                    "No input files provided")));
            }
            
            if (inputs.size() == 1 && !JellyfishHelper::isSequenceFile(inputs[0].extension())) {
                // No need to jellyfish count, input is already a jellyfish hash
                return inputs[0];                
            }
            else {
                
                for (path p : inputs) {
                    if (!JellyfishHelper::isSequenceFile(p.extension())) {
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
                
                executeJellyfishCount(jellyfishCountCmd(inputs, output, merLen, hashSize, threads, canonical), verbose);            
            
                return output;
            }
        }
        
    protected:
        
        static void executeJellyfishCount(const string& cmd, bool verbose) {
            
            if (verbose) {
                auto_cpu_timer timer(1, "Kmer counting total runtime: %ws\n\n");
            
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
    };
}
