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
using std::cerr;
using std::endl;
using std::shared_ptr;
using std::make_shared;

#include <boost/filesystem/path.hpp>
#include <boost/timer/timer.hpp>
using boost::filesystem::path;
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
                cerr << "KAT does not currently support bloom counted kmer hashes.  Please create a binary hash with jellyfish and use that instead.";
                throw;
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
                cerr << "Processing a text format hash will be painfully slow, so we don't support it.  Please create a binary hash with jellyfish and use that instead.";
                throw;
            } else {
                cerr << "Unknown format '" << header.format() << "'";
                throw;
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
        
        static string jellyfishCountCmd(path input, string output, uint8_t merLen, uint16_t threads) {
            
            vector<path> paths;
            paths.push_back(input);
            return jellyfishCountCmd(paths, output, merLen, threads);
        }
        
        static string jellyfishCountCmd(vector<path>& input, bool canonical, string output, uint8_t merLen, uint64_t hashSize, uint16_t threads) {
            
            string i;
            for (path p : input) {
                i += p.string();                
                i += " ";
            }
            
            return string("jellyfish count ") 
                    + (canonical ? "-C " : "") 
                    + "-m " + merLen + 
                    " -s " + hashSize + 
                    " -t " + threads + 
                    " -o " + output + 
                    " " + i;
        }
        
        static void jellyfishCount(path& input, bool canonical, string output, uint8_t merLen, uint64_t hashSize, uint16_t threads) {
            
            executeJellyfishCount(jellyfishCountCmd(input, canonical, output, merLen, hashSize, threads));
        }
        
        static void jellyfishCount(vector<path>& input, bool canonical, string output, uint8_t merLen, uint64_t hashSize, uint16_t threads) {
            
            executeJellyfishCount(jellyfishCountCmd(input, canonical, output, merLen, hashSize, threads));            
        }
        
    protected:
        
        void executeJellyfishCount(const string& cmd, bool verbose) {
            
            if (verbose) {
                auto_cpu_timer timer(1, "Kmer counting total runtime: %ws\n\n");
            
                cout << "Counting kmers...";
                cout.flush();
            }

            int res = system(cmd);
            
            if (res != 0) {
                BOOST_THROW_EXCEPTION(CompException() << CompErrorInfo(string(
                        "Could not find third jellyfish hash file at: ") + args.input3 + "; please check the path and try again."));
            }

            if (verbose)
                cout << " done" << endl;
        }
    };
}
