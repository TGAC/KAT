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

using std::ifstream;
using std::ostream;
using std::string;
using std::cerr;
using std::endl;

#include <boost/shared_ptr.hpp>

using boost::shared_ptr;

#include <jellyfish/err.hpp>
#include <jellyfish/file_header.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/jellyfish.hpp>
#include <fstream_default.hpp>

using jellyfish::mer_dna;
using jellyfish::file_header;
using jellyfish::mapped_file;

namespace kat {

    enum AccessMethod {
        SEQUENTIAL,
        RANDOM
    };
    
    class JellyfishHelper {
    private:
        string jfHashPath;
        AccessMethod accessMethod;
        shared_ptr<ifstream> in;
        file_header header;
        shared_ptr<binary_reader> reader;
        shared_ptr<mapped_file> map;
        shared_ptr<binary_query> query;
        ostream* out;

    public:

        JellyfishHelper(string _jfHashPath, AccessMethod _accessMethod) :
        jfHashPath(_jfHashPath), accessMethod(_accessMethod) {

            in = shared_ptr<ifstream>(new ifstream(jfHashPath, std::ios::in | std::ios::binary));
            header = file_header(*in);

            if (!in->good())
                die << "Failed to parse header of file '" << jfHashPath << "'";

            mer_dna::k(header.key_len() / 2);

            // Output jellyfish has details if requested
            if (out) {
                header.write(*out);
            }

            if (header.format() == "bloomcounter") {
                cerr << "KAT does not currently support bloom counted kmer hashes.  Please create a binary hash with jellyfish and use that instead.";
                throw;
            } else if (header.format() == binary_dumper::format) {

                // Create a binary reader for the input file, configured using the header properties
                reader = shared_ptr<binary_reader>(new binary_reader(*in, &header));

                // Create a binary map for the input file
                map = shared_ptr<mapped_file>(new mapped_file(jfHashPath.c_str()));

                if (accessMethod == SEQUENTIAL) {
                    map->sequential();
                }
                else {
                    map->random();
                }
                
                query = shared_ptr<binary_query>(new binary_query(map->base() + header.offset(), 
                        header.key_len(), header.counter_len(), header.matrix(),
                        header.size() - 1, map->length() - header.offset()));

                
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

        uint64_t getCount(mer_dna& kmer) {
            if (header.canonical())
                kmer.canonicalize();
            return (*query)[kmer];
        }
        
        void setOut(ostream* out) {
            this->out = out;
        }

        


    };
}
