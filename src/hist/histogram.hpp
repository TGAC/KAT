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

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <thread>
using std::shared_ptr;
using std::make_shared;
using std::thread;

#include <matrix/matrix_metadata_extractor.hpp>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish_helper.hpp>

#include "hist_args.hpp"

namespace kat {

    class Histogram : public jellyfish::thread_exec {
    private:
        // Arguments from user
        HistArgs args;

        // Jellyfish mapped file hash vars
        shared_ptr<JellyfishHelper> jfh;
        
        // Internal vars
        uint64_t base, ceil, inc, nb_buckets, nb_slices;
        uint64_t *data;
        uint64_t slice_id;

    public:

        Histogram(HistArgs& _args) : args(_args) {
            
            // Some validation first
            if (args.high < args.low)
                throw "High count value must be >= to low count value";

            if (args.verbose)
                cerr << "Setting up histo tool..." << endl;

            // Setup handles to load hashes
            jfh = make_shared<JellyfishHelper>(args.db_path, AccessMethod::SEQUENTIAL);

            // Calculate other vars required for this run
            base = args.calcBase();
            ceil = args.calcCeil();
            nb_buckets = ceil + 1 - base;
            nb_slices = args.threads * 100;

            data = new uint64_t[args.threads * nb_buckets];
            memset(data, '\0', args.threads * nb_buckets * sizeof (uint64_t));

            if (args.verbose)
                cerr << "Histo tool setup successfully." << endl;
        }

        virtual ~Histogram() {
            if (data)
                delete [] data;
        }

        void execute() {
            if (args.verbose) {
                cerr << "Loading hash..." << endl;
            }

            std::ostream* out_stream = args.verbose ? &cerr : (std::ostream*)0;

            if (args.verbose)
                cerr << endl
                    << "Hash loaded successfully." << endl
                    << "Starting threads...";

            startAndJoinThreads(args.threads);

            if (args.verbose)
                cerr << "done." << endl;
        }
        
        void print(std::ostream &out) {
            // Output header
            out << mme::KEY_TITLE << "K-mer spectra for: " << args.db_path << endl;
            out << mme::KEY_X_LABEL << "K" << jfh->getKeyLen() << " multiplicity: " << args.db_path << endl;
            out << mme::KEY_Y_LABEL << "Number of distinct K" << jfh->getKeyLen() << " mers" << endl;
            out << mme::MX_META_END << endl;

            uint64_t col = base;
            for (uint64_t i = 0; i < nb_buckets; i++, col += inc) {
                uint64_t count = 0;

                for (uint16_t j = 0; j < args.threads; j++)
                    count += data[j * nb_buckets + i];

                out << col << " " << count << "\n";
            }
        }

    protected:
        
        void startAndJoinThreads(const uint16_t nbThreads) {
            
            thread t[nbThreads];
            
            for(int i = 0; i < nbThreads; i++) {
                t[i] = thread(&Histogram::start, this, i);
            }
            
            for(int i = 0; i < nbThreads; i++){
                t[i].join();
            }
        }
         
        void start(int th_id) {
            uint64_t *hist = &data[th_id * nb_buckets];

            for (size_t i = slice_id++; i < nb_slices; i = slice_id++) {
                lha::region_iterator it = jfh->getSlice(i, nb_slices);
                while (it.next()) {
                    if (it.val() < base)
                        ++hist[0];
                    else if (it.val() > ceil)
                        ++hist[nb_buckets - 1];
                    else
                        ++hist[(it.val() - base) / inc];
                }
            }
        }

        
    };
}
