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

#include <matrix/matrix_metadata_extractor.hpp>

#include <jellyfish/hash.hpp>
#include <jellyfish/counter.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/jellyfish_helper.hpp>

#include "hist_args.hpp"

namespace kat
{

    template<typename hash_t>
    class Histogram : public thread_exec
    {
    private:
        // Arguments from user
        HistArgs *args;

        // Jellyfish mapped file hash vars
        JellyfishHelper         *jfh;
        hash_t                  *hash;		// Shortcut to Jellyfish hash (don't try to delete)

        // Internal vars
        uint64_t  base, ceil, inc, nb_buckets, nb_slices;
        uint64_t        *data;
        counter_t       slice_id;

    public:
        Histogram(HistArgs* _args) : args(_args)
        {
            // Some validation first
            if(args->high < args->low)
                throw "High count value must be >= to low count value";

            if (args->verbose)
                cerr << "Setting up histo tool..." << endl;

            // Setup handles to load hashes
            jfh = new JellyfishHelper(args->db_path);

            // Ensure hash is null at this stage (we'll load them later)
            hash = NULL;

            // Calculate other vars required for this run
            base = args->calcBase();
            ceil = args->calcCeil();
            nb_buckets = ceil + 1 - base;
            nb_slices = args->threads * 100;

            data = new uint64_t[args->threads * nb_buckets];
            memset(data, '\0', args->threads * nb_buckets * sizeof(uint64_t));

            if (args->verbose)
                cerr << "Histo tool setup successfully." << endl;
        }

        ~Histogram()
        {
            if(data)
                delete [] data;

            if (jfh)
                delete jfh;

            jfh = NULL;
        }

        void do_it()
        {
            if (args->verbose)
            {
                cerr << "Loading hash..." << endl;
            }

            std::ostream* out_stream = args->verbose ? &cerr : (std::ostream*)0;

            // Load the hashes
            hash = jfh->loadHash(true, out_stream);

            // Whether to treat this hash as double stranded or not.
            // Ideally it would be nice to determine this directly from the hash but I'm
            // not sure how to do that at the moment... it might not be possible
            if (args->both_strands)
            {
                hash->set_canonical(true);
            }

            if (args->verbose)
                cerr << endl
                     << "Hash loaded successfully." << endl
                     << "Starting threads...";

            exec_join(args->threads);

            if (args->verbose)
                cerr << "done." << endl;
        }

        void start(int th_id)
        {
            uint64_t *hist = &data[th_id * nb_buckets];

            for(size_t i = slice_id++; i < nb_slices; i = slice_id++)
            {
                typename hash_t::iterator it = hash->iterator_slice(i, nb_slices);
                while(it.next())
                {
                    if(it.get_val() < base)
                        ++hist[0];
                    else if(it.get_val() > ceil)
                        ++hist[nb_buckets - 1];
                    else
                        ++hist[(it.get_val() - base) / inc];
                }
            }
        }


        void print(std::ostream &out, bool full)
        {
            // Output header
            out << mme::KEY_TITLE << "K-mer spectra for: " << args->db_path << endl;
            out << mme::KEY_X_LABEL << "K" << hash->get_mer_len() << " multiplicity: " << args->db_path << endl;
            out << mme::KEY_Y_LABEL << "Number of distinct K" << hash->get_mer_len() << " mers" << endl;
            out << mme::MX_META_END << endl;

            uint64_t col = base;
            for(uint64_t i = 0; i < nb_buckets; i++, col += inc)
            {
                uint64_t count = 0;

                for(uint_t j = 0; j < args->threads; j++)
                    count += data[j * nb_buckets + i];

                out << col << " " << count << "\n";
            }
        }
    };
}
