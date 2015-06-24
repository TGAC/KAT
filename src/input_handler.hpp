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

#include <memory>
using std::shared_ptr;

#include <jellyfish_helper.hpp>
using kat::JellyfishHelper;

namespace kat {
    
    class InputHandler {
        
    public:
        
        enum InputMode {
            LOAD,
            COUNT
        };
        
        uint16_t index;
        vector<path> input;
        InputMode mode = InputMode::COUNT;
        bool canonical = false;
        uint64_t hashSize = DEFAULT_HASH_SIZE;
        bool dumpHash = false;
        bool disableHashGrow = false;
        HashCounterPtr hashCounter = nullptr;
        shared_ptr<HashLoader> hashLoader = nullptr;
        LargeHashArrayPtr hash = nullptr;
        shared_ptr<file_header> header;         // Only applicable if loaded

        void setSingleInput(path p) { input.clear(); input.push_back(p); }
        path getSingleInput() { return input[0]; }
        string pathString();
        void validateInput();   // Throws if input is not present.  Sets input mode.
        void loadHeader();
        void validateMerLen(uint16_t merLen);   // Throws if incorrect merlen
        void count(uint16_t merLen, uint16_t threads);   // Uses the jellyfish library to count kmers in the input
        void loadHash() { loadHash(false); }
        void loadHash(bool verbose);
        void dump(const path& outputPath, uint16_t threads, bool verbose);
    };
    
}