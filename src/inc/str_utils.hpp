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

#include <stdint.h>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

using std::string;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::vector;
using std::endl;

namespace kat
{
    // TODO Should write some templates for these methods
    // Note that none of these methods are particularly efficient so don't use them
    // for lots of data.

    static uint16_t strToInt16(string s)
    {
        istringstream str_val(s);
        uint16_t int_val;
        str_val >> int_val;
        return int_val;
    }

    static uint32_t strToInt32(string s)
    {
        istringstream str_val(s);
        uint32_t int_val;
        str_val >> int_val;
        return int_val;
    }

    static uint64_t strToInt64(string s)
    {
        istringstream str_val(s);
        uint64_t int_val;
        str_val >> int_val;
        return int_val;
    }

    static double strToDouble(string s)
    {
        istringstream str_val(s);
        double double_val;
        str_val >> double_val;
        return double_val;
    }


    static vector<uint32_t> &splitUInt32(const string& s, const char delim, vector<uint32_t> &elems)
    {
        stringstream ss(s);
        string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(strToInt32(item));
        }
        return elems;
    }

    static vector<uint32_t> splitUInt32(const string &s, char delim) {
        vector<uint32_t> elems;
        splitUInt32(s, delim, elems);
        return elems;
    }

    static vector<uint64_t> &splitUInt64(const string& s, const char delim, vector<uint64_t> &elems)
    {
        stringstream ss(s);
        string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(strToInt64(item));
        }
        return elems;
    }

    static vector<uint64_t> splitUInt64(const string &s, char delim) {
        vector<uint64_t> elems;
        splitUInt64(s, delim, elems);
        return elems;
    }


    static vector<string> &splitString(const string &s, const char delim, vector<string> &elems) {
        stringstream ss(s);
        string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
        return elems;
    }


    static vector<string> splitString(const string &s, char delim) {
        vector<string> elems;
        splitString(s, delim, elems);
        return elems;
    }


    static string lineBreakString(string &s, const uint16_t line_length, const string line_prefix)
    {
        istringstream ss(s);
        string word;
        ostringstream out_str;

        out_str << line_prefix;

        uint16_t char_count = 0;
        while (std::getline(ss, word, ' ')) {

            if (word.compare("</br>") == 0)
            {
                out_str << endl << endl << line_prefix;
                char_count = 0;
            }
            else
            {
                char_count += word.length();

                if (char_count > line_length)
                {
                    out_str << endl << line_prefix;
                    char_count = word.length();
                }

                out_str << word << " ";
            }
        }

        return out_str.str();
    }
    
    static uint32_t gcCount(const string& seq) {
    
        uint32_t g_or_c = 0;

        for (const auto& c : seq) {
            if (c == 'G' || c == 'g' || c == 'C' || c == 'c')
                g_or_c++;
        }
        
        return g_or_c;
    }
}
