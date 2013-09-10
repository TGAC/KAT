//  ********************************************************************
//  This file is part of KAT - the Kmer Analysis Toolkit.
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

using std::string;
using std::istringstream;
using std::vector;

namespace kat
{
    static uint32_t strToInt(string s)
    {
        istringstream str_val(s);
        uint32_t int_val;
        str_val >> int_val;
        return int_val;
    }

    // This is horribly inefficient, and annoying to use! :(  Fix later
    static void split(const string& txt, vector<uint32_t> &strs, const char ch)
    {
        strs.clear();
        istringstream iss(txt);
        do
        {
            string sub;
            iss >> sub;
            uint32_t intVal = strToInt(sub);
            strs.push_back(intVal);
        } while (iss);
    }
}
