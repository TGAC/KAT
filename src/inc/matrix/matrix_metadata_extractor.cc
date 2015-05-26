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

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>

#include <matrix/matrix_metadata_extractor.hpp>

using std::ifstream;
using std::string;

void mme::trim(string& str) {
    string::size_type pos = str.find_last_not_of(' ');
    if (pos != string::npos) {
        str.erase(pos + 1);
        pos = str.find_first_not_of(' ');
        if (pos != string::npos)
            str.erase(0, pos);
    } else
        str.erase(str.begin(), str.end());
}

int mme::getNumeric(const path& path, const string& key) {
    string line;
    ifstream infile;
    infile.open(path.c_str());
    int32_t val = -1;
    while (!infile.eof() && line.compare(mme::MX_META_END) != 0) {
        getline(infile, line);
        size_t pos = line.find(key);

        if (pos != string::npos) {
            size_t start = pos + key.length();
            string str_val = line.substr(start, string::npos);
            trim(str_val);
            val = atoi(str_val.c_str());
        }
    }
    infile.close();

    return val;
}

string mme::getString(const path& path, const string& key) {
    string line;
    ifstream infile;
    infile.open(path.c_str());
    string val = "";
    while (!infile.eof() && line.compare(mme::MX_META_END) != 0) {
        getline(infile, line);
        size_t pos = line.find(key);

        if (pos != string::npos) {
            size_t start = pos + key.length();
            string strVal = line.substr(start, string::npos);
            trim(strVal);
            val = strVal;
        }
    }
    infile.close();

    return val;
}
