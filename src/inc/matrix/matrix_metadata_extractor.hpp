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

#include <string.h>
#include <stdint.h>

using std::string;


namespace mme
{
    const string KEY_NB_COLUMNS = "# Columns:";
    const string KEY_NB_ROWS    = "# Rows:";
    const string KEY_X_LABEL    = "# XLabel:";
    const string KEY_Y_LABEL    = "# YLabel:";
    const string KEY_Z_LABEL    = "# ZLabel:";
    const string KEY_TITLE      = "# Title:";
    const string KEY_MAX_VAL    = "# MaxVal:";
    const string MX_META_END    = "###";

    void trim(string& str);
    int getNumeric(const string& path, const string& key);
    string getString(const string& path, const string& key);
}

