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

