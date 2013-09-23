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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <iostream>

#include "kat_args.hpp"
#include "sect/sect_main.hpp"
#include "comp/comp_main.hpp"
#include "gcp/gcp_main.hpp"
#include "hist/hist_main.hpp"
#include "plot/plot_main.hpp"

using std::string;

using kat::KatArgs;

/**
 * Start point for KAT.  Processes the start of the command line and then
 * delegates the rest of the command line to the child tool.
 */
int main(int argc, char *argv[])
{
    // Parse args
    KatArgs args(argc, argv);

    // Shortcut to mode
    string mode = args.getMode();

    // Pass remaining args to relevant child tool
    if (mode.compare(kat::KAT_SECT_ID) == 0)
    {
        kat::sectStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (mode.compare(kat::KAT_COMP_ID) == 0)
    {
        kat::compStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (mode.compare(kat::KAT_GCP_ID) == 0)
    {
        kat::gcpStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (mode.compare(kat::KAT_HIST_ID) == 0)
    {
        kat::histStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (mode.compare(kat::KAT_PLOT_ID) == 0)
    {
        kat::plotStart(args.getModeArgC(), args.getModeArgV());
    }

    return 0;
}
