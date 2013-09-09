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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "kat_args.hpp"
#include "sect/sect_main.hpp"
#include "comp/comp_main.hpp"
#include "gcp/gcp_main.hpp"
#include "histo/histo_main.hpp"
#include "plot/plot_main.hpp"

/**
 * Start point for KAT.  Processes the start of the command line and then
 * delegates the rest of the command line to the child tool.
 */
int main(int argc, char *argv[])
{
    // Parse args
    KatArgs args(argc, argv);

    // Pass remaining args to relevant child tool
    if (args.getMode().compare("sect") == 0)
    {
        sectStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("comp") == 0)
    {
        compStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("gcp") == 0)
    {
        gcpStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("histo") == 0)
    {
        kat::histoStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("plot") == 0)
    {
        plotStart(args.getModeArgC(), args.getModeArgV());
    }

    return 0;
}
