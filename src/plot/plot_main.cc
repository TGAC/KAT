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

#include <string.h>

#include "flame/flame_plot_main.hpp"
#include "asm/asm_plot_main.hpp"
#include "sect/sect_plot_main.hpp"
#include "spectra/spectra_plot_main.hpp"

#include "plot_args.hpp"
#include "plot_main.hpp"

using std::string;

using kat::PlotArgs;

// Start point
int kat::plotStart(int argc, char *argv[])
{
    // Parse args
    PlotArgs args(argc, argv);

    // Pass remaining args to relevant child tool
    if (args.getMode().compare("flame") == 0)
    {
        kat::flamePlotStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("asm") == 0)
    {
        kat::asmPlotStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("sect") == 0)
    {
        kat::sectPlotStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (args.getMode().compare("spectra") == 0)
    {
        kat::spectraPlotStart(args.getModeArgC(), args.getModeArgV());
    }

    return 0;
}
