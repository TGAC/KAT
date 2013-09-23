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

#include "density/density_plot_main.hpp"
#include "profile/profile_plot_main.hpp"
#include "spectra-cn/spectra_cn_plot_main.hpp"
#include "spectra-hist/spectra_hist_plot_main.hpp"

#include "plot_args.hpp"
#include "plot_main.hpp"

using std::string;

using kat::PlotArgs;

// Start point
int kat::plotStart(int argc, char *argv[])
{
    // Parse args
    PlotArgs args(argc, argv);

    // Shortcut to mode
    string mode = args.getMode();

    // Pass remaining args to relevant child tool
    if (mode.compare(KAT_PLOT_DENSITY_ID) == 0)
    {
        kat::densityPlotStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (mode.compare(KAT_PLOT_PROFILE_ID) == 0)
    {
        kat::profilePlotStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (mode.compare(KAT_PLOT_SPECTRA_CN_ID) == 0)
    {
        kat::spectraCnPlotStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (mode.compare(KAT_PLOT_SPECTRA_HIST_ID) == 0)
    {
        kat::spectraHistPlotStart(args.getModeArgC(), args.getModeArgV());
    }

    return 0;
}
