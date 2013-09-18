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

// OTHER INCLUDES

// MODIFY FOR YOUR TOOL
#include "template_args.hpp"
#include "template_main.hpp"


// Start point
int templateStart(int argc, char *argv[])
{
    // Parse args (MODIFY FOR YOUR TOOL)
    TemplateArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    // ADD TOOL LOGIC HERE



    return 0;
}
