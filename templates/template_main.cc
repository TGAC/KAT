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
