#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>

#include <gnuplot/gnuplot_i.hpp>

#include <matrix/matrix_metadata_extractor.hpp>

#include "flame_plot_args.hpp"
#include "flame_plot_main.hpp"

using std::string;
using std::ifstream;
using std::istringstream;

void configureFlamePlot(Gnuplot* plot, string* type, const char* output_path,
	uint canvas_width, uint canvas_height)
{
    
    std::ostringstream term_str;

    if (type->compare("png") == 0)
    {
        term_str << "set terminal png";
    }
    else if (type->compare("ps") == 0)
    {
        term_str << "set terminal postscript color";
    }
    else if (type->compare("pdf") == 0)
    {
        term_str << "set terminal pdf color";
    }
    else
    {
        std::cerr << "Unknown file type, assuming PNG\n";
        term_str << "set terminal png";
    }

    term_str << " large";
    term_str << " size " << canvas_width << "," << canvas_height;

    plot->cmd(term_str.str());

    std::ostringstream output_str;
    output_str << "set output \"" << output_path << "\"";
    plot->cmd(output_str.str());
}


// Start point
int flamePlotStart(int argc, char *argv[])
{
    // Parse args
    FlamePlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();


    // Get plotting properties, either from file, or user.  User args have precedence.
    uint16_t x_range = args.x_max != DEFAULT_X_MAX ? args.x_max : mme::getNumeric(*args.mx_arg, mme::KEY_NB_COLUMNS);
    uint16_t y_range = args.y_max != DEFAULT_Y_MAX ? args.y_max : mme::getNumeric(*args.mx_arg, mme::KEY_NB_ROWS);
    uint32_t z_cap = args.z_cap != DEFAULT_Z_CAP ? args.z_cap : mme::getNumeric(*args.mx_arg, mme::KEY_MAX_VAL);

    string x_label = args.x_label.compare(DEFAULT_X_LABEL) != 0 ? args.x_label : mme::getString(*args.mx_arg, mme::KEY_X_LABEL);
    string y_label = args.y_label.compare(DEFAULT_Y_LABEL) != 0 ? args.y_label : mme::getString(*args.mx_arg, mme::KEY_Y_LABEL);

    string title = args.title.compare(DEFAULT_TITLE) != 0 ? args.title : mme::getString(*args.mx_arg, mme::KEY_TITLE);



    // If neither the user or the data file contain any ideas of what values to use then use defaults
    x_range = x_range == -1 ? 1001 : x_range;
    y_range = y_range == -1 ? 1001 : y_range;
    z_cap = z_cap == -1 ? 10000 : z_cap / 2;    // Saturate the hot spots a bit to show more detail around the edges

    x_label = x_label.empty() ? DEFAULT_X_LABEL : x_label;
    y_label = y_label.empty() ? DEFAULT_Y_LABEL : y_label;

    title = title.empty() ? DEFAULT_TITLE : title;

    if (args.verbose)
    {
        cerr << "Actual variables used to create plot:" << endl;
        cerr << "X Range: " << x_range << endl;
        cerr << "Y Range: " << y_range << endl;
        cerr << "Z cap: " << z_cap << endl;
        cerr << "X Label: " << x_label << endl;
        cerr << "Y Label: " << y_label << endl;
        cerr << "Title: " << title << endl;
    }


    // Start defining the plot
    Gnuplot* flame = new Gnuplot("lines");

    configureFlamePlot(flame, args.output_type, args.output_path->c_str(), args.width, args.height);

    flame->set_title(title);
    flame->set_xlabel(x_label);
    flame->set_ylabel(y_label);

    flame->set_xrange(0, x_range);
    flame->set_yrange(0, y_range);

    //flame->set_xlogscale();
    //flame->set_ylogscale();
    //flame->set_zlogscale();

    flame->cmd("set palette rgb 21,22,23");
    flame->cmd("set size ratio 1");

    std::ostringstream rangestr;
    rangestr << "set cbrange [0:" << z_cap << "]";
    flame->cmd(rangestr.str());

    std::ostringstream plotstr;
    plotstr << "plot '" << args.mx_arg->c_str() << "' matrix with image";

    flame->cmd(plotstr.str());

    delete flame;

    return 0;
}
