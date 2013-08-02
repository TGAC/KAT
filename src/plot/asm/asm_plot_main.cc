#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>

#include <gnuplot/gnuplot_i.hpp>

#include "asm_plot_args.hpp"
#include "asm_plot_main.hpp"

using std::string;
using std::ostringstream;
using std::vector;


string createLineStyleStr(uint16_t i, const char* colour)
{
    ostringstream line_style_str;
    line_style_str << "set style line " << i << " lc rgb \"" << colour << "\"";
    return line_style_str.str();
}

string createSinglePlotString(const char* data_file, uint16_t idx)
{
    uint16_t col = idx+1;

    ostringstream plot_str;
    plot_str << "'" << data_file << "' u " << col << " t \"" << idx << "x\"";
    return plot_str.str();
}


vector<uint16_t>* getStandardCols(AsmPlotArgs* args)
{
    vector<uint16_t>* col_defs = new vector<uint16_t>();

    if (!args->ignore_absent)
    {
        col_defs->push_back(0);
    }

    for(uint16_t i = 1; i <= args->max_duplication; i++)
    {
        col_defs->push_back(i);
    }

    return col_defs;
}


vector<uint16_t>* getUserDefinedCols(AsmPlotArgs* args)
{
    vector<uint16_t>* col_defs = new vector<uint16_t>();

    string delimiter = ",";
    string s(args->columns->c_str());

    size_t pos = 0;
    string token;
    while ((pos = s.find(delimiter)) != string::npos) {
        token = s.substr(0, pos);

        col_defs->push_back(atoi(s.c_str()));

        s.erase(0, pos + delimiter.length());
    }

    col_defs->push_back(atoi(s.c_str()));

    return col_defs;
}




// Start point
int asmPlotStart(int argc, char *argv[])
{
    // Parse args
    AsmPlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();




    vector<uint16_t>* plot_cols = !args.columns ? getStandardCols(&args) : getUserDefinedCols(&args);

    if (!plot_cols->empty())
    {
        // Determine configuration
        bool request_absent = (*plot_cols)[0] == 0 ? true : false;
        uint16_t level_count = request_absent ? plot_cols->size() - 1 : plot_cols->size();


        if (args.verbose)
        {
            cerr << "Request plot fo absent kmers" << endl
                 << level_count << " levels of present kmers requested for plotting" << endl << endl;
        }

        // Initialise gnuplot
        Gnuplot* asm_plot = new Gnuplot("lines");


        // Work out the output path to use (either user specified or auto generated)
        string output_path = args.determineOutputPath();


        asm_plot->configurePlot(*(args.output_type), output_path, args.width, args.height);

        asm_plot->set_title(args.title);
        asm_plot->set_xlabel(args.x_label);
        asm_plot->set_ylabel(args.y_label);


        // Get plot strings
        ostringstream plot_str;



        bool first = true;

        if (request_absent)
        {
            plot_str << createSinglePlotString(args.mx_arg->c_str(), 0) << " lt rgb \"black\"";
            first = false;
        }

        for(uint16_t i = 0; i < level_count; i++)
        {
            if (first)
                first = false;
            else
                plot_str << ", ";


            double col_frac = 1.0 - ((double)i / (double)(level_count - 1));
            uint16_t index = request_absent ? i+1 : i;

            plot_str << createSinglePlotString(args.mx_arg->c_str(), (*plot_cols)[index]) << " lt palette frac " << std::fixed << col_frac;
        }

        asm_plot->cmd("set palette rgb 33,13,10");
        asm_plot->cmd("unset colorbox");

        asm_plot->cmd("set style fill solid 1 noborder");
        asm_plot->cmd("set style histogram rowstacked");
        asm_plot->cmd("set style data histograms");

        asm_plot->set_xrange(0, args.x_max);
        asm_plot->set_yrange(0, args.y_max);

        ostringstream plot_cmd;
        plot_cmd << "plot " << plot_str.str();
        asm_plot->cmd(plot_cmd.str());

        if (args.verbose)
            cerr << "Gnuplot command: " << plot_cmd.str() << endl;

        delete asm_plot;
    }


    delete plot_cols;

    return 0;
}
