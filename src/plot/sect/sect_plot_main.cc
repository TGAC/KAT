#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <gnuplot/gnuplot_i.hpp>

#include "sect_plot_args.hpp"
#include "sect_plot_main.hpp"

using std::string;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::vector;

using seqan::SequenceStream;
using seqan::CharString;
using seqan::Dna5String;

// Finds a particular fasta header in a fasta file and returns the associated sequence
const string* getEntryFromFasta(const string *fasta_path, const string header)
{
    // Setup stream to sequence file and check all is well
    SequenceStream seqStream(fasta_path->c_str(), SequenceStream::READ);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the sequence file: " << *fasta_path << endl;
        return NULL;
    }

    CharString id;
    Dna5String seq;

    while(!atEnd(seqStream))
    {
        // Read record
        if (readRecord(id, seq, seqStream) != 0)
        {
            std::cerr << "ERROR: Could not read record!" << endl;
            return NULL;
        }

        stringstream ssHeader;
        ssHeader << id;

        if (header.compare(ssHeader.str()) == 0)
        {
            stringstream ssSeq;
            ssSeq << seq;

            return new string(ssSeq.str());
        }
    }

    return NULL;
}

// Finds the nth entry from the fasta file and returns the associated sequence
const string* getEntryFromFasta(const string *fasta_path, uint32_t n)
{
    // Setup stream to sequence file and check all is well
    SequenceStream seqStream(fasta_path->c_str(), SequenceStream::READ);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the sequence file: " << *fasta_path << endl;
        return NULL;
    }

    CharString header;
    Dna5String seq;
    uint32_t i = 1;

    while(!atEnd(seqStream))
    {
        // Read record
        if (readRecord(header, seq, seqStream) != 0)
        {
            std::cerr << "ERROR: Could not read record!" << endl;
            return NULL;
        }

        if (i == n)
        {
            stringstream ssSeq;
            ssSeq << seq;
            return new string(ssSeq.str());
        }

        i++;
    }

    return NULL;
}



void configureSectPlot(Gnuplot* plot, string* type, const char* output_path,
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

uint32_t strToInt(string s)
{
    istringstream str_val(s);
    uint32_t int_val;
    str_val >> int_val;
    return int_val;
}

unsigned int split(const string *txt, vector<uint32_t> *strs, const char ch)
{
    unsigned int pos = txt->find( ch );
    unsigned int initialPos = 0;
    strs->clear();

    // Decompose statement
    while( pos != string::npos ) {

        strs->push_back( strToInt(txt->substr( initialPos, pos - initialPos + 1 )) );
        initialPos = pos + 1;

        pos = txt->find( ch, initialPos );
    }

    // Add the last one
    strs->push_back( strToInt(txt->substr( initialPos, std::min( pos, (unsigned int)txt->size() ) - initialPos + 1 )) );

    return strs->size();
}


// Start point
int sectPlotStart(int argc, char *argv[])
{
    // Parse args
    SectPlotArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    const string* coverages = NULL;

    if (args.fasta_header)
    {
        string header(args.fasta_header);
        coverages = getEntryFromFasta(args.sect_file_arg, header);
    }
    else if (args.fasta_index > 0)
    {
        coverages = getEntryFromFasta(args.sect_file_arg, args.fasta_index);
    }

    if (!coverages)
    {
        cerr << "Could not find requested fasta header in sect coverages fasta file" << endl;
    }
    else
    {
        if (args.verbose)
            cerr << "Found requested sequence" << endl;

        // Split coverages
        vector<uint32_t> cvs(coverages->length() / 2);
        split(coverages, &cvs, ' ');

        if (args.verbose)
            cerr << "Acquired kmer counts" << endl;

        // Initialise gnuplot
        Gnuplot* sect_plot = new Gnuplot("lines");

        configureSectPlot(sect_plot, args.output_type, args.output_arg->c_str(), args.width, args.height);

        sect_plot->set_title(args.title);
        sect_plot->set_xlabel(args.x_label);
        sect_plot->set_ylabel(args.y_label);

        sect_plot->cmd("set style data linespoints");

        std::ostringstream data_str;

        for(uint32_t i = 0; i < cvs.size(); i++)
        {
            uint32_t index = i+1;
            data_str << index << " " << cvs[i] << "\n";
        }

        std::ostringstream plot_str;

        plot_str << "plot '-'\n" << data_str.str() << "e\n";

        sect_plot->cmd(plot_str.str());

        delete sect_plot;
        delete coverages;

        if (args.verbose)
            cerr << "Plotted data" << endl;
    }

    return 0;
}
