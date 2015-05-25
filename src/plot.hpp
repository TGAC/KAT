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

#pragma once

#include <string>
#include <vector>
using std::string;
using std::vector;

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;

namespace kat {
    
    typedef boost::error_info<struct KatPlotError,string> KatPlotErrorInfo;
    struct KatPlotException: virtual boost::exception, virtual std::exception { };
    
    class BasePlotArgs {
    private:
        bool title_mod;
        bool x_label_mod;
        bool y_label_mod;


        const string currentStatus() const {
            ostringstream status;

            status << "Output type: " << output_type.c_str() << endl;
            status << "Output file specified: " << output_arg.c_str() << endl;
            status << "Plot title: " << title << endl;
            status << "X Label: " << x_label << endl;
            status << "Y Label: " << y_label << endl;
            status << "Width: " << width << endl;
            status << "Height: " << height << endl;

            return status.str();
        }
    public:
        string output_type;
        string output_arg;
        string title;
        string x_label;
        string y_label;
        uint16_t width;
        uint16_t height;

        /**
         * @brief BasePlotArgs Constructor.  Requires the number of trailing arguments to be specified for this plotting tool.
         * @param min_args
         */
        BasePlotArgs(uint16_t min_args) : BaseArgs(min_args), output_type(DEFAULT_OUTPUT_TYPE), output_arg("") {
            title_mod = false;
            x_label_mod = false;
            y_label_mod = false;
        }

        /**
         * @brief ~BasePlotArgs Virtual destructor makes this class abstract
         */
        virtual ~BasePlotArgs() {
        }


        // **************************************************
        // Default values for properties which must be defined
        // by the child class

        virtual const string defaultOutputPrefix() const = 0;
        virtual const string defaultTitle() const = 0;
        virtual const string defaultXLabel() const = 0;
        virtual const string defaultYLabel() const = 0;
        virtual const uint16_t defaultWidth() const = 0;
        virtual const uint16_t defaultHeight() const = 0;

        /**
         * @brief determineOutputPath Work out the output path to use (either user specified or auto generated)
         * @return
         */
        string determineOutputPath() {
            std::ostringstream output_str;
            output_str << defaultOutputPrefix() << "." << output_type;
            return output_arg.empty() ? output_str.str() : output_arg;
        }

        const bool titleModified() {
            return title_mod;
        }

        const bool xLabelModified() {
            return x_label_mod;
        }

        const bool yLabelModified() {
            return y_label_mod;
        }

    };
    
    class Plot {
    
    protected:
        
        enum PlotMode {
            DENSITY,
            PROFILE,
            SPECTRA_CN,
            SPECTRA_HIST,
            SPECTRA_MX
        };
        

        PlotMode parseMode(string mode) {

            string upperMode = boost::to_upper_copy(mode);

            if (upperMode == string("DENSITY")) {
                return DENSITY;                
            }
            else if (upperMode == string("PROFILE")) {
                return PROFILE;
            }
            else if (upperMode == string("SPECTRA_CN")) {
                return SPECTRA_CN;
            }
            else if (upperMode == string("SPECTRA_HIST")) {
                return SPECTRA_HIST;
            }
            else if (upperMode == string("SPECTRA_MX")) {
                return SPECTRA_MX;
            }
            else {
                BOOST_THROW_EXCEPTION(KatPlotException() << KatPlotErrorInfo(string(
                            "Could not recognise mode string: ") + mode));
            }
        }
        
        const string helpMessage() const {
            return string("Usage: kat plot <mode>\n\n") +
                    "Create K-mer Plots\n\n" +
                    "First argument should be the plot mode you wish to use:\n" \
                    "  - density:         Creates a density plot from a matrix created with the \"comp\" tool.  Typically this is\n" \
                    "                     used to compare two K-mer hashes produced by different NGS reads.\n" \
                    "  - profile:         Creates a K-mer coverage plot for a single sequence.  Takes in fasta coverage output\n" \
                    "                     coverage from the \"sect\" tool\n" \
                    "  - spectra-cn:      Creates a stacked histogram using a matrix created with the \"comp\" tool.  Typically\n" \
                    "                     this is used to compare a jellyfish hash produced from a read set to a jellyfish hash\n" \
                    "                     produced from an assembly. The plot shows the amount of distinct K-mers absent, as well\n" \
                    "                     as the copy number variation present within the assembly.\n" \
                    "  - spectra-hist:    Creates a K-mer spectra plot for a set of K-mer histograms produced either by jellyfish-\n" \
                    "                     histo or kat-histo.\n" \
                    "  - spectra-mx:      Creates a K-mer spectra plot for a set of K-mer histograms that are derived from\n" \
                    "                     selected rows or columns in a matrix produced by the \"comp\".";
        }
        
    public:
        static int main(int argc, char *argv[]) {
            
            string modeStr;
            vector<string> others;
            bool help;
        
            // Declare the supported options.
            po::options_description generic_options(Plot::helpMessage());
            generic_options.add_options()
                    ("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
                    ;

            // Hidden options, will be allowed both on command line and
            // in config file, but will not be shown to the user.
            po::options_description hidden_options("Hidden options");
            hidden_options.add_options()
                    ("mode", po::value<string>(&modeStr), "KAT PLOT mode.")
                    ("others", po::value< std::vector<string> >(&others), "Other options.")
                    ;

            // Positional option for the input bam file
            po::positional_options_description p;
            p.add("mode", 1);
            p.add("others", 100);

            // Parse command line
            po::variables_map vm;
            po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).allow_unregistered().run(), vm);
            po::notify(vm);

            // Output help information the exit if requested
            if (argc == 1 || argc == 2 && help) {
                cout << generic_options << endl;
                return 1;
            }
        
            Mode mode = parseMode(modeStr);
        
            const int modeArgC = argc-1;
            char** modeArgV = argv+1;
        
            switch (mode) {
            case DENSITY:
                PlotDensity::main(modeArgC, modeArgV);
                break;
            case PROFILE:
                PlotProfile::main(modeArgC, modeArgV);
                break;
            case SPECTRA_CN:
                PlotSpectraCn::main(modeArgC, modeArgV);
                break;
            case SPECTRA_HIST:
                PlotSpectraHist::main(modeArgC, modeArgV);
                break;
            case SPECTRA_MX
                PlotSpectraMx::main(modeArgC, modeArgV);            
                break;
            default:
                BOOST_THROW_EXCEPTION(KatPlotException() << KatPlotErrorInfo(string(
                        "Unrecognised KAT PLOT mode: ") + modeStr));
            }

            return 0;
        }
    };
    
}

