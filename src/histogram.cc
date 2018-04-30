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

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <thread>
#include <sys/ioctl.h>
using std::shared_ptr;
using std::make_shared;
using std::thread;
using std::ofstream;

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
using bfs::path;
using boost::lexical_cast;

#include <jellyfish/mer_dna.hpp>

#include <kat/matrix_metadata_extractor.hpp>
#include <kat/jellyfish_helper.hpp>

#include "plot.hpp"
using kat::Plot;

#include "histogram.hpp"

kat::Histogram::Histogram(vector<path> _inputs, uint64_t _low, uint64_t _high, uint64_t _inc) {

	input.setMultipleInputs(_inputs);
	input.index = 1;
	outputPrefix = "kat-hist";
	low = _low;
	high = _high;
	inc = _inc;
	threads = 1;

	// Calculate other vars required for this run
	base = calcBase();
	ceil = calcCeil();
	nb_buckets = ceil + 1 - base;
}

void kat::Histogram::execute() {

	// Some validation first
	if (high < low) {
		BOOST_THROW_EXCEPTION(HistogramException() << HistogramErrorInfo(string(
			"High count value must be >= to low count value.  High: ") + lexical_cast<string>(high) +
			"; Low: " + lexical_cast<string>(low)));
	}

	// Validate input
	input.validateInput();

	// Create output directory
	path parentDir = bfs::absolute(outputPrefix).parent_path();
	KatFS::ensureDirectoryExists(parentDir);

	// Either count or load input
	if (input.mode == InputHandler::InputHandler::InputMode::COUNT) {
		input.count(threads);
	} else {
		input.loadHeader();
		input.loadHash();
	}

	data = vector<uint64_t>(nb_buckets, 0);
	threadedData = vector<shared_ptr<vector < uint64_t>>>();

	// Do the work
	bin();

	// Dump any hashes that were previously counted to disk if requested
	// NOTE: MUST BE DONE AFTER COMPARISON AS THIS CLEARS ENTRIES FROM HASH ARRAY!
	if (input.dumpHash) {
		path outputPath(outputPrefix.string() + "-hash.jf" + lexical_cast<string>(input.merLen));
		input.dump(outputPath, threads);
	}
	// Merge results
	merge();


}

void kat::Histogram::save() {

	auto_cpu_timer timer(1, "  Time taken: %ws\n\n");

	cout << "Saving results to disk ...";
	cout.flush();

	// Send main matrix to output file
	ofstream main_hist_out_stream(outputPrefix.c_str());
	print(main_hist_out_stream);
	main_hist_out_stream.close();

	cout << " done.";
	cout.flush();
}

void kat::Histogram::print(std::ostream &out) {
	// Output header
	out << mme::KEY_TITLE << input.merLen << "-mer spectra for: " << input.fileName() << endl;
	out << mme::KEY_X_LABEL << input.merLen << "-mer frequency" << endl;
	out << mme::KEY_Y_LABEL << "# distinct " << input.merLen << "-mers" << endl;
	out << mme::KEY_KMER << input.merLen << endl;
	out << mme::KEY_INPUT_1 << input.pathString() << endl;
	out << mme::MX_META_END << endl;

	uint64_t col = base;
	for (uint64_t i = 0; i < nb_buckets; i++, col += inc) {
		out << col << " " << data[i] << "\n";
	}
}

void kat::Histogram::merge() {
	auto_cpu_timer timer(1, "  Time taken: %ws\n\n");

	cout << "Merging counts ...";
	cout.flush();

	for (size_t i = 0; i < nb_buckets; i++) {
		for (size_t j = 0; j < threads; j++) {
			data[i] += threadedData[j]->at(i);
		}
	}

	cout << " done.";
	cout.flush();
}

void kat::Histogram::bin() {

	auto_cpu_timer timer(1, "  Time taken: %ws\n\n");

	cout << "Bining kmers ...";
	cout.flush();

	vector<thread> t(threads);

	for (uint16_t i = 0; i < threads; i++) {
		t[i] = thread(&Histogram::binSlice, this, i);
	}

	for (uint16_t i = 0; i < threads; i++) {
		t[i].join();
	}

	cout << " done.";
	cout.flush();
}

void kat::Histogram::binSlice(int th_id) {

	shared_ptr<vector < uint64_t>> hist = make_shared<vector < uint64_t >> (nb_buckets);

	LargeHashArray::region_iterator it = input.hash->region_slice(th_id, threads);
	while (it.next()) {
		uint64_t val = it.val();
		if (val < base)
			++(*hist)[0];
		else if (val > ceil)
			++(*hist)[nb_buckets - 1];
		else
			++(*hist)[(val - base) / inc];
	}

	threadedData.push_back(hist);
}

void kat::Histogram::analysePeaks() {
#ifdef HAVE_PYTHON
	cout << "Analysing peaks" << endl
		 << "---------------" << endl;
	
	path dascript = path("kat") / "distanalysis.py";

	vector<string> args;
	args.push_back(dascript.string());
	if (verbose) {
		args.push_back("--verbose");
	}
	args.push_back("--from_kat");
	args.push_back("--output_prefix=" + outputPrefix.string());
	args.push_back(outputPrefix.string());

	char* char_args[50];

	for (size_t i = 0; i < args.size(); i++) {
		char_args[i] = strdup(args[i].c_str());
	}

	PyHelper::getInstance().execute(dascript.string(), (int) args.size(), char_args);

	for (size_t i = 0; i < args.size(); i++) {
		free(char_args[i]);
	}

	cout << endl;
#endif
}

void kat::Histogram::plot(const string& output_type) {

#ifdef HAVE_PYTHON
	auto_cpu_timer timer(1, "  Time taken: %ws\n\n");

	cout << "Creating plot ...";
	cout.flush();

	path outputFile1 = path(outputPrefix.string() + "." + output_type);

	vector<string> args;
	args.push_back("kat/plot/spectra-hist.py");
	args.push_back(string("--output=") + outputFile1.string());
	if (verbose) {
		args.push_back("--verbose");
	}
	args.push_back(outputPrefix.string());
	Plot::executePythonPlot(Plot::PlotMode::SPECTRA_HIST, args);

	cout << " done.";
	cout.flush();
#endif
}

int kat::Histogram::main(int argc, char *argv[]) {

	vector<path> inputs;
	path output_prefix;
	uint16_t threads;
	uint64_t low;
	uint64_t high;
	uint64_t inc;
	string trim5p;
	bool non_canonical;
	uint16_t mer_len;
	uint64_t hash_size;
	bool dump_hash;
	string plot_output_type;
	bool verbose;
	bool help;

	struct winsize w;
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);


	// Declare the supported options.
	po::options_description generic_options(Histogram::helpMessage(), w.ws_col);
	generic_options.add_options()
		("output_prefix,o", po::value<path>(&output_prefix)->default_value("kat.hist"),
		"Path prefix for files generated by this program.")
		("threads,t", po::value<uint16_t>(&threads)->default_value(1),
		"The number of threads to use")
		("low,l", po::value<uint64_t>(&low)->default_value(1),
		"Low count value of histogram")
		("high,h", po::value<uint64_t>(&high)->default_value(10000),
		"High count value of histogram")
		("inc,i", po::value<uint64_t>(&inc)->default_value(1),
		"Increment for each bin")
		("5ptrim", po::value<string>(&trim5p)->default_value("0"),
		"Ignore the first X bases from reads.  If more that one file is provided you can specify different values for each file by seperating with commas.")
		("non_canonical,N", po::bool_switch(&non_canonical)->default_value(false),
		"If counting fast(a/q), store explicit kmer as found.  By default, we store 'canonical' k-mers, which means we count both strands.")
		("mer_len,m", po::value<uint16_t>(&mer_len)->default_value(DEFAULT_MER_LEN),
		"The kmer length to use in the kmer hashes.  Larger values will provide more discriminating power between kmers but at the expense of additional memory and lower coverage.")
		("hash_size,H", po::value<uint64_t>(&hash_size)->default_value(DEFAULT_HASH_SIZE),
		"If kmer counting is required for the input, then use this value as the hash size.  If this hash size is not large enough for your dataset then the default behaviour is to double the size of the hash and recount, which will increase runtime and memory usage.")
		("dump_hash,d", po::bool_switch(&dump_hash)->default_value(false),
		"Dumps any jellyfish hashes to disk that were produced during this run. Normally, this is not recommended, and will likely consume a significant amount of disk space.")
		("output_type,p", po::value<string>(&plot_output_type)->default_value(DEFAULT_HIST_PLOT_OUTPUT_TYPE),
		"The plot file type to create: png, ps, pdf.")
		("verbose,v", po::bool_switch(&verbose)->default_value(false),
		"Print extra information.")
		("help", po::bool_switch(&help)->default_value(false), "Produce help message.")
		;

	// Hidden options, will be allowed both on command line and
	// in config file, but will not be shown to the user.
	po::options_description hidden_options("Hidden options");
	hidden_options.add_options()
		("inputs", po::value<std::vector < path >> (&inputs), "Path to the input file(s) to process.")
		;

	// Positional option for the input bam file
	po::positional_options_description p;
	p.add("inputs", -1);

	// Combine non-positional options
	po::options_description cmdline_options;
	cmdline_options.add(generic_options).add(hidden_options);

	// Parse command line
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	// Output help information the exit if requested
	if (help || argc <= 1) {
		cout << generic_options << endl;
		return 1;
	}

	vector<string> d1_5ptrim_strs;
	vector<uint16_t> d1_5ptrim_vals;
	boost::split(d1_5ptrim_strs, trim5p, boost::is_any_of(","));
	for (auto& v : d1_5ptrim_strs) d1_5ptrim_vals.push_back(boost::lexical_cast<uint16_t>(v));


	auto_cpu_timer timer(1, "KAT HIST completed.\nTotal runtime: %ws\n\n");

	cout << "Running KAT in HIST mode" << endl
		<< "------------------------" << endl << endl;

	// Create the sequence coverage object
	Histogram histo(inputs, low, high, inc);
	histo.setOutputPrefix(output_prefix);
	histo.setThreads(threads);
	histo.setTrim(d1_5ptrim_vals);
	histo.setCanonical(!non_canonical);
	histo.setMerLen(mer_len);
	histo.setHashSize(hash_size);
	histo.setDumpHash(dump_hash);
	histo.setVerbose(verbose);

	// Do the work
	histo.execute();

	// Save results
	histo.save();

	// Plot
#ifdef HAVE_PYTHON
	histo.plot(plot_output_type);
	histo.analysePeaks();
#endif

	return 0;
}
