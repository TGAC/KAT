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

#include <glob.h>

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/timer/timer.hpp>
#include <boost/algorithm/string.hpp>
namespace bfs = boost::filesystem;
using bfs::path;
using boost::timer::auto_cpu_timer;
using boost::split;

#include <kat/jellyfish_helper.hpp>
using kat::JellyfishHelper;

#include <kat/input_handler.hpp>

void kat::InputHandler::setMultipleInputs(const vector<path>& inputs) {
    for(auto& p : inputs) {
        input.push_back(p);      
    }
}

void kat::InputHandler::validateInput() {
    
    bool start = true;
    
    // Check input file(s) exists
    for(auto& rp : input) {
    
        path p(rp);
        
        if (bfs::is_symlink(p)) {
            if (bfs::symbolic_link_exists(p)) {
                p = bfs::canonical(p);
            }
            else {
                BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
                    "Could not find input file at: ") + p.string() + "; please check the path and try again."));
            }
        }
        
        if (!bfs::exists(p)) {
            BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
                    "Could not find input file at: ") + p.string() + "; please check the path and try again."));
        }
        
        InputMode m = JellyfishHelper::isSequenceFile(p) ? InputMode::COUNT : InputMode::LOAD;
        
        if (start) {
            mode = m;
        }
        else {
            if (m != mode) {
                BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
                    "Cannot mix sequence files and jellyfish hashes.  Input: ") + p.string()));
            }
        }
    }
    
    
}

void kat::InputHandler::loadHeader() {
    if (mode == InputMode::LOAD) {
        header = JellyfishHelper::loadHashHeader(input[0]);
    }    
}

void kat::InputHandler::validateMerLen(const uint16_t merLen) {    
    
    if (mode == InputMode::LOAD) {
        if (header->key_len() != merLen * 2) {

            BOOST_THROW_EXCEPTION(JellyfishException() << JellyfishErrorInfo(string(
                "Cannot process hashes that were created with different K-mer lengths.  Expected: ") +
                lexical_cast<string>(merLen) + 
                ".  Key length was " + 
                lexical_cast<string>(header->key_len() / 2) + 
                " for : " + input[0].string()));
        }
    }
}

string kat::InputHandler::pathString() {
    
    string s;
    for(auto& p : input) {
        s += p.string() + " ";
    }
    return boost::trim_right_copy(s);
}

void kat::InputHandler::count(const uint16_t threads) {
    
    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");      
    
    hashCounter = make_shared<HashCounter>(hashSize, merLen * 2, 7, threads);
    hashCounter->do_size_doubling(!disableHashGrow);
        
    cout << "Input is a sequence file.  Counting kmers for " << pathString() << "...";
    cout.flush();

    hash = JellyfishHelper::countSeqFile(input, *hashCounter, canonical, threads);
    
    // Create header for newly counted hash
    header = make_shared<file_header>();
    header->fill_standard();
    header->update_from_ary(*hash);
    header->counter_len(4);  // Hard code for now.
    header->canonical(canonical);
    header->format(binary_dumper::format);
    
    cout << " done.";
    cout.flush();    
}

void kat::InputHandler::loadHash() {
    
    auto_cpu_timer timer(1, "  Time taken: %ws\n\n");        

    cout << "Loading hashes into memory...";
    cout.flush();  
    
    hashLoader = make_shared<HashLoader>();
    hashLoader->loadHash(input[0], false); 
    hash = hashLoader->getHash();
    canonical = hashLoader->getCanonical();
    merLen = hashLoader->getMerLen();
    
    cout << " done.";
    cout.flush();    
}

void kat::InputHandler::dump(const path& outputPath, const uint16_t threads) {
    
    // Remove anything that exists at the target location
    if (bfs::is_symlink(outputPath) || bfs::exists(outputPath)) {
        bfs::remove(outputPath.c_str());
    }

    // Either dump or symlink as appropriate
    if (mode == InputHandler::InputHandler::InputMode::COUNT) {
    
        auto_cpu_timer timer(1, "  Time taken: %ws\n\n"); 
        cout << "Dumping hash to " << outputPath.string() << " ...";
        cout.flush();

        JellyfishHelper::dumpHash(hash, *header, threads, outputPath);

        cout << " done.";
        cout.flush();
    }
    else {
        bfs::create_symlink(getSingleInput(), outputPath);
    }
}

void kat::InputHandler::globFiles(const string& input, vector<path>& globbed) {

    vector<string> inputvec;
    boost::split(inputvec, input, boost::is_any_of(" "));    
    
    vector<path> pathvec;
    for(auto& s : inputvec) {
        pathvec.push_back(s);
    }
    
    globFiles(pathvec, globbed);
}

void kat::InputHandler::globFiles(const vector<path>& input, vector<path>& globbed) {

    glob_t globbuf;

    // Translate glob patterns into real paths
    int i = 0;
    for(auto& g : input) {           
        glob(g.c_str(), i > 0 ? GLOB_TILDE | GLOB_APPEND : GLOB_TILDE, NULL, &globbuf);
        i++;
    }

    for( size_t i = 0; i < globbuf.gl_pathc; ++i )
        globbed.push_back( path(globbuf.gl_pathv[i]) );

    if( globbuf.gl_pathc > 0 )
        globfree( &globbuf );

    // Check for content.  If there isn't any then probably the input file doesn't
    // exist.  But we add the basic input regardless, the user will have to check
    // for file non-existence later.
    if (globbed.empty()) {
        globbed.push_back(input[0]);
    }    
}
