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

#include <iostream>
using std::endl;
using std::string;
using std::cout;

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/exception/all.hpp>
namespace fs = boost::filesystem;
using fs::exists;
using fs::path;


namespace kat {
    
    typedef boost::error_info<struct FileSystemError,string> FileSystemErrorInfo;
    struct FileSystemException: virtual boost::exception, virtual std::exception { };

 
    class KatFS {
    private:
        
        path exe;
        
        bool isAbsolute;
        bool isRelative;
        bool isOnPath;
        
        // Executables
        path canonicalExe;
        
        // Directories
        path scriptsDir;
        
    public:
       
        /**
         * Assume on PATH by default
         */
        KatFS() : KatFS("kat") {}
        
        /**
         * 
         * @param exe Full path to the exe, probably derived from argv0.
         */
        KatFS(const char* argv) {
            
            isAbsolute = false;
            isRelative = false;
            isOnPath = false;
            
            exe = argv;
            
            
            if(exe.is_absolute()) {
                
                // Absolute path provided...  Easy job... nothing special to do, 
                // resolve symlink then take two levels up
                canonicalExe = fs::canonical(exe);
                isAbsolute = true;
            }
            else if (exe.string().find('/') != string::npos) {
                
                // Relative with some parent paths... get absolute path, resolving 
                // symlinks then take two levels up
                canonicalExe = fs::canonical(fs::system_complete(exe));
                isRelative = true;
            }
            else {

                // Only name provided
                // In this case we just have to assume everything is properly installed 
                canonicalExe = exe;
                isOnPath = true;     
            }
            
                            
            // We assume scripts are on the path if exe was on the path
            if (isAbsolute || isRelative) {
                
                // Check to see if scripts are adjacent to exe first
                path kda(canonicalExe.parent_path());
                kda /= "kat_distanalysis.py";
                if (exists(kda)) {
                    scriptsDir = canonicalExe.parent_path();
                }
                else {
                
                    // Not 100% sure how far back we need to go (depends on whether using KAT exe or tests) 
                    // so try 2, 3 and 4 levels.
                    scriptsDir = canonicalExe.parent_path().parent_path();
                    scriptsDir /= "scripts";                 

                    if (!exists(scriptsDir)) {
                        scriptsDir = canonicalExe.parent_path().parent_path().parent_path();
                        scriptsDir /= "scripts";       
                        
                        if (!exists(scriptsDir)) {
                            scriptsDir = canonicalExe.parent_path().parent_path().parent_path().parent_path();
                            scriptsDir /= "scripts";        

                            if (!exists(scriptsDir)) {
                                BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                                    "Could not find suitable directory containing KAT scripts relative to provided exe: ") + canonicalExe.c_str()));
                            }
                        }
                    }

                    // Also double check the kat_distanalysis.py script exists
                    kda = scriptsDir;
                    kda /= "kat_distanalysis.py";

                    if (!exists(kda)) {
                        BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                            "Found the scripts directory where expected") + scriptsDir.string() + 
                                ". However, could not find the \"kat_distanalysis.py\" script inside."));
                    }
                }
            }
        }
        
        
        // **** Destructor ****
        virtual ~KatFS() { }
        
        path GetCanonicalExe() const {
            return canonicalExe;
        }

        path GetScriptsDir() const {
            return scriptsDir;
        }
        
        
        friend std::ostream& operator<<(std::ostream &strm, const KatFS& pfs) {
            
            return strm << "KAT paths: "<< endl
                        << " - argv: " << pfs.exe << endl
                        << "  - type: " << (pfs.isAbsolute ? "absolute" : pfs.isRelative ? "relative" : "on PATH") << endl
                        << " - Canonical path: " << pfs.canonicalExe << endl
                        << " - Scripts dir: " << (pfs.scriptsDir.empty() ? "assuming scripts on PATH" : pfs.scriptsDir) << endl << endl;                   
        }     
    };
    
       
    // Make available everywhere
    extern KatFS katFileSystem;
}

