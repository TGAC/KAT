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

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
namespace fs = boost::filesystem;
using fs::exists;
using fs::path;

namespace kat {
    
    typedef boost::error_info<struct FileSystemError,string> FileSystemErrorInfo;
    struct FileSystemException: virtual boost::exception, virtual std::exception { };

 
    class KatFS {
    private:
        
        // Executables
        path katExe;
        
        // Directories
        path binDir;
        path rootDir;
        path libsDir;
        path scriptsDir;
        
        
        string exec(const char* cmd) {
            FILE* pipe = popen(cmd, "r");
            if (!pipe) return "ERROR";
            char buffer[512];
            string result = "";
            while(!feof(pipe)) {
                if(fgets(buffer, 512, pipe) != NULL)
                        result += buffer;
            }
            pclose(pipe);
            return result;
        }
    
    public:
       
        KatFS() {}
        
        /**
         * 
         * @param exe Full path to the exe, probably derived from argv0.
         */
        KatFS(const char* argv) {
            
            path exe(argv);
            
            //cout << exe << endl;
            
            if(exe.is_absolute()) {
                
                //cout << "Absolute" << endl;
                
                // Easy job... nothing special to do, resolve symlink then take two levels up
                katExe = fs::canonical(exe);
                rootDir = katExe.parent_path().parent_path();
            }
            else if (exe.string().find('/') != string::npos) {
                
                //cout << "Relative" << endl;
                
                // Relative with some parent paths... get absolute path, resolving symlinks then take two levels up
                katExe = fs::canonical(fs::system_complete(exe));
                rootDir = katExe.parent_path().parent_path();
            }
            else {

                //cout << "name only" << endl;
                
                // Tricky one... just exe name, no indication of where if comes from. Now we have to resort to using which.
                string cmd = string("which ") + exe.string();
                string res = exec(cmd.c_str());
                string fullpath = res.substr(0, res.length() - 1);

                //cout << "fullpath" << fullpath << endl;
                katExe = fs::canonical(path(fullpath));
                rootDir = katExe.parent_path().parent_path();
            }
            
            binDir = path(rootDir);
            binDir /= "bin";
            
            libsDir = path(rootDir);
            libsDir /= ".libs";
            
            path srcDir = path(rootDir);
            srcDir /= "src";
            
            path testDir = path(rootDir);
            testDir /= "tests";
                            
            if (katExe.parent_path() == srcDir || katExe.parent_path() == libsDir || katExe.parent_path() == testDir) {
                scriptsDir = path(rootDir);
                scriptsDir /= "scripts";                 
            }
            else {
                scriptsDir = path(binDir);
            }
            
            if (!exists(scriptsDir)) {
                BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                    "Could not find suitable directory containing KAT scripts at: ") + scriptsDir.c_str()));
            }
            
        }
        
        
        
        
        // **** Destructor ****
        virtual ~KatFS() { }
        
        path GetBinDir() const {
            return binDir;
        }

        path GetKatExe() const {
            return katExe;
        }

        path GetRootDir() const {
            return rootDir;
        }

        path GetScriptsDir() const {
            return scriptsDir;
        }

        
        friend std::ostream& operator<<(std::ostream &strm, const KatFS& pfs) {
            
            return strm << "Directories: "<< endl
                        << " - Root: " << pfs.rootDir << endl
                        << " - Bin: " << pfs.binDir << endl
                        << " - Scripts: " << pfs.scriptsDir << endl << endl                    
                        << "Executables: " << endl
                        << " - kat: " << pfs.katExe << endl;                        
        }     
    };
    
       
    
}

