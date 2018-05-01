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

#include <unistd.h>
#include <limits.h>
#include <iostream>
using std::endl;
using std::string;
using std::cout;

#ifdef OS_MAC
#include <mach-o/dyld.h>
#endif

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/exception/exception.hpp>
#include <boost/exception/info.hpp>
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
        KatFS() {}

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
#ifdef OS_LINUX
                canonicalExe = do_readlink();
#elif OS_MAC
                canonicalExe = get_mac_exe();
#else
                canonicalExe = do_readlink(); // Assume linux
#endif
                isOnPath = true;
            }


            // Check to see if scripts are adjacent to exe first
            path kda(canonicalExe.parent_path());
            if (kda.stem().string() == "bin") {
#ifdef HAVE_PYTHON
                // Ok, so we are in a installed location.  Wind back the eprefix to get to root (this may not be '/' if an alternate root is specified.)
				path kep(KAT_EPREFIX);
				path root = kep;
				path altroot = kda.parent_path();
				while (root.has_parent_path()) {
					root = root.parent_path();
					altroot = altroot.parent_path();					
				}
				
				// Looks like we are running from an installed location.  Don't try to use the
                // scripts from here.  We will try the KAT_SCRIPTS path instead.
                this->scriptsDir = altroot / KAT_SCRIPTS;
                if (!exists(this->scriptsDir)) {
	    	    	BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
		                        "Could not find KAT scripts at the expected installed location: ") + this->scriptsDir.c_str()));
                }
#else
                this->scriptsDir = canonicalExe.parent_path();
#endif
            }
            else {
                path kcc(canonicalExe.parent_path());
		kcc /= "kat.cc";
                if (exists(kcc)) {
                    // If we are here then we are not running from an installed location,
                    // we are running from the source tree.
                    // Not 100% sure how far back we need to go (depends on whether using KAT exe or tests)
                    // so try 2, 3 and 4 levels.
                    this->scriptsDir = canonicalExe.parent_path().parent_path();
                    this->scriptsDir /= "scripts";

                    if (!exists(this->scriptsDir)) {
                        this->scriptsDir = canonicalExe.parent_path().parent_path().parent_path();
                        this->scriptsDir /= "scripts";

                        if (!exists(this->scriptsDir)) {
                            this->scriptsDir = canonicalExe.parent_path().parent_path().parent_path().parent_path();
                            this->scriptsDir /= "scripts";

                            if (!exists(this->scriptsDir)) {
                                BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                                    "Could not find suitable directory containing KAT scripts relative to provided exe: ") + canonicalExe.c_str()));
                            }

                        }
                    }
                    kda = this->scriptsDir;
                    kda /= "setup.py";
                    if (!exists(kda)) {
                        BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                             "Could not find suitable directory containing KAT scripts derived from relative path of executable")));
                    }
                }
            }

        }


        std::string do_readlink() {
            char buff[PATH_MAX];
            ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
            if (len != -1) {
                buff[len] = '\0';
                return std::string(buff);
            }
            BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                     "Could not find locations of executable from /proc/self/exe")));

        }

#ifdef OS_MAC
        std::string get_mac_exe() {
            char path[1024];
            uint32_t size = sizeof(path);
            _NSGetExecutablePath(path, &size);
            return path;
        }
#endif


        // **** Destructor ****
        virtual ~KatFS() { }

        path GetCanonicalExe() const {
            return canonicalExe;
        }

        path GetScriptsDir() const {
            return scriptsDir;
        }

        bool IsAbsolute() const {
            return isAbsolute;
        }

        bool IsOnPath() const {
            return isOnPath;
        }

        bool IsRelative() const {
            return isRelative;
        }



        friend std::ostream& operator<<(std::ostream &strm, const KatFS& pfs) {

            return strm << "KAT paths: "<< endl
                        << " - argv: " << pfs.exe << endl
                        << "  - type: " << (pfs.isAbsolute ? "absolute" : pfs.isRelative ? "relative" : "on PATH") << endl
                        << " - Canonical path: " << pfs.canonicalExe << endl
                        << " - Scripts dir: " << (pfs.scriptsDir.empty() ? "assuming scripts on PATH" : pfs.scriptsDir) << endl;
        }

        /**
         * Ensures a directory exists
         * @param dir
         */
        static void ensureDirectoryExists(const path& dir) {

            path canDir = fs::absolute(dir);
            if (!fs::exists(canDir) || !fs::is_directory(canDir)) {
                if (!fs::create_directories(canDir)) {
                    if (!fs::exists(canDir) || !fs::is_directory(canDir)) { // Check again before throwing
                        BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                                "Could not create output directory: ") + canDir.string()));
                    }
                }
            }
        }
    };


    // Make available everywhere
    extern KatFS katFileSystem;
}
