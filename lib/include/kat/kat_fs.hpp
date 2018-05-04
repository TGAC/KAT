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
        } else if (exe.string().find('/') != string::npos) {

            // Relative with some parent paths... get absolute path, resolving
            // symlinks then take two levels up
            canonicalExe = fs::canonical(fs::system_complete(exe));
            isRelative = true;
        } else {

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

#ifdef KAT_SCRIPTS
        path kat_scripts(KAT_SCRIPTS);
#else
        path kat_scripts("");
#endif


#ifdef HAVE_PYTHON

        // If python is installed we need to figure out where the scripts are located relative to the
        // running executable.  This can be in various different places depending on how everything is
        // setup: installed kat, running compiled binary from source directory, or running unit tests.

        // First get the executable directory
        path exe_dir(canonicalExe.parent_path());

        // If the KAT_SCRIPTS variable is defined then we are either in an installed location, or running unit tests
        if (exe_dir.leaf().string() == "bin") {
            // Ok, so we are in a installed location.  Figuring out the scripts directory isn't as straight
            // forward as it may seem because we might have installed to a alternate root.  So wind back the
            // exec_prefix to get to root (or alternate root) directory.
            path kep(KAT_EXECPREFIX);
            path root = kep;
            path altroot = exe_dir.parent_path();
            while (root.has_parent_path()) {
                root = root.parent_path();
                altroot = altroot.parent_path();
            }
            this->scriptsDir = altroot / kat_scripts;
        } else if (exe_dir.leaf().string() == ".libs" && exists(exe_dir.parent_path() / "kat.cc")) {
            // If we are here then we are running the kat executable from the source directory but linked dynamically
            this->scriptsDir = exe_dir.parent_path().parent_path() / "scripts";
		} else if (exe_dir.leaf().string() == "src" && exists(exe_dir / "kat.cc")) {
			// If we are here then we are running the kat executable from the source directory but linked statically
			this->scriptsDir = exe_dir.parent_path() / "scripts";
		} else if (exe_dir.leaf().string() == "tests" && exists(exe_dir / "check_main.cc")) {
            // Presumably if we are here then we are running unit tests
            this->scriptsDir = exe_dir.parent_path() / "scripts";
        } else {
            // So if we got here then I'm not sure what config we are in.  So just use whatever scripts directory was provided.
            this->scriptsDir = kat_scripts;
        }

        // Validate the existence of the scripts directory and scripts file.
        if (!exists(this->scriptsDir)) {
            BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                                      "Could not find suitable directory containing KAT scripts at the expected location: ") + this->scriptsDir.string()));
        }

        path dascript = this->scriptsDir / "kat" / "distanalysis.py";
        if (!exists(dascript)) {
            BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                                      "Found a suitable KAT scripts directory but could not find distribution analysis script at: ") + dascript.string()));
        }

#endif	// HAVE_PYTHON

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
