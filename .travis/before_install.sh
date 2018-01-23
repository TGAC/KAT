#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    brew update
    brew install autoconf automake libtool
else
    sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
    sudo apt-get update -qq
    sudo apt-get install -qq libc6-dev
    if [ "$COMPILER" == "GCC5" ]; then 
        sudo apt-get install -qq gcc-5 g++-5
        sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 100
        sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 100
        export CXX="g++-5"
        export CC="gcc-5"
    else 
        sudo apt-get install -qq gcc-4.9 g++-4.9
        sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.9 100
        sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 100
        export CXX="g++-4.9"
        export CC="gcc-4.9"
    fi
    gcc --version
    g++ --version
    gcc --print-search-dirs
fi
