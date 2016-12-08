#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    brew update
else
    sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
    sudo apt-get update -qq
    sudo apt-get install -qq libc6-dev
    if [ "$COMPILER" == "GCC5" ]; then 
        sudo apt-get install -qq gcc-5.3 g++-5.3
        sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5.3 100
        sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5.3 100
        export CXX="g++-5.3"
        export CC="gcc-5.3"
    else 
        sudo apt-get install -qq gcc-4.9 g++-4.9
        sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.9 100
        sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 100
        export CXX="g++-4.9"
        export CC="gcc-4.9"
    fi
    gcc --version
fi
