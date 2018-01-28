#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    brew update
    brew install autoconf automake libtool

    # install anaconda
    if [[ "$PYTHON" == "YES" ]]; then
        wget https://repo.continuum.io/miniconda/Miniconda3-4.3.31-MacOSX-x86_64.sh -O miniconda.sh;
    fi

else
    sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
    sudo apt-get update -qq
    
    # Install Autotools
    sudo apt-get install -qq autoconf automake libtool
    
    # Install zlib
    sudo apt-get install -qq zlib1g zlib1g-dev
    
    if [ "$COMPILER" == "GCC5" ]; then
        sudo apt-get install -qq gcc-5 g++-5
        sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 100
        sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 100
        export CXX="g++-5"
        export CC="gcc-5"
    elif [ "$COMPILER" == "GCC6" ]; then
        sudo apt-get install -qq gcc-6 g++-6
        sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-6 100
        sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 100
        export CXX="g++-6"
        export CC="gcc-6"
    elif [ "$COMPILER" == "GCC7" ]; then
        sudo apt-get install -qq gcc-7 g++-7
        sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-7 100
        sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 100
        export CXX="g++-7"
        export CC="gcc-7"
    else
        sudo apt-get install -qq gcc-4.9 g++-4.9
        sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.9 100
        sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 100
        export CXX="g++-4.9"
        export CC="gcc-4.9"
    fi
    gcc --version
    g++ --version


    # Plotting installation
    if [[ "$PYTHON" == "YES" ]]; then
        wget https://repo.continuum.io/miniconda/Miniconda3-4.3.31-Linux-x86_64.sh -O miniconda.sh;
    fi

fi


if [[ "$PYTHON" == "YES" ]]; then

    bash miniconda.sh -b -p $HOME/miniconda;
    export PATH="$HOME/miniconda/bin:$PATH";
    hash -r;
    conda config --set always_yes yes --set changeps1 no;
    conda info -a
    conda create -q -n test-environment python=3.6 numpy scipy matplotlib sphinx;
fi
