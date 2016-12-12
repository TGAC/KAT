#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    # Install boost -- looks like boost 1.55.0_2 is already installed with brew
    #brew install boost


    # install anaconda
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
    bash miniconda.sh -b -p $HOME/miniconda;
    export PATH="$HOME/miniconda/bin:$PATH";
    hash -r; 
    conda config --set always_yes yes --set changeps1 no; 
    conda update -q conda;  
    conda info -a
    conda create -q -n test-environment python=3.5 anaconda; 

else

    # Boost installation
    wget -q http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.gz/download
    mv download boost.tar.gz
    tar -xf boost.tar.gz
    cd boost_1_59_0
    sudo ./bootstrap.sh --with-libraries=chrono,timer,program_options,filesystem,system
    if [[ "$COMPILER" == "GCC5" ]]; then 
        export CXX="g++-5"
        export CC="gcc-5"
        sudo ./b2 -d0 --toolset=gcc-5 install; 
    else 
        export CXX="g++-4.9"
        export CC="gcc-4.9"
        sudo ./b2 -d0 --toolset=gcc-4.9 install; 
    fi
    cd ..


    # Plotting installation
    if [[ "$PLOT" == "python" ]]; then 
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
        bash miniconda.sh -b -p $HOME/miniconda;
        export PATH="$HOME/miniconda/bin:$PATH"; 
        hash -r; 
        conda config --set always_yes yes --set changeps1 no; 
        conda update -q conda;	
        conda info -a
        conda create -q -n test-environment python=3.5 anaconda;	
    elif [ "$PLOT" == "gnuplot" ]; then 
        sudo apt-get install gnuplot
        gnuplot --version; 
    fi

fi

