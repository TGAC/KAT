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
