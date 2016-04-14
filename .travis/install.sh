#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

	# install anaconda
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
	bash miniconda.sh -b -p $HOME/miniconda;
	export PATH="$HOME/miniconda/bin:$PATH";
	hash -r; 
	conda config --set always_yes yes --set changeps1 no; 
	conda update -q conda;  
	conda info -a; conda create -q -n test-environment python=3.5 anaconda; 
	source activate test-environment; 

	# Install boost
	brew install boost

else

	# Boost installation
	wget -q http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.gz/download
	mv download boost.tar.gz
	tar -xf boost.tar.gz
	cd boost_1_59_0
	sudo ./bootstrap.sh --with-libraries=chrono,timer,program_options,filesystem,system
	if [ "$COMPILER" == "GCC5" ]; then 
		sudo ./b2 -d0 --toolset=gcc-5 install; 
	else 
		sudo ./b2 -d0 --toolset=gcc-4.9 install; 
	fi
	cd ..


	# Plotting installation
	if [ "$PLOT" == "python" ]; then 
		wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
		bash miniconda.sh -b -p $HOME/miniconda;
		export PATH="$HOME/miniconda/bin:$PATH"; 
		hash -r; 
		conda config --set always_yes yes --set changeps1 no; 
		conda update -q conda;	
		conda info -a; conda create -q -n test-environment python=3.5 anaconda;	
		source activate test-environment; 
	elif [ "$PLOT" == "gnuplot" ]; then 
		sudo apt-get install gnuplot; gnuplot --version; 
	fi

fi

