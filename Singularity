Bootstrap: docker
From: ubuntu:latest

%post

	apt-get update

	# Install Common Dependencies via apt
	apt-get install -qq autoconf automake libtool zlib1g zlib1g-dev gcc-7 g++-7 wget git make python3.6 python3-pip

	# Configure GCC
        update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-7 100
        update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 100
        export CXX="g++-7"
        export CC="gcc-7"

	# Install python and required packages using miniconda
	pip3 install numpy scipy matplotlib sphinx tabulate

	# Install KAT
	git clone https://github.com/TGAC/KAT.git
	cd KAT
	./build_boost.sh
	./autogen.sh
	./configure --with-sse
	make -j4 V=1
	make install

	# Test
	kat --version

