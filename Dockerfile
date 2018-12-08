FROM ubuntu:18.04

# Install Dependencies via apt and pip
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y autoconf automake libtool zlib1g zlib1g-dev gcc-7 g++-7 wget git make python3-dev python3-matplotlib python3-numpy python3-scipy python3-sphinx python3-setuptools python3-tabulate && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-7 100 && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 100 && export CXX="g++-7" && export CC="gcc-7"

# Install KAT
RUN git clone https://github.com/TGAC/KAT.git && cd KAT && ./build_boost.sh && ./autogen.sh && ./configure --with-sse && make -j4 V=1 && make install

