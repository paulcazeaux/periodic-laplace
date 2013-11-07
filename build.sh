mkdir build
cd build

# cmake parameters
 export CC=/usr/local/bin/gcc-4.8
 export CXX=/usr/local/bin/g++-4.8

cmake ..
make install
