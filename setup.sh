#!/bin/bash

cd $(dirname $(readlink -f "$0")) &&
rm -rf build && cmake -E make_directory build &&
cd build && CXX=g++-10 cmake ..
