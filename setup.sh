#!/bin/bash

cd $(dirname $(readlink -f "$0")) && rm -rf build && mkdir -p build && cd build && cmake ..
