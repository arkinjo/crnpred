#!/bin/sh

cd src
make clean
make NDIM=5000 install
make clean
make NDIM=2000 install
make NDIM=2000 clean
cd ..
