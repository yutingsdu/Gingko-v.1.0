#!/bin/bash


## Please first dowload the demo data from https://sourceforge.net/projects/transassembly/files/Gingko/DemoData/

# !!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!

export LD_LIBRARY_PATH=/your/boost/dir/lib:$LD_LIBRARY_PATH

#export LD_LIBRARY_PATH=/yuting/local/boost/lib:$LD_LIBRARY_PATH

../Gingko -B bamlist -s first -o gingko_oudir -p 2
