#!/bin/bash

# this only needs to be run when using a new .dat file for use as an input file
# delete the mean and the covariance matrix bit if it exists in the input file
cat $1 |awk {'print $1 " " $2'} >> temp$1.dat
mv temp$1.dat $1
