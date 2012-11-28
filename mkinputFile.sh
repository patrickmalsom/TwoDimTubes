#!/bin/bash

# delete the mean and the covariance matrix bit if it exists in the input file
cat inputPos.txt |awk {'print $1 " " $2'} >> TempinputPos.txt
rm inputPos.txt
mv TempinputPos.txt inputPos.txt
