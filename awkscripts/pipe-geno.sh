#!/bin/bash

# cat a binary .geno file through this to get the posterior genotype probabilities
# with number of individuals the argument

NINDS=$1

od --format=fD -w$((8*3*$NINDS)) -v - |  cut -f 2- -d ' '
