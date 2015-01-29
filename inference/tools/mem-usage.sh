#!/bin/bash
# from https://stat.ethz.ch/pipermail/r-sig-hpc/2013-June/001713.html

RBIN=$(which R) ## wherever the R you are timing lives

SCRIPT=$1
POST=$(date +"%H-%M_%m-%d-%Y")

rm free.RAM.txt

## Decrease the value after sleep if you are concerned about
## missing a peak

while true; do free | grep 'buffers/cache' | awk '{print $4}' >> free.RAM.txt; sleep 0.5; done &
FREE_RAM_PID=$!

$RBIN --vanilla < $SCRIPT > $SCRIPT.$POST.Rout 

kill $FREE_RAM_PID

TOT_RAM_USAGE=$(Rscript --vanilla -e 'tmp <- scan("free.RAM.txt"); usage <- (max(tmp) - min(tmp))/(1024^2); cat(usage)')
mv free.RAM.txt free.RAM.txt.$SCRIPT.$POST

echo
echo
echo $TOT_RAM_USAGE

echo "Total RAM usage = " $TOT_RAM_USAGE >> $SCRIPT.$POST.summary
