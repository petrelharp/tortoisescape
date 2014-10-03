#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:dl165:ppn=24
#PBS -l walltime=200:00:00
#PBS -l mem=48000mb
#PBS -l vmem=48000mb
#PBS -l pmem=2000mb

set -eu
set -o pipefail

if [ -e /home/rcf-40/pralph/cmb/lib/openblas-usc-R-setup.sh ]
then
    source /home/rcf-40/pralph/cmb/lib/openblas-usc-R-setup.sh
    cd $PBS_O_WORKDIR
fi

# get initial hitting times on 500x grid
#   makes 500x/six-raster-list-hitting-times.tsv
Rscript make-resistance-distances.R ../geolayers/TIFF/500x/500x_ 500x ../inference/six-raster-list simple-init-params-six-raster-list.tsv analytic

# push these up to 100x grid
Rscript disaggregate-ht.R ../geolayers/TIFF/500x/500x_ ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_ 500x 100x 500x/six-raster-list-hitting-times.tsv 5

# now use those to find hitting times on 100x grid
Rscript make-resistance-distances.R ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_ 100x ../inference/six-raster-list simple-init-params-six-raster-list.tsv CG 500x/six-raster-list-hitting-times.tsv 3600 100x/six-raster-list-hitting-times-iter01.tsv

