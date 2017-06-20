#!/bin/sh

FILE_LIST=$1
OUTDIR=$2

cd $OUTDIR

while read FILENAME; do
        AFTERSLASHNAME=${FILENAME##*/}
        BASENAME=${AFTERSLASHNAME%.*}
        sbatch /lustre/nyx/hades/user/$USER/macro/MPD/start_flow.sh ${FILENAME} ${OUTDIR}/${BASENAME}_res_out.root     
done <$FILE_LIST


