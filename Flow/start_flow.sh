#!/bin/bash

INFILE=$1
OUTFILE=$2

#SBATCH --time=0:15:00
#SBATCH -D /tmp

PROJECT_DIR=/lustre/nyx/hades/user/$USER/macro/

. /lustre/nyx/hades/user/$USER/Soft/MPDRoot/build/config.sh

cd $PROJECT_DIR

OUTDIR=${OUTFILE%/*}
AFTERSLASHNAME=${INFILE##*/}
BASENAME=${AFTERSLASHNAME%.*}

root -l -b -q "MPD/start_flow.C(\"${INFILE}\",\"${OUTFILE}\")" 1>> ${OUTDIR}/${BASENAME}_rec_calc.OUT 2>> ${OUTDIR}/${BASENAME}_rec_calc.ERR


mv BASENAME_res_out.root OUTDIR/.

