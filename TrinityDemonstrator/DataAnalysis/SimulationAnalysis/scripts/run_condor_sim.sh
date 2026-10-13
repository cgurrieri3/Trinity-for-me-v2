#!/bin/bash
# run_condor_sim.sh
#
# One condor job = one run folder of a simulation campaign (e.g. job_000004 of
# 20260916_nathanwmycare). A campaign holds tens of those, each with its own
# GROPT/CARE files and its own ISOG file, so they are completely independent and
# a whole campaign no longer has to be walked on a single node.
#
# Condor drops the transferred inputs flat in the job's scratch directory, so the
# tree SimRootFiles expects is rebuilt here before it is run:
#
#   <scratch>/simdata/sofiadata/<DATE>_<SIMRUN>/<RUNFOLDER>/{GROPT|CARE,ISOG}/
#   <scratch>/DataAnalysis/SimulationAnalysis/scripts/SimRootFiles
#   <scratch>/DataAnalysis/SimulationAnalysis/data/<DATE>/<SIMTYPE>/   <- output

source /opt/root/bin/thisroot.sh

# Setup ExACT
export EXACT_DIR=/exact
export LD_LIBRARY_PATH=/exact/dict:/usr/lib/oracle/21/client64/lib:/usr/local/lib:$LD_LIBRARY_PATH

DATE=$1
SIMRUN=$2
SIMTYPE=$3
RUNFOLDER=$4

if [ -z "$RUNFOLDER" ]; then
    echo "Usage: run_condor_sim.sh <YYYYMMDD> <simrun> <GROPT|CARE> <runFolder>"
    exit 1
fi

# SimRootFiles takes the mount as a prefix and appends simdata/... and
# DataAnalysis/... to it, so it has to end in a slash.
WORKDIR=$(pwd)
echo "Running job for $DATE $SIMRUN ($SIMTYPE), run folder $RUNFOLDER"
echo "Scratch directory: $WORKDIR"
ls -lh

RUNPATH=simdata/sofiadata/${DATE}_${SIMRUN}/${RUNFOLDER}
mkdir -p ${RUNPATH}
mkdir -p DataAnalysis/SimulationAnalysis/scripts
mkdir -p DataAnalysis/SimulationAnalysis/data

# the stage folder and the ISOG folder arrive at the top of the scratch dir
mv ${SIMTYPE} ${RUNPATH}/
mv ISOG ${RUNPATH}/
mv SimRootFiles DataAnalysis/SimulationAnalysis/scripts/

echo "Rebuilt input tree:"
ls -lhR ${RUNPATH} | head -40

cd DataAnalysis/SimulationAnalysis/scripts/
chmod +x SimRootFiles
echo "Running SimRootFiles $DATE $SIMRUN ${WORKDIR}/ $SIMTYPE $RUNFOLDER"
./SimRootFiles "$DATE" "$SIMRUN" "${WORKDIR}/" "$SIMTYPE" "$RUNFOLDER"
STATUS=$?
echo "SimRootFiles exited with $STATUS"

cd "$WORKDIR"
echo "Output produced:"
ls -lhR DataAnalysis/SimulationAnalysis/data/${DATE}

# Every file is named with the run folder in it (datafiles_<date>_<simrun>_<runFolder>.root
# and so on), so the per job output directories merge into the shared
# data/<DATE>/<SIMTYPE>/ without any job overwriting another's results.
chmod -R 774 DataAnalysis/SimulationAnalysis/data/${DATE} 2>/dev/null

exit $STATUS
