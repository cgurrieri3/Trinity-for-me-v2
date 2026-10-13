#!/bin/bash
# run_script.sh
#
#   ./run_script.sh <YYYYMMDD> <simrun> [GROPT|CARE] [runFolder]
#
# Runs the whole campaign locally. The optional 4th argument restricts the pass
# to one run folder (e.g. job_000004), which is what submit_sim.sh hands to each
# condor job - useful here for reproducing a single failed job by hand.

DATE=$1
SIMRUN=$2
SIMTYPE=${3:-GROPT}   # GROPT (default) or CARE
RUNFOLDER=${4:-all}   # one run folder, or "all" for the whole campaign
echo "Running for $DATE ($SIMTYPE, run folder $RUNFOLDER)"
apptainer exec --bind /projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator:/mnt /projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator/DataAnalysis/containers/rootandexact.sif /mnt/DataAnalysis/SimulationAnalysis/scripts/SimRootFiles $DATE $SIMRUN y $SIMTYPE $RUNFOLDER
