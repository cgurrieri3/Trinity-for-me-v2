#!/bin/bash
# run_script.sh

echo "Opening container to compile code"
value=$1
# ExACT is baked read-only into the .sif; bind the git checkout over it so builds pick up
# the library rebuilt from source there instead of the image's older copy.
apptainer exec --bind /data/TrinityLabComputer/Trinity/exact:/exact,/projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator:/mnt /projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator/DataAnalysis/containers/rootandexact.sif make $value