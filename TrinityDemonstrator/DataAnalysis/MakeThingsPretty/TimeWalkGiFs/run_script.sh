#!/bin/bash
# run_script.sh

# ./run_script.sh /DataAnalysis/MakeThingsPretty/TimeWalkGiFs/TimeWalkGifs.txt
# The list path is inside TrinityDemonstrator, and the File column inside the list uses /cedar/ paths.
LIST=$1
echo "Making time-walk gifs from $LIST"
# Same ExACT bind as compile_script.sh: the binary loads libExACT.so from /exact/dict at runtime.
apptainer exec --bind /data/TrinityLabComputer/Trinity/exact:/exact,/projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator:/cedar /projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator/DataAnalysis/containers/rootandexact.sif /cedar/DataAnalysis/MakeThingsPretty/TimeWalkGiFs/MakeTimeWalkGif $LIST /cedar/
