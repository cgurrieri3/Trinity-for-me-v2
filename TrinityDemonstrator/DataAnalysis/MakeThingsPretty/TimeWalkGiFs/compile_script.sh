#!/bin/bash
# compile_script.sh

echo "Opening container to compile code"
value=$1
# Same ExACT bind as event_cleaning/ClusterCleaning/compile_script.sh: build against the live ExACT
# checkout instead of the older copy baked into the .sif.
apptainer exec --bind /data/TrinityLabComputer/Trinity/exact:/exact,/projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator:/cedar /projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator/DataAnalysis/containers/rootandexact.sif make $value
