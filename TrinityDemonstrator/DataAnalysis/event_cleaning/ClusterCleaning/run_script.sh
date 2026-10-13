#!/bin/bash
# run_script.sh

DATE=$1
SIMRUN=$2
shift $(( $# < 2 ? $# : 2 ))   # anything left over is forwarded to the binary, e.g. "event 18" to focus the diagnostics
# RISE_FRACTION=<percent> in front picks the rise fraction for that run (default 50), e.g.
#   RISE_FRACTION=90 ./run_script.sh sim20260827_taudavid n
# apptainer passes it through to the binary, which prints the value it used.
echo "Running for $DATE"
# Same ExACT bind as compile_script.sh: the binary resolves libExACT.so from /exact/dict at
# runtime, so it has to see the same rebuilt library it was compiled against.
apptainer exec --bind /data/TrinityLabComputer/Trinity/exact:/exact,/data/TrinityLabComputer/TrinityDemonstrator:/mnt,/projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator:/cedar /projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator/DataAnalysis/containers/rootandexact.sif /cedar/DataAnalysis/event_cleaning/ClusterCleaning/EventCleaning $DATE /cedar/ $SIMRUN "$@"
