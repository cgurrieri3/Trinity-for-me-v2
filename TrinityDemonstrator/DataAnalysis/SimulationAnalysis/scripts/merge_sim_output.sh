#!/bin/bash
# merge_sim_output.sh
#
#   ./merge_sim_output.sh <YYYYMMDD> <simrun> [GROPT|CARE]
#
# A condor pass leaves one datafiles_<date>_<simrun>_<runFolder>.root per run
# folder in data/<date>/<simtype>/. Downstream analysis wants the single tree a
# local run used to produce, so this hadds them back together into
#   datafiles_<date>_<simrun>.root
# the same name a whole campaign pass writes. The per run folder files are left
# in place; they are the only way back to which folder an event came from.

DATE=$1
SIMRUN=$2
SIMTYPE=${3:-GROPT}

if [ -z "$DATE" ] || [ -z "$SIMRUN" ]; then
  echo "Usage: $0 <YYYYMMDD> <simrun> [GROPT|CARE]"
  exit 1
fi

BASE=${TDEM_BASE:-/storage/osg-otte1/shared/TrinityDemonstrator}
CONTAINER=${BASE}/DataAnalysis/containers/rootandexact.sif
DATA_DIR=${BASE}/DataAnalysis/SimulationAnalysis/data/${DATE}/${SIMTYPE}

if [ ! -d "$DATA_DIR" ]; then
  echo "No output directory: $DATA_DIR"
  exit 1
fi

# only the per run folder files, never a previous merge of them. Kept as an
# array and joined with spaces: the list is handed to a shell inside the
# container, so a newline separated list would be read as separate commands.
PARTS=("${DATA_DIR}"/datafiles_${DATE}_${SIMRUN}_*.root)
if [ ! -e "${PARTS[0]}" ]; then
  echo "No per run folder files matching datafiles_${DATE}_${SIMRUN}_*.root in $DATA_DIR"
  exit 1
fi
NPARTS=${#PARTS[@]}
PARTS_ARG="${PARTS[*]}"

MERGED=${DATA_DIR}/datafiles_${DATE}_${SIMRUN}.root
echo "Merging $NPARTS file(s) into $MERGED"

apptainer exec --bind "${BASE}":/mnt "$CONTAINER" bash -c \
  "export LD_LIBRARY_PATH=/exact/dict:\$LD_LIBRARY_PATH; hadd -f '$MERGED' $PARTS_ARG"
STATUS=$?

if [ $STATUS -eq 0 ]; then
  chmod 774 "$MERGED" 2>/dev/null
  echo "Merged output: $MERGED"
else
  echo "hadd failed with $STATUS"
fi
exit $STATUS
