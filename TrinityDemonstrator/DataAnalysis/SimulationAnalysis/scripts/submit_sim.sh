#!/bin/bash
# submit_sim.sh
#
# Submits one condor job per run folder of a simulation campaign.
#
#   ./submit_sim.sh <YYYYMMDD> <simrun> [GROPT|CARE]
#   ./submit_sim.sh 20260916 nathanwmycare CARE
#
# The campaign simdata/sofiadata/<date>_<simrun>/ holds one run folder per ISOG
# batch (job_000004, job_000005, ...). Each is self contained, so this writes the
# list of them to condor_lists/ and lets condor run them side by side instead of
# walking tens of thousands of showers on one node.
#
# Set DRYRUN=1 to write the list and print the submit command without submitting.

DATE=$1
SIMRUN=$2
SIMTYPE=${3:-GROPT}

if [ -z "$DATE" ] || [ -z "$SIMRUN" ]; then
  echo "Usage: $0 <YYYYMMDD | date_file> <simrun> [GROPT|CARE]"
  exit 1
fi

if [ "$SIMTYPE" != "GROPT" ] && [ "$SIMTYPE" != "CARE" ]; then
  echo "Error: sim type must be GROPT or CARE, got '$SIMTYPE'"
  exit 1
fi

# The submit node sees the campaign under /storage; override for a different
# mount (the analysis box has the same tree under /projects).
BASE=${TDEM_BASE:-/projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator}
SUBMIT_TEMPLATE="plotcamera.submit"

SIM_DIR="${BASE}/simdata/sofiadata/${DATE}_${SIMRUN}"
LIST_DIR="${BASE}/DataAnalysis/SimulationAnalysis/condor_lists"
LOG_DIR="${BASE}/DataAnalysis/SimulationAnalysis/.logs/${USER}"

if [ ! -d "$SIM_DIR" ]; then
  echo "Campaign not found: $SIM_DIR"
  exit 1
fi

echo "Campaign:     $SIM_DIR"
echo "Sim type:     $SIMTYPE"

echo "Ensuring log directory exists: $LOG_DIR"
mkdir -p "$LOG_DIR"
chmod 775 "$LOG_DIR" 2>/dev/null
mkdir -p "$LIST_DIR"
chmod 775 "$LIST_DIR" 2>/dev/null

OUTLIST="${LIST_DIR}/run_folders_${DATE}_${SIMRUN}_${SIMTYPE}.txt"
rm -f "$OUTLIST"

# A run folder only earns a job if it has both halves of the work: the stage
# files to read and the ISOG file describing the same showers. Half a campaign
# is usually still running when the finished half is analysed, and a folder
# without one of the two would only be skipped on the node after paying for the
# file transfer.
NSKIP=0
for dir in "$SIM_DIR"/*/; do
  RUNFOLDER=$(basename "$dir")

  STAGE_TILT=$(find "${dir}${SIMTYPE}" -maxdepth 1 -type d -name 'Tilt_*' 2>/dev/null | head -1)
  if [ -z "$STAGE_TILT" ]; then
    echo "  $RUNFOLDER: no $SIMTYPE/Tilt_* folder, skipping"
    NSKIP=$((NSKIP+1))
    continue
  fi
  if [ -z "$(ls -A "$STAGE_TILT"/*.root 2>/dev/null)" ]; then
    echo "  $RUNFOLDER: $SIMTYPE folder has no .root files, skipping"
    NSKIP=$((NSKIP+1))
    continue
  fi
  if [ -z "$(find "${dir}ISOG" -name '*.txt' 2>/dev/null | head -1)" ]; then
    echo "  $RUNFOLDER: no ISOG .txt file, skipping"
    NSKIP=$((NSKIP+1))
    continue
  fi

  echo "$RUNFOLDER" >> "$OUTLIST"
done

if [ ! -s "$OUTLIST" ]; then
  echo "No complete run folders found for $DATE $SIMRUN $SIMTYPE - nothing to submit."
  exit 1
fi

NJOBS=$(wc -l < "$OUTLIST")
chmod 774 "$OUTLIST" 2>/dev/null
echo "---------------------------------------------"
echo "Wrote $OUTLIST"
echo "  $NJOBS run folder(s) to submit, $NSKIP skipped"

if [ "$DRYRUN" = "1" ]; then
  echo "DRYRUN=1, not submitting. Would run:"
  echo "  condor_submit Date=$DATE SimRun=$SIMRUN SimType=$SIMTYPE $SUBMIT_TEMPLATE"
  exit 0
fi

echo "Submitting $NJOBS jobs..."
condor_submit Date="$DATE" SimRun="$SIMRUN" SimType="$SIMTYPE" "$SUBMIT_TEMPLATE"

echo "All submissions complete."
