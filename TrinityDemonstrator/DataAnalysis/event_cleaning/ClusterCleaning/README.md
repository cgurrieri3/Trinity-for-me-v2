# ClusterCleaning — Trinity Demonstrator event cleaning

Cleans merged demonstrator data, background samples, muon runs and CARE simulations. Each event
is cut down to the pixels of one shower image, parametrised (Hillas/PCA), and written to ROOT and
PDF files. Every event gets a flag that says which cut removed it, or that it was kept.

---

## How to run

All commands are run from this folder
(`/projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator/DataAnalysis/event_cleaning/ClusterCleaning`).
Both scripts run inside the `rootandexact.sif` Apptainer container, so nothing needs to be
installed on the host.

### 1. Build

```bash
# a) ExACT: only needed when exact/ has changed (it has in this upload).
#    The scripts bind this copy over the older ExACT baked into the container.
apptainer exec --bind /data/TrinityLabComputer/Trinity/exact:/exact --pwd /exact \
  /projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator/DataAnalysis/containers/rootandexact.sif make

# b) EventInfo dictionary: needed whenever EventInfo.h changes (it has in this upload)
./compile_script.sh rootcling

# c) The cleaning itself
touch EventCleaning.cpp      # the Makefile does not notice header-only edits (e.g. a cut in EventCleaning.h)
./compile_script.sh
```

Stop any night that is still running before you rebuild, because it is using the binary.

### 2. Run

```bash
./run_script.sh <target> <file> [verbose | event <N>]
```

| `<target>` | What gets cleaned | Input location (under `TrinityDemonstrator/`) | Output folder (`DataAnalysis/event_cleaning/`) |
|---|---|---|---|
| `YYYYMMDD` | one night of data | `DataAnalysis/MergedData/Output/YYYYMMDD/` | `Output/` |
| `bkg<name>` | background sample from a file list | `DataAnalysis/DataTxtFiles/BackgroundSamples/RandomSampling/bkg<name>.txt` | `OutputBkg/` |
| `muon<name>` | muon run from a file list | `DataAnalysis/DataTxtFiles/Muon/muon<name>.txt` | `OutputSim/` |
| `sim<folder>` | a CARE simulation production | `simdata/sofiadata/<folder>/<run>/CARE/Tilt_*/*_tdem.root` | `OutputSim/` |

- **`<file>`**: `n` cleans every file. For data you can instead give one merged file name from
  that night's folder, e.g. `Merged_CoBo0_AsAd0_2025-05-17T02:32:07.828_0000.root`. Background,
  muon and sim runs always clean the whole list or production.
- **Diagnostics (optional)**: a plain run prints only a line per file and the totals.
  - `verbose` prints a per-event block for every event: sim truth, pedestal, rail levels, the
    timing walk, and a `cut passed` / `CUT FAILED` line for each cut, with the measured value and
    what the cut requires.
  - `event <N>` prints the same block for one event only. `N` is the `E#` shown in the plot titles:
    the sim event number for sims, and the Test-tree entry for data.
  - Neither option changes what is cleaned or written. Redirect the output to a file for a whole
    night, because it is large.
- **Rise fraction (optional)**: put `RISE_FRACTION=<1-100>` in front of the command to set which
  point on each pixel's rising edge is used as its time (default 50). A value that isn't a whole
  number from 1 to 100 stops the run.

Examples:

```bash
./run_script.sh 20250517 n                           # whole night, quiet
./run_script.sh 20250517 n event 21                  # whole night, diagnostics for event 21 only
./run_script.sh bkg10hours n                         # background list bkg10hours.txt
./run_script.sh muon20241011 n                       # muon list muon20241011.txt
./run_script.sh sim20260827_taudavid n verbose       # every sim event, full diagnostics
RISE_FRACTION=90 ./run_script.sh sim20260827_taudavid n
```

**Before cleaning data**, the night has to be merged (`DataAnalysis/MergedData`) and have a flasher
calibration in `DataAnalysis/flasher_calibration/Output/`. A night whose calibration is all NaN
gives NaN amplitudes for every pixel.

**On another machine**, change the bind paths in `compile_script.sh` and `run_script.sh`. They
point at the `cos-lab-aotte6` cluster: `/projects/...` is mounted as `/cedar`, and
`/data/TrinityLabComputer/Trinity/exact` is mounted as `/exact`.

### 3. Outputs

Every file name carries the cut settings, so runs with different settings never overwrite each other:

```
EventCleanedCluster<target>_TC_200_TB_239_NP_200_s_256_FA_800_mp_3_er_100_tr_1_ts_1_md_4_sp_0_rf_50.{root,pdf}
DataFiles/Data_EventCleaning<target>_<same tags>.root
```

| Tag | Setting (in `EventCleaning.h`) | Tag | Setting |
|---|---|---|---|
| `TC` | triggered-pixel amplitude cut (ADC) | `er` | width/length cut ×100 |
| `TB` | expected pulse bin: 239 data, 184 sims | `tr` | remove top row (1 = keep) |
| `NP` | panel-2 pixel amplitude cut (ADC) | `ts` | time-walk step allowance (bins) |
| `s` | saturated-pixel cut | `md` | max drift from the triggered pixel (bins) |
| `FA` | flasher average-amplitude cut | `sp` | saturated peak choice (0 = first railed sample) |
| `mp` | minimum surviving pixels | `rf` | rise fraction (%) |

- **`EventCleanedCluster…root`**: event displays for kept events. `N…` is the 2×2 canvas, `PreN…TH2F`
  the calibrated image, `PostN…TH2F` the cleaned image, and `TimeN…` / `TimePanel2N…` the timing
  maps (rise bin minus the triggered pixel's rise bin).
- **`EventCleanedCluster…pdf`**: the same displays plus the summary plots.
- **`DataFiles/Data_EventCleaning…root`**: tree `EventCleaning`, branch `Cleaned` (`EventInfo`), one
  entry per event including the removed ones. Sims also get tree `Sim` (branch `GrOptics`), which is
  the input truth tree copied as-is.

Event flags stored in `EventInfo`:

| Flag | Meaning |
|---|---|
| 1 | too many saturated pixels, triggered pixel too dim, or a known hot spot |
| 2 | too few pixels above the panel-2 amplitude cut |
| 3 | panel 3: configure event, too few pixels after the time walk, too small, or a core pixel was cleaned away |
| 4 | too round (width/length above the cut) |
| 5 | **kept**: passed every cut |
| 9 | evenly lit camera (flasher/HLED) |

---

## Changes in this upload (2026-09-27)

These are relative to GitHub commit `baec777` (2026-09-13), the last commit that touched this folder.
The work was done between 2026-09-07 and 2026-09-23. Outputs made with the previous version are not
comparable to new ones. The file-name tags changed, so the old outputs are not overwritten.

### Timing: showers take several ns to cross the camera

- **Per-pixel peak search.** Every pixel's peak is searched in a window around the expected pulse:
  [230, 250) for data, HLED and muons, and [175, 195) for sims.
- **Amplitude snapshot.** The triggered pixel is picked from the calibrated per-pixel peak
  amplitudes. If the reported MUSIC's mate holds the brighter pixel, the mate is used instead. The
  whole camera is then read at *that pixel's* peak bin. Before, data was read at fixed bin 239. Muon
  amplitudes are still each pixel's own peak.
- **Pixel time is the rise bin.** A pixel's time is the first sample of its rising edge at or above
  `RISE_FRACTION`% (default 50) of its maximum, looking back at most 20 samples. The rising edge comes
  before saturation and crosstalk, so it is more reliable than the peak.
- **Panel-3 time walk replaces the fixed 237–240 window.** The cluster grows outward from the triggered
  pixel. A neighbour is admitted when its time is within ±1 bin of the pixel it was reached from
  (`TimeStepAllowance`) and within ±4 bins of the triggered pixel (`MaxTimeDrift`, -1 turns it off).
  Core pixels look at all 8 neighbours; other pixels look broadside only.

### Saturation

- Railed samples are found in the same peak-search window. A saturated pixel's peak bin is its
  **first** railed sample (`SaturationPeakMode = kSatFirst`), the last sample that still follows the signal.
- Data rails at the ADC floor (≤ 1 ADC). A sim pixel rails at its own level, marked by a flat bottom
  of at least 2 samples that is at least 2000 ADC below its pedestal.
- Pixels railed at the snapshot bin always count as **core pixels**. Their amplitude is only the rail
  value, so it can't be fairly compared with the triggered pixel's.

### Simulations

- Sims are now read from `simdata/sofiadata/<folder>/<run>/CARE/Tilt_<angle>/*_tdem.root`. Before,
  they were read from `SimulationAnalysis/data/<date>/CARE/`.
- Sim traces are no longer shifted by 55 bins (`ShiftTrace` removed). The peak search is centred on
  where the CARE pulse actually sits (bin 184).
- Sim pedestals come from CARE samples 157–166. Samples 0–99, which data uses, are constant filler in sims.
- Truth matching handles both layouts of the Test tree: "triggered showers only" maps to the k-th
  triggered Sim entry, and "every thrown shower" maps one-to-one. An event with no match gets event
  number -1 and run `noSimTruth` instead of silently reusing the previous event's truth. A warning is
  printed when several events share a (run, event number).
- The flasher cut (flag 9) now applies to sims too.
- Sim event displays draw the horizon and a sky arrow, using the tilt from the SEvent or the `Tilt_` folder.

### Cut values

| Setting | Before | Now |
|---|---|---|
| `TriggeredChannelAmpCutOff` | 481 | 200 |
| `CorePixelAmpCutOff` | 481 | 200 |
| `CorePixelCutOff` | 1.5 | 3 |

### Bug fixes (in ExACT `CEvent`)

- **Second core pixel not next to the triggered pixel.** The crosstalk step zeroed any second core
  pixel that was not next to the triggered pixel. That in turn zeroed both core amplitudes and
  removed the event with flag 3. Core pixels are now always kept in that step.
- **`SecondCoreAmp` was an `int`.** It is now a `float`. `CoreRatio` had been computed with a truncated
  second-core amplitude.

### New outputs

- `EventInfo` has a new field, `OnCameraEdge`: true when the cleaned image touches the outer ring
  (row/column 0 or 15). It is only filled for kept events (flag 5). `ClassDef` went from 5 to 6.
- Timing maps (`TimeN…`, `TimePanel2N…`) are added to the display ROOT file.
- Panel 1 shows the image time bin (and the CARE sample number for sims). Panel 4 shows `Edge:`.
- File names gain the `_ts_ _md_ _sp_ _rf_` tags, and `_TB_` is now the expected pulse bin.

### Diagnostics and scripts

- The run is quiet by default. Use `verbose` or `event <N>` (see *How to run*). The PCA matrix dumps
  follow the same switch.
- `run_script.sh` passes extra arguments on to the binary. Both scripts bind the rebuilt ExACT as `/exact`.

---

## Files in this upload

**This folder** (`TrinityDemonstrator/DataAnalysis/event_cleaning/ClusterCleaning/`):

| File | Status |
|---|---|
| `EventCleaning.cpp` | modified |
| `EventCleaning.h` | modified |
| `EventInfo.cpp`, `EventInfo.h` | modified (`OnCameraEdge`, ClassDef 6) |
| `PCA.cpp`, `PCA.h` | modified (verbose switch) |
| `compile_script.sh`, `run_script.sh` | modified |
| `EventCleaning` | rebuilt binary |
| `README.md` | new (this file); replaces the old `README` |
| `EventInfoDict.cpp`, `PlotHelp.cpp`, `PlotHelp.h`, `Makefile`, `LinkDef.h`, `neighbors/` | unchanged |

**ExACT** (`exact/`). The cleaning will not compile against the version of ExACT on GitHub without these:

| File | Change |
|---|---|
| `exact/dict/include/CEvent.h`, `exact/dict/src/CEvent.cpp` | time walk, saturated cores, crosstalk fix, `SecondCoreAmp` float, verbose output |
| `exact/dict/include/IUtilities.h`, `exact/dict/src/IUtilities.cpp` | `SatSpan`, `GetSaturationSpans`, `kSatADCFloor` |
| `exact/dict/libExACT.so`, `exact/dict/exactDict.cpp`, `exact/dict/exactDict_rdict.pcm`, `exact/ExACT`, `exact/include/libIExACT.a`, `exact/dict/obj/*.o`, `exact/obj/Display.o` | rebuilt build outputs |

**Do not upload**: `*.o` in this folder (git-ignored), `EventInfoDict_rdict.pcm` (made by
`./compile_script.sh rootcling`), and the `.nfs*` files (leftover NFS handles from old binaries).
