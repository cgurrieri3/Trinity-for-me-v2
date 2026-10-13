# Event cleaning: core pixels, saturation and rise-time timing — handoff

Started 2026-09-16, picking up from the multi-ns timing change made earlier that day (per-pixel peak
search, amplitude snapshot at the triggered pixel's peak bin, panel-3 time walk). "Changes, in order"
and "Results" below are that day's work, built and run on night 20250517.

A second batch was added **2026-09-22** — runtime diagnostics, the sim pedestal and sim saturation.
It has its own section near the end and was built and run on sim20260827.

## Where the code is

| What | Path |
|---|---|
| Cleaning driver | `/projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator/DataAnalysis/event_cleaning/ClusterCleaning/EventCleaning.{cpp,h}` |
| Panel-3 logic (`CEvent`) | `/data/TrinityLabComputer/Trinity/exact/dict/{src/CEvent.cpp,include/CEvent.h}` (the ExACT the scripts bind as `/exact`) |

Build (order matters, ExACT first):
```
apptainer exec --bind /data/TrinityLabComputer/Trinity/exact:/exact --pwd /exact \
  /projects/cos-lab-aotte6/cos-4a17118/Trinity/TrinityDemonstrator/DataAnalysis/containers/rootandexact.sif make
cd .../event_cleaning/ClusterCleaning && touch EventCleaning.cpp && ./compile_script.sh
./run_script.sh 20250517 n                 # quiet: per-file lines and the totals only
./run_script.sh sim20260827 n verbose      # per-event diagnostics for every event
./run_script.sh sim20260827 n event 18     # per-event diagnostics for one event
```
**Gotcha:** the ClusterCleaning Makefile does not rebuild when only `EventCleaning.h` changes (e.g. a
setting). `touch EventCleaning.cpp` first, or the old binary runs silently. Stop a running night
before rebuilding.

## Changes, in order

### 1. Second core pixel was deleted by the crosstalk step (bug fix)
- **Symptom:** 20250517 04:25 file, event 21: 9 pixels survived the time walk but the event was
  removed ("Removed Panel 3 Pixel amount", `0 , 0`).
- **Cause:** `CEvent::CheckCloseNeigborsOnly` keeps pixels that are neighbours of the triggered or
  core pixels, but a pixel is never in its own neighbour list. A second core not touching the
  triggered pixel (66 and 68) was set to 0, the time walk skipped it silently, and the final loop
  in `SetPanel3` then set **both** core amplitudes to 0. This predates the time walk.
- **Fix:** core pixels are always kept. The neighbour lists are now read once instead of once per pixel.

### 2. `SecondCoreAmp` was an int
`SetSecondCoreAmp`/`GetSecondCoreAmp` took and returned `int`, so the CoreRatio denominator was
truncated (event 21: 1.01938 instead of 1.01843). Now `float`. Outputs made before today carry this.

### 3. Clearer printouts (all under `PrintTimingDiagnostics`)
- Crosstalk step: which core pixels were kept, and which pixels were removed for not being next to a core.
- "Core pixel N is not in the panel-3 cluster -> ... set to 0" when a core is dropped.
- Removal messages now say why: `Removed Panel 3: too few pixels (1 < 3)`,
  `Removed Panel 3: size ...`, `Removed Panel 3: a core pixel was cleaned (MaxAmplitude ..., SecondCoreAmp ...)`.
- Removed the `YOU`, `0 , 0` and repeated `5` debug lines.

### 4. Triggered pixel marked on the timing map
The `TimeN...` canvas in the display ROOT file draws a grey box with "T" over the triggered pixel
(`DrawTimingHeatMaps`). The bare `TimeN...TH2F` histogram has no marker.
Since 2026-09-23 the canvas holds two maps side by side: panel-2 survivors (left) and panel-3
survivors (right), with the bare histograms saved as `TimePanel2N...TH2F` and `TimeN...TH2F`. Both
share one colour range, +/-MaxTimeDrift, widened for an event whose panel-2 pixels sit further out
(`TimingMapRange`), since COLZ leaves pixels below the minimum unpainted. `kTimeMapSource` is gone.

### 5. Saturated pixels count as cores
- **Symptom:** 04:59 file, event 3 had been saved before but now had "Second core pixel was not found".
- **Cause:** a railed pixel's calibrated amplitude is only its rail value times its calibration, and that
  differs pixel to pixel (137–164 PE here). A saturated trigger (23 at 163.6 PE) set a core cutoff of
  144.4 PE that the other railed pixels (89, 90, 91 at 139.6, 137.4, 143.9) could not reach.
  The old fixed bin 239 happened to catch 23 coming off the rail (155.6 PE), which lowered the cutoff.
- **Rule now:** a pixel whose trace is at or below `IUtilities::kSatADCFloor` at the snapshot bin is a
  core, whatever its amplitude (`CEvent::SetSaturatedAtSnapshot`, filled in EventCleaning.cpp).
  Printout: `Saturated at snapshot bin ... (count as cores)` and `(below cutoff, core because saturated at snapshot bin)`.

### 6. Saturation uses the first railed bin
`SaturationPeakMode = kSatFirst`, since the first railed sample is the last one that still follows the signal.
With the last railed bin, neighbouring saturated pixels got fake 2-bin gaps (90: 237–240 → 240 vs
91: 237–238 → 238), and the time walk rejected the second core. The panel amplitudes are read at the
triggered pixel's first railed bin.

### 7. Time walk uses 50% rise bins
- **Settings** (EventCleaning.h):
  ```cpp
  int RiseFractionPercent = 50;        // 50, 75, 90; 100 = peak bin (old behaviour)
  const int MaxRiseLookBack     = 20;  // max samples the rise search walks back from the peak
  ```
  Since 2026-09-23 the fraction can be picked per run without recompiling:
  `RISE_FRACTION=90 ./run_script.sh <date> <file>` (whole percent 1-100; anything else stops the run).
  The run prints `Rise fraction: 90% (from RISE_FRACTION)`, and the outputs carry `_rf_90`, so 50/75/90
  runs of the same night sit side by side.
- **Method:** `FindRiseTimeBins` starts at each pixel's peak (230–250 search, pedestal from `GetPedestal()`).
  It walks back while `pedestal - sample >= RiseFractionPercent% × max amplitude` and keeps the
  earliest such sample. Bins are whole samples.
- **What uses them:** the time walk, where each pixel's rise bin is compared with its neighbour's and
  with the triggered pixel's rise bin (`MaxTimeDrift`), and the timing heat map (rise bin − triggered rise bin).
- **Unchanged:** the triggered pixel is still picked by peak amplitude, and every panel still reads all
  pixels at the triggered pixel's peak bin (first railed bin if saturated).
- **Printout:** `Rise time (50% of max): triggered pixel ... -> rise bin ... (peak bin ...)`, plus
  `Rise bins of pixels passing the panel-2 cut (ID: rise/peak)`.
- Output names gained `_rf_50`.

### 8. Core pixels step to all 8 neighbours in the time walk (2026-09-23)
- **Rule:** in `CEvent::TimePropagatedNeighborhoodChecker` (ExACT), a core pixel (the triggered pixel or
  any of `GetCoreIDs()`) looks at all 8 neighbours from `neighbor.csv`, the same file the crosstalk step
  uses to define the boundary pixels around each core. Every other pixel still steps only to its 4
  broadside neighbours (`neighborbroadside.csv`). Step and drift tests are unchanged.
- **Why:** the crosstalk step deliberately keeps a core's diagonal neighbours, but the walk could only
  reach them through a broadside path, so they (or a diagonal second core) could be lost.
- **Printout:** time-walk lines now say `from core N` when the parent is a core, `from parent N` otherwise.
- **Checked on sim20260827, first 10 files:** no flag changed (17 kept, 9 flag 3). Three kept events gained
  pixels: 11414705 E33 9 -> 10, 11414707 E196 6 -> 10, 11414708 E145 4 -> 6. In E196, 151 was rejected
  on a +2 step from cores 155/150 and joins diagonally from core 216. The flag-3 core losses left in
  those files are all step > 1 rejections, not connectivity.
- **Not changed:** the connectivity check inside `isSecondCorePixel` and the `CountSameTimeBinSurvivors`
  diagnostic still use broadside neighbours (the first one's result is discarded anyway).

## Results on 20250517 (all 25,867 events)

| Run | Removed panel 2 | Panel 3 too few pixels | Panel 3 core cleaned | Reach PCA |
|---|---|---|---|---|
| Last railed bin, fixes 1–4 (`_sp_2`) | 5,033 | 2,809 | 0 | — |
| First railed bin + saturated cores (`_sp_0`) | 5,118 | 2,722 | 0 | 12 |
| + 50% rise bins (`_sp_0_rf_50`) | 5,118 | 2,723 | 1 | 10 |

The two events lost with rise bins:
- **04:41 E11:** a whole camera row (224–239) lit up together, which is not a shower.
- **05:19 E4:** saturated trigger 74 has rise bin 235, but unsaturated neighbour 73 has 237. That +2 step
  isolates the trigger. See open issue 1.

Checked events:
- **04:25 E21:** survives, but with 4 pixels (66, 65, 69, 68). It had 10 with the last railed bin.
- **04:59 E3:** survives with 8 pixels, all time-walk steps within ±1.

## Open issues / next steps
1. **Saturated pixels rise "early".** Their max is the rail, not the true peak, so 50% of the rail is reached
   before the true 50% point. Next to unsaturated neighbours this makes +2 steps (05:19 E4). Options:
   exclude saturated pixels from the step test against unsaturated ones, or allow ±2 when one side is saturated.
2. **Late neighbours are cut by amplitude, not timing.** 04:25 E21: pixels 64, 67, 70, 71, 73 peak at
   241–242 but are read at 239 (trigger's first railed bin), fall under 8 PE and are removed in panel 2.
   Rise bins cannot help, since timing is only used after panel 2.
3. ~~**Stale diagnostic.**~~ Done 2026-09-22: `Old fixed window [237,240] would have kept: N` is gone
   (with `SetOldTimeWindow` and `CountOldTimeWindowSurvivors`). In its place `CountSameTimeBinSurvivors`
   prints `No time steps: N pixels sit in the triggered pixel's own time bin (B), M of them connect to it`,
   which uses the same rise bins the walk runs on. `PrintRiseTimes` adds the peak-bin comparison:
   on sim20260827 E18 all 6 panel-2 pixels share peak bin 194 but only the trigger sits in rise bin 191,
   so every time step the walk took came from the rise-bin spread, not from the shower crossing the camera.
4. **Panel-2 cutoff rounded down.** `CEvent::SetPanel2` takes an `int`, so the 8.29875 PE cutoff is applied as 8 PE (pre-existing).
5. **Whole-row events.** A full row (e.g. 112–127, 224–239) lights up with a flat offset whose "peak" sits at
   the window start (230). All 17 look-back-limit hits in the rise search come from these.
6. **Muon runs** still use each pixel's own peak amplitude, not the snapshot (unchanged, not discussed further).
7. **Earlier fixes not yet compared at 100%.** The 100% check against the `_sp_0` output still needs to be run.
8. **Replace the per-pixel `Pulse` with a local peak search (speed only, DO NOT do this yet).**
   Deliberately deferred 2026-09-22: the sim pedestal and sim saturation changes made that day have to be
   shown to behave as expected first. Cutting time before the physics is settled would mean debugging two
   things at once. Details, so this does not have to be worked out again:

   `EventCleaning.cpp` builds one `Pulse` per pixel per event just to get the peak bin:
   ```cpp
   pulse = new Pulse(traces[k], PeakSearchStart, PeakSearchEnd, 512, true);
   Pedestals.push_back(PedestalForPixel(traces[k], pulse->GetPedestal()));
   PeakTimeBin.push_back(pulse->GetTimePeak());
   ```
   That constructor (`Pulse::Pulse(vector<Int_t> samples, Int_t tStart, Int_t tEnd, int, bool)`,
   `/data/TrinityLabComputer/Trinity/exact/dict/src/Pulse.cpp:75`) costs, per pixel per event:
   - **two full 512-int copies of the trace** — `samples` is taken BY VALUE and then `trace = samples`;
   - `CalcPedestal(100)`: two passes over samples 0-99, one for the mean and one for the RMS;
   - `FindPeakTime(tStart,tEnd)`: 20 samples;
   - `CalcCharge(tStart,tEnd,3,1)`: `isRevolving` is 1, so it runs `FindPeakTime` a SECOND time, then sums 3 samples.

   Roughly 1200 operations, of which only `GetTimePeak()` is still used:
   - for sims `GetPedestal()` is discarded (`PedestalForPixel` overrides it with samples 157-166, since
     samples 0-99 are `create_root_file_care` filler),
   - `GetCharge()` is never read at all,
   - `GetAmplitude()` was already replaced by `AmplitudeAtBin`.

   `GetTimePeak()` is just the deepest sample in `[PeakSearchStart, PeakSearchEnd)` — the same 20-sample
   window `GetRailLevels` already walks, which costs about 40 operations. So the replacement is a small
   local helper returning `argmin(trace, PeakSearchStart, PeakSearchEnd)`, matching `Pulse::FindPeakTime`'s
   tie-breaking (it keeps the FIRST deepest sample, and starts from `pkTmp = tWinSt` with `aCalc = 5000`).

   Do it for sims first and leave data on `Pulse` until a night is run both ways and the output trees
   diff clean; `PedestalForPixel`/`AmplitudeAtBin` already isolate the two paths. Expected saving is over
   an order of magnitude of the per-pixel work, and it also removes the 2x512-int copy, which is the part
   that shows up as I/O-shaped wall time rather than CPU.

## Slowdown (investigated, not fixed)
- **Run times on the same night:** 7.5 min before the rise bins, 20 min with them. But two runs of identical
  code also took 7.5 and 13 min, so the shared filesystem adds a lot of noise.
- **One-file timing (04:59, 127 events):** 31 s wall but only 4.1 s CPU, so most of the time is spent
  waiting on file I/O, not computing.
- **Output volume doubled** (51 MB → 94 MB for the night). Fixed 2026-09-22: `PrintTimingDiagnostics`
  and `PrintCutDiagnostics` are now runtime flags, both off by default, so a production run prints only
  per-file lines and the final totals (200 sim events -> 44 lines). Turn them back on per run with
  `./run_script.sh <date> <file> verbose`, or `... event <N>` for one event. The 256-number amplitude
  dump is gone, replaced by the ten brightest pixels.
- **Suspect I/O:** the flasher calibration ROOT file is now opened twice per event (`GetCalibrationFactors` plus
  `StoreCalibratedAmplitudes`, both through `CEvent::LoadFlasherCalibration`). Before the timing change it was once.
  Loading the factors once per night instead of per event would remove this. Not yet measured.

## 2026-09-22: diagnostics, sim pedestal, sim saturation

Built and run on `sim20260827` (200 thrown showers, 5 triggered). Nothing here changes a cut for
data: `PedestalForPixel`, `AmplitudeAtBin` and `GetRailLevels` all return the previous behaviour
unless `whatData == "sim"`, and the diagnostics are printing only.

### Diagnostics are now runtime, and off by default
`PrintTimingDiagnostics` and `PrintCutDiagnostics` were `const bool = true`; they are now plain
`bool`, both `false`, set by optional trailing arguments that `run_script.sh` forwards
(`ParseDiagnosticArguments`): `verbose` for every event, `event <N>` for one event (`FocusEventNumber`,
matched against the E# in the plot titles). A quiet run of the 200 sim events prints 44 lines.
This also closes the output-volume problem noted under "Slowdown" above.

To keep a quiet run actually quiet, prints that ignored the switch were gated too: the core search and
cluster dumps in `CEvent::isSecondCorePixel` / `SetPanel3` (both now take and honour `verbose`), the
PCA matrix dumps (`PCA::SetVerbose`, called once from `main`), and the cross-point and RMS lines in
`CheckCrossPoints`. The 256-number amplitude dump was replaced by `PrintBrightestPixels` (top ten).

### Every cut says why
`PrintCutVerdict` prints one line per cut — measured value, requirement, and the flag if it removed
the event — at all ten cut sites. Before this, the RMS-ratio, saturated-pixel, hot-spot, configure and
W/L cuts printed nothing at all, so an event could vanish with no explanation. `PrintSimTruth`,
`PrintTraceSummary`, `PrintPedestalSummary` and `PrintRailLevels` give the per-event context.

### Sim pedestal and sim saturation (EventCleaning only, no ExACT change)
- **Pedestal.** Sims took theirs from `Pulse::CalcPedestal(100)` = samples 0-99, which for a sim are
  the constant filler `create_root_file_care` writes outside the 200-sample CARE trace (zero noise).
  `PedestalForPixel` now averages samples `SimPedestalStart`=157 .. 166 (`SimPedestalSamples`=10)
  instead, and `AmplitudeAtBin` measures peak, snapshot and rise-time amplitudes against it, so all
  three stay consistent. Data still uses Pulse's own pedestal, byte for byte. E152 triggered pixel:
  3783 (filler) -> 3780.8. The window is deliberately short: at 20 samples it reached 176, already
  ~100 ADC into the rising edge, and gave 3774.4.
- **Saturation.** `IUtilities::kSatADCFloor`=1 never fires for sims, because the flip/shift is by each
  pixel's own CARE pedestal and the cap lands somewhere different per pixel (E152: -15, 35, 37 ADC).
  `GetRailLevels` marks a sim pixel railed when its window minimum repeats over 2+ consecutive samples
  AND sits `SimSaturationMinDepth`=2000 ADC below the pedestal. Both tests are needed: the plateau
  alone flagged 47 of 256 pixels on E152 (empty pixels repeat their noise minimum by chance) and fed
  eight zero-amplitude pixels to the core search. E152 now reports exactly the three capped pixels:
  `0: -15 ADC -> 180-181, 1: 35 ADC -> 178-184, 5: 37 ADC -> 179-184`. Raw traces for reference:
  pixel 1 holds 35 for seven consecutive samples and pixel 5 holds 37 for six, while pixel 6 (a real,
  uncapped peak) sits at its minimum of 408 for one sample only.

### Not verified yet
- Data is unchanged **by construction**, not by measurement: the three sim-only helpers return the old
  behaviour for any non-sim `whatData`, so `FindSaturationSpans` reduces to `IUtilities::GetSaturationSpans`
  and `AmplitudeAtBin` to `Pulse::CalcAmplitude`. A night run both ways, diffing the output trees, is
  still worth doing before this goes into production.
- `CorePixelCutOff` was changed from 1.5 to 3 outside this work. That alone moved sim20260827 from
  3 panel-3 removals / 2 survivors to 2 / 3, so do not read run-to-run differences as coming from the
  changes above.

## Backups (originals before each step)

2026-09-22, `event_cleaning/Claude/backups_20260922/`
- `before_diagnostics/`: EventCleaning.{cpp,h}, run_script.sh and CEvent.{cpp,h} as they were before
  the runtime diagnostics, the cut verdicts and the no-time-steps diagnostic
- `before_sim_pedestal/`: EventCleaning.{cpp,h} after those, before the sim pedestal and sim saturation

2026-09-16, `event_cleaning/Claude/backups_20260916/`
- `backup_20260916_corepixel/`: CEvent.{cpp,h} and EventCleaning.{cpp,h} before step 1, plus `EventCleaning.*.before_trigmarker`
- `backup_before_satcores/`: before steps 5 and 6
- `backup_before_risetime/`: EventCleaning.{cpp,h} before step 7

ExACT changes can also be seen with `git diff` in `/data/TrinityLabComputer/Trinity`. The run logs are
only in temporary scratch space (`/tmp/claude-3180259/.../scratchpad/run_20250517*.log`) and may be gone.
