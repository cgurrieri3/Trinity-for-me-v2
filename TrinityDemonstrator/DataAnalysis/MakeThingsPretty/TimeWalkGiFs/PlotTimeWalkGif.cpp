#include "PlotTimeWalkGif.h"

std::vector<std::vector<std::string>> splitCSVLines(const std::vector<std::string>& lines) {
    std::vector<std::vector<std::string>> result;

    for (const auto& line : lines) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string item;

        while (std::getline(ss, item, ',')) {
            row.push_back(item);
        }

        result.push_back(row);
    }

    return result;
}

// First 8-digit run starting with "20" in the text. On a merged data path that is the night folder,
// which is also the night the flasher calibration is named after.
int extractYYYYMMDD(const std::string& text){
    for (size_t pos = 0; pos + 8 <= text.size(); pos++){
        if (text.compare(pos, 2, "20") != 0) continue;

        bool allDigits = true;
        for (size_t digit = pos; digit < pos + 8; digit++){
            if (!std::isdigit(static_cast<unsigned char>(text[digit]))){
                allDigits = false;
                break;
            }
        }

        if (allDigits){
            return std::atoi(text.substr(pos, 8).c_str());
        }
    }

    return 0;
}

// Test entries carry no event identifier, so each one is matched to its Sim (truth) entry from the
// layout of the file, the same rule EventCleaning uses:
//   - Test holds only the showers that triggered: Test entry k is the k-th triggered Sim entry.
//   - Test holds every thrown shower: the two trees are parallel, Test entry k is Sim entry k.
// Returns an empty vector when neither layout fits.
std::vector<Long64_t> MatchTestToSimEntries(TTree* testTree, TTree* simTree, SEvent* simEvent){
    std::vector<Long64_t> triggeredSimEntries;
    for (Long64_t simEntry = 0; simEntry < simTree->GetEntries(); simEntry++){
        simTree->GetEntry(simEntry);
        if (simEvent->GetTriggered()) triggeredSimEntries.push_back(simEntry);
    }

    const Long64_t nTest = testTree->GetEntries();
    std::vector<Long64_t> simEntryForTestEntry;
    if ((Long64_t)triggeredSimEntries.size() == nTest){
        simEntryForTestEntry = triggeredSimEntries;
    } else if (simTree->GetEntries() == nTest){
        for (Long64_t testEntry = 0; testEntry < nTest; testEntry++) simEntryForTestEntry.push_back(testEntry);
    } else {
        cout << "WARNING: Test entries match neither Sim layout; the gif gets no run name or sim event number" << endl;
    }
    return simEntryForTestEntry;
}


// ---- The cleaning's signal rules (copied from TimeWalkCleaning/EventCleaning.h, keep in step) ----

// Sample the pulse is expected near for the data type being cleaned.
int ExpectedPulseBin() { return (whatData == "sim") ? SimPulseBin : TimeBinAll; }

// Data: Pulse's pedestal (mean of samples 0-99). Sims: samples 157-166, since 0-155 are filler.
float PedestalForPixel(const std::vector<Int_t>& trace, float pulsePedestal)
{
    if (whatData != "sim") return pulsePedestal;

    const int start = std::max(0, SimPedestalStart);
    const int end   = std::min((int)trace.size(), SimPedestalStart + SimPedestalSamples);
    if (start >= end) return pulsePedestal;

    double sum = 0.0;
    for (int sample = start; sample < end; sample++) sum += trace[sample];
    return (float)(sum / (end - start));
}

// The pulse goes down from the pedestal, and a sample above the pedestal reads as zero.
float AmplitudeAtBin(const std::vector<Int_t>& trace, float pedestal, int bin)
{
    if (bin < 0 || bin >= (int)trace.size()) return 0.0;
    const float amplitude = pedestal - trace[bin];
    return (amplitude < 0.0) ? 0.0 : amplitude;
}

// Data rails at the fixed ADC floor; a sim pixel rails at its own deep plateau.
std::vector<float> GetRailLevels(const std::vector<std::vector<Int_t>>& traces, const std::vector<float>& pedestals,
                                 int windowStart, int windowEnd)
{
    std::vector<float> railLevels(traces.size(), (float)IUtilities::kSatADCFloor);
    if (whatData != "sim") return railLevels;

    for (std::size_t pixelID = 0; pixelID < traces.size(); pixelID++) {
        const std::vector<Int_t>& trace = traces[pixelID];
        const int start = std::max(0, windowStart);
        const int end   = std::min((int)trace.size(), windowEnd);
        if (start >= end) continue;

        const Int_t lowest = *std::min_element(trace.begin() + start, trace.begin() + end);

        // longest run of consecutive samples sitting on exactly that lowest value
        int longestRun = 0, run = 0;
        for (int sample = start; sample < end; sample++) {
            run = (trace[sample] == lowest) ? run + 1 : 0;
            if (run > longestRun) longestRun = run;
        }

        const bool hasPlateau = longestRun >= SimSaturationPlateauSamples;
        const bool isDeep     = (pedestals[pixelID] - lowest) >= SimSaturationMinDepth;
        railLevels[pixelID] = (hasPlateau && isDeep) ? (float)lowest : (float)IUtilities::kSatADCFloor;
    }
    return railLevels;
}

// The samples of each pixel at or below its rail level, inside the peak search window.
std::vector<IUtilities::SatSpan> FindSaturationSpans(const std::vector<std::vector<Int_t>>& traces,
                                                     const std::vector<float>& railLevels, int windowStart, int windowEnd)
{
    std::vector<IUtilities::SatSpan> spans;
    for (std::size_t pixelID = 0; pixelID < traces.size(); pixelID++) {
        const std::vector<Int_t>& trace = traces[pixelID];
        const int start = std::max(0, windowStart);
        const int end   = std::min((int)trace.size(), windowEnd);

        IUtilities::SatSpan span = {-1, -1, 0};
        for (int sample = start; sample < end; sample++) {
            if (trace[sample] <= railLevels[pixelID]) {
                if (span.count == 0) span.first = sample;
                span.last = sample;
                span.count++;
            }
        }
        spans.push_back(span);
    }
    return spans;
}

// A railed pixel's peak bin becomes the railed sample SaturationPeakMode picks.
void ApplySaturationPeakMode(std::vector<int>& peakTimeBin, const std::vector<IUtilities::SatSpan>& spans)
{
    for (std::size_t pixelID = 0; pixelID < peakTimeBin.size() && pixelID < spans.size(); pixelID++) {
        const IUtilities::SatSpan& span = spans[pixelID];
        if (span.count == 0) continue;

        switch (SaturationPeakMode) {
            case kSatFirst: peakTimeBin[pixelID] = span.first;                    break;
            case kSatAvg:   peakTimeBin[pixelID] = (span.first + span.last) / 2;  break;
            case kSatLast:  peakTimeBin[pixelID] = span.last;                     break;
        }
    }
}

// Earliest sample of the rising edge that is still at or above RiseFractionPercent of the pixel's max.
std::vector<int> FindRiseTimeBins(const std::vector<std::vector<Int_t>>& traces, const std::vector<float>& pedestals,
                                  const std::vector<float>& peakAmplitudes, const std::vector<int>& peakTimeBin)
{
    std::vector<int> riseTimeBin = peakTimeBin;
    if (RiseFractionPercent >= 100) return riseTimeBin;

    for (std::size_t pixelID = 0; pixelID < traces.size(); pixelID++) {
        if (peakAmplitudes[pixelID] <= 0) continue;

        const std::vector<Int_t>& trace = traces[pixelID];
        const float threshold   = peakAmplitudes[pixelID] * RiseFractionPercent / 100.0;
        const int   peakBin     = peakTimeBin[pixelID];
        const int   earliestBin = std::max(0, peakBin - MaxRiseLookBack);

        int riseBin = peakBin;
        while (riseBin - 1 >= earliestBin && pedestals[pixelID] - trace[riseBin - 1] >= threshold) {
            riseBin--;
        }
        riseTimeBin[pixelID] = riseBin;
    }
    return riseTimeBin;
}

// The same product as CEvent::SetAmplitudeValuesTimeBin, 1 / (flasher factor x SiPM gain), for flasher
// files CEvent cannot read. Only the 2025 single-CamFlasher files come through here.
std::vector<float> FactorsFromOldFlasherHistogram(TH2F* flasher, const std::vector<double>& gain)
{
    std::vector<float> factors(MaxNofChannels, 0.0);
    for (int pixelID = 0; pixelID < MaxNofChannels && pixelID < (int)gain.size(); pixelID++) {
        int nx, ny;
        IPlotTools::FindBin(pixelID, &nx, &ny);
        factors[pixelID] = 1.0 / (flasher->GetBinContent(nx + 1, ny + 1) * gain[pixelID]);
    }
    return factors;
}

// ADC -> PE factor per pixel, through the same ExACT call the cleaning uses (CEvent::SetAmplitudeValuesTimeBin):
// data divides by the night's flasher factor and the SiPM gain; sims use the "all1new" flasher file
// (-1 everywhere, hence the -1 input) and 24.1 ADC per PE. Empty when the flasher file is unusable.
std::vector<float> GetCalibrationFactors(int night, const std::vector<double>& gain)
{
    const std::string flasherName = (whatData == "sim") ? "all1new" : std::to_string(night);
    const std::string flasherPath = CalibrationFactorDir + flasherName + "_FlasherCalibration_Factor.root";

    // CEvent crashes on a missing file or histogram, so check both first and say which one it is.
    TFile flasherFile(flasherPath.c_str(), "READ");
    if (flasherFile.IsZombie()) {
        cerr << "Error: could not open flasher calibration " << flasherPath << endl;
        return std::vector<float>();
    }
    if (!flasherFile.Get("CamFlasher44")) {
        // 2025 nights such as 20250125 hold one histogram named CamFlasher. The current cleaning asks
        // ExACT for CamFlasher44 and stops on these nights, so the gif shows a walk it never ran.
        TH2F* oldFlasher = (TH2F*)flasherFile.Get("CamFlasher");
        if (whatData == "sim" || !oldFlasher) {
            cerr << "Error: no CamFlasher44 or CamFlasher calibration in " << flasherPath << endl;
            return std::vector<float>();
        }
        cout << "WARNING: " << flasherPath << " has only the old CamFlasher histogram. The gif uses it, but"
             << " TimeWalkCleaning stops on this night (ExACT asks for CamFlasher44), so there is no cleaning"
             << " output to compare this walk with." << endl;
        return FactorsFromOldFlasherHistogram(oldFlasher, gain);
    }
    flasherFile.Close();

    CEvent* calibration = new CEvent();
    if (whatData == "sim") {
        calibration->SetAmplitudeValuesTimeBin(std::vector<float>(MaxNofChannels, -1.0), CalibrationFactorDir,
                                               flasherName, std::vector<double>(MaxNofChannels, 24.1), 44.0);
    } else {
        calibration->SetAmplitudeValuesTimeBin(std::vector<float>(MaxNofChannels, 1.0), CalibrationFactorDir,
                                               flasherName, gain, 44.0);
    }
    std::vector<float> factors = calibration->GetAmplitudeValuesTimeBin();
    delete calibration;
    return factors;
}

// The trigger reports one MUSIC, but its mate can hold the brighter pixel. Use whichever is brighter.
int FindTriggeredMUSIC(const std::vector<float>& calibratedAmplitudes, int reportedMUSICID)
{
    int mateMUSICID = util->GetMUSICMate(reportedMUSICID, true);
    if (util->GetMaximum(calibratedAmplitudes, reportedMUSICID) < util->GetMaximum(calibratedAmplitudes, mateMUSICID)) {
        return mateMUSICID;
    }
    return reportedMUSICID;
}

// ---- end of the copied rules ----


float EventSignals::PEAtBin(int pixelID, int bin) const
{
    return AmplitudeAtBin(traces[pixelID], pedestals[pixelID], bin) * calibrationFactors[pixelID];
}

std::vector<float> EventSignals::CameraAtBin(int bin) const
{
    std::vector<float> pePerPixel(MaxNofChannels, 0.0);
    for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++) pePerPixel[pixelID] = PEAtBin(pixelID, bin);
    return pePerPixel;
}

// Pass 1 of the cleaning's event loop: every pixel's pedestal, peak bin and rise bin, the triggered
// pixel and the snapshot bin. calibrationFactors comes back empty when the flasher file is unusable.
EventSignals MeasureEvent(IEvent* ev, int night, const std::vector<double>& gain)
{
    EventSignals signals;
    signals.traces.resize(MaxNofChannels);
    const int peakSearchStart = ExpectedPulseBin() - PeakSearchBinsBefore;
    const int peakSearchEnd   = ExpectedPulseBin() + PeakSearchBinsAfter;

    for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++) {
        signals.traces[pixelID] = ev->GetSignalValue(pixelID);
        Pulse* pulse = new Pulse(signals.traces[pixelID], peakSearchStart, peakSearchEnd, 512, true);
        signals.pedestals.push_back(PedestalForPixel(signals.traces[pixelID], pulse->GetPedestal()));
        signals.peakTimeBin.push_back(pulse->GetTimePeak());
        delete pulse;
    }
    // the max used for the rise fraction is taken before the saturation rule moves the peak bin, as in the cleaning
    std::vector<float> peakAmplitudes;
    for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++) {
        peakAmplitudes.push_back(AmplitudeAtBin(signals.traces[pixelID], signals.pedestals[pixelID], signals.peakTimeBin[pixelID]));
    }

    std::vector<float> railLevels = GetRailLevels(signals.traces, signals.pedestals, peakSearchStart, peakSearchEnd);
    ApplySaturationPeakMode(signals.peakTimeBin, FindSaturationSpans(signals.traces, railLevels, peakSearchStart, peakSearchEnd));
    signals.riseTimeBin = FindRiseTimeBins(signals.traces, signals.pedestals, peakAmplitudes, signals.peakTimeBin);

    signals.calibrationFactors = GetCalibrationFactors(night, gain);
    if (signals.calibrationFactors.empty()) return signals;

    std::vector<float> calibratedPeakAmplitudes(MaxNofChannels);
    for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++) {
        calibratedPeakAmplitudes[pixelID] = peakAmplitudes[pixelID] * signals.calibrationFactors[pixelID];
    }
    signals.triggeredMUSICID = FindTriggeredMUSIC(calibratedPeakAmplitudes, ev->GetROIMusicID()[0]);
    signals.triggeredPixelID = util->GetMaximumPixelID(calibratedPeakAmplitudes, signals.triggeredMUSICID);
    signals.snapshotBin      = signals.peakTimeBin[signals.triggeredPixelID];
    return signals;
}

// The cleaning removes some events before the walk runs. The gif still shows the walk, but says so.
void WarnIfCleaningStopsBeforeWalk(const EventSignals& signals)
{
    const std::vector<float> snapshotPE = signals.CameraAtBin(signals.snapshotBin);

    double sumOfSquares = 0.0;
    for (float pe : snapshotPE) sumOfSquares += pe * pe;
    const float rms = sqrt(sumOfSquares) / MaxNofChannels;  // CEvent::SetRMS
    const float averagePE = IUtilities::GetEventAverageAmplitude(snapshotPE);
    if (rms / averagePE < 0.15) {
        cout << "WARNING: RMS / average amplitude " << Form("%.3f", rms / averagePE)
             << " < 0.15: the cleaning removes this event as a flasher (flag 9) before the walk" << endl;
    }

    const float triggeredPE = util->GetMaximum(snapshotPE, signals.triggeredMUSICID);
    const float triggeredCutPE = TriggeredChannelAmpCutOff / IUtilities::GetADCtoPEratio();
    if (triggeredPE < triggeredCutPE) {
        cout << "WARNING: triggered MUSIC maximum " << Form("%.2f", triggeredPE) << " PE < " << Form("%.2f", triggeredCutPE)
             << " PE: the cleaning removes this event (flag 1) before the walk" << endl;
    }
}


// ---- The time walk (GrowClusterInTime / RunTimeWalk of the cleaning), recording every step ----

// Neighbour lists from one of the neighbors/*.csv files, indexed by pixel ID.
const std::vector<std::vector<int>>& NeighborTable(const std::string& csvPath)
{
    static std::map<std::string, std::vector<std::vector<int>>> tablesRead;
    auto alreadyRead = tablesRead.find(csvPath);
    if (alreadyRead != tablesRead.end()) return alreadyRead->second;

    std::vector<std::vector<int>>& table = tablesRead[csvPath];
    table.assign(MaxNofChannels, std::vector<int>());
    std::ifstream csv(csvPath);
    if (!csv.is_open()) {
        cout << "WARNING: could not open neighbour file " << csvPath
             << ", the time walk will keep only the triggered pixel" << endl;
        return table;
    }
    std::string line;
    std::getline(csv, line); // header: Index,Neighbors
    while (std::getline(csv, line)) {
        // e.g.  3,"2,7,64"  ->  pixel 3 has neighbours 2, 7 and 64
        std::replace(line.begin(), line.end(), '"', ' ');
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int pixelID, neighborID;
        if (!(fields >> pixelID) || pixelID < 0 || pixelID >= MaxNofChannels) continue;
        while (fields >> neighborID) table[pixelID].push_back(neighborID);
    }
    return table;
}

// Parent time bin +/- TimeWalkWindowHalfWidth, cut down to time bin zero +/- MaxTimeDrift and to the trace.
BinWindow TimeWalkSearchWindow(int parentTimeBin, int timeZeroBin, int traceLength)
{
    BinWindow window = {parentTimeBin - TimeWalkWindowHalfWidth, parentTimeBin + TimeWalkWindowHalfWidth};
    if (MaxTimeDrift >= 0) {
        window.first = std::max(window.first, timeZeroBin - MaxTimeDrift);
        window.last  = std::min(window.last,  timeZeroBin + MaxTimeDrift);
    }
    window.first = std::max(window.first, 0);
    window.last  = std::min(window.last, traceLength - 1);
    return window;
}

// One stage of the walk, step for step as the cleaning takes it: every parent in pixelsToSearchFrom
// checks its neighbours that are not in the cluster yet, and a neighbour joins when its signal reaches
// the stage's threshold in any bin of the parent's window. Each parent is recorded as one WalkStep.
// The printout is the cleaning's own "event <N>" format, so the two can be diffed.
std::vector<int> GrowClusterInTime(std::vector<int>& clusterStage, std::queue<int> pixelsToSearchFrom, WalkStage stage,
                                   bool keepSearching, const EventSignals& signals, WalkRecord& record)
{
    const std::vector<std::vector<int>>& neighbors =
        NeighborTable(neighborDir + (TimeWalkUsesDiagonals ? "neighbor.csv" : "neighborbroadside.csv"));
    const int thresholdPE = (stage == kSeedStage) ? SeedThresholdPE : ExtensionThresholdPE;
    const int triggeredPixelID = signals.triggeredPixelID;
    const std::vector<int>& timeBin = signals.riseTimeBin;
    const int timeZeroBin = timeBin[triggeredPixelID];
    std::vector<int> addedPixels;
    std::vector<bool> rejectedOnTime(MaxNofChannels, false);

    while (!pixelsToSearchFrom.empty()) {
        WalkStep step;
        step.stage = stage;
        step.parentPixelID = pixelsToSearchFrom.front();
        pixelsToSearchFrom.pop();
        step.window = TimeWalkSearchWindow(timeBin[step.parentPixelID], timeZeroBin,
                                           (int)signals.traces[step.parentPixelID].size());
        step.clusterStageBefore = clusterStage;
        const BinWindow& window = step.window;

        for (int neighborID : neighbors[step.parentPixelID]) {
            if (neighborID < 0 || neighborID >= MaxNofChannels || clusterStage[neighborID] != 0) continue;

            // the neighbour's brightest bin inside the window
            NeighborCheck check = {neighborID, -1, 0.0, false};
            for (int bin = window.first; bin <= window.last; bin++) {
                const float signalPE = signals.PEAtBin(neighborID, bin);
                if (check.bestBin < 0 || signalPE > check.bestPE) {
                    check.bestBin = bin;
                    check.bestPE  = signalPE;
                }
            }
            check.joins = (check.bestBin >= 0 && check.bestPE >= thresholdPE);
            if (check.joins) {
                clusterStage[neighborID] = stage;
                addedPixels.push_back(neighborID);
                if (keepSearching) pixelsToSearchFrom.push(neighborID);
            }
            step.neighbors.push_back(check);

            const float ownPeakPE = signals.PEAtBin(neighborID, signals.peakTimeBin[neighborID]);
            const bool brightEnoughAtOwnPeak = ownPeakPE >= thresholdPE;
            if (!check.joins && brightEnoughAtOwnPeak) rejectedOnTime[neighborID] = true;
            if (!PrintWalkDiagnostics || !(check.joins || brightEnoughAtOwnPeak)) continue;

            cout << "    " << (check.joins ? "admitted" : "rejected") << "  pixel " << std::setw(3) << neighborID
                 << "  from " << (step.parentPixelID == triggeredPixelID ? "trigger " : "pixel ") << step.parentPixelID
                 << " (time bin " << timeBin[step.parentPixelID] << ")";
            if (window.first > window.last) {
                cout << "  no bins left inside max drift +/-" << MaxTimeDrift << endl;
                continue;
            }
            cout << "  window " << window.first << "-" << window.last
                 << "  best bin " << check.bestBin << " (" << std::showpos << check.bestBin - timeZeroBin << std::noshowpos << ")"
                 << "  " << Form("%.1f", check.bestPE) << " PE";
            if (!check.joins) {
                cout << " < " << thresholdPE << " PE, its own peak is " << Form("%.1f", ownPeakPE)
                     << " PE at bin " << signals.peakTimeBin[neighborID];
            }
            cout << "  own time bin " << timeBin[neighborID] << endl;
        }
        record.steps.push_back(step);
    }

    if (PrintWalkDiagnostics) {
        int rejectedCount = 0;
        for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++) {
            if (rejectedOnTime[pixelID] && clusterStage[pixelID] == 0) rejectedCount++;
        }
        cout << "    stage added " << addedPixels.size() << " pixels:";
        if (addedPixels.empty()) cout << " none";
        for (int pixelID : addedPixels) cout << " " << pixelID;
        cout << "   (" << rejectedCount << " pixels reach " << thresholdPE << " PE at their own peak but not in any window)" << endl;
    }
    return addedPixels;
}

// Both stages: the 8 PE seed walk from the triggered pixel, then the 4 PE extension walk from the
// triggered pixel and every seed pixel.
WalkRecord RunTimeWalk(const EventSignals& signals)
{
    WalkRecord record;
    std::vector<int> clusterStage(MaxNofChannels, 0);
    clusterStage[signals.triggeredPixelID] = kSeedStage;
    const int triggeredPixelID = signals.triggeredPixelID;

    if (PrintWalkDiagnostics) {
        cout << "--- Time walk (panel 2: " << SeedThresholdPE << " PE seed stage, panel 3: + "
             << ExtensionThresholdPE << " PE extension stage) ---" << endl;
        cout << "  Time bin zero: triggered pixel " << triggeredPixelID << " " << RiseFractionPercent
             << "% rise bin " << signals.riseTimeBin[triggeredPixelID] << " (peak bin " << signals.peakTimeBin[triggeredPixelID]
             << ", " << Form("%.1f", signals.PEAtBin(triggeredPixelID, signals.peakTimeBin[triggeredPixelID])) << " PE)" << endl;
        cout << "  Seed stage (>= " << SeedThresholdPE << " PE in the window, every new pixel searches on):" << endl;
    }
    std::queue<int> seedStart;
    seedStart.push(triggeredPixelID);
    std::vector<int> seedPixels = GrowClusterInTime(clusterStage, seedStart, kSeedStage, true, signals, record);

    if (PrintWalkDiagnostics) {
        cout << "  Extension stage (>= " << ExtensionThresholdPE << " PE in the window, "
             << (ExtensionStageIsRecursive ? "every new pixel searches on" : "neighbours of seed pixels only") << "):" << endl;
    }
    std::queue<int> extensionStart;
    extensionStart.push(triggeredPixelID);
    for (int pixelID : seedPixels) extensionStart.push(pixelID);
    GrowClusterInTime(clusterStage, extensionStart, kExtensionStage, ExtensionStageIsRecursive, signals, record);

    record.clusterStage = clusterStage;
    return record;
}


// ---- Horizon overlay (sims only), copied from TimeWalkCleaning/EventCleaning.h ----
const double kCameraFOVdeg = 3.87;
const double kCameraNPix   = 16.0;
const double kDegPerPixel  = kCameraFOVdeg / kCameraNPix;
const double kCameraCentre = (kCameraNPix - 1.0) / 2.0;   // row 7.5

// Tilt from the Tilt_<angle> folder in a sim file's path; 0 when there is none.
double TiltFromPath(const std::string& path)
{
    std::string::size_type tiltPos = path.find("Tilt_");
    if (tiltPos == std::string::npos) return 0.0;
    return atof(path.c_str() + tiltPos + 5);
}

double HorizonRow(double tiltDeg)
{
    return kCameraCentre - (90.0 - tiltDeg)/kDegPerPixel;
}

// Dashed horizon line plus a small arrow and "S" pointing to the sky.
void DrawHorizonOverlay(TH2* hist, double tiltDeg, Color_t color)
{
    const double xLow  = hist->GetXaxis()->GetXmin();
    const double xHigh = hist->GetXaxis()->GetXmax();
    const double yLow  = hist->GetYaxis()->GetXmin();
    const double yHigh = hist->GetYaxis()->GetXmax();
    const double ux = (xHigh - xLow)/kCameraNPix;
    const double uy = (yHigh - yLow)/kCameraNPix;

    double yHorizon = yLow + (HorizonRow(tiltDeg) + 0.5)*uy;
    bool inView = (yHorizon > yLow && yHorizon < yHigh);

    double yAnchor   = std::max(yLow, std::min(yHorizon, yHigh));
    double arrowX    = xLow + 0.7*ux;
    double arrowLen  = 1.2*uy;
    double arrowEnd  = yHigh - 0.2*uy;
    double skyFrom   = yAnchor + 0.2*uy;
    double skyTo     = skyFrom + arrowLen;
    if (skyTo > arrowEnd){
        skyTo   = arrowEnd;
        skyFrom = arrowEnd - arrowLen;
    }

    if (inView){
        TLine line;
        line.SetLineColor(color);
        line.SetLineWidth(3);
        line.SetLineStyle(2);
        line.DrawLine(xLow, yHorizon, xHigh, yHorizon);
    }

    TArrow arrow;
    arrow.SetLineColor(color);
    arrow.SetFillColor(color);
    arrow.SetLineWidth(2);
    arrow.DrawArrow(arrowX, skyFrom, arrowX, skyTo, 0.015, "|>");

    TLatex label;
    label.SetTextColor(color);
    label.SetTextSize(0.035);
    label.SetTextAlign(12);
    label.DrawLatex(arrowX + 0.45*ux, 0.5*(skyFrom + skyTo) + 0.25*uy, "S");
}


// ---- The gif ----

// What one frame shows on top of the camera image at timeBin.
struct FrameContent {
    int timeBin = 0;
    std::vector<int> clusterStage;           // pixels already in the cluster, outlined by stage
    int parentPixelID = -1;                  // the pixel whose neighbours are searched
    std::vector<int> searchedPixels;         // neighbours being searched, under the threshold in this bin
    std::vector<int> passingPixels;          // neighbours being searched, at or above the threshold in this bin
    std::vector<int> rejectedPixels;         // neighbours that failed (result frames)
    std::map<int, std::string> pixelText;    // number written on a pixel
    std::string stageLine;
    std::string stepLine;
    std::string timeLine;
    int delayCs = SweepFrameDelay;
};

// "Window  185  186  [187]  188  189   (time zero: bin 187)", with the bin shown highlighted.
std::string WindowTimeLine(const BinWindow& window, int shownBin, int timeZeroBin)
{
    std::string line = "Window bins: ";
    for (int bin = window.first; bin <= window.last; bin++) {
        if (bin == shownBin) line += Form("  #color[%d]{#bf{[%d]}}", (int)ParentPixelColor, bin);
        else                 line += Form("  %d", bin);
    }
    line += Form("      (shown: %+d from time zero, bin %d)", shownBin - timeZeroBin, timeZeroBin);
    return line;
}

// The pixel IDs of a list, comma separated, for the step lines.
std::string PixelList(const std::vector<int>& pixelIDs)
{
    if (pixelIDs.empty()) return "none";
    std::string list;
    for (std::size_t i = 0; i < pixelIDs.size(); i++) list += Form("%s%d", i ? ", " : "", pixelIDs[i]);
    return list;
}

// PE written on a pixel: one decimal near the thresholds, whole PE from 100 up so it fits in the pixel.
std::string PixelPEText(float pe)
{
    return (pe >= 100.0) ? Form("%.0f", pe) : Form("%.1f", pe);
}

int CountInCluster(const std::vector<int>& clusterStage, int stage)
{
    return (int)std::count(clusterStage.begin(), clusterStage.end(), stage);
}

// Turns the recorded walk into frames: time bin zero, then for each stage its title frame and, for
// every parent with neighbours left to search, one frame per window bin and a result frame, and
// last the finished cluster at the snapshot bin.
std::vector<FrameContent> PlanFrames(const EventSignals& signals, const WalkRecord& record)
{
    std::vector<FrameContent> frames;
    const int triggeredPixelID = signals.triggeredPixelID;
    const int timeZeroBin = signals.riseTimeBin[triggeredPixelID];

    FrameContent intro;
    intro.timeBin = timeZeroBin;
    intro.clusterStage.assign(MaxNofChannels, 0);
    intro.clusterStage[triggeredPixelID] = kSeedStage;
    intro.stageLine = "Time bin zero";
    intro.stepLine = Form("Triggered pixel %d: %d%% rise bin %d (peak bin %d, %.1f PE). The walk starts here.",
                          triggeredPixelID, RiseFractionPercent, timeZeroBin, signals.peakTimeBin[triggeredPixelID],
                          signals.PEAtBin(triggeredPixelID, signals.peakTimeBin[triggeredPixelID]));
    intro.timeLine = Form("Time bin %d (time zero)", timeZeroBin);
    intro.delayCs = IntroFrameDelay;
    frames.push_back(intro);

    int currentStage = 0;
    for (const WalkStep& step : record.steps) {
        const int thresholdPE = (step.stage == kSeedStage) ? SeedThresholdPE : ExtensionThresholdPE;
        const std::string stageName = (step.stage == kSeedStage)
            ? Form("Stage 1: %d PE seed walk (panel 2)", SeedThresholdPE)
            : Form("Stage 2: %d PE extension walk (panel 3)", ExtensionThresholdPE);

        if (step.stage != currentStage) {
            currentStage = step.stage;
            FrameContent title;
            title.timeBin = timeZeroBin;
            title.clusterStage = step.clusterStageBefore;
            title.stageLine = stageName;
            title.stepLine = (step.stage == kSeedStage)
                ? Form("From the triggered pixel; a neighbour joins at >= %d PE in any bin of its parent's time bin +/-%d",
                       SeedThresholdPE, TimeWalkWindowHalfWidth)
                : Form("Again from the triggered pixel and the %d seed pixels, now at >= %d PE",
                       CountInCluster(step.clusterStageBefore, kSeedStage) - 1, ExtensionThresholdPE);
            title.timeLine = Form("Time bin %d (time zero)", timeZeroBin);
            title.delayCs = IntroFrameDelay;
            frames.push_back(title);
        }

        // every neighbour of this parent is already in the cluster: nothing to show
        if (step.neighbors.empty()) continue;

        const int parentTimeBin = signals.riseTimeBin[step.parentPixelID];
        const bool windowIsEmpty = step.window.first > step.window.last;
        const std::string stepStart = Form("Pixel %d (time bin %d, %+d) searches %d neighbours",
                                           step.parentPixelID, parentTimeBin, parentTimeBin - timeZeroBin,
                                           (int)step.neighbors.size());

        if (ShowEveryWindowBin && !windowIsEmpty) {
            for (int bin = step.window.first; bin <= step.window.last; bin++) {
                FrameContent sweep;
                sweep.timeBin = bin;
                sweep.clusterStage = step.clusterStageBefore;
                sweep.parentPixelID = step.parentPixelID;
                for (const NeighborCheck& check : step.neighbors) {
                    const float pe = signals.PEAtBin(check.pixelID, bin);
                    if (pe >= thresholdPE) sweep.passingPixels.push_back(check.pixelID);
                    else                   sweep.searchedPixels.push_back(check.pixelID);
                    sweep.pixelText[check.pixelID] = PixelPEText(pe);
                }
                sweep.stageLine = stageName;
                sweep.stepLine = stepStart + Form(" for >= %d PE in bins %d-%d", thresholdPE, step.window.first, step.window.last);
                sweep.timeLine = WindowTimeLine(step.window, bin, timeZeroBin);
                sweep.delayCs = SweepFrameDelay;
                frames.push_back(sweep);
            }
        }

        FrameContent result;
        result.timeBin = windowIsEmpty ? parentTimeBin : std::max(step.window.first, std::min(parentTimeBin, step.window.last));
        result.clusterStage = step.clusterStageBefore;
        result.parentPixelID = step.parentPixelID;
        std::vector<int> joined;
        for (const NeighborCheck& check : step.neighbors) {
            if (check.joins) {
                result.clusterStage[check.pixelID] = step.stage;
                joined.push_back(check.pixelID);
            } else {
                result.rejectedPixels.push_back(check.pixelID);
            }
            if (!windowIsEmpty) result.pixelText[check.pixelID] = PixelPEText(check.bestPE);
        }
        result.stageLine = stageName;
        if (windowIsEmpty) {
            result.stepLine = stepStart + Form(": no bins left inside time zero +/-%d, none can join", MaxTimeDrift);
            result.timeLine = Form("Time bin %d (parent's time bin)", result.timeBin);
        } else {
            result.stepLine = stepStart + Form(" -> joined (>= %d PE): %s", thresholdPE, PixelList(joined).c_str());
            result.timeLine = Form("Result: numbers = brightest PE in bins %d-%d   (camera at bin %d)",
                                   step.window.first, step.window.last, result.timeBin);
        }
        result.delayCs = ResultFrameDelay;
        frames.push_back(result);
    }

    FrameContent finish;
    finish.timeBin = signals.snapshotBin;
    finish.clusterStage = record.clusterStage;
    const int seedCount = CountInCluster(record.clusterStage, kSeedStage);
    const int extensionCount = CountInCluster(record.clusterStage, kExtensionStage);
    finish.stageLine = "Walk finished";
    finish.stepLine = Form("Panel 2 (%d PE): %d pixels.   Panel 3: %d pixels (+%d at %d PE).",
                           SeedThresholdPE, seedCount, seedCount + extensionCount, extensionCount, ExtensionThresholdPE);
    finish.timeLine = Form("Time bin %d: the triggered pixel's peak, where the cleaning reads the images", signals.snapshotBin);
    finish.delayCs = FinalFrameDelay;
    frames.push_back(finish);
    return frames;
}


// Outline inside one pixel: a black line under the coloured one, so it reads on any palette colour.
void DrawPixelOutline(int pixelID, Color_t color, int lineWidth, int lineStyle = 1)
{
    int nx, ny;
    IPlotTools::FindBin(pixelID, &nx, &ny);
    const double inset = 0.08; // keeps neighbouring outlines from sitting on each other
    const double x1 = nx - 0.5 + inset, x2 = nx + 0.5 - inset;
    const double y1 = ny - 0.5 + inset, y2 = ny + 0.5 - inset;

    TBox outline;
    outline.SetFillStyle(0);
    outline.SetLineColor(kBlack);
    outline.SetLineWidth(lineWidth + 2);
    outline.DrawBox(x1, y1, x2, y2);
    outline.SetLineColor(color);
    outline.SetLineWidth(lineWidth);
    outline.SetLineStyle(lineStyle);
    outline.DrawBox(x1, y1, x2, y2);
}

// Text centred on a pixel, black on the bright end of the palette and white on the dark end.
void DrawPixelText(int pixelID, const std::string& text, float pixelPE, float maxPE, float textSize, bool bold)
{
    int nx, ny;
    IPlotTools::FindBin(pixelID, &nx, &ny);
    TLatex label;
    label.SetTextAlign(22);
    label.SetTextSize(textSize);
    label.SetTextFont(bold ? 62 : 42);
    label.SetTextColor(pixelPE > 0.55 * maxPE ? kBlack : kWhite);
    label.DrawLatex(nx, ny, text.c_str());
}

// One legend entry under the camera: a coloured square with a black edge, then its label (NDC).
void DrawLegendEntry(double ndcX, double ndcY, Color_t color, const std::string& text)
{
    const double padX1 = gPad->GetX1(), padX2 = gPad->GetX2();
    const double padY1 = gPad->GetY1(), padY2 = gPad->GetY2();
    auto userX = [&](double x) { return padX1 + x * (padX2 - padX1); };
    auto userY = [&](double y) { return padY1 + y * (padY2 - padY1); };

    TBox square;
    square.SetFillColor(color);
    square.SetFillStyle(1001);
    square.SetLineColor(kBlack);
    square.SetLineWidth(2);
    square.DrawBox(userX(ndcX), userY(ndcY - 0.009), userX(ndcX + 0.018), userY(ndcY + 0.009));
    TBox edge;
    edge.SetFillStyle(0);
    edge.SetLineColor(kBlack);
    edge.SetLineWidth(2);
    edge.DrawBox(userX(ndcX), userY(ndcY - 0.009), userX(ndcX + 0.018), userY(ndcY + 0.009));

    TLatex label;
    label.SetNDC();
    label.SetTextAlign(12);
    label.SetTextSize(0.019);
    label.DrawLatex(ndcX + 0.025, ndcY, text.c_str());
}

// Draws one frame on the canvas: camera at frame.timeBin, the outlines and numbers, the text lines.
void DrawFrame(TH2F& camera, const FrameContent& frame, const EventSignals& signals, float maxPE,
               const std::string& eventLabelText, double horizonTiltDeg)
{
    const std::vector<float> cameraPE = signals.CameraAtBin(frame.timeBin);
    for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++) {
        int nx, ny;
        IPlotTools::FindBin(pixelID, &nx, &ny);
        camera.SetBinContent(nx + 1, ny + 1, cameraPE[pixelID]);
    }

    canvas->Clear();
    canvas->cd();
    camera.Draw("COLZ");
    canvas->Update(); // paints the palette and sets the pad range the legend squares are placed in
    TPaletteAxis* palette = (TPaletteAxis*)camera.GetListOfFunctions()->FindObject("palette");
    if (palette) {
        palette->SetX1NDC(0.86);
        palette->SetX2NDC(0.89);
        palette->SetY1NDC(0.17);
        palette->SetY2NDC(0.86);
    }
    plottools->DrawMUSICBoundaries();
    if (whatData == "sim" && horizonTiltDeg != 0.0) DrawHorizonOverlay(&camera, horizonTiltDeg, kWhite);

    // cluster so far, by the stage each pixel joined in
    for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++) {
        if (frame.clusterStage[pixelID] == kSeedStage)      DrawPixelOutline(pixelID, SeedPixelColor, 3);
        if (frame.clusterStage[pixelID] == kExtensionStage) DrawPixelOutline(pixelID, ExtensionPixelColor, 3);
    }
    DrawPixelText(signals.triggeredPixelID, "T", cameraPE[signals.triggeredPixelID], maxPE, 0.03, true);

    // the pixels being worked on
    for (int pixelID : frame.searchedPixels) DrawPixelOutline(pixelID, SearchedPixelColor, 3);
    for (int pixelID : frame.passingPixels)  DrawPixelOutline(pixelID, PassingPixelColor, 4);
    for (int pixelID : frame.rejectedPixels) DrawPixelOutline(pixelID, RejectedPixelColor, 3, 2);
    if (frame.parentPixelID >= 0) DrawPixelOutline(frame.parentPixelID, ParentPixelColor, 5);
    for (const auto& text : frame.pixelText) {
        const bool passing = std::find(frame.passingPixels.begin(), frame.passingPixels.end(), text.first) != frame.passingPixels.end();
        DrawPixelText(text.first, text.second, cameraPE[text.first], maxPE, 0.018, passing);
    }

    TLatex line;
    line.SetNDC();
    line.SetTextAlign(22);
    line.SetTextSize(0.019);
    line.DrawLatex(0.5, 0.978, eventLabelText.c_str());
    line.SetTextSize(0.028);
    line.SetTextFont(62);
    line.DrawLatex(0.5, 0.948, frame.stageLine.c_str());
    line.SetTextFont(42);
    line.SetTextSize(0.021);
    line.DrawLatex(0.5, 0.916, frame.stepLine.c_str());
    line.DrawLatex(0.5, 0.886, frame.timeLine.c_str());

    DrawLegendEntry(0.05, 0.060, SeedPixelColor,      Form("%d PE stage (T = triggered)", SeedThresholdPE));
    DrawLegendEntry(0.37, 0.060, ExtensionPixelColor, Form("%d PE stage", ExtensionThresholdPE));
    DrawLegendEntry(0.58, 0.060, ParentPixelColor,    "pixel searching");
    DrawLegendEntry(0.05, 0.025, SearchedPixelColor,  "neighbour searched (PE in this bin)");
    DrawLegendEntry(0.45, 0.025, PassingPixelColor,   "reaches the threshold");
    DrawLegendEntry(0.75, 0.025, RejectedPixelColor,  "did not join");

    canvas->Modified();
    canvas->Update();
}

void MakeTimeWalkGif(const EventSignals& signals, const WalkRecord& record, const std::string& gifPath,
                     const std::string& eventLabelText, double horizonTiltDeg)
{
    const std::vector<FrameContent> frames = PlanFrames(signals, record);

    // brightest pixel in any frame sets the top of the colour scale, so frames compare
    float maxPE = 1.0;
    int totalDelayCs = 0;
    for (const FrameContent& frame : frames) {
        for (float pe : signals.CameraAtBin(frame.timeBin)) maxPE = std::max(maxPE, pe);
        totalDelayCs += frame.delayCs;
    }
    cout << "Gif: " << frames.size() << " frames, " << Form("%.1f", totalDelayCs / 100.0) << " s, colour scale up to "
         << Form("%.1f", maxPE) << " PE" << endl;

    // same axes, titles and colour scale floor as the cleaning's camera panels
    TH2F camera("walkCamera", ";Pixel Bin along x axis of Camera;Pixel Bin along y axis of Camera;PE",
                16, -0.5, 15.5, 16, -0.5, 15.5);
    camera.SetDirectory(0);
    camera.SetStats(0);
    camera.GetXaxis()->SetLabelSize(0.025);
    camera.GetYaxis()->SetLabelSize(0.025);
    camera.GetZaxis()->SetLabelSize(0.022);
    camera.GetXaxis()->SetTitleSize(0.028);
    camera.GetYaxis()->SetTitleSize(0.028);
    camera.GetZaxis()->SetTitleSize(0.028);
    camera.GetXaxis()->SetTitleOffset(1.2);
    camera.GetYaxis()->SetTitleOffset(1.4);
    camera.GetZaxis()->SetTitleOffset(1.3);
    camera.SetMinimum(-1);
    camera.SetMaximum(maxPE);

    canvas->cd();
    canvas->SetLeftMargin(0.11);
    canvas->SetRightMargin(0.15);
    canvas->SetTopMargin(0.14);
    canvas->SetBottomMargin(0.17);

    gSystem->Unlink(gifPath.c_str()); // a rerun would otherwise add its frames to the old gif
    for (std::size_t f = 0; f < frames.size(); f++) {
        DrawFrame(camera, frames[f], signals, maxPE, eventLabelText, horizonTiltDeg);
        // "gif+N" adds a frame shown for N cs; "gif++N" adds the last frame and makes the gif loop
        const bool lastFrame = (f == frames.size() - 1);
        canvas->Print(Form("%s%s%d", gifPath.c_str(), lastFrame ? "++" : "+", frames[f].delayCs));
    }
    canvas->Clear(); // the camera goes out of scope here, so take it off the canvas first
    cout << "Wrote " << gifPath << endl;
}


// Reads one event of one file, measures it the way the cleaning does, reruns the walk and makes the gif.
void MakeGifForRequest(const GifRequest& request)
{
    cout << "file path: " << request.filePath << endl;
    TFile *f = TFile::Open(request.filePath.c_str());
    if (!f || f->IsZombie()){
        cerr << "Error opening file: " << request.filePath << endl;
        return;
    }
    TTree *testTree = (TTree*)f->Get("Test");
    if (!testTree){
        cerr << "Error: no Test tree in " << request.filePath << endl;
        f->Close();
        return;
    }
    if (request.testEntry < 0 || request.testEntry >= testTree->GetEntries()){
        cerr << "Error: event " << request.testEntry << " is not in the Test tree (" << testTree->GetEntries() << " entries)" << endl;
        f->Close();
        return;
    }

    TTree *simTree = (TTree*)f->Get("Sim");
    whatData = simTree ? "sim" : "data";
    cout << "Data type: " << (simTree ? "simulation (Sim tree found)" : "demonstrator data (no Sim tree)") << endl;

    IEvent *ev = new IEvent();
    testTree->SetBranchAddress("Events", &ev);
    ISiPM *sipmInfo = new ISiPM();
    if (!simTree) testTree->SetBranchAddress("SiPM", &sipmInfo);
    testTree->GetEntry(request.testEntry);

    const int night = extractYYYYMMDD(request.filePath);
    const std::string fileName = request.filePath.substr(request.filePath.find_last_of('/') + 1);
    std::string eventLabelText;
    double horizonTiltDeg = 0.0;
    if (simTree){
        // the run name and sim event number tie the gif to the cleaning PDF, whose E# is the sim event number
        SEvent *simEvent = new SEvent();
        simTree->SetBranchAddress("GrOptics", &simEvent);
        std::vector<Long64_t> simEntryForTestEntry = MatchTestToSimEntries(testTree, simTree, simEvent);
        std::string run = "unknown run";
        int simEventNumber = -1;
        if (request.testEntry < (Long64_t)simEntryForTestEntry.size()){
            simTree->GetEntry(simEntryForTestEntry[request.testEntry]);
            run = simEvent->GetSumRun();
            simEventNumber = simEvent->GetSimEventNumber();
            horizonTiltDeg = simEvent->GetTilt();
        }
        if (horizonTiltDeg == 0.0) horizonTiltDeg = TiltFromPath(request.filePath);
        eventLabelText = Form("Run %s   E# %d (Test entry %lld)   tilt %.1f", run.c_str(), simEventNumber,
                              request.testEntry, horizonTiltDeg);
    } else {
        // the same N# / F# / E# the cleaning PDF titles carry
        std::string fileTime = fileName.substr(fileName.find("T") + 1, 5);
        eventLabelText = Form("N# %d   F# %s   E# %lld      %s", night, fileTime.c_str(), request.testEntry, fileName.c_str());
    }
    cout << eventLabelText << endl;

    const EventSignals signals = MeasureEvent(ev, night, simTree ? std::vector<double>() : sipmInfo->GetGain());
    if (signals.calibrationFactors.empty()){
        f->Close();
        return;
    }
    WarnIfCleaningStopsBeforeWalk(signals);

    const WalkRecord record = RunTimeWalk(signals);

    const std::string gifPath = Form("%s%s_E%lld_tw_%d_rf_%d.gif", outDir.c_str(), request.output.c_str(),
                                     request.testEntry, TimeWalkWindowHalfWidth, RiseFractionPercent);
    MakeTimeWalkGif(signals, record, gifPath, eventLabelText, horizonTiltDeg);
    f->Close();
}

// RISE_FRACTION=<percent> and TIME_WINDOW=<n> work as they do for the cleaning's run_script.sh.
bool ReadWholeNumberFromEnvironment(const std::string& variableName, int& value, int lowest, int highest, const std::string& label)
{
    const char* environmentValue = getenv(variableName.c_str());
    if (environmentValue == nullptr || *environmentValue == '\0') {
        cout << label << ": " << value << " (default; set " << variableName << " to change it)" << endl;
        return true;
    }
    char* end = nullptr;
    long requested = strtol(environmentValue, &end, 10);
    if (*end != '\0' || requested < lowest || requested > highest) {
        cout << variableName << "=" << environmentValue << " is not a whole number from " << lowest
             << " to " << highest << ", stopping" << endl;
        return false;
    }
    value = (int)requested;
    cout << label << ": " << value << " (from " << variableName << ")" << endl;
    return true;
}

int main(int argc, char **argv){
    if(argc < 3){
        cout << "Too few arguments; please include the gif list file and mount flag" << endl;
        cout << "  ./MakeTimeWalkGif /DataAnalysis/MakeThingsPretty/TimeWalkGiFs/TimeWalkGifs.txt /cedar/" << endl;
        return 1;
    }

    std::string ListPath = argv[1];
    std::string mount = argv[2];
    if (mount == "y"){
        mnt = "/mnt/";
    } else if (mount != "n"){
        mnt = mount;
    }
    outDir               = Form("%sDataAnalysis/MakeThingsPretty/TimeWalkGiFs/Gifs/", mnt.c_str());
    CalibrationFactorDir = Form("%sDataAnalysis/flasher_calibration/Output/", mnt.c_str());
    neighborDir          = Form("%sDataAnalysis/event_cleaning/TimeWalkCleaning/neighbors/", mnt.c_str());
    gSystem->mkdir(outDir.c_str(), true);
    cout << "Gifs go to: " << outDir << endl;

    if (!ReadWholeNumberFromEnvironment("RISE_FRACTION", RiseFractionPercent, 1, 100, "Rise fraction (%)")) return 1;
    if (!ReadWholeNumberFromEnvironment("TIME_WINDOW", TimeWalkWindowHalfWidth, 0, 20, "Time-walk window half-width (bins)")) return 1;

    util = new IUtilities();
    plottools = new IPlotTools();

    std::vector<std::string> listLines = util->readFileToVectorString(Form("%s%s", mnt.c_str(), ListPath.c_str()));
    if (listLines.size() < 2){
        cout << "The gif list needs a header line and at least one row: " << mnt << ListPath << endl;
        return 1;
    }
    auto rows = splitCSVLines(listLines);

    gErrorIgnoreLevel = kWarning; // ROOT prints an Info line for every gif frame otherwise
    canvas = new TCanvas("WalkDisplay", "Time walk", 900, 1000);

    for (int i = 1; i < (int)rows.size(); i++){ // row 0 is the header
        if (rows[i].size() < 3){
            if (!listLines[i].empty()) cerr << "Error: row needs Output,File,Event: " << listLines[i] << endl;
            continue;
        }
        GifRequest request;
        request.output    = rows[i][0];
        request.filePath  = rows[i][1];
        request.testEntry = std::atoll(rows[i][2].c_str());

        cout << "======== " << request.output << " ========" << endl;
        MakeGifForRequest(request);
    }

    cout << "Completed!" << endl;
}
