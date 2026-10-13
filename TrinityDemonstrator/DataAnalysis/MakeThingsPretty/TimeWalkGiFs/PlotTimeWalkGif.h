#include <TH2.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TBox.h>
#include <TLine.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TPaletteAxis.h>
#include <IUtilities.h>
#include <IPlotTools.h>
#include <IEvent.h>
#include <ISiPM.h>
#include <CEvent.h>
#include <Pulse.h>
#include <SEvent.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

// Makes an animated gif of the event cleaning's time walk on one event: the 8 PE seed stage, then the
// 4 PE extension stage. Each step of the walk is one pixel (the parent) searching its neighbours for
// signal inside a window of time bins around its own time bin. The gif sweeps through that window one
// time bin per frame: the whole camera shows the PE at the bin, the parent and the neighbours it is
// searching are outlined, and a result frame then shows which neighbours joined.
// The camera is drawn the way the cleaning PDF panels draw it (same orientation, MUSIC boundaries,
// horizon for sims), so pixels sit where they sit in panels 2 and 3.


IUtilities *util;
IPlotTools *plottools;
TCanvas *canvas = 0;

std::string mnt = "/cedar/";
std::string outDir = "/cedar/DataAnalysis/MakeThingsPretty/TimeWalkGiFs/Gifs/";
std::string CalibrationFactorDir = "/cedar/DataAnalysis/flasher_calibration/Output/";
std::string neighborDir = "/cedar/DataAnalysis/event_cleaning/TimeWalkCleaning/neighbors/";

int MaxNofChannels = 256;

// "data" or "sim", decided per file by whether it has a Sim tree. The copied cleaning rules below
// branch on it the same way EventCleaning does.
std::string whatData = "data";


// ---- Copied from event_cleaning/TimeWalkCleaning/EventCleaning.h: keep in step with it ----
// The gif reruns the cleaning's walk rather than reading it back from the cleaning's output, so every
// rule the walk depends on is repeated here with the same name. If one of these changes there, change
// it here too, or the gif shows a different walk from the one the cleaning did.
// RISE_FRACTION=<percent> and TIME_WINDOW=<n> in the environment override them the same way.
int TimeWalkWindowHalfWidth = 2;               // n: a neighbour is searched in parent time bin +/- n
const int SeedThresholdPE      = 8;            // stage 1 threshold, PE (panel 2)
const int ExtensionThresholdPE = 4;            // stage 2 threshold, PE (panel 3)
const bool ExtensionStageIsRecursive = true;   // extension pixels search their own neighbours as well
const bool TimeWalkUsesDiagonals     = true;   // all 8 neighbours (neighbor.csv), not the 4 broadside ones
const int MaxTimeDrift = 4;                    // windows never reach further than this from time bin zero; -1 = no cap
int RiseFractionPercent = 50;                  // a pixel's time bin is where its rising edge reaches this % of its max
const int MaxRiseLookBack = 20;                // furthest the rise search walks back from the peak
enum SaturationPeakChoice { kSatFirst = 0, kSatAvg = 1, kSatLast = 2 };
const int SaturationPeakMode = kSatFirst;      // which railed sample counts as a saturated pixel's peak
const int TriggeredChannelAmpCutOff = 200;     // ADC; events with a dimmer triggered pixel never reach the walk
const int TimeBinAll = 239;                    // data: expected pulse bin
const int SimCareTraceStart = 156;             // sims: CARE sample 0 sits at this bin of the 512-sample trace
const int SimPulseBin = SimCareTraceStart + 28;// sims: expected pulse bin (184)
const int PeakSearchBinsBefore = 9;            // peak search window [expected - before, expected + after)
const int PeakSearchBinsAfter  = 11;
const int SimPedestalStart   = 157;            // sims: pedestal = mean of samples 157-166
const int SimPedestalSamples = 10;
const int SimSaturationPlateauSamples = 2;     // sims: a railed pixel sits this many samples on its minimum
const int SimSaturationMinDepth       = 2000;  // ... at least this many ADC below its pedestal
// ---- end of the copied settings ----


// ---- Gif settings ----
// Frame delays in centiseconds (GIF delays are whole cs).
const int IntroFrameDelay  = 250;  // the first frame (time bin zero) and the stage title frames
const int SweepFrameDelay  = 35;   // one time bin of a parent's search window
const int ResultFrameDelay = 120;  // which neighbours joined, after each parent's sweep
const int FinalFrameDelay  = 400;  // the finished panel-2 and panel-3 pixels at the snapshot bin

// true: one frame per time bin of every search window, then the result frame.
// false: only the result frame of each step, for a shorter gif of a big event.
const bool ShowEveryWindowBin = true;

// Outline colours. MUSIC boundaries are red, so no outline is red.
const Color_t SeedPixelColor      = kWhite;    // in the cluster from the 8 PE stage (and the triggered pixel)
const Color_t ExtensionPixelColor = kMagenta;  // in the cluster from the 4 PE stage
const Color_t ParentPixelColor    = kOrange+1; // the pixel whose neighbours are being searched
const Color_t SearchedPixelColor  = kCyan;     // a neighbour being searched
const Color_t PassingPixelColor   = kGreen;    // a neighbour at or above the threshold in the bin shown
const Color_t RejectedPixelColor  = kGray+1;   // a neighbour that never reached the threshold in the window

// Prints the walk's decisions in the same format as the cleaning's "event <N>" printout, so the two
// can be compared line by line.
const bool PrintWalkDiagnostics = true;


// One row of the gif list: Output,File,Event
struct GifRequest {
    std::string output;    // gif name prefix, the event number is added to it
    std::string filePath;  // full path to a merged data file or a sim _tdem.root file
    Long64_t testEntry;    // entry in the Test tree (the E# of the cleaning PDF for data)
};

// Everything the walk needs from one event, measured the way the cleaning measures it.
struct EventSignals {
    std::vector<std::vector<Int_t>> traces;
    std::vector<float> pedestals;
    std::vector<float> calibrationFactors;  // ADC -> PE, per pixel
    std::vector<int> peakTimeBin;           // own peak bin, first railed bin if saturated
    std::vector<int> riseTimeBin;           // the time bin the walk uses
    int triggeredPixelID = -1;
    int triggeredMUSICID = -1;
    int snapshotBin = -1;                   // triggered pixel's peak bin, where the cleaning's images are read

    // One pixel's calibrated signal at one sample, in PE
    float PEAtBin(int pixelID, int bin) const;
    // Every pixel's PE at one time bin
    std::vector<float> CameraAtBin(int bin) const;
};

// Search window of a step: [first, last]; first > last means nothing was left to search
struct BinWindow { int first; int last; };

// One neighbour checked in one step of the walk
struct NeighborCheck {
    int pixelID;
    int bestBin;   // brightest bin inside the window, -1 when the window is empty
    float bestPE;
    bool joins;
};

// One step of the walk: a parent pixel searching all its neighbours that are not in the cluster yet
enum WalkStage { kSeedStage = 1, kExtensionStage = 2 };
struct WalkStep {
    WalkStage stage;
    int parentPixelID;
    BinWindow window;
    std::vector<NeighborCheck> neighbors;
    std::vector<int> clusterStageBefore;  // per pixel: 0 not in the cluster, else the WalkStage it joined in
};

// The recorded walk: every step in the order the cleaning takes them
struct WalkRecord {
    std::vector<WalkStep> steps;
    std::vector<int> clusterStage;        // final membership, per pixel (0 / kSeedStage / kExtensionStage)
};


// Event loading, with the cleaning's pedestal, calibration, peak, saturation and rise rules
std::vector<std::vector<std::string>> splitCSVLines(const std::vector<std::string>& lines);
int extractYYYYMMDD(const std::string& text);
std::vector<Long64_t> MatchTestToSimEntries(TTree* testTree, TTree* simTree, SEvent* simEvent);
int ExpectedPulseBin();
float PedestalForPixel(const std::vector<Int_t>& trace, float pulsePedestal);
float AmplitudeAtBin(const std::vector<Int_t>& trace, float pedestal, int bin);
std::vector<float> GetRailLevels(const std::vector<std::vector<Int_t>>& traces, const std::vector<float>& pedestals,
                                 int windowStart, int windowEnd);
std::vector<IUtilities::SatSpan> FindSaturationSpans(const std::vector<std::vector<Int_t>>& traces,
                                                     const std::vector<float>& railLevels, int windowStart, int windowEnd);
void ApplySaturationPeakMode(std::vector<int>& peakTimeBin, const std::vector<IUtilities::SatSpan>& spans);
std::vector<int> FindRiseTimeBins(const std::vector<std::vector<Int_t>>& traces, const std::vector<float>& pedestals,
                                  const std::vector<float>& peakAmplitudes, const std::vector<int>& peakTimeBin);
std::vector<float> FactorsFromOldFlasherHistogram(TH2F* flasher, const std::vector<double>& gain);
std::vector<float> GetCalibrationFactors(int night, const std::vector<double>& gain);
int FindTriggeredMUSIC(const std::vector<float>& calibratedAmplitudes, int reportedMUSICID);
EventSignals MeasureEvent(IEvent* ev, int night, const std::vector<double>& gain);
void WarnIfCleaningStopsBeforeWalk(const EventSignals& signals);

// The walk itself, recorded step by step
const std::vector<std::vector<int>>& NeighborTable(const std::string& csvPath);
BinWindow TimeWalkSearchWindow(int parentTimeBin, int timeZeroBin, int traceLength);
std::vector<int> GrowClusterInTime(std::vector<int>& clusterStage, std::queue<int> pixelsToSearchFrom, WalkStage stage,
                                   bool keepSearching, const EventSignals& signals, WalkRecord& record);
WalkRecord RunTimeWalk(const EventSignals& signals);

// Drawing
double TiltFromPath(const std::string& path);
void DrawHorizonOverlay(TH2* hist, double tiltDeg, Color_t color);
void MakeTimeWalkGif(const EventSignals& signals, const WalkRecord& record, const std::string& gifPath,
                     const std::string& eventLabelText, double horizonTiltDeg);
