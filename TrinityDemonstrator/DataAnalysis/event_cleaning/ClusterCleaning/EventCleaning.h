//R__LOAD_LIBRARY(libExACT.so)
//#include "libExACT.so"
#include "PCA.h"
#include "PlotHelp.h"
#include <set>
#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <typeinfo>
#include "TLine.h"
#include <fstream>
#include <vector>
#include <sstream>
#include <TPaletteAxis.h>
#include <algorithm>
#include <TEllipse.h>
#include <IUtilities.h>
#include <TObjString.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <sstream>
#include <iomanip>
#include <string>
#include <iostream>
#include <dirent.h>
#include <TArrow.h>
#include <sys/types.h>
#include <vector>
#include <TLatex.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <Event.h>
#include <BiFocal.h>
#include <Pulse.h>
#include <TFile.h>
#include <TBox.h>
#include <TMarker.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <Getline.h>
#include <TRandom.h>
#include <TTimer.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <TPaveStats.h>
#include <TList.h>
#include <TFile.h>
#include <TText.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TPrincipal.h>
#include "Math/Integrator.h"
#include <TPaveText.h>
#include <random>
#include <utility>
#include <CEvent.h>
#include <IEvent.h>
#include <IPlotTools.h>
#include <ISiPM.h>
#include <execution>
#include "EventInfo.h"
#include <SEvent.h>
#include <TExec.h>
#include <TColor.h>


// Define variables (global)
std::vector<double> MaxAmplitudePixel;
std::vector<int> MaxPixelID;
std::vector<int> MaxPeakTimePixel;
std::vector<int> MaxMUSICID;
std::vector<double> MaxAmplitudePixelTimeBin;
std::vector<int> MaxPixelIDTimeBin;
std::vector<double> AvgAmplitudePerEvent;
std::vector<double> CameraAmplitudeAtTimeBin;

IUtilities *util;
IPlotTools *plottools;
CEvent *cev;
// Saturated pixel IDs as measured on the raw trace, before cleaning has run. Held here
// so CompletePanel4 can intersect them with the surviving pixels once those are known.
std::vector<Int_t> SaturatedPixelIDsRaw;
PlotHelp *plothelp;
EventInfo *eventInfo=0;

TTree *tree = 0;
TTree *treeHLED = 0;
IEvent *ev;
IEvent *evHLED;
ISiPM *sipmInfo;


TH2F* hcam_panel1=0;
TH2F* hcam_panel2=0;
TH2F* hcam_panel3=0;
TH2F* hcam_panel4=0;
TH2F* COGgraph=0;
TH2F* COGgraphweighted=0;
TH2F* pixeldist=0;

TFile *fO;
TFile *file;

TCanvas *c_cleaned =0;
TCanvas *c_PCA=0;
int MaxNofChannels = 256;
std::string OutputFileRoot = "";
std::string OutputFilePDF = "";
std::string OutputFilePDFOpen = "";
std::string OutputFilePDFClose = "";
std::string OutputFileEventCleaningDataRoot = "";
std::string whatData = "Data";

std::string mnt = "/cedar/";
std::string dataDir = "/cedar/DataAnalysis/MergedData/Output/";
std::string outDir = "/cedar/DataAnalysis/event_cleaning/Output/";


std::string neighborDir = "/cedar/DataAnalysis/event_cleaning/ClusterCleaning/neighbors/";
std::string CalibrationFactorDir = "/cedar/DataAnalysis/flasher_calibration/Output/";

// FIX THESE 2025-09-23
std::string muonDir = "/cedar/DataAnalysis/DataTxtFiles/Muon/";
std::string bkgDir = "/cedar/DataAnalysis/DataTxtFiles/BackgroundSamples/RandomSampling/";
// Simulation root. The old productions lived flat under DataAnalysis/SimulationAnalysis/data/;
// the sofia/taudavid productions live under simdata/ with one folder per thrown run.
// std::string simDir = "/cedar/DataAnalysis/SimulationAnalysis/data/";
std::string simDir = "/cedar/simdata/sofiadata/";

std::string outDirBkg = "/cedar/DataAnalysis/event_cleaning/OutputBkg/";
std::string outDirMuon = "/cedar/DataAnalysis/event_cleaning/OutputSim/";
std::string outDirSim = "/cedar/DataAnalysis/event_cleaning/OutputSim/";


const int TriggeredChannelAmpCutOff = 200; // (200 ADC/8 PE) Cut off for the triggered music channel
// Sample the demonstrator pulse is expected near in data. It centres the peak search window (see
// ExpectedPulseBin). It no longer sets where amplitudes are read: that is the triggered pixel's own
// peak bin (see TakeAmplitudeSnapshot).
const int TimeBinAll = 239;
const int CorePixelAmpCutOff = 200; // (200ish ADC/8PE) // Cut off for the pixels surrounding the triggered music channel
const int SaturatedPixelCutoff = 256; // 0 removed saturated pixels max channels from being cleaned. 1 allows them to be cleaned (changed from 0 on 2/26/2025)
int FlasherEventsCutOff = 800; // (350 ish ADC/15 PE ) Average amplitude across the camera ~1200 is  Flasher event
int PixelSurviveCutOff = 3; //How many pixels need to survive cleaning to plot
const int WLRatioCutOff = 100; // length width How elliptical you require the events to be after the Principle Compomnent Analysis (helps removed crosstalk events.)
const int rmTopRow = 1; // removed the top row of the camera (sky) 0 = remove, 1 = dont remove
int sizeCutoff = 0; // in PE so 100, or 200 pe removed the top row of the camera (sky) 0 = remove, 1 = dont remove
int NumberOfCoresCutoff = 2; // 1, 2,3 cores removed the top row of the camera (sky) 0 = remove, 1 = dont remove
float CrosstalkAllowanceCutoff = 0.3; // typically 0.3 removed the top row of the camera (sky) 0 = remove, 1 = dont remove
float CorePixelCutOff = 3;

// ---- Multi-ns shower timing ----
// Showers cross the camera over several ns, so pixels far from the triggered pixel peak in later
// FADC bins than the triggered pixel does. These settings control how the cleaning follows that.
const int TimeStepAllowance = 1; // +/- time bins allowed between neighbouring pixels in the panel-3 time walk
const int MaxTimeDrift      = 4; // cap on |pixel time bin - triggered pixel time bin| in the time walk; -1 disables the cap

// A pixel's time bin (time walk and timing map) is its rise bin: the first sample of the rising edge
// that reaches RiseFractionPercent of that pixel's max amplitude. The rising edge comes before the
// pulse saturates and before crosstalk builds up, so it can be trusted where the peak cannot.
// 100 uses the peak bin itself (the old behaviour). Amplitudes are still read at the triggered pixel's peak bin.
// 50 is the default; RISE_FRACTION=<percent> in the environment overrides it for one run, e.g.
//     RISE_FRACTION=90 ./run_script.sh <date> <file>
// The value goes into the output names (_rf_<percent>), so runs at different fractions sit side by side.
int RiseFractionPercent = 50;
const int MaxRiseLookBack     = 20; // furthest the rise search walks back from a pixel's peak, in samples

// A railed (saturated) pulse has a flat bottom, so its peak bin is ambiguous. Pick which railed
// sample counts as the peak: the first, the middle of the span, or the last.
// First is used: it is the last sample that still follows the signal; after it the ADC sits on the
// rail and no longer tracks the pulse.
enum SaturationPeakChoice { kSatFirst = 0, kSatAvg = 1, kSatLast = 2 };
const int SaturationPeakMode = kSatFirst;

// ---- Per-event diagnostics ----
// OFF by default, so a production run over a whole night prints only per-file lines and the final
// totals. Turned on from the command line, which run_script.sh forwards:
//     ./run_script.sh <date> <file> verbose         every event
//     ./run_script.sh <date> <file> event 18        only event 18 (implies verbose)
// Nothing about the cleaning itself changes; only how much is printed.

// The per-event timing block: snapshot bin, saturation spans, rise times, every time-walk decision.
bool PrintTimingDiagnostics = false;

// A labelled verdict for every cut an event passes through: what was measured, what the cut asks
// for, and which flag removed the event. Answers "why was this event thrown away?" from stdout
// alone, without reopening the output files.
bool PrintCutDiagnostics = false;

// Restrict the per-event printout to a single event, so one event can be read without wading
// through the whole run. -1 prints every event.
// N is the event number the plot titles show (E#): the sim event number for sims, the entry
// number in the Test tree for data.
int FocusEventNumber = -1;

// Event number of the event being processed, so the print helpers below need not all be handed it.
// Set once per event by PrintEventBanner.
int DiagnosedEventNumber = -1;

// True when the event being processed is the one the printout is restricted to.
bool DiagnoseThisEvent() { return FocusEventNumber < 0 || DiagnosedEventNumber == FocusEventNumber; }

// True when this event's per-event lines should be printed at all. Every diagnostic print in this
// file and the ones handed to CEvent go through it, so one run stays quiet and the next is verbose.
bool VerboseThisEvent() { return PrintCutDiagnostics && DiagnoseThisEvent(); }

// Reads the optional trailing command-line arguments ("verbose", "event <N>") and sets the
// switches above. Returns nothing; unrecognised arguments are ignored.
void ParseDiagnosticArguments(int argc, char** argv);

// Sims are cleaned exactly like data; only where the pulse sits in the trace differs. The 200-sample
// CARE trace is centred in the 512-sample window at offset 156 (PlotCamera::create_root_file_care)
// and the CARE pulse is at kCareTimeBin = 28 within it (PlotCamera.h). Measured on sim20260827:
// bright-pixel peaks have mode 183, mean 184, 90% within 179-190.
const int SimCareTraceStart = 156;              // CARE sample 0 sits at this bin of the 512-sample trace
const int SimPulseBin = SimCareTraceStart + 28; // = 184

// Peak search window around the expected pulse: [expected - before, expected + after).
// For data that is [230,250); for sims [175,195).
const int PeakSearchBinsBefore = 9;
const int PeakSearchBinsAfter  = 11;

// ---- Sim pedestal ----
// Pulse averages samples 0-99 of the 512-sample trace whatever window the peak is searched in
// (Pulse::CalcPedestal, called as CalcPedestal(100) by every constructor used here). For the
// demonstrator those samples are real sky. For sims they are not sampled at all:
// PlotCamera::create_root_file_care fills the 512-sample trace with carePedestal+3400 and drops the
// 200-sample CARE trace in at offset 156, so samples 0-99 are constant filler carrying no noise.
// Sims therefore take their pedestal from inside the CARE trace instead, ahead of the pulse at ~180.
// Averaged here rather than in Pulse so ExACT keeps one pedestal definition.
// Kept short because the pulse starts earlier than it looks: on sim20260827 event 152 sample 176 is
// already 100 ADC into the rising edge, so a 20-sample window pulled the pedestal ~9 ADC low.
const int SimPedestalStart   = 157; // first sample averaged; the CARE trace starts at sample 156
const int SimPedestalSamples = 10;  // so the window is samples 157-166

// ---- Sim saturation ----
// The demonstrator rails at the bottom of the ADC range, so one fixed level (IUtilities::kSatADCFloor)
// says when any pixel has stopped following the signal. A CARE pixel caps at its own value instead:
// create_root_file_care flips and shifts each trace by that pixel's own CARE pedestal, so the cap
// lands somewhere different for every pixel (sim20260827 event 152: -15 ADC for pixel 0, 35 for
// pixel 1, 37 for pixel 5, and the fixed floor catches only pixel 0). What capped pixels have in
// common is the shape: consecutive samples sitting on exactly the same value at the bottom of the
// pulse, where no more signal is recoverable. That plateau is what marks a sim pixel railed.
// A plateau on its own is not enough: an empty pixel's baseline noise repeats its minimum value by
// chance every few samples, which on sim20260827 event 152 marked 47 of the 256 pixels railed and
// handed eight zero-amplitude pixels to the core search. A real cap is also deep, so the plateau
// must sit far below the pedestal as well. Capped pixels there are ~3740 ADC deep, the accidental
// ones 20-40, so anything between the two works and the exact value does not matter.
const int SimSaturationPlateauSamples = 2;    // consecutive samples on the minimum
const int SimSaturationMinDepth       = 2000; // ADC the plateau must sit below the pedestal

// Sample the pulse is expected near for the data type being cleaned.
int ExpectedPulseBin() { return (whatData == "sim") ? SimPulseBin : TimeBinAll; }


void ParseDiagnosticArguments(int argc, char** argv)
{
    for (int arg = 4; arg < argc; arg++) {
        std::string option = argv[arg];
        if (option == "verbose" || option == "-v") {
            PrintCutDiagnostics   = true;
            PrintTimingDiagnostics = true;
            std::cout << "Verbose diagnostics on for every event" << std::endl;
        } else if (option == "event" && arg + 1 < argc) {
            FocusEventNumber      = std::atoi(argv[arg + 1]);
            PrintCutDiagnostics   = true;
            PrintTimingDiagnostics = true;
            arg++;
            std::cout << "Verbose diagnostics on for event " << FocusEventNumber << " only" << std::endl;
        }
    }
}


// Sets RiseFractionPercent from the RISE_FRACTION environment variable when it is set, and says which
// value this run uses. Returns false for a value that is not a whole number from 1 to 100, so a typo
// stops the run instead of quietly running at the default.
bool ReadRiseFractionFromEnvironment()
{
    const char* riseFractionEnv = getenv("RISE_FRACTION");
    if (riseFractionEnv == nullptr || *riseFractionEnv == '\0') {
        cout << "Rise fraction: " << RiseFractionPercent << "% (default; set RISE_FRACTION to change it)" << endl;
        return true;
    }
    char* end = nullptr;
    long requested = strtol(riseFractionEnv, &end, 10);
    if (*end != '\0' || requested < 1 || requested > 100) {
        cout << "RISE_FRACTION=" << riseFractionEnv << " is not a whole percent from 1 to 100, stopping" << endl;
        return false;
    }
    RiseFractionPercent = (int)requested;
    cout << "Rise fraction: " << RiseFractionPercent << "% (from RISE_FRACTION)" << endl;
    return true;
}


// Collect the CARE files of a sofia-style simulation production.
// Those productions give every thrown run its own folder, so the files sit at
//   <date>/<run>/CARE/Tilt_<angle>/<shower>_tdem.root
// e.g. 20260827_taudavid/20260827_11415796_daro/CARE/Tilt_90.000000/nutausim_..._tdem.root
// Walk one level of run folders, then one level of Tilt folders, and keep only the
// _tdem.root files (the plain .root and .vbf siblings are the pre-Trinity CARE output).
// The names come back relative to <date>/ so the caller can still build FolderPath + name
// and so the filename stored on each event says which run it came from.
std::vector<std::string> GetSimTdemFiles(const std::string& dateFolder)
{
    std::vector<std::string> simFiles;
    const std::string tdemSuffix = "_tdem.root";

    DIR* runFolders = opendir(dateFolder.c_str());
    if (runFolders == nullptr) {
        std::cout << "Could not open simulation folder: " << dateFolder << std::endl;
        return simFiles;
    }

    while (struct dirent* run = readdir(runFolders)) {
        std::string runName = run->d_name;
        if (runName == "." || runName == "..") continue;

        std::string careFolder = dateFolder + runName + "/CARE/";
        DIR* tiltFolders = opendir(careFolder.c_str());
        if (tiltFolders == nullptr) continue; // not a run folder, or the run never reached CARE

        while (struct dirent* tilt = readdir(tiltFolders)) {
            std::string tiltName = tilt->d_name;
            if (tiltName.find("Tilt") != 0) continue; // skips ".", "..", care.log

            std::string tiltFolder = careFolder + tiltName + "/";
            DIR* careFiles = opendir(tiltFolder.c_str());
            if (careFiles == nullptr) continue;

            while (struct dirent* careFile = readdir(careFiles)) {
                std::string fileName = careFile->d_name;
                if (fileName.size() <= tdemSuffix.size()) continue;
                if (fileName.compare(fileName.size() - tdemSuffix.size(),
                                     tdemSuffix.size(), tdemSuffix) != 0) continue;
                simFiles.push_back(runName + "/CARE/" + tiltName + "/" + fileName);
            }
            closedir(careFiles);
        }
        closedir(tiltFolders);
    }
    closedir(runFolders);

    // readdir hands back whatever order the filesystem likes; sort so runs are processed
    // in run-number order and reruns of the same date are reproducible.
    std::sort(simFiles.begin(), simFiles.end());
    return simFiles;
}


void SetBranches(IEvent *evD);
void SetBranchesHLED(IEvent *evD);
void LoadEvents(std::string filename, std::string treeString);
void LoadEventsHLED(string NameofFile, std::string treeString);
void saveEventInfo(EventInfo* evI, TTree* treeSims);

void LoadDataPCA(PCA& pca, TH2F* hist, int totalAmp);
std::vector<double> CreateWLRatio(PCA& pca, TVectorD& eigenVals, TMatrixD& eigenVecs);
std::vector<int> CheckCrossPoints(TArrow* arrow, TH2F* hist, EventInfo* eventInfo);
void DrawSaturationMarkers(const std::vector<Int_t>& pixelIDs);
void CompletePanel4(PCA& pca, TH2F* hcam_panel4, CEvent* cev, std::vector<double> EllipicRatio, TVectorD& eigenVals, TMatrixD& eigenVecs, EventInfo* eventInfo);
bool projectOnMajorAxis(double xcog, double ycog, double anglerad, std::vector<int> sur_pix,
                       std::vector<float> amps, std::vector<double>& d, std::vector<double>& q);
double getM3Long(double xcog, double ycog, double anglerad, std::vector<int> sur_pix, std::vector<float> amps);
double getM3LongPow3(double xcog, double ycog, double anglerad, std::vector<int> sur_pix, std::vector<float> amps);
std::vector<double> generateRandomNumbers();
void CreateFileName(std::string filename, bool bkg);

// ---- Multi-ns shower timing helpers ----
// Moves each saturated pixel's peak bin onto the railed sample chosen by SaturationPeakMode.
void ApplySaturationPeakMode(std::vector<int>& peakTimeBin, const std::vector<IUtilities::SatSpan>& spans);
// Each pixel's rise bin: first sample of the rising edge at or above RiseFractionPercent of its max amplitude.
std::vector<int> FindRiseTimeBins(const std::vector<std::vector<Int_t>>& traces, const std::vector<float>& pedestals,
                                  const std::vector<float>& peakAmplitudes, const std::vector<int>& peakTimeBin);
// Prints the triggered pixel's rise-time numbers and the rise/peak bins of pixels that pass the panel-2 cut.
void PrintRiseTimes(const std::vector<float>& pedestals, const std::vector<float>& peakAmplitudes,
                    const std::vector<int>& peakTimeBin, const std::vector<int>& riseTimeBin,
                    const std::vector<float>& calibratedAmplitudes, int triggeredPixelID);
// Applies the flasher and absolute gain calibration to raw amplitudes and stores the result in cev.
void StoreCalibratedAmplitudes(const std::vector<float>& rawAmplitudes);
// Per-pixel calibration factor (flasher x absolute gain), so amplitudes can be calibrated locally.
std::vector<float> GetCalibrationFactors();
// Multiplies each pixel's amplitude by its calibration factor.
std::vector<float> ApplyCalibration(const std::vector<float>& amplitudes, const std::vector<float>& factors);
// Picks the triggered MUSIC: the reported one, or its mate if the mate holds the brighter pixel.
int FindTriggeredMUSIC(const std::vector<float>& calibratedAmplitudes, int reportedMUSICID);
// Reads every pixel's amplitude at one common bin (one snapshot in time for the whole camera).
std::vector<float> TakeAmplitudeSnapshot(const std::vector<std::vector<Int_t>>& traces,
                                         const std::vector<float>& pedestals, int snapshotBin);
// The pedestal this pixel is measured against: Pulse's own for data, the CARE window for sims.
float PedestalForPixel(const std::vector<Int_t>& trace, float pulsePedestal);
// One pixel's amplitude at one sample, measured against the pedestal actually in use.
float AmplitudeAtBin(const std::vector<Int_t>& trace, float pedestal, int bin);
// The level at or below which each pixel counts as railed: the fixed ADC floor for data, each sim
// pixel's own plateau value.
std::vector<float> GetRailLevels(const std::vector<std::vector<Int_t>>& traces, const std::vector<float>& pedestals,
                                 int windowStart, int windowEnd);
// Where each pixel sat on its own rail level inside the window.
std::vector<IUtilities::SatSpan> FindSaturationSpans(const std::vector<std::vector<Int_t>>& traces,
                                                     const std::vector<float>& railLevels,
                                                     int windowStart, int windowEnd);
// Prints the per-event timing block: triggered pixel, snapshot bin, saturation spans.
void PrintTimingSummary(const std::vector<std::vector<Int_t>>& traces, const std::vector<float>& calibrationFactors,
                        const std::vector<float>& pedestals,
                        const std::vector<int>& peakTimeBin, const std::vector<IUtilities::SatSpan>& saturationSpans,
                        int reportedMUSICID, int triggeredMUSICID, int triggeredPixelID, int snapshotBin, bool usesSnapshot);
// Which pixels count as railed and at what level, so a sim cap that the fixed ADC floor misses is visible.
void PrintRailLevels(const std::vector<float>& railLevels, const std::vector<IUtilities::SatSpan>& saturationSpans);
// ---- Per-event cut diagnostics ----
// Opens an event's printout block and sets DiagnosedEventNumber, so every line below can be tied
// to one event.
void PrintEventBanner(int entryNumber, int eventNumber, int totalEvents);
// One cut's verdict: what was measured, what the cut requires, and the flag if it removed the event.
void PrintCutVerdict(const std::string& cutName, bool removed, const std::string& measured,
                     const std::string& requirement, int eventFlag);
// Closing line of an event's block, saying which cut ended it or that it survived everything.
void PrintEventOutcome(const std::string& outcome, int eventFlag);
// The brightest pixels after calibration. For sims this is where an empty or mispositioned CARE
// pulse shows up first: the whole camera reads ~0 PE.
void PrintBrightestPixels(const std::vector<float>& calibratedAmplitudes, int howMany);
// How much signal the raw traces hold, so an all-zero CARE file is told apart from a real dim event.
void PrintTraceSummary(const std::vector<std::vector<Int_t>>& traces, int triggeredPixelID, int snapshotBin);
// Sims only: the truth record behind this event (run, trigger flag, energy, geometry).
void PrintSimTruth(SEvent* simTruth, Long64_t simEntry, int entryNumber);
// Where each pixel's pedestal came from, and how flat the samples it was averaged over are.
void PrintPedestalSummary(const std::vector<float>& pedestals, const std::vector<std::vector<Int_t>>& traces,
                          int triggeredPixelID);

// Half-width in bins of the colour scale both timing maps share: MaxTimeDrift, widened to reach every panel-2 pixel.
int TimingMapRange(CEvent* cev, int triggeredTimeBin);
// Builds the 16x16 map of (pixel time bin - triggered pixel time bin) for the given pixels; panelLabel heads the title.
TH2F* FillTimingHeatMap(CEvent* cev, int triggeredTimeBin, const std::string& filenameTitle,
                        const std::vector<int>& pixelsToShow, const std::string& panelLabel, int colourRange);
// Draws the panel-2 and panel-3 timing maps side by side on c_time with their own red-black-green palette, so
// they open in colour from the ROOT file. The triggered pixel is covered by a grey box labelled "T".
void DrawTimingHeatMaps(TH2F* hcam_time_panel2, TH2F* hcam_time_panel3, int triggeredPixelID);

void CreateFileName(std::string filename, std::string dataType) {
    if (dataType == "bkg") {
        outDir = outDirBkg;
    }
    if (dataType == "muon") {
        outDir = outDirMuon;
    }
    if (dataType == "sim") {
        outDir = outDirSim;
    }

    // Every cut goes into the file name so outputs made with different settings never overwrite
    // each other and can be sorted by setting.
    // _TC_ = TriggeredChannelAmpCutOff   _TB_ = ExpectedPulseBin()   _NP_ = CorePixelAmpCutOff
    // _s_  = SaturatedPixelCutoff        _FA_ = FlasherEventsCutOff  _mp_ = PixelSurviveCutOff
    // _er_ = WLRatioCutOff               _tr_ = rmTopRow
    // _ts_ = TimeStepAllowance           _md_ = MaxTimeDrift         _sp_ = SaturationPeakMode
    // _rf_ = RiseFractionPercent
    std::string cutTag = Form("_TC_%i_TB_%i_NP_%i_s_%i_FA_%i_mp_%i_er_%i_tr_%i_ts_%i_md_%i_sp_%i_rf_%i",
        TriggeredChannelAmpCutOff,
        ExpectedPulseBin(),
        CorePixelAmpCutOff,
        SaturatedPixelCutoff,
        FlasherEventsCutOff,
        PixelSurviveCutOff,
        WLRatioCutOff,
        rmTopRow,
        TimeStepAllowance,
        MaxTimeDrift,
        SaturationPeakMode,
        RiseFractionPercent
    );
    std::string displayFileBase = Form("%sEventCleanedCluster%s%s", outDir.c_str(), filename.c_str(), cutTag.c_str());

    OutputFileRoot = displayFileBase + ".root";
    OutputFileEventCleaningDataRoot = Form("%sDataFiles/Data_EventCleaning%s%s.root", outDir.c_str(), filename.c_str(), cutTag.c_str());

    int check = mkdir(Form("%sDataFiles",outDir.c_str()),0777);

    // check if directory is created or not
    if (!check)
        printf("Directory created\n");
    else {
        printf("Unable to create directory\n");
    }

    // "[" opens a multi-page PDF and "]" closes it
    OutputFilePDF      = displayFileBase + ".pdf";
    OutputFilePDFOpen  = displayFileBase + ".pdf[";
    OutputFilePDFClose = displayFileBase + ".pdf]";
}


void LoadEvents(string NameofFile, std::string treeString)
{
    NameofFile = NameofFile;
    cout << "Loading file: " << NameofFile << endl;
    TFile *fO = new TFile(NameofFile.c_str(), "READ");
    tree = (TTree*)fO->Get(treeString.c_str());
}

void LoadEventsHLED(string NameofFile, std::string treeString)
{
    NameofFile = NameofFile;
    //cout << "Loading file for HLED: " << NameofFile << endl;
    TFile *fO = new TFile(NameofFile.c_str(), "READ");
    treeHLED = (TTree*)fO->Get(treeString.c_str());
}


// A railed pulse sits on the ADC floor for several samples, and Pulse reports the earliest of them
// as the peak. Replace that with the sample chosen by SaturationPeakMode. Unsaturated pixels keep
// the peak bin Pulse found.
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


// Walks back from each pixel's peak along the rising edge while the samples stay at or above
// RiseFractionPercent of the pixel's max amplitude, and returns the earliest such sample. Pulses are
// negative-going, so the amplitude of a sample is pedestal - sample (raw ADC, same units as the max).
// Pixels with no pulse (max 0) keep their peak bin; the walk never goes back more than MaxRiseLookBack.
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


void PrintRiseTimes(const std::vector<float>& pedestals, const std::vector<float>& peakAmplitudes,
                    const std::vector<int>& peakTimeBin, const std::vector<int>& riseTimeBin,
                    const std::vector<float>& calibratedAmplitudes, int triggeredPixelID)
{
    const float triggeredThreshold = peakAmplitudes[triggeredPixelID] * RiseFractionPercent / 100.0;
    cout << "Rise time (" << RiseFractionPercent << "% of max): triggered pixel " << triggeredPixelID
         << " pedestal " << pedestals[triggeredPixelID] << " ADC, max " << peakAmplitudes[triggeredPixelID]
         << " ADC, threshold " << triggeredThreshold << " ADC -> rise bin " << riseTimeBin[triggeredPixelID]
         << " (peak bin " << peakTimeBin[triggeredPixelID] << ")" << endl;

    // same cutoff SetPanel2 applies (it takes an int, so the PE value is truncated there too)
    const int panel2CutoffPE = (int)(CorePixelAmpCutOff / util->GetADCtoPEratio());
    cout << "Rise bins of pixels passing the panel-2 cut (ID: rise/peak):";
    int printed = 0;
    for (std::size_t pixelID = 0; pixelID < riseTimeBin.size(); pixelID++) {
        if (calibratedAmplitudes[pixelID] < panel2CutoffPE) continue;
        cout << (printed == 0 ? " " : ", ") << pixelID << ": " << riseTimeBin[pixelID] << "/" << peakTimeBin[pixelID];
        if (riseTimeBin[pixelID] == peakTimeBin[pixelID] - MaxRiseLookBack) cout << " (hit look-back limit)";
        printed++;
    }
    if (printed == 0) cout << " none";
    cout << endl;

    // How much of the image is genuinely in step, on each of the two clocks. The time walk runs on
    // the rise bin, but sim pulses all have the same shape and so tend to share one peak bin: if the
    // peak line says "all" and the rise line says "1", every time step the walk had to take came
    // from the rise-bin spread, not from the shower crossing the camera.
    int aboveCut = 0, sharePeakBin = 0, shareRiseBin = 0;
    for (std::size_t pixelID = 0; pixelID < riseTimeBin.size(); pixelID++) {
        if (calibratedAmplitudes[pixelID] < panel2CutoffPE) continue;
        aboveCut++;
        if (peakTimeBin[pixelID] == peakTimeBin[triggeredPixelID]) sharePeakBin++;
        if (riseTimeBin[pixelID] == riseTimeBin[triggeredPixelID]) shareRiseBin++;
    }
    cout << "Pixels above the panel-2 cut sharing the triggered pixel's peak bin ("
         << peakTimeBin[triggeredPixelID] << "): " << sharePeakBin << " of " << aboveCut
         << "   rise bin (" << riseTimeBin[triggeredPixelID] << "): " << shareRiseBin << " of " << aboveCut << endl;
}


// The only place sims and data are calibrated differently. Sims have no flasher run and no SiPM gain
// record, so they use the "all1new" flasher file and a fixed gain.
void StoreCalibratedAmplitudes(const std::vector<float>& rawAmplitudes)
{
    if (whatData == "sim") {
        // PixelSurviveCutOff=1;
        // sizeCutoff = 0; 
        // NumberOfCoresCutoff = 0; 
        // PixelSurviveCutOff=5;
        // sizeCutoff = 200; 
        // NumberOfCoresCutoff = 3;

        // CorePixelCutOff = 100;

        // all1new holds -1 in every pixel, so the amplitudes are negated here to come out positive.
        std::vector<float> signCorrected(rawAmplitudes.size());
        for (std::size_t pixelID = 0; pixelID < rawAmplitudes.size(); pixelID++) {
            signCorrected[pixelID] = -rawAmplitudes[pixelID];
        }
        // for simulations the photon detection efficiency was applied in a prior step, so only the ADC to PE conversion is needed here
        std::vector<double> fakeGain(256, 24.1);
        cev->SetAmplitudeValuesTimeBin(signCorrected, CalibrationFactorDir, "all1new", fakeGain, 44.0);
    } else {
        cev->SetAmplitudeValuesTimeBin(rawAmplitudes, CalibrationFactorDir, std::to_string(cev->GetEventDate()), sipmInfo->GetGain(), 44.0);
    }
}


// Calibration (flasher correction and absolute gain) is one multiplicative factor per pixel, so
// calibrating a vector of ones returns exactly those factors. Having them lets the peak amplitudes
// and the snapshot amplitudes be calibrated without reopening the flasher file for each one.
std::vector<float> GetCalibrationFactors()
{
    StoreCalibratedAmplitudes(std::vector<float>(MaxNofChannels, 1.0));
    return cev->GetAmplitudeValuesTimeBin();
}


std::vector<float> ApplyCalibration(const std::vector<float>& amplitudes, const std::vector<float>& factors)
{
    std::vector<float> calibrated(amplitudes.size(), 0.0);
    for (std::size_t pixelID = 0; pixelID < amplitudes.size() && pixelID < factors.size(); pixelID++) {
        calibrated[pixelID] = amplitudes[pixelID] * factors[pixelID];
    }
    return calibrated;
}


// The trigger reports one MUSIC, but its mate can hold the brighter pixel. Use whichever of the two
// has the higher maximum amplitude.
int FindTriggeredMUSIC(const std::vector<float>& calibratedAmplitudes, int reportedMUSICID)
{
    int mateMUSICID = util->GetMUSICMate(reportedMUSICID, true);
    if (util->GetMaximum(calibratedAmplitudes, reportedMUSICID) < util->GetMaximum(calibratedAmplitudes, mateMUSICID)) {
        return mateMUSICID;
    }
    return reportedMUSICID;
}


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


// Same convention as Pulse::CalcAmplitude: the pulse goes down from the pedestal, and a sample above
// the pedestal reads as zero rather than as a negative amplitude.
float AmplitudeAtBin(const std::vector<Int_t>& trace, float pedestal, int bin)
{
    if (bin < 0 || bin >= (int)trace.size()) return 0.0;
    const float amplitude = pedestal - trace[bin];
    return (amplitude < 0.0) ? 0.0 : amplitude;
}


std::vector<float> GetRailLevels(const std::vector<std::vector<Int_t>>& traces, const std::vector<float>& pedestals,
                                 int windowStart, int windowEnd)
{
    std::vector<float> railLevels(traces.size(), (float)IUtilities::kSatADCFloor);
    if (whatData != "sim") return railLevels; // the demonstrator rails at the ADC floor

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
        // without a deep plateau the pixel falls back to the fixed floor, which it only reaches if
        // its trace really did run off the bottom of the range
        railLevels[pixelID] = (hasPlateau && isDeep) ? (float)lowest : (float)IUtilities::kSatADCFloor;
    }
    return railLevels;
}


// Same span bookkeeping as IUtilities::GetSaturationSpans, except each pixel is compared with its
// own rail level instead of one level shared by the whole camera.
std::vector<IUtilities::SatSpan> FindSaturationSpans(const std::vector<std::vector<Int_t>>& traces,
                                                     const std::vector<float>& railLevels,
                                                     int windowStart, int windowEnd)
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


// One snapshot in time for the whole camera: every pixel's amplitude is read at snapshotBin,
// against the pedestal that pixel is measured with.
std::vector<float> TakeAmplitudeSnapshot(const std::vector<std::vector<Int_t>>& traces,
                                         const std::vector<float>& pedestals, int snapshotBin)
{
    std::vector<float> snapshot;
    for (std::size_t pixelID = 0; pixelID < traces.size(); pixelID++) {
        snapshot.push_back(AmplitudeAtBin(traces[pixelID], pedestals[pixelID], snapshotBin));
    }
    return snapshot;
}


void PrintTimingSummary(const std::vector<std::vector<Int_t>>& traces, const std::vector<float>& calibrationFactors,
                        const std::vector<float>& pedestals,
                        const std::vector<int>& peakTimeBin, const std::vector<IUtilities::SatSpan>& saturationSpans,
                        int reportedMUSICID, int triggeredMUSICID, int triggeredPixelID, int snapshotBin, bool usesSnapshot)
{
    const char* satModeNames[] = {"first", "middle", "last"};

    cout << "--- Timing diagnostics, event " << cev->GetEventNumber() << " ---" << endl;

    cout << "Triggered MUSIC: reported " << reportedMUSICID << ", used " << triggeredMUSICID;
    if (triggeredMUSICID != reportedMUSICID) cout << " (swapped to mate, brighter peak amplitude)";
    cout << endl;

    if (usesSnapshot) {
        cout << "Snapshot bin (triggered pixel " << triggeredPixelID << " peak): " << snapshotBin
             << "   [expected pulse bin: " << ExpectedPulseBin() << "]" << endl;

        // Would reading every pixel at the fixed expected pulse bin have picked a different triggered pixel?
        std::vector<float> fixedBinAmplitudes = ApplyCalibration(TakeAmplitudeSnapshot(traces, pedestals, ExpectedPulseBin()), calibrationFactors);
        int fixedBinTriggeredPixelID = util->GetMaximumPixelID(fixedBinAmplitudes, FindTriggeredMUSIC(fixedBinAmplitudes, reportedMUSICID));
        if (fixedBinTriggeredPixelID != triggeredPixelID) {
            cout << "WARNING: triggered pixel changed: fixed bin " << ExpectedPulseBin() << " picked pixel " << fixedBinTriggeredPixelID
                 << " (peak bin " << peakTimeBin[fixedBinTriggeredPixelID] << "), peak amplitudes pick pixel " << triggeredPixelID << endl;
        }
    } else {
        // muon amplitudes are each pixel's own peak amplitude, so there is no common snapshot
        cout << "Triggered pixel " << triggeredPixelID << " peak bin (time walk reference): " << snapshotBin
             << "   [amplitudes are per-pixel peaks for " << whatData << "]" << endl;
    }

    const IUtilities::SatSpan& triggeredSpan = saturationSpans[triggeredPixelID];
    if (triggeredSpan.count > 0) {
        cout << "Triggered pixel saturated: yes, railed bins " << triggeredSpan.first << "-" << triggeredSpan.last
             << ", mode=" << satModeNames[SaturationPeakMode] << " -> peak bin " << peakTimeBin[triggeredPixelID];
        if (usesSnapshot) cout << " (this moves the snapshot for the whole camera)";
        cout << endl;
    } else {
        cout << "Triggered pixel saturated: no" << endl;
    }

    cout << "Saturated pixels (ID: first-last -> used):";
    int saturatedCount = 0;
    for (std::size_t pixelID = 0; pixelID < saturationSpans.size(); pixelID++) {
        const IUtilities::SatSpan& span = saturationSpans[pixelID];
        if (span.count == 0) continue;
        cout << (saturatedCount == 0 ? " " : ", ") << pixelID << ": " << span.first << "-" << span.last << " -> " << peakTimeBin[pixelID];
        saturatedCount++;
    }
    if (saturatedCount == 0) cout << " none";
    cout << endl;
}


// ---- Per-event cut diagnostics ----
// All of these print nothing unless PrintCutDiagnostics is on and the event is the one
// FocusEventNumber selects, so a focused run reads as one event's story.

void PrintEventBanner(int entryNumber, int eventNumber, int totalEvents)
{
    DiagnosedEventNumber = eventNumber; // every helper below reads this
    if (!VerboseThisEvent()) return;
    cout << "\n================ " << whatData << " event " << eventNumber
         << "  (Test entry " << entryNumber << " of " << totalEvents << ") ================" << endl;
}


void PrintCutVerdict(const std::string& cutName, bool removed, const std::string& measured,
                     const std::string& requirement, int eventFlag)
{
    if (!VerboseThisEvent()) return;
    cout << (removed ? "  CUT FAILED  " : "  cut passed  ") << std::left << std::setw(26) << cutName
         << std::right << "  measured " << std::setw(12) << measured << "   requires " << requirement;
    if (removed) cout << "   -> REMOVED, flag " << eventFlag;
    cout << endl;
}


void PrintEventOutcome(const std::string& outcome, int eventFlag)
{
    if (!VerboseThisEvent()) return;
    cout << "  EVENT " << DiagnosedEventNumber << ": " << outcome << " (flag " << eventFlag << ")" << endl;
}


void PrintBrightestPixels(const std::vector<float>& calibratedAmplitudes, int howMany)
{
    if (!VerboseThisEvent()) return;

    // sort pixel IDs by amplitude without disturbing the caller's vector
    std::vector<int> pixelOrder(calibratedAmplitudes.size());
    for (std::size_t pixelID = 0; pixelID < pixelOrder.size(); pixelID++) pixelOrder[pixelID] = (int)pixelID;
    std::sort(pixelOrder.begin(), pixelOrder.end(),
              [&calibratedAmplitudes](int a, int b) { return calibratedAmplitudes[a] > calibratedAmplitudes[b]; });

    const int panel2CutoffPE = (int)(CorePixelAmpCutOff / util->GetADCtoPEratio());
    int pixelsAboveCut = 0;
    for (float amplitude : calibratedAmplitudes) {
        if (amplitude >= panel2CutoffPE) pixelsAboveCut++;
    }

    cout << "  Brightest pixels (ID: PE at the snapshot bin):";
    for (int rank = 0; rank < howMany && rank < (int)pixelOrder.size(); rank++) {
        cout << (rank == 0 ? " " : ", ") << pixelOrder[rank] << ": " << calibratedAmplitudes[pixelOrder[rank]];
    }
    cout << endl;
    cout << "  Pixels at or above the panel-2 cut (" << panel2CutoffPE << " PE): " << pixelsAboveCut << endl;
}


void PrintTraceSummary(const std::vector<std::vector<Int_t>>& traces, int triggeredPixelID, int snapshotBin)
{
    if (!VerboseThisEvent()) return;

    // A CARE file whose FADC traces were never filled is flat everywhere; say so plainly rather
    // than letting it look like a dim shower.
    int pixelsWithSignal = 0;
    for (const std::vector<Int_t>& trace : traces) {
        if (trace.empty()) continue;
        Int_t minimum = *std::min_element(trace.begin(), trace.end());
        Int_t maximum = *std::max_element(trace.begin(), trace.end());
        if (minimum != maximum) pixelsWithSignal++;
    }
    cout << "  Traces: " << pixelsWithSignal << " of " << traces.size() << " pixels are not flat";
    if (pixelsWithSignal == 0) cout << "  <-- the FADC traces in this file are empty, nothing can trigger";
    cout << endl;

    if (triggeredPixelID >= 0 && triggeredPixelID < (int)traces.size() && !traces[triggeredPixelID].empty()) {
        const std::vector<Int_t>& trace = traces[triggeredPixelID];
        int deepestBin = (int)(std::min_element(trace.begin(), trace.end()) - trace.begin());
        cout << "  Triggered pixel " << triggeredPixelID << " trace: deepest sample " << trace[deepestBin]
             << " ADC at bin " << deepestBin << ", search window ["
             << ExpectedPulseBin() - PeakSearchBinsBefore << "," << ExpectedPulseBin() + PeakSearchBinsAfter
             << "), snapshot bin " << snapshotBin << " holds " << trace[snapshotBin] << " ADC" << endl;
        if (deepestBin < ExpectedPulseBin() - PeakSearchBinsBefore || deepestBin >= ExpectedPulseBin() + PeakSearchBinsAfter) {
            cout << "  WARNING: the deepest sample sits outside the peak search window, so the pulse this"
                    " event is cleaned on is not the real one" << endl;
        }
    }
}


void PrintSimTruth(SEvent* simTruth, Long64_t simEntry, int entryNumber)
{
    if (!VerboseThisEvent()) return;
    if (!simTruth) {
        cout << "  Sim truth: Test entry " << entryNumber << " has no matching Sim entry"
                " (event number -1, run noSimTruth)" << endl;
        return;
    }
    cout << "  Sim truth: Test entry " << entryNumber << " -> Sim entry " << simEntry
         << ", run " << simTruth->GetSumRun() << ", sim event " << simTruth->GetSimEventNumber()
         << ", triggered flag " << (simTruth->GetTriggered() ? "true" : "false") << endl;
    cout << "  Sim truth: neutrino energy " << simTruth->GetNeutrinoEnergy() << " GeV, tau energy "
         << simTruth->GetTauEnergy() << " GeV, decay length " << simTruth->GetDecayLength()
         << " m, emergence angle " << util->GetEmergenceAngle(simTruth->GetAzimuthAngle()) << " deg, "
         << "distance to emergence point "
         << util->GetEmergenceDistance(simTruth->GetTelescope_Xpos(), simTruth->GetTelescope_Ypos(),
                                       simTruth->GetTelescope_Zpos()) << " m" << endl;
}


// Samples Pulse averages for the demonstrator pedestal (Pulse::CalcPedestal, called as
// CalcPedestal(100) by every constructor EventCleaning uses). Sims override this window entirely,
// see SimPedestalStart / SimPedestalSamples.
const int kPedestalSamples = 100;

void PrintPedestalSummary(const std::vector<float>& pedestals, const std::vector<std::vector<Int_t>>& traces,
                          int triggeredPixelID)
{
    if (!VerboseThisEvent()) return;
    if (triggeredPixelID < 0 || triggeredPixelID >= (int)traces.size() || traces[triggeredPixelID].empty()) return;

    const std::vector<Int_t>& trace = traces[triggeredPixelID];
    const bool  isSim      = (whatData == "sim");
    const int   firstSample = isSim ? SimPedestalStart : 0;
    const int   lastSample  = std::min(isSim ? SimPedestalStart + SimPedestalSamples : kPedestalSamples,
                                       (int)trace.size()) - 1;
    if (firstSample > lastSample) return;

    Int_t lowest  = *std::min_element(trace.begin() + firstSample, trace.begin() + lastSample + 1);
    Int_t highest = *std::max_element(trace.begin() + firstSample, trace.begin() + lastSample + 1);

    cout << "  Pedestal: triggered pixel " << triggeredPixelID << " -> " << pedestals[triggeredPixelID]
         << " ADC, averaged over samples " << firstSample << "-" << lastSample
         << " (spread " << lowest << "-" << highest << " ADC)" << endl;
    if (isSim) {
        // Pulse's own window would be samples 0-99, which for a sim hold only the constant filler
        // create_root_file_care writes outside the CARE trace. Print both so the shift is visible.
        Pulse fillerPulse(trace, ExpectedPulseBin(), ExpectedPulseBin() + 1, 512, true);
        cout << "  Pedestal: Pulse's own samples 0-" << kPedestalSamples - 1 << " would give "
             << fillerPulse.GetPedestal() << " ADC, but they are the flat filler outside the 200-sample"
                " CARE trace (which starts at sample " << SimPulseBin - 28 << ")" << endl;
    }
}


void PrintRailLevels(const std::vector<float>& railLevels, const std::vector<IUtilities::SatSpan>& saturationSpans)
{
    if (!VerboseThisEvent()) return;

    cout << "  Railed pixels (ID: level -> railed bins):";
    int railedCount = 0;
    for (std::size_t pixelID = 0; pixelID < saturationSpans.size(); pixelID++) {
        if (saturationSpans[pixelID].count == 0) continue;
        cout << (railedCount == 0 ? " " : ", ") << pixelID << ": " << railLevels[pixelID] << " ADC -> "
             << saturationSpans[pixelID].first << "-" << saturationSpans[pixelID].last;
        railedCount++;
    }
    if (railedCount == 0) cout << " none";
    cout << endl;
    if (whatData == "sim") {
        cout << "  Railed pixels: each sim pixel caps at its own level, so a plateau of "
             << SimSaturationPlateauSamples << "+ samples on its minimum marks it railed, not the fixed "
             << IUtilities::kSatADCFloor << " ADC floor" << endl;
    }
}


// The panel-3 time walk already caps every survivor at MaxTimeDrift, but panel-2 pixels were never
// checked, and COLZ leaves anything below the minimum unpainted, so an early pixel would look absent.
// Both maps use the same range so a colour means the same offset on each.
int TimingMapRange(CEvent* cev, int triggeredTimeBin)
{
    const int baseRange = (MaxTimeDrift >= 0) ? MaxTimeDrift : 10;
    int colourRange = baseRange;
    int furthestPixelID = -1;
    std::vector<int> timeBin = cev->GetPeakTimeBin(); // rise bins, see FindRiseTimeBins
    for (int pixelID : cev->GetSurvivingPixelPanel2()) {
        int offset = std::abs(timeBin[pixelID] - triggeredTimeBin);
        if (offset > colourRange) {
            colourRange = offset;
            furthestPixelID = pixelID;
        }
    }
    if (PrintTimingDiagnostics && DiagnoseThisEvent()) {
        cout << "Timing maps: colour range +/-" << colourRange << " bins";
        if (furthestPixelID >= 0) {
            cout << " (widened from +/-" << baseRange << " for panel-2 pixel " << furthestPixelID << " at "
                 << timeBin[furthestPixelID] - triggeredTimeBin << " bins)";
        }
        cout << endl;
    }
    return colourRange;
}


TH2F* FillTimingHeatMap(CEvent* cev, int triggeredTimeBin, const std::string& filenameTitle,
                        const std::vector<int>& pixelsToShow, const std::string& panelLabel, int colourRange)
{
    TString histName = Form("hcam_time_%s", panelLabel.c_str());
    histName.ReplaceAll(" ", ""); // "Panel 2" -> hcam_time_Panel2
    TH2F* hcam_time = new TH2F(histName,
        Form("%s: %i%% rise bin - triggered rise bin --N# %i --F# %s-- E# %i ", panelLabel.c_str(), RiseFractionPercent,
             cev->GetEventDate(), filenameTitle.c_str(), cev->GetEventNumber()),
        16, -0.5, 15.5, 16, -0.5, 15.5);
    hcam_time->SetDirectory(0);
    hcam_time->SetStats(0);
    hcam_time->SetXTitle("Pixel Bin along x axis of Camera");
    hcam_time->SetYTitle("Pixel Bin along y axis of Camera");
    hcam_time->SetZTitle(Form("%i%% rise bin - triggered rise bin", RiseFractionPercent));

    // Symmetric colour range centred on the triggered pixel's bin (delta 0)
    double zMin = -colourRange - 0.5;
    double zMax =  colourRange + 0.5;
    hcam_time->SetMinimum(zMin);
    hcam_time->SetMaximum(zMax);
    // With ROOT's default 20 colour levels, 0 sits on a level boundary and is drawn with the "later"
    // colour. 255 levels keep 0 inside its own level, so it is drawn black.
    hcam_time->SetContour(255);

    // A content of 0 means "same bin as the triggered pixel", so it cannot also mean "not in the map".
    // Pixels left out get a value below the axis minimum instead, which COLZ leaves unpainted (white).
    for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++) {
        int nx, ny;
        plottools->FindBin(pixelID, &nx, &ny);
        hcam_time->SetBinContent(nx + 1, ny + 1, zMin - 1.0);
    }

    std::vector<int> timeBin = cev->GetPeakTimeBin(); // rise bins, see FindRiseTimeBins
    for (int pixelID : pixelsToShow) {
        int nx, ny;
        plottools->FindBin(pixelID, &nx, &ny);
        hcam_time->SetBinContent(nx + 1, ny + 1, timeBin[pixelID] - triggeredTimeBin);
    }
    return hcam_time;
}


// A TH2F does not store a palette, and changing the global palette here would recolour the PDF panels.
// So each map is drawn on its own pad between two TExec objects: the first switches to
// red (earlier) -> black (same bin as trigger) -> green (later) when the pad is painted, the second
// switches back afterwards. Nothing is painted during this program; the canvas is only written to file.
TCanvas* c_time = 0;
TExec* setTimingPalette = 0;
TExec* restoreDefaultPalette = 0;

// One timing map on the current pad, with the triggered pixel marked.
void DrawTimingMapPad(TH2F* hcam_time, int triggeredPixelID)
{
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.18);
    hcam_time->Draw("colz");
    gPad->GetListOfPrimitives()->AddFirst(setTimingPalette); // palette must be set before the map is painted
    restoreDefaultPalette->Draw();
    // No MUSIC boundaries here: they are drawn in red, which would read as "earlier" on this map.

    // The triggered pixel is always delta 0 (black), like every pixel that peaks in the same bin, so
    // cover it with a grey box and a "T" to tell it apart. kCanDelete lets c_time->Clear("D") free them.
    int nx, ny;
    plottools->FindBin(triggeredPixelID, &nx, &ny);
    TBox* triggeredBox = new TBox(nx - 0.5, ny - 0.5, nx + 0.5, ny + 0.5);
    triggeredBox->SetFillColor(kGray + 1);
    triggeredBox->SetBit(kCanDelete);
    triggeredBox->Draw();
    TLatex* triggeredLabel = new TLatex(nx, ny, "T");
    triggeredLabel->SetTextAlign(22); // centred on the pixel
    triggeredLabel->SetTextSize(0.035);
    triggeredLabel->SetBit(kCanDelete);
    triggeredLabel->Draw();
}

void DrawTimingHeatMaps(TH2F* hcam_time_panel2, TH2F* hcam_time_panel3, int triggeredPixelID)
{
    if (!c_time) {
        // Divided once here; Clear("D") empties the two pads each event but keeps them.
        c_time = new TCanvas("TimingDisplay", "Rise time relative to triggered pixel", 1200, 550);
        c_time->Divide(2, 1);
        // Black is kept to a narrow band around 0 and the colour ramps up steeply just outside it,
        // so a +/-1 bin offset is clearly red or green rather than nearly black.
        setTimingPalette = new TExec("setTimingPalette",
            "{ Double_t stops[6] = {0.0, 0.44, 0.48, 0.52, 0.56, 1.0};"
            "  Double_t red[6]   = {1.0, 0.55, 0.0,  0.0,  0.0,  0.0};"
            "  Double_t green[6] = {0.0, 0.0,  0.0,  0.0,  0.55, 0.9};"
            "  Double_t blue[6]  = {0.0, 0.0,  0.0,  0.0,  0.0,  0.0};"
            "  TColor::CreateGradientColorTable(6, stops, red, green, blue, 255); }");
        restoreDefaultPalette = new TExec("restoreDefaultPalette", "gStyle->SetPalette(kBird);");
    }
    c_time->Clear("D");
    c_time->cd(1);
    DrawTimingMapPad(hcam_time_panel2, triggeredPixelID);
    c_time->cd(2);
    DrawTimingMapPad(hcam_time_panel3, triggeredPixelID);
}


void SetBranches(IEvent *evD)
{
    tree->SetBranchAddress("Events", &ev);
    // tree->SetBranchAddress("Event", &ev);
    // tree->Scan("Events.Event");
}

void SetBranchesHLED(IEvent *evD)
{
    treeHLED->SetBranchAddress("Events", &evHLED);
}

void saveEventInfo(EventInfo* evI, TTree* treeSims) {
    // TFile *fileOutput = new TFile("Test.root", "UPDATE");
    // treeSims->Branch("Cleaned","EventInfo",evI);
    treeSims->Fill();
    // fileOutput->Write();
    // fileOutput->Close();
}



// /***
// Start of Functions

// ***/



void LoadDataPCA(PCA& pca, TH2F* hist, int totalAmp){
    


    // (*COG)[0] = 0.0;
    // (*COG)[1] = 0.0;
   

    // Loop over the bins of the histogram and extract the bin content (data points)
    for (int bin = 0; bin <= MaxNofChannels; bin++) {
       
        int nx, ny;
        plottools->FindBin(bin, &nx, &ny);
        // Get the bin center for X and Y
        double weight = hist->GetBinContent(nx+1, ny+1);
         // Check if the value is infinity
        if (std::isinf(weight)) {
            weight = 0.0; // Set to 0
        }
        
        // Add the data point to the PCA object weighted by the bin content
        if (weight > 0) {  // Only add rows for non-zero bins

            std::vector<double> randomNumbers;
            for (double w = 0; w < weight; w++){
                randomNumbers = generateRandomNumbers();
                // cout << "nx: " << nx <<  " ny: " << ny << endl;
                // cout << "rand 1: " << nx + randomNumbers[0] << " rand 2: " << ny + randomNumbers[1] << endl;
                std::vector<double> data = {nx + randomNumbers[0], ny+randomNumbers[1]};
                // cout << "COG " << weight << "total" << totalAmp << endl;
                
                pca.AddRow(data);  
            }
            if (nx == 0){
                nx = 0.1;
            }
            if (ny == 0){
                ny = 0.1;
            }

            // (*COG)[0] += (nx * (weight/totalAmp));
            // (*COG)[1] += (ny * (weight/totalAmp));
        }   
    }
}

std::vector<double> generateRandomNumbers() {
    // Define the random number generator and distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-0.49, 0.49);

    // Generate two random numbers and return them as a pair
    return {dist(gen), dist(gen)};
}

std::vector<double> CreateWLRatio(PCA& pca, TVectorD& eigenVals, TMatrixD& eigenVecs) {
    eigenVals.ResizeTo(2);
    eigenVecs.ResizeTo(2,2);
    std::vector<double> sigmas;

    pca.ComputeEigenDecomposition();
    eigenVals = pca.GetEigenValues();
    eigenVecs = pca.GetEigenVectors();
    sigmas = pca.GetSigmas();                            
    
    double r1 = sqrt(eigenVals[0]);
    double r2 = sqrt(eigenVals[1]);
    std::vector<double> EllipicRatio = {r1,r2,r2/r1};  // length, width, width/length

    return EllipicRatio;
}

std::vector<int> CheckCrossPoints(TArrow* arrow, TH2F* hist,EventInfo* eventInfo) {
    // get the arrow slope
    if (VerboseThisEvent()) cout << "Checking cross points" << endl;
    double x1 = arrow->GetX1();
    double y1 = arrow->GetY1();
    double x2 = arrow->GetX2();
    double y2 = arrow->GetY2();

    // Line vector components
    double dx = x2 - x1;
    double dy = y2 - y1;
    double len = std::sqrt(dx*dx + dy*dy);

    // Degenerate line guard
    if (len == 0.0) return {0, 0};

    int cLineIntersectPixels    = 0;
    int cLineNonIntersectPixels = 0;

    int nBinsX = hist->GetNbinsX();
    int nBinsY = hist->GetNbinsY();

    double RMSMajorAxis = 0.0;
    double WeightedRMSMajorAxis = 0.0;
    double totalWeight = 0.0;

    for (int ix = 1; ix <= nBinsX; ++ix) {
        for (int iy = 1; iy <= nBinsY; ++iy) {

            if (hist->GetBinContent(ix, iy) <= 0) continue;

            // Bin center in axis coordinates
            double cx = hist->GetXaxis()->GetBinCenter(ix);
            double cy = hist->GetYaxis()->GetBinCenter(iy);

            // Perpendicular distance from bin center to the infinite line
            // d = |(P - P1) x d_hat|  (2D cross product gives scalar)
            double dist = std::abs((cy - y1)*dx - (cx - x1)*dy) / len;
            plothelp->AddtoDistance2MajorAxis(dist);
            eventInfo->SetDistance2MajorAxis(dist);
            RMSMajorAxis += dist*dist;
            WeightedRMSMajorAxis += dist*dist*hist->GetBinContent(ix, iy);
            totalWeight += hist->GetBinContent(ix, iy);
            // Half-diagonal of the bin as intersection threshold
            double hw = 0.5 * hist->GetXaxis()->GetBinWidth(ix);
            double hh = 0.5 * hist->GetYaxis()->GetBinWidth(iy);
            double threshold = 0.905*std::sqrt(hw*hw + hh*hh);
            if (VerboseThisEvent()) cout << threshold << endl;
            if (dist <= threshold) {
                ++cLineIntersectPixels;
                if (VerboseThisEvent()) cout << "Intersecting pixel at (" << cx << ", " << cy << ") with content " << hist->GetBinContent(ix, iy) << " distance " << dist << endl;

            }
            else {
                ++cLineNonIntersectPixels;
                if (VerboseThisEvent()) cout << "Non-intersecting pixel at (" << cx << ", " << cy << ") with content " << hist->GetBinContent(ix, iy) << " distance " << dist << endl;
            }
        }
    }
    RMSMajorAxis = std::sqrt(RMSMajorAxis/(cLineIntersectPixels+cLineNonIntersectPixels));
    WeightedRMSMajorAxis = std::sqrt(WeightedRMSMajorAxis/(totalWeight));
    if (VerboseThisEvent()) {
        cout << "RMS Major Axis: " << RMSMajorAxis << endl;
        cout << "Weighted RMS Major Axis: " << WeightedRMSMajorAxis << endl;
    }

    plothelp->AddtoPixelsonMajorAxis(cLineIntersectPixels);
    plothelp->AddtoPixelsoffMajorAxis(cLineNonIntersectPixels);
    plothelp->AddtoRatioPixelsMajorAxis((double)cLineIntersectPixels/(cLineIntersectPixels+cLineNonIntersectPixels));
    plothelp->AddtoRatioPixelsMajorAxisPerPixel((double)cLineIntersectPixels/(cLineIntersectPixels+cLineNonIntersectPixels),
                                                cLineIntersectPixels+cLineNonIntersectPixels);
    plothelp->AddtoRMSMajorAxis(RMSMajorAxis);
    plothelp->AddtoWeightedRMSMajorAxis(WeightedRMSMajorAxis);

    eventInfo->SetRMSMajorAxis(RMSMajorAxis);
    eventInfo->SetWeightedRMSMajorAxis(WeightedRMSMajorAxis);
    eventInfo->SetPixelsonMajorAxis(cLineIntersectPixels);
    eventInfo->SetPixelsoffMajorAxis(cLineNonIntersectPixels);
    eventInfo->SetRatioPixelsMajorAxis((double)cLineIntersectPixels/(cLineIntersectPixels+cLineNonIntersectPixels));


    return {cLineIntersectPixels, cLineNonIntersectPixels};

}

// Dots the given pixels on whatever pad is current, so saturation is visible on the camera
// image itself and not only in the header text. Call it AFTER the panel histogram has been drawn:
// TH1::Draw() without "same" clears the pad, so the next event's panel takes these markers with it
// rather than letting them pile up, which is the same lifetime the arrows and ellipse below rely on.
void DrawSaturationMarkers(const std::vector<Int_t>& pixelIDs)
{
    for (std::vector<Int_t>::size_type s = 0; s < pixelIDs.size(); s++) {
        int nx, ny;
        plottools->FindBin((int)pixelIDs[s], &nx, &ny);
        // camera bins are one unit wide and centred on the integer bin index, so the bin index
        // is the pixel centre the dot goes on
        TMarker* satDot = new TMarker(nx, ny, kFullCircle);
        satDot->SetMarkerColor(kBlack);
        satDot->SetMarkerSize(0.7); // small enough to leave the amplitude colour around it readable
        satDot->Draw();
    }
}


// Panel-1 line saying which sample the camera images are read at. Data, HLED and sims are one
// snapshot at the triggered pixel's peak bin; sims also give it as the CARE sample, the numbering
// PlotCamera's "PE at timebin" uses. Muon images are each pixel's own peak, so there is no one bin.
std::string ImageTimeBinLabel(bool usesSnapshot, int snapshotBin)
{
    if (!usesSnapshot) return "Image time bin: each pixel's own peak";
    if (whatData == "sim") return Form("Image time bin: %i (CARE sample %i)", snapshotBin, snapshotBin - SimCareTraceStart);
    return Form("Image time bin: %i", snapshotBin);
}


// ---- Horizon overlay (sims only) ----
// The same geometry PlotCamera (SimulationAnalysis/scripts/PlotCamera.h) and
// TCameraPlotEvents_Trinity_Fancy draw; keep the three in step. The camera sees 3.87 x 3.87 deg
// through 16 x 16 pixels, so one pixel is 0.242 deg and the horizon is only in the picture within
// about two degrees of tilt 90. Tilt is the zenith angle of the optical axis: 90 puts the horizon
// on the camera centre, less than 90 points (90 - tilt) deg above it.
const double kCameraFOVdeg = 3.87;
const double kCameraNPix   = 16.0;
const double kDegPerPixel  = kCameraFOVdeg / kCameraNPix;
const double kCameraCentre = (kCameraNPix - 1.0) / 2.0;   // row 7.5

// Tilt the horizon is drawn for on the current event's panels. Set per event for sims; stays 0
// for data and for a sim whose tilt could not be found, and 0 means no overlay is drawn.
double HorizonTiltDeg = 0.0;

// Tilt from the Tilt_<angle> folder in a sim file's path, e.g. .../CARE/Tilt_90.000000/x_tdem.root
// -> 90. Returns 0 when the path has no Tilt_ folder.
double TiltFromPath(const std::string& path)
{
    std::string::size_type tiltPos = path.find("Tilt_");
    if (tiltPos == std::string::npos) return 0.0;
    return atof(path.c_str() + tiltPos + 5); // atof stops at the '/' after the number
}

// Row of the camera the horizon falls on, counting pixel rows 0..15 with the frame running
// -0.5..15.5. The axis looks (90 - tilt) deg above the horizon, so the horizon sits that far below
// it. The optics invert the view but the image is flipped back before it is written, so it lands
// BELOW the centre and the sky is at the top.
double HorizonRow(double tiltDeg)
{
    return kCameraCentre - (90.0 - tiltDeg)/kDegPerPixel;
}

// Dashed horizon line plus a small arrow and "S" pointing to the sky, on the pad hist has just
// been drawn into. Works in whatever units hist is binned in by scaling one camera pixel to hist's
// axis range. Everything is made with Draw* so it carries kCanDelete and goes when the next
// event's panel is drawn. Only the pad is touched, never the histogram or the data files.
void DrawHorizonOverlay(TH2* hist, double tiltDeg, Color_t color)
{
    const double xLow  = hist->GetXaxis()->GetXmin();
    const double xHigh = hist->GetXaxis()->GetXmax();
    const double yLow  = hist->GetYaxis()->GetXmin();
    const double yHigh = hist->GetYaxis()->GetXmax();
    const double ux = (xHigh - xLow)/kCameraNPix;   // one pixel, in hist units
    const double uy = (yHigh - yLow)/kCameraNPix;

    double yHorizon = yLow + (HorizonRow(tiltDeg) + 0.5)*uy;
    bool inView = (yHorizon > yLow && yHorizon < yHigh);

    // The arrow hangs off the horizon, or off the bottom of the frame once the horizon has run
    // below it. With no room above (an all-ground image) it sits against the top edge, still
    // pointing out of frame towards the sky.
    double yAnchor   = std::max(yLow, std::min(yHorizon, yHigh));
    double arrowX    = xLow + 0.7*ux;       // inside the leftmost column
    double arrowLen  = 1.2*uy;
    double arrowEnd  = yHigh - 0.2*uy;      // highest the tip may go
    double skyFrom   = yAnchor + 0.2*uy;    // clear of the line itself
    double skyTo     = skyFrom + arrowLen;
    if (skyTo > arrowEnd){
        skyTo   = arrowEnd;
        skyFrom = arrowEnd - arrowLen;
    }

    // Dashed so it is not confused with the MUSIC boundaries.
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
    label.SetTextAlign(12);                 // left aligned, vertically centred
    label.DrawLatex(arrowX + 0.45*ux, 0.5*(skyFrom + skyTo) + 0.25*uy, "S");  // S for sky
}

// Draws the horizon on the current panel for sims with a known tilt; does nothing for data.
// White suits panel 1, where every pixel is painted; the cleaned panels leave most pixels empty
// (white), so they pass kBlack.
void DrawSimHorizon(TH2* hist, Color_t color)
{
    if (whatData != "sim" || HorizonTiltDeg == 0.0) return;
    DrawHorizonOverlay(hist, HorizonTiltDeg, color);
}

// Where the horizon lands for this tilt, and where the tilt came from, in the diagnostics format.
void PrintHorizonTilt(double tiltDeg, const std::string& tiltSource)
{
    if (!VerboseThisEvent()) return;
    if (tiltDeg == 0.0) {
        cout << "  Horizon: no tilt in the SEvent or the file path, no horizon drawn" << endl;
        return;
    }
    double horizonRow = HorizonRow(tiltDeg);
    cout << "  Horizon: tilt " << tiltDeg << " deg (from " << tiltSource << "), row "
         << Form("%.2f", horizonRow) << " of 0-15 (centre " << kCameraCentre << ")";
    if (horizonRow <= -0.5)      cout << ", below the frame: whole image is sky";
    else if (horizonRow >= 15.5) cout << ", above the frame: whole image is ground";
    cout << endl;
}


// ---- Camera edge indicator ----
// True when any of the given pixels sits on the camera's outer ring (column or row 0 or 15), so
// the image may run off the camera and its parameters describe only the part that is on it. False
// when the image is fully contained. edgePixelIDs gets the pixels that touch the ring.
bool TouchesCameraEdge(const std::vector<int>& pixelIDs, std::vector<int>& edgePixelIDs)
{
    const int lastRowOrColumn = (int)kCameraNPix - 1;
    edgePixelIDs.clear();
    for (int pixelID : pixelIDs) {
        int nx, ny;
        plottools->FindBin(pixelID, &nx, &ny);
        if (nx == 0 || ny == 0 || nx == lastRowOrColumn || ny == lastRowOrColumn) edgePixelIDs.push_back(pixelID);
    }
    return !edgePixelIDs.empty();
}

void PrintCameraEdge(bool onCameraEdge, const std::vector<int>& edgePixelIDs, int survivingPixelCount)
{
    if (!VerboseThisEvent()) return;
    if (!onCameraEdge) {
        cout << "  Camera edge: false, all " << survivingPixelCount
             << " surviving pixels are inside rows and columns 1-14" << endl;
        return;
    }
    cout << "  Camera edge: true, " << edgePixelIDs.size() << " of " << survivingPixelCount
         << " surviving pixels on the outer ring (ID (x,y)):";
    for (int pixelID : edgePixelIDs) {
        int nx, ny;
        plottools->FindBin(pixelID, &nx, &ny);
        cout << " " << pixelID << " (" << nx << "," << ny << ")";
    }
    cout << endl;
}

void CompletePanel4(PCA& pca, TH2F* hcam_panel4, CEvent* cev, std::vector<double> EllipicRatio, TVectorD& eigenVals, TMatrixD& eigenVecs,  EventInfo* eventInfo) {
    
    
    
    // The parameter text under the title needs two lines now, so pull the frame down to make
    // room. Without this the second line lands inside the frame, on top of the y axis label.
    gPad->SetTopMargin(0.20);
    hcam_panel4->Draw("colz");
    plottools->DrawMUSICBoundaries();
    DrawSimHorizon(hcam_panel4, kBlack);

    double meanx = hcam_panel4->GetMean(1);
    double meany = hcam_panel4->GetMean(2);
    // cout << "MEAN " << meanx << "," << meany << endl;
    
    double angledeg = atan(eigenVecs[1][0]/eigenVecs[0][0])*(180.0/3.141592653589793238463);
    double anglerad = atan(eigenVecs[1][0]/eigenVecs[0][0]);
    
    
    double r1 = sqrt(eigenVals[0]);
    double r2 = sqrt(eigenVals[1]);
    TArrow* arrow = new TArrow(meanx - 2*(r1*cos(anglerad)), meany - 2*(r1*sin(anglerad)), meanx + 2*(r1*cos(anglerad)), meany + 2*(r1*sin(anglerad)), 0.01, "|"); // "|>" option gives an arrowhead
    TArrow* arrow1 = new TArrow(meanx - 2*(r2*cos(anglerad+1.5708)), meany - 2*(r2*sin(anglerad+1.5708)), meanx + 2*(r2*cos(anglerad+1.5708)), meany + 2*(r2*sin(anglerad+1.5708)), 0.01, "|"); // "|>" option gives an arrowhead
    
    std::vector<int> CrossPoints = CheckCrossPoints(arrow, hcam_panel4,eventInfo);
    if (VerboseThisEvent()) {
        cout << "# Intersecting pixels with major axis: " << CrossPoints[0] << endl;
        cout << "# Non-intersecting pixels with major axis: " << CrossPoints[1] << endl;
    }
    arrow->SetLineColor(kRed); // Optional: Set the color of the arrow
    arrow->Draw("SAME");             // Draw the arrow on the same canvas
    arrow1->SetLineColor(kBlue); // Optional: Set the color of the arrow
    arrow1->Draw("SAME");             // Draw the arrow on the same canvas
    
    TEllipse* ell = new TEllipse(meanx, meany,r1,r2,0,360,angledeg);
    ell->SetFillColorAlpha(kGreen,0.00);
    ell->Draw("SAME");
    // Keep only the saturated pixels that also survived cleaning, so the stored IDs describe
    // the cleaned image. Events that never reach here keep the empty list set by Clear().
    std::vector<int> survivingPixels = cev->GetSurvivingPixelPanel3();
    std::vector<Int_t> SaturatedSurvivingIDs;
    for (std::vector<Int_t>::size_type s = 0; s < SaturatedPixelIDsRaw.size(); s++) {
        if (std::find(survivingPixels.begin(), survivingPixels.end(),
                      (int)SaturatedPixelIDsRaw[s]) != survivingPixels.end()) {
            SaturatedSurvivingIDs.push_back(SaturatedPixelIDsRaw[s]);
        }
    }

    // Show where the saturation actually is, not just how much of it there is.
    DrawSaturationMarkers(SaturatedSurvivingIDs);

    TLatex* title = new TLatex();
    title->SetNDC(); // Set to Normalized Device Coordinates (NDC)
    title->SetTextSize(0.03);
    
    float areaEllipse = (3.8/16) * (3.8/16) * (cev->GetSurvivingPixelPanel3().size()); // area of an ellipse
    double conc = (cev->GetTotalCoreAmp())/cev->GetSurvivingPixelTotalAmpPanel3();
    
    
    // Two longitudinal 3rd moments, both projected on the major axis. They do NOT agree in sign
    // and are not interchangeable -- see the header comments on each function.
    double M3Long     = getM3Long(meanx, meany, anglerad, cev->GetSurvivingPixelPanel3(),cev->GetAmplitudeValuesTimeBin());
    double M3LongPow3 = getM3LongPow3(meanx, meany, anglerad, cev->GetSurvivingPixelPanel3(),cev->GetAmplitudeValuesTimeBin());
    // Up-down symmetry: the major axis dotted with the vertical, oriented head-to-tail so the
    // value says which way the image actually points. The eigenvector alone cannot: atan()
    // forces cos(angle) > 0, so it always points toward +x and carries no up/down information.
    // M3Long supplies the missing direction. M3Long > 0 means the faint tail runs toward +x, so
    // the bright head sits toward -x and the head direction is -sign(M3Long)*(cos,sin). Dotting
    // that with the vertical gives +1 for a head pointing straight up, -1 straight down, and 0
    // for an image lying along the camera x axis.
    //
    // getM3Long returns exactly 0.0 only from its guard paths (no usable charge, or an axis too
    // close to vertical to project onto), never from a real image, so that case has no head-tail
    // direction to orient by. It is recorded as NaN rather than flattened onto 0, which would
    // otherwise read as a perfectly horizontal image. PlothUpDownSym skips those entries.
    double upDownSym = std::numeric_limits<double>::quiet_NaN();
    if (M3Long != 0.0) {
        upDownSym = (M3Long > 0.0 ? -1.0 : 1.0) * sin(anglerad);
    }

    // Whether the cleaned image reaches the camera's outer ring (true) or is fully contained (false)
    std::vector<int> edgePixelIDs;
    const bool onCameraEdge = TouchesCameraEdge(survivingPixels, edgePixelIDs);
    PrintCameraEdge(onCameraEdge, edgePixelIDs, (int)survivingPixels.size());

    // cev->GetSurvivingPixelTotalAmpPanel3()/util->GetADCtoPEratio(),
    // Split over two lines: one line no longer fits the pad width now that it carries both
    // moments. Sat is the saturated pixels that SURVIVED cleaning, matching Pixels and the
    // SaturatedPixelIDs stored on EventInfo, not the raw count seen on the trace.
    title->DrawLatex(0.1, 0.915, Form("WL:%.2f Pixels:%i Sat:%i Area:%.2f Size:%.2f Conc:%.2f",
        EllipicRatio[2],
        static_cast<int>((cev->GetSurvivingPixelPanel3()).size()),
        static_cast<int>(SaturatedSurvivingIDs.size()),
        areaEllipse,
        cev->GetSurvivingPixelTotalAmpPanel3(),
        conc));
    title->DrawLatex(0.1, 0.870, Form("M3Long:%.2f M3Pow3:%.3g UpDownSym:%.2f Edge:%s",
        M3Long,
        M3LongPow3,
        upDownSym,
        onCameraEdge ? "true" : "false"));
    delete title;

    // add values  for creating the plots
    plothelp->AddtoWL(EllipicRatio[2]);
    plothelp->AddtoL(EllipicRatio[0]);
    plothelp->AddtoW(EllipicRatio[1]);
    plothelp->AddtoSize(cev->GetSurvivingPixelTotalAmpPanel3());
    plothelp->AddtoSurvivingPixelCount(cev->GetSurvivingPixelPanel3().size());
    plothelp->AddtoConcentation(conc);
    plothelp->AddtoCOGx(meanx);
    plothelp->AddtoCOGy(meany);
    plothelp->AddtoTriggeredPixelsID(cev->GetMaxAmplitudePixelID());
    plothelp->AddtoCoreRatio(cev->GetCoreRatio());
    plothelp->AddtoNumberOfCores(cev->GetNumberofCorePixels());
    plothelp->AddtoAngle(anglerad);
    plothelp->AddtoUpDownSym(upDownSym);
    eventInfo->SetSaturatedPixels(SaturatedSurvivingIDs);
    plothelp->AddtoSaturatedPixels((int)SaturatedSurvivingIDs.size());
    // camera position of every saturated pixel that survived, one entry per pixel
    for (std::vector<Int_t>::size_type s = 0; s < SaturatedSurvivingIDs.size(); s++) {
        int satx, saty;
        plottools->FindBin((int)SaturatedSurvivingIDs[s], &satx, &saty);
        plothelp->AddtoSaturatedPixelPosition(satx, saty);
    }
    // camera position of every pixel that survived cleaning, one entry per pixel
    for (int survivingPixelID : cev->GetSurvivingPixelPanel3()) {
        int spx, spy;
        plottools->FindBin(survivingPixelID, &spx, &spy);
        plothelp->AddtoSurvivingPixelPosition(spx, spy);
    }
    
    // add info for all the events
    eventInfo->SetHPanel4(hcam_panel4);
    eventInfo->SetL(EllipicRatio[0]);
    eventInfo->SetW(EllipicRatio[1]);
    eventInfo->SetWLRatio(EllipicRatio[2]);
    eventInfo->SetAngle(anglerad);
    eventInfo->Setr1(r1);
    eventInfo->Setr2(r2);
    eventInfo->SetSize(cev->GetSurvivingPixelTotalAmpPanel3());
    eventInfo->SetSurvivngPixels(cev->GetSurvivingPixelPanel3().size());
    eventInfo->SetConc(conc);
    eventInfo->SetCOGx(meanx);
    eventInfo->SetCOGy(meany);
    eventInfo->SetTriggeredPixelID(cev->GetMaxAmplitudePixelID());
    eventInfo->SetCoreRatio(cev->GetCoreRatio());
    eventInfo->SetNumberOfCores(cev->GetNumberofCorePixels());
    eventInfo->SetArea(areaEllipse);
    eventInfo->SetM3Long(M3Long);
    eventInfo->SetM3LongPow3(M3LongPow3);
    eventInfo->SetUpDownSym(upDownSym);
    eventInfo->SetOnCameraEdge(onCameraEdge);



            
}

// Projects every surviving pixel onto the major axis and returns the signed distance d_i from the
// charge-weighted centroid together with the charge q_i. Pixels with q <= 0 are dropped: in the
// q^1 sum a negative charge flips the weighting, and in the q^3 sum it survives the cube with the
// wrong sign. Returns false when there is nothing usable to form a moment from.
//
// The major axis angle comes from atan(), so it lies in (-90,+90) and cos(theta) > 0: the unit
// vector always points toward +x, which is what makes the sign comparable across events. A
// perfectly vertical axis leaves the direction riding on the arbitrary sign of the eigenvector,
// so it is refused rather than reported as a coin flip.
bool projectOnMajorAxis(double xcog, double ycog, double anglerad, std::vector<int> sur_pix,
                        std::vector<float> amps, std::vector<double>& d, std::vector<double>& q){
    d.clear();
    q.clear();
    double ct = cos(anglerad);
    double st = sin(anglerad);
    if (std::fabs(ct) < 1e-9) {
        return false;
    }
    for (std::vector<int>::size_type p = 0; p < sur_pix.size(); p++) {
        int nx, ny;
        plottools->FindBin(sur_pix[p], &nx, &ny);
        double charge = amps[sur_pix[p]];
        if (charge <= 0) continue;
        d.push_back((nx*1.0 - xcog)*ct + (ny*1.0 - ycog)*st);
        q.push_back(charge);
    }
    return !d.empty();
}

// Standard Hillas longitudinal third moment, weighted by q^1 and normalised by the total charge:
//     M3Long = cbrt( sum(q_i * d_i^3) / sum(q_i) )
//
// The sign gives the direction of the SKEW, not of the bright end. (xcog,ycog) is the
// charge-weighted centroid, so piling charge at one end drags the centroid there and leaves a long
// faint tail at the other; the cubed distance makes that tail dominate the sum. So M3Long > 0
// means the tail runs toward +x and the brighter, denser end sits toward -x.
// Verified against a hand calculation and synthetic blobs in Claude/test_m3long.cpp.
//
// Note on normalisation: the thesis this analysis follows normalises by the pixel count N rather
// than by sum(q). Both denominators are positive and cbrt preserves sign, so that choice cannot
// change the sign and therefore cannot change any head-tail decision. It only rescales the value.
//
// cbrt is used rather than pow(x,1/3.), which returns NaN for a negative base. The sign is the
// entire point of the quantity.
double getM3Long(double xcog, double ycog, double anglerad, std::vector<int> sur_pix, std::vector<float> amps){
    std::vector<double> d, q;
    if (!projectOnMajorAxis(xcog, ycog, anglerad, sur_pix, amps, d, q)) {
        return 0.0;
    }
    double num = 0.0, qsum = 0.0;
    for (std::vector<double>::size_type i = 0; i < d.size(); i++) {
        num  += q[i] * d[i] * d[i] * d[i];
        qsum += q[i];
    }
    if (qsum <= 0) {
        return 0.0; // no charge, no moment
    }
    return std::cbrt(num / qsum);
}

// M3LongPow3: the thesis variant that weights the third moment by the THIRD POWER of the charge
// and normalises by the pixel count N:
//     M3LongPow3 = cbrt( (1/N) * sum(q_i^3 * d_i^3) )
//
// Its purpose is low-energy images, which go roundish and lose their head-tail asymmetry in the
// noise; giving bright pixels far more weight is meant to recover it.
//
// IMPORTANT: this is NOT a rescaled M3Long and its sign is not the same. Cubing the charge lets
// the bright compact head dominate the sum instead of the faint extended tail, so in practice it
// tends to come out with the OPPOSITE sign to M3Long on the same image. Decide which variable a
// head-tail cut is written against and stay with it; do not mix the two.
double getM3LongPow3(double xcog, double ycog, double anglerad, std::vector<int> sur_pix, std::vector<float> amps){
    std::vector<double> d, q;
    if (!projectOnMajorAxis(xcog, ycog, anglerad, sur_pix, amps, d, q)) {
        return 0.0;
    }
    double num = 0.0;
    for (std::vector<double>::size_type i = 0; i < d.size(); i++) {
        num += q[i] * q[i] * q[i] * d[i] * d[i] * d[i];
    }
    // N is the number of pixels actually summed, i.e. surviving pixels carrying charge.
    return std::cbrt(num / (double)d.size());
}