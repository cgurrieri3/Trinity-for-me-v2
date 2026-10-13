#include <TH2.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLine.h>
#include <TText.h>
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
#include <SEvent.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Makes an animated gif of one event: one camera image per time bin, at FramesPerSecond.
// Simulation or demonstrator data is decided by the file itself: only sim files have a "Sim" tree.



IUtilities *util;
IPlotTools *plottools;
TCanvas *canvas = 0;

std::string mnt = "/cedar/";
std::string outDir = "/cedar/DataAnalysis/MakeThingsPretty/EventGiFs/Gifs/";
std::string CalibrationFactorDir = "/cedar/DataAnalysis/flasher_calibration/Output/";

int MaxNofChannels = 256;

// Time bins shown in the gif, one frame per bin. A row of the gif list can override them.
// Data: the demonstrator pulse sits around bin 239.
const int DataFirstTimeBin = 225;
const int DataLastTimeBin  = 265;
// Sims: the 200-sample CARE trace is pasted into the 512-sample trace starting at bin 156, so the
// frames start where the CARE signal starts. The pulses themselves sit around bins 176-200.
const int SimFirstTimeBin = 156;
const int SimLastTimeBin  = 230;

const int FramesPerSecond = 6;

// Pedestals, the same ones EventCleaning measures amplitudes against:
// data uses Pulse::CalcPedestal(100), the mean of samples 0-99;
// sims use samples 157-166, because samples 0-155 are constant filler with no CARE data in them.
const int DataPedestalStart   = 0;
const int DataPedestalSamples = 100;
const int SimPedestalStart    = 157;
const int SimPedestalSamples  = 10;

// Calibration, the same as EventCleaning's StoreCalibratedAmplitudes:
// data is divided by the night's flasher factor and each pixel's SiPM gain (44 V bias);
// sims have no flasher run and use a fixed 24.1 ADC per PE.
const int   FlasherBiasVoltage = 44;
const float SimADCtoPE         = 24.1;

// Prints one line per frame (brightest pixel and camera sum) so each gif can be checked from stdout.
const bool PrintFrameDiagnostics = true;

// Sim horizon, the same geometry EventCleaning (DrawSimHorizon) and PlotCamera draw. The camera sees
// 3.87 x 3.87 deg through 16 x 16 pixels, so one pixel is 0.242 deg. Tilt is the zenith angle of the
// optical axis: 90 puts the horizon on the camera centre (row 7.5), less than 90 points the camera
// (90 - tilt) deg above it, so the horizon drops below the centre.
const double kCameraFOVdeg = 3.87;
const double kCameraNPix   = 16.0;
const double kDegPerPixel  = kCameraFOVdeg / kCameraNPix;
const double kCameraCentre = (kCameraNPix - 1.0) / 2.0;


// One row of the gif list: Output,File,Event[,FirstTimeBin,LastTimeBin]
struct GifRequest {
    std::string output;      // gif name prefix, the event number is added to it
    std::string filePath;    // full path to a merged data file or a sim _tdem.root file
    std::string event;       // Test tree entry, or "all" for every triggered event in the file
    int firstTimeBin = -1;   // -1 = use the default for the data type
    int lastTimeBin  = -1;
};


std::vector<std::vector<std::string>> splitCSVLines(const std::vector<std::string>& lines);
int extractYYYYMMDD(const std::string& text);
bool FileHasSimTree(TFile* f);
std::vector<Long64_t> MatchTestToSimEntries(TTree* testTree, TTree* simTree, SEvent* simEvent);
std::vector<Long64_t> EventsToAnimate(const GifRequest& request, TTree* testTree, TTree* simTree,
                                      SEvent* simEvent, const std::vector<Long64_t>& simEntryForTestEntry);
TH2F* LoadFlasherCalibration(int night);
std::vector<float> DataCalibrationFactors(TH2F* flasher, const std::vector<double>& gain);
float MeanOfSamples(const std::vector<Int_t>& trace, int firstSample, int nSamples);
std::vector<float> CameraAtTimeBin(const std::vector<std::vector<Int_t>>& traces, const std::vector<float>& pedestals,
                                   const std::vector<float>& calibrationFactors, int timeBin);
void FillCamera(TH2F* camera, const std::vector<float>& pePerPixel, bool flipBothAxes);
double TiltFromPath(const std::string& path);
double HorizonRow(double tiltDeg);
double DataSkyGroundLineY(int night);
void DrawSkyGroundLine(double skyGroundLineY);
int FrameDelayCentiseconds(int frame);
