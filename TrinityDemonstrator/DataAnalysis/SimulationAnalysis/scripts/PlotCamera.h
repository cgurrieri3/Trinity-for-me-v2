#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TBox.h>
#include <set>
#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <typeinfo>
#include "TLine.h"
#include <fstream>
#include <sstream>
#include <TPaletteAxis.h>
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
#include <sys/stat.h>
#include <vector>
#include <TLatex.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <Event.h>
#include <BiFocal.h>
#include <Pulse.h>
#include <TBox.h>
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
#include <Event.h>
#include <IEvent.h>
#include <IPlotTools.h>
#include <ISiPM.h>
#include <SEvent.h>




IUtilities *util;
IPlotTools *plottools;
SEvent *simRun=0;
IEvent *ev=0;
CEvent *cev;
TFile *file;
TFile *fileOutput;
TFile *fileOutputPlots;
TTree *treeSims;
TTree *treeNew;

TH1F *hPhotonX=0;
TH1F *hPhotonY=0;
TH1F *hPhotonDcosX=0;
TH1F *hPhotonDcosY=0;
TH1F *hPhotonZ=0;
TH2F* hPhotonDist1x=0;
TH2F* hPhotonDist1x6mm=0;
TH2F* hPhotonDist4x=0;

TCanvas *canvas=0;



std::vector<float>* vPhotonX = 0;
std::vector<float>* vPhotonY = 0;
std::vector<float>* vPhotonDcosX = 0;
std::vector<float>* vPhotonDcosY = 0;
int nResMult = 100;

int CameraStep = 16;
float CameraMin = -0.5;
float CameraMax = 15.5;
float RadtoDeg = 57.2958;
float PtoPEconversion = 0.4*0.83;
// 0.4? photon detection effiency 
// 0.83 is the average mirror reflectivity 
// simulations accounts for atmopsheric absorbtion

// ---------------------------------------------------------------------------
// Simulation source. Selected with the 4th argument, GROPT or CARE.
//   GROPT -> ray traced photon positions, tree T1, photonX/photonY/...
//   CARE  -> simulated FADC traces,       tree Events/T0, vFADCTraces0..255
// The two produce different histogram contents but share every label, the
// ISOG derived shower information and the output file/tree layout.
// ---------------------------------------------------------------------------
std::string simType = "GROPT";

// CARE readout constants (taken from the Sofia_CARE_Simanalysis branch)
const int   kNPixels     = 256;   // one vFADCTraces<N> branch per pixel
const int   kCareTimeBin = 38;    // fallback sample, used only on a flat camera
const int   kCarePedBins = 10;    // leading samples averaged for the pedestal
const float kCareADCtoPE = 24.1;  // ADC counts per photo electron

std::vector<int>* vFADCTraces[kNPixels] = {};

// Per pixel CARE readout dump, one line a pixel and 256 lines an event, so it
// is off by default and turned on with the CAREVERBOSE environment variable
// when the pedestals of a single run have to be inspected by hand.
bool careVerbose = false;

// Pedestal of one CARE trace: the average of its leading samples.
//
// This stays the RAW pedestal. The demonstrator's ~3400 ADC baseline is added
// where a demonstrator style trace is built (the corrected overlay in the event
// loop and create_root_file_care), not here: the camera image needs the signal
// measured above the raw baseline, and folding 3400 in here drives every pixel
// negative so nothing ever triggers.
int GetCarePedestal(const std::vector<int>& trace)
{
        int nPed = std::min(kCarePedBins, (int)trace.size());
        if (nPed <= 0) return 0;
        int pedestal = 0;
        for (int t = 0; t < nPed; t++) pedestal += trace[t];
        return pedestal / nPed;
}

// Where in the trace the camera image is taken.
//
// A camera image is one FADC sample of every pixel, so which sample is read out
// decides what the image shows. CARE puts the pulse wherever the shower's
// arrival time puts it, and a sample fixed in advance therefore lands on the
// pulse for some events and on bare baseline for the rest - those come out
// empty and are written as untriggered even though the shower is in the trace.
// The sample is found per event instead: the brightest sample of the brightest
// pixel, which is the pixel that fires the trigger. Every pixel is then read out
// at that one sample, so the image stays a simultaneous snapshot of the camera
// and does not become a per pixel peak map.
struct CareSnapshot {
        int   bin;        // sample every pixel is read out at
        int   pixel;      // pixel that peaks there, -1 if the camera is flat
        float amplitude;  // that pixel's pedestal subtracted amplitude, in ADC
};

CareSnapshot FindCareSnapshot(std::vector<int>* traces[], const int pedestals[])
{
        CareSnapshot snapshot;
        snapshot.bin       = kCareTimeBin;
        snapshot.pixel     = -1;
        snapshot.amplitude = 0.0f;

        for (int pix = 0; pix < kNPixels; pix++){
                if (!traces[pix]) continue;
                const std::vector<int>& trace = *traces[pix];
                // The pedestal samples are skipped: a pulse sitting inside them
                // has already spoiled that pixel's pedestal, so a peak found
                // there is measured against itself and means nothing.
                for (int t = kCarePedBins; t < (int)trace.size(); t++){
                        float amplitude = trace[t] - pedestals[pix];
                        if (amplitude > snapshot.amplitude){
                                snapshot.amplitude = amplitude;
                                snapshot.bin       = t;
                                snapshot.pixel     = pix;
                        }
                }
        }
        return snapshot;
}

// Size of the event label drawn over the displays. The label is split across
// two lines, so each half fits across a PDF page at this size.
const float kLabelTextSize = 0.025;

// Pad margins. The top one is opened up so the two line label sits above the
// frame instead of over the top of the image. The histogram title still draws
// above the label, hence the gap left between kLabelY2 and 1.
const float kPadLeft   = 0.10;
const float kPadRight  = 0.15;   // room for the colour palette
const float kPadBottom = 0.10;
const float kPadTop    = 0.18;   // header band for the label
const float kLabelY1   = 0.83;   // just above the frame at 1 - kPadTop = 0.82
const float kLabelY2   = 0.92;   // just below the histogram title

// Histogram title size. ROOT auto-shrinks long titles to fit the pad, which
// made them unreadable, so it is pinned here instead.
const float kTitleFontSize = 0.030;

// Gap between the colour palette's numbers and its axis title. At the ROOT
// default of 1.0 the title sits on top of the numbers.
const float kZTitleOffset = 1.4;

// Canvas geometry. Square plot area needs W*(1-kPadLeft-kPadRight) to equal
// H*(1-kPadBottom-kPadTop), i.e. W = H*0.72/0.75 = 0.96*H per pad. Height is
// set first so each pad stays wide enough for the title.
const int kPadPixelsH  = 1300;                 // canvas height
const int kPadPixelsW  = 1248;                 // one square pad (0.96 * 1300)

// Usable drawing area on the A4 sheet, in cm. The canvas is fitted inside this
// at its own aspect so it fills the page without running off the edge.
const double kPdfAreaLongCm  = 28.0;
const double kPdfAreaShortCm = 19.0;

// Trigger threshold in PE on the brightest pixel, shared by both pipelines.
const float kTriggerPE = 20.0;

// Tilt sub-directory. The angle changes from campaign to campaign (91.560000
// under daviddata, 90.000000 under sofiadata) and the number of decimals is not
// fixed either, so the folder is found by its "Tilt_" prefix rather than being
// matched against a list of spellings.
const char* kTiltPrefix    = "Tilt_";
const size_t kTiltPrefixLen = 5;


// std::string dataDir = "/storage/hive/project/phy-otte/shared/Trinity/Simulations/TDemSims/";
std::string mnt = "/storage/osg-otte1/shared/TrinityDemonstrator/";
std::string dataDir = "/storage/osg-otte1/shared/TrinityDemonstrator/simdata/roots/diffuse/";
std::string outDir = "/storage/osg-otte1/shared/TrinityDemonstrator/DataAnalysis/SimulationAnalysis/data/";

// std::string OutputFileRoot=Form("%sSimulations20250529.root",outDir.c_str());

// Define a struct to represent one simulation event
// struct SimEventData {
//     float tb_time;
//     float cobo_time;
//     float unix_time;
//     float eventType;
//     float signalValue;
//     float roiMusicID;
//     float roiPixelID;
// };

// struct SimRunData{
//     float energy;
//     float time;
//     float posX;
//     float posY;
//     // std::vector<float> float posZ;
//     float zenith;
//     float azimuth;
//     float bsm;
//     float tilt;
// };

// struct SimAuxData {
//     Float_t IcorrectionFactor;
//     Float_t sunAzimuth;
//     Float_t sunElevation;
//     Float_t moonAzimuth;
//     Float_t moon_Elevation;
//     Float_t moon_Illumination;
//     Float_t CameraRA;
//     Float_t CameraDEC;
//     Float_t TrigEvent;
//     Float_t hv;
//     Float_t hvc;
//     Float_t ucTemp;
//     Float_t siabMPWR;
//     Float_t hvSW;
//     Float_t sipmTemp;
//     Float_t siabCurr;
//     Float_t ASADcurr;
//     Float_t TBcurr;
//     Float_t outsideTempature;
//     Float_t humidity;
// };
// One piece of work: a single <run>/<stage>/<tilt>/ folder paired with the ISOG
// text file describing the same showers. A flat run folder yields one of these,
// a campaign folder yields one per run (and per tilt) underneath it.
struct RunUnit {
    std::string name;       // run folder name, e.g. 20260901_10455802_ss
    std::string tilt;       // tilt folder name, e.g. Tilt_90.000000
    std::string stagePath;  // <run>/<stage>/<tilt>/, holds the .root files
    std::string isogPath;   // <run>/ISOG/<tilt>/<shower>.txt
    int         isogRows;   // shower rows in that ISOG file, -1 until counted

    RunUnit() : isogRows(-1) {}
};

// One .root file to read, together with the run it belongs to. The campaign is
// flattened into a list of these so the event loop stays a single pass over all
// of the runs instead of a loop nested inside a loop.
struct FileJob {
    size_t      unit;   // index into the RunUnit list
    std::string file;   // file name inside that unit's stagePath
};

void DrawMUSICBoundaries();
void create_root_file(TH2F* hist, TTree* treeNew);
void create_root_file_care(TH2F* hist, TTree* treeNew, std::vector<int>* traces[]);
void SetPdfPageToCanvas(TCanvas *c);
std::ifstream& GotoLine(std::ifstream& file, unsigned int num);
bool IsDirectory(const std::string &path);
std::vector<std::string> ListSubDirs(const std::string &dir, const std::string &prefix);
std::vector<std::string> ListFilesWithExt(const std::string &dir, const std::string &ext);
std::vector<std::string> ListTiltDirs(const std::string &runDir, const std::string &stage);
std::string ResolveTiltDir(const std::string &runDir, const std::string &stage);
std::string FindISOGFile(const std::string &runDir, const std::string &tilt);
std::vector<RunUnit> CollectRunUnits(const std::string &topDir, const std::string &stage);
bool LoadISOG(SEvent *simRun, const std::string &ISOGFilePath, int showerIndex);
bool CheckISOGHeader(const std::string &ISOGFilePath);
int  CountISOGRows(const std::string &ISOGFilePath);


bool IsDirectory(const std::string &path){
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    return S_ISDIR(st.st_mode);
}

// Sub-directories of <dir> whose name starts with <prefix>, sorted so a campaign
// is walked in the same order every time. dirent's d_type is not filled in on
// every filesystem, so each candidate is stat'ed rather than trusted.
std::vector<std::string> ListSubDirs(const std::string &dir, const std::string &prefix){
    std::vector<std::string> names;
    DIR *d = opendir(dir.c_str());
    if (!d) return names;
    struct dirent *entry;
    while ((entry = readdir(d)) != NULL){
        std::string name = entry->d_name;
        if (name == "." || name == "..") continue;
        if (!prefix.empty() && name.compare(0, prefix.size(), prefix) != 0) continue;
        if (!IsDirectory(dir + name)) continue;
        names.push_back(name);
    }
    closedir(d);
    std::sort(names.begin(), names.end());
    return names;
}

// Files in <dir> ending in <ext>, sorted. Sorting matters now that one job walks
// hundreds of run folders: the event order has to be reproducible.
std::vector<std::string> ListFilesWithExt(const std::string &dir, const std::string &ext){
    std::vector<std::string> names;
    DIR *d = opendir(dir.c_str());
    if (!d) return names;
    struct dirent *entry;
    while ((entry = readdir(d)) != NULL){
        std::string name = entry->d_name;
        if (name.size() <= ext.size()) continue;
        if (name.compare(name.size()-ext.size(), ext.size(), ext) != 0) continue;
        if (IsDirectory(dir + name)) continue;
        names.push_back(name);
    }
    closedir(d);
    std::sort(names.begin(), names.end());
    return names;
}

// Every Tilt_* folder under <runDir>/<stage>/.
std::vector<std::string> ListTiltDirs(const std::string &runDir, const std::string &stage){
    return ListSubDirs(runDir + stage + "/", kTiltPrefix);
}

// The first Tilt_* folder under <runDir>/<stage>/, or an empty string when the
// stage folder has none. Kept for callers that only ever want one tilt.
std::string ResolveTiltDir(const std::string &runDir, const std::string &stage){
    std::vector<std::string> tilts = ListTiltDirs(runDir, stage);
    if (tilts.empty()) return std::string();
    return runDir + stage + "/" + tilts[0] + "/";
}

// The ISOG text file for <runDir> at <tilt>. ISOG is normally written under the
// same tilt as the stage being read; if that folder is missing or empty, any
// other ISOG tilt folder in the same run is used instead.
std::string FindISOGFile(const std::string &runDir, const std::string &tilt){
    std::vector<std::string> candidates;
    if (!tilt.empty()) candidates.push_back(runDir + "ISOG/" + tilt + "/");
    std::vector<std::string> isogTilts = ListTiltDirs(runDir, "ISOG");
    for (size_t i = 0; i < isogTilts.size(); i++)
        candidates.push_back(runDir + "ISOG/" + isogTilts[i] + "/");

    for (size_t i = 0; i < candidates.size(); i++){
        std::vector<std::string> txt = ListFilesWithExt(candidates[i], ".txt");
        if (!txt.empty()) return candidates[i] + txt[0];
    }
    return std::string();
}

// Builds the list of folders one invocation has to process under <topDir>.
//
// Two layouts are handled:
//   flat      <topDir>/<stage>/Tilt_*/         one run, the original layout
//   campaign  <topDir>/<run>/<stage>/Tilt_*/   many runs side by side
// The campaign form is what simdata/sofiadata/<date>_<name>/ looks like now, so
// every run folder underneath it is picked up automatically instead of being
// named on the command line. Runs with no <stage> folder - the ones that never
// finished - simply contribute nothing; runs that have one but no ISOG file are
// reported by the caller as it walks the list.
std::vector<RunUnit> CollectRunUnits(const std::string &topDir, const std::string &stage){
    std::vector<RunUnit> units;
    if (!IsDirectory(topDir)) return units;

    // flat layout: the stage folder sits directly under topDir
    std::vector<std::string> tilts = ListTiltDirs(topDir, stage);
    if (!tilts.empty()){
        std::string name = topDir;
        if (!name.empty() && name[name.size()-1] == '/') name.erase(name.size()-1);
        size_t slash = name.rfind('/');
        if (slash != std::string::npos) name = name.substr(slash+1);
        for (size_t t = 0; t < tilts.size(); t++){
            RunUnit u;
            u.name      = name;
            u.tilt      = tilts[t];
            u.stagePath = topDir + stage + "/" + tilts[t] + "/";
            u.isogPath  = FindISOGFile(topDir, tilts[t]);
            units.push_back(u);
        }
        return units;
    }

    // campaign layout: one level of run folders in between
    std::vector<std::string> runs = ListSubDirs(topDir, "");
    for (size_t r = 0; r < runs.size(); r++){
        std::string runDir = topDir + runs[r] + "/";
        std::vector<std::string> runTilts = ListTiltDirs(runDir, stage);
        for (size_t t = 0; t < runTilts.size(); t++){
            RunUnit u;
            u.name      = runs[r];
            u.tilt      = runTilts[t];
            u.stagePath = runDir + stage + "/" + runTilts[t] + "/";
            u.isogPath  = FindISOGFile(runDir, runTilts[t]);
            units.push_back(u);
        }
    }
    return units;
}


// ---------------------------------------------------------------------------
// ISOG row layout.
//
// Column order is fixed by ISOG's own second header line:
//   CHASM: PolarAngle(rad) AzimuthAngle(rad) XTel(m) YTel(m) ZTel(m)
//          TelescopeRadius(m) XMax(g/cm^-2) NMax(Nparticles) X0(g/cm^-2)
//          Lambda NeutrinoEnergy(GeV) DecayLength(m) ShowerStartX(m)
//          ShowerStartY(m) ShowerStartZ(m) TauEmergenceProb EmergencePointX(m)
//          EmergencePointY(m) EmergencePointZ(m) GlobeTheta(rad) GlobePhi(rad)
//          TauEnergy(GeV) PExit LastSurfaceCode RunCHASMFlag MasterSeed
//
// The names are spelled out here so the mapping into the SEvent can be checked
// against the file instead of being taken on trust; CheckISOGHeader does that
// once per run. Only the first kISOGNColumns are stored: SEvent has no field
// for PExit, LastSurfaceCode, RunCHASMFlag or MasterSeed, and MasterSeed is a
// 128 bit integer that would not survive a Float_t anyway.
// ---------------------------------------------------------------------------
enum EISOGColumn {
    kISOGPolarAngle = 0, kISOGAzimuthAngle,   kISOGXTel,             kISOGYTel,
    kISOGZTel,           kISOGTelRadius,      kISOGXMax,             kISOGNMax,
    kISOGX0,             kISOGLambda,         kISOGNeutrinoEnergy,   kISOGDecayLength,
    kISOGShowerStartX,   kISOGShowerStartY,   kISOGShowerStartZ,     kISOGTauEmergenceProb,
    kISOGEmergencePointX,kISOGEmergencePointY,kISOGEmergencePointZ,  kISOGGlobeTheta,
    kISOGGlobePhi,       kISOGTauEnergy,
    kISOGNColumns
};

const char* kISOGColumnNames[kISOGNColumns] = {
    "PolarAngle",    "AzimuthAngle",    "XTel",            "YTel",
    "ZTel",          "TelescopeRadius", "XMax",            "NMax",
    "X0",            "Lambda",          "NeutrinoEnergy",  "DecayLength",
    "ShowerStartX",  "ShowerStartY",    "ShowerStartZ",    "TauEmergenceProb",
    "EmergencePointX","EmergencePointY","EmergencePointZ", "GlobeTheta",
    "GlobePhi",      "TauEnergy"
};

// Number of header lines above the first shower row.
const int kISOGHeaderLines = 2;

// Per shower ISOG dump. One shower prints 22 lines, so a 40000 shower campaign
// would bury the condor log; off by default and turned on with the ISOGVERBOSE
// environment variable when a single run needs to be inspected by hand.
bool isogVerbose = false;

// Normalises an ISOG column name for comparison: drops the unit in parentheses
// and any spaces, and lowercases the rest. "Telescope Radius(m)" and
// "TelescopeRadius(m)" are the same column written by two ISOG versions, and
// both have to compare equal to "TelescopeRadius".
std::string NormaliseISOGName(const std::string &raw){
    std::string out;
    for (size_t i = 0; i < raw.size(); i++){
        if (raw[i] == '(') break;          // unit, not part of the name
        if (isspace((unsigned char)raw[i])) continue;
        out += tolower((unsigned char)raw[i]);
    }
    return out;
}

// Splits an ISOG header line into column names. ISOG writes the header with
// U+2423 (open box, bytes e2 90 a3) between the names rather than a space, so
// that sequence is treated as a separator alongside ordinary whitespace.
std::vector<std::string> SplitISOGHeader(const std::string &line){
    std::string flat = line;
    const std::string openBox = "\xe2\x90\xa3";
    size_t pos = 0;
    while ((pos = flat.find(openBox, pos)) != std::string::npos)
        flat.replace(pos, openBox.size(), " ");

    std::vector<std::string> names;
    std::istringstream iss(flat);
    std::string tok;
    while (iss >> tok) names.push_back(tok);
    // the line opens with the literal "CHASM:", which is not a column
    if (!names.empty() && NormaliseISOGName(names[0]) == "chasm:") names.erase(names.begin());
    return names;
}

// Confirms that the columns this code reads really are the columns the file
// carries, by name, before a single shower is filled. Checking the first
// kISOGNColumns only: past TauEnergy the header has a name with a space inside
// it ("...nan=no topography") which would split into two tokens, and those
// columns are not stored anyway. Returns false on a mismatch so the caller can
// refuse the run instead of writing a tree full of misassigned energies.
bool CheckISOGHeader(const std::string &ISOGFilePath){
    std::ifstream f(ISOGFilePath.c_str());
    if (!f.is_open()){
        cout << "Could not open ISOG file to check its header: " << ISOGFilePath << endl;
        return false;
    }
    std::string line;
    std::getline(f, line);   // "The First 11 Arguments are for CHASM..."
    std::getline(f, line);   // the column names
    std::vector<std::string> names = SplitISOGHeader(line);

    if ((int)names.size() < kISOGNColumns){
        cout << "ISOG header lists only " << names.size() << " columns, expected at least "
             << kISOGNColumns << ": " << ISOGFilePath << endl;
        return false;
    }
    bool ok = true;
    for (int c = 0; c < kISOGNColumns; c++){
        if (NormaliseISOGName(names[c]) != NormaliseISOGName(kISOGColumnNames[c])){
            cout << "ISOG column " << c << " is '" << names[c] << "', expected '"
                 << kISOGColumnNames[c] << "'" << endl;
            ok = false;
        }
    }
    if (ok) cout << "ISOG header check passed (" << names.size() << " columns, "
                 << kISOGNColumns << " stored)" << endl;
    return ok;
}

// Number of shower rows in an ISOG file: everything below the header. The event
// loop indexes this file by tree entry, so the caller compares this against the
// tree and drops the run if the two disagree rather than shifting every shower's
// energy onto the wrong event.
int CountISOGRows(const std::string &ISOGFilePath){
    std::ifstream f(ISOGFilePath.c_str());
    if (!f.is_open()) return -1;
    int nLines = 0;
    std::string line;
    while (std::getline(f, line)){
        // a trailing newline leaves one empty line that is not a shower
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
        nLines++;
    }
    return nLines - kISOGHeaderLines;
}

// One ISOG field. ISOG writes "nan" wherever a quantity does not apply to a
// shower - ShowerStartX/Y/Z on anything that never decayed, LastSurfaceCode
// with no topography - and reading the row with "iss >> Float_t" stops dead on
// that token: the stream sets failbit, the extraction loop ends, and every
// column after it is left unread. That is why tau energy, column 21, came back
// missing on 20260916: the first nan sits at column 12, so only 12 of the 26
// columns were ever parsed. Converting one token at a time keeps a nan confined
// to its own column, and anything unparsable becomes NaN rather than ending the
// row.
Float_t ParseISOGValue(const std::string &token){
    const char* start = token.c_str();
    char* end = 0;
    double value = strtod(start, &end);
    if (end == start) return std::numeric_limits<Float_t>::quiet_NaN();
    return (Float_t)value;
}

// Reads one shower row out of the ISOG text file into the SEvent. The file has
// two header lines, so event <showerIndex> lives on line <showerIndex>+3.
// Shared by GROPT and CARE so both carry identical shower information.
//
// Every field is written on every call, with NaN standing in for a column the
// row does not have. The SEvent is reused across events and has no reset, so a
// field left unset would silently keep the previous shower's value - which is
// how a short row turns into wrong numbers rather than missing ones.
bool LoadISOG(SEvent *simRun, const std::string &ISOGFilePath, int showerIndex){
    const Float_t kNaN = std::numeric_limits<Float_t>::quiet_NaN();
    std::vector<Float_t> v(kISOGNColumns, kNaN);

    std::ifstream ISOGfile(ISOGFilePath.c_str());
    bool ok = true;
    if(!ISOGfile.is_open()){
        cout << "Could not open ISOG file: " << ISOGFilePath << endl;
        ok = false;
    }

    int nParsed = 0;
    if (ok){
        std::ifstream& file = GotoLine(ISOGfile, showerIndex + kISOGHeaderLines + 1);
        std::string line;
        std::getline(file, line);
        if (isogVerbose) cout << "ISOG line " << showerIndex + kISOGHeaderLines + 1 << ": " << line << endl;

        std::istringstream iss(line);
        std::string token;
        while (iss >> token){
            if (nParsed < kISOGNColumns) v[nParsed] = ParseISOGValue(token);
            nParsed++;   // keep counting past the stored columns so the total is honest
        }
        if (nParsed == 0){
            cout << "ISOG row empty for shower " << showerIndex << " in " << ISOGFilePath << endl;
            ok = false;
        } else if (nParsed < kISOGNColumns){
            cout << "ISOG row " << showerIndex << " has only " << nParsed << " columns, expected "
                 << kISOGNColumns << "; the rest are set to NaN" << endl;
        }
    }

    simRun->SetPolarAngle       (v[kISOGPolarAngle]);
    simRun->SetAzimuthAngle     (v[kISOGAzimuthAngle]);
    simRun->SetTelescope_Xpos   (v[kISOGXTel]);
    simRun->SetTelescope_Ypos   (v[kISOGYTel]);
    simRun->SetTelescope_Zpos   (v[kISOGZTel]);
    simRun->SetTelescope_Radius (v[kISOGTelRadius]);
    simRun->SetXMax             (v[kISOGXMax]);
    simRun->SetNMax             (v[kISOGNMax]);
    simRun->SetX0Max            (v[kISOGX0]);
    simRun->SetLambda           (v[kISOGLambda]);
    simRun->SetNeutrinoEnergy   (v[kISOGNeutrinoEnergy]);
    simRun->SetDecayLength      (v[kISOGDecayLength]);
    simRun->SetShowerStartX     (v[kISOGShowerStartX]);
    simRun->SetShowerStartY     (v[kISOGShowerStartY]);
    simRun->SetShowerStartZ     (v[kISOGShowerStartZ]);
    simRun->SetTauEmergenceProb (v[kISOGTauEmergenceProb]);
    simRun->SetEmergencePointX  (v[kISOGEmergencePointX]);
    simRun->SetEmergencePointY  (v[kISOGEmergencePointY]);
    simRun->SetEmergencePointZ  (v[kISOGEmergencePointZ]);
    simRun->SetGlobeTheta       (v[kISOGGlobeTheta]);
    simRun->SetGlobePhi         (v[kISOGGlobePhi]);
    simRun->SetTauEnergy        (v[kISOGTauEnergy]);

    if (isogVerbose){
        cout << "  ISOG shower " << showerIndex << " (" << nParsed << " columns read)" << endl;
        for (int c = 0; c < kISOGNColumns; c++)
            cout << "    " << kISOGColumnNames[c] << ": " << v[c] << endl;
    }

    return ok;
}

std::ifstream& GotoLine(std::ifstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}


void DrawMUSICBoundaries()
{
        TBox *b = new TBox(-48.5,-48.5,-36.5,-24.5);
        b->SetFillStyle(0);
        b->SetLineColor(kRed);
        for(int i=1; i<32; i++)
        {
                TBox *bn = (TBox*)b->Clone();
                bn->SetX1((i%8)*12-48.5);
                bn->SetX2((i%8)*12-36.5);
                bn->SetY1((i/8)*24-48.5);
                bn->SetY2((i/8)*24-24.5);
                bn->Draw();
        }
}