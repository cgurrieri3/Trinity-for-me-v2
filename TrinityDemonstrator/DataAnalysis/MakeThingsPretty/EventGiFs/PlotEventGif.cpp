#include "PlotEventGif.h"

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

// Only simulation files carry the "Sim" (GrOptics / SEvent) truth tree.
bool FileHasSimTree(TFile* f){
    TTree *simTree = (TTree*)f->Get("Sim");
    return simTree != nullptr;
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
    cout << "Sim entries: " << simTree->GetEntries()
         << "  triggered: " << triggeredSimEntries.size()
         << "  Test entries: " << nTest << endl;

    std::vector<Long64_t> simEntryForTestEntry;
    if ((Long64_t)triggeredSimEntries.size() == nTest){
        cout << "Sim layout: Test holds triggered showers only -> Test entry k = k-th triggered Sim entry" << endl;
        simEntryForTestEntry = triggeredSimEntries;
    } else if (simTree->GetEntries() == nTest){
        cout << "Sim layout: Test holds every thrown shower -> Test entry k = Sim entry k" << endl;
        for (Long64_t testEntry = 0; testEntry < nTest; testEntry++) simEntryForTestEntry.push_back(testEntry);
    } else {
        cout << "WARNING: Test entries match neither Sim layout; the gifs get no run name or trigger flag" << endl;
    }
    return simEntryForTestEntry;
}

// The Test entries to make gifs of. A number is used as given. "all" means every event that triggered:
// every Test entry for data, and for sims only the entries whose Sim truth is flagged triggered
// (a parallel-layout sim file also holds the showers that never triggered, and their traces are flat).
std::vector<Long64_t> EventsToAnimate(const GifRequest& request, TTree* testTree, TTree* simTree,
                                      SEvent* simEvent, const std::vector<Long64_t>& simEntryForTestEntry){
    std::vector<Long64_t> events;
    if (request.event != "all"){
        events.push_back(std::atoll(request.event.c_str()));
        return events;
    }

    for (Long64_t testEntry = 0; testEntry < testTree->GetEntries(); testEntry++){
        if (simTree && testEntry < (Long64_t)simEntryForTestEntry.size()){
            simTree->GetEntry(simEntryForTestEntry[testEntry]);
            if (!simEvent->GetTriggered()) continue;
        }
        events.push_back(testEntry);
    }
    return events;
}

// The night's flasher calibration. Newer files name the histogram by bias voltage (CamFlasher44),
// 2025 files have a single CamFlasher. Returns nullptr when the night has no usable file.
TH2F* LoadFlasherCalibration(int night){
    std::string flasherPath = Form("%s%d_FlasherCalibration_Factor.root", CalibrationFactorDir.c_str(), night);
    TFile *flasherFile = TFile::Open(flasherPath.c_str());
    if (!flasherFile || flasherFile->IsZombie()){
        cerr << "Error: could not open flasher calibration " << flasherPath << endl;
        return nullptr;
    }

    std::string histName = Form("CamFlasher%d", FlasherBiasVoltage);
    TH2F *flasher = (TH2F*)flasherFile->Get(histName.c_str());
    if (!flasher){
        histName = "CamFlasher";
        flasher = (TH2F*)flasherFile->Get(histName.c_str());
    }
    if (!flasher){
        cerr << "Error: no CamFlasher" << FlasherBiasVoltage << " or CamFlasher histogram in " << flasherPath << endl;
        flasherFile->Close();
        return nullptr;
    }

    flasher->SetDirectory(0);
    flasherFile->Close();
    cout << "Flasher calibration: " << flasherPath << " (" << histName << ")" << endl;
    return flasher;
}

// One factor per pixel that turns ADC into PE: 1 / (flasher factor x SiPM gain), the same as
// CEvent::SetAmplitudeValuesTimeBin. Pixels without a finite factor (NaN flasher nights, zero gain) get 0.
std::vector<float> DataCalibrationFactors(TH2F* flasher, const std::vector<double>& gain){
    std::vector<float> factors(MaxNofChannels, 0.0);
    if ((int)gain.size() < MaxNofChannels){
        cerr << "Error: the SiPM branch has " << gain.size() << " gains, expected " << MaxNofChannels << endl;
        return factors;
    }

    int pixelsWithoutCalibration = 0;
    for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++){
        int nx, ny;
        IPlotTools::FindBin(pixelID, &nx, &ny);
        const double factor = 1.0 / (flasher->GetBinContent(nx + 1, ny + 1) * gain[pixelID]);
        if (!std::isfinite(factor)){
            pixelsWithoutCalibration++;
            continue;
        }
        factors[pixelID] = factor;
    }

    if (pixelsWithoutCalibration > 0){
        cout << "WARNING: " << pixelsWithoutCalibration << " pixels have no finite calibration factor"
             << " (NaN flasher bin or zero gain) and are drawn as 0 PE" << endl;
    }
    return factors;
}

float MeanOfSamples(const std::vector<Int_t>& trace, int firstSample, int nSamples){
    double sum = 0.0;
    for (int sample = firstSample; sample < firstSample + nSamples; sample++){
        sum += trace[sample];
    }
    return sum / nSamples;
}

// Every pixel's PE at one time bin. The pulse goes down from the pedestal, and a sample above the
// pedestal reads as 0, the same convention as EventCleaning's AmplitudeAtBin.
std::vector<float> CameraAtTimeBin(const std::vector<std::vector<Int_t>>& traces, const std::vector<float>& pedestals,
                                   const std::vector<float>& calibrationFactors, int timeBin){
    std::vector<float> pePerPixel(MaxNofChannels, 0.0);
    for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++){
        const float amplitudeADC = pedestals[pixelID] - traces[pixelID][timeBin];
        if (amplitudeADC > 0) pePerPixel[pixelID] = amplitudeADC * calibrationFactors[pixelID];
    }
    return pePerPixel;
}

// Puts every pixel's PE in its camera bin. The unflipped bin is (nx + 1, ny + 1), the layout
// EventCleaning draws. Data is flipped in both camera axes to match the CreateCameraImageFromHist
// stills (makeSkyTopHistogram there), which puts the sky on top. Sims are written sky-up already,
// so they go in unflipped.
void FillCamera(TH2F* camera, const std::vector<float>& pePerPixel, bool flipBothAxes){
    const int nBins = camera->GetNbinsX();
    for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++){
        int nx, ny;
        IPlotTools::FindBin(pixelID, &nx, &ny);
        if (flipBothAxes){
            // the mirror of bin b is nBins + 1 - b
            camera->SetBinContent(nBins - nx, nBins - ny, pePerPixel[pixelID]);
        } else {
            camera->SetBinContent(nx + 1, ny + 1, pePerPixel[pixelID]);
        }
    }
}

// Tilt from the Tilt_<angle> folder in a sim file's path, e.g. .../CARE/Tilt_90.000000/x_tdem.root
// -> 90. Returns 0 when the path has no Tilt_ folder.
double TiltFromPath(const std::string& path){
    std::string::size_type tiltPos = path.find("Tilt_");
    if (tiltPos == std::string::npos) return 0.0;
    return std::atof(path.c_str() + tiltPos + 5); // atof stops at the '/' after the number
}

// Camera y (rows 0-15, frame -0.5 to 15.5) the sim horizon falls on, the same as EventCleaning's
// HorizonRow. The axis looks (90 - tilt) deg above the horizon, so the horizon sits that far below it.
double HorizonRow(double tiltDeg){
    return kCameraCentre - (90.0 - tiltDeg) / kDegPerPixel;
}

// Sky/ground boundary for data, as drawn in CreateCameraImageFromHist. Sofia's commissioning-paper
// note gives the horizon by camera-row ratios:
// - Oct. 2024 through Sept. 2025: 12 ground rows and 4 sky rows -> line at y=11.5
// - Oct. 2025 onwards:             9 ground rows and 7 sky rows -> line at y=8.5
double DataSkyGroundLineY(int night){
    const int newHorizonRatioDate = 20251001;
    const double yMin = -0.5;
    const double oldGroundSkyBoundary = yMin + 12.0;
    const double newGroundSkyBoundary = yMin + 9.0;
    return (night >= newHorizonRatioDate) ? newGroundSkyBoundary : oldGroundSkyBoundary;
}

// The CreateCameraImageFromHist Sky/Ground line at camera y = skyGroundLineY, for data and sims alike.
// A horizon outside the frame (more than ~2 deg from tilt 90) draws nothing, and a label that would
// fall off the frame is left out.
void DrawSkyGroundLine(double skyGroundLineY){
    const double frameLow = -0.5, frameHigh = 15.5;
    if (skyGroundLineY <= frameLow || skyGroundLineY >= frameHigh) return;

    TLine *l = new TLine(frameLow, skyGroundLineY, frameHigh, skyGroundLineY);
    l->SetLineWidth(6);
    l->SetLineColor(kGray);
    l->Draw();

    const double skyLabelY    = skyGroundLineY + 0.75;
    const double groundLabelY = skyGroundLineY - 0.6;

    if (skyLabelY < frameHigh){
        TText *t1 = new TText(13.1, skyLabelY, "Sky");
        t1->SetTextAlign(22);
        t1->SetTextSize(0.05);
        t1->SetTextColor(kGray);
        t1->Draw();
    }

    if (groundLabelY > frameLow){
        TText *t2 = new TText(13.1, groundLabelY, "Ground");
        t2->SetTextAlign(22);
        t2->SetTextSize(0.05);
        t2->SetTextColor(kGray);
        t2->Draw();
    }
}

// GIF delays are whole centiseconds, and 24 fps is 4.17 cs per frame. Each frame gets the delay that
// keeps the running total on the 24 fps clock (4, 4, 4, 5, 4, ...), so the gif plays at 24 fps on average.
int FrameDelayCentiseconds(int frame){
    const int frameStart = TMath::Nint(100.0 * frame / FramesPerSecond);
    const int frameEnd   = TMath::Nint(100.0 * (frame + 1) / FramesPerSecond);
    return frameEnd - frameStart;
}

// Draws one event once, then steps through the time bins, refilling the camera and the time bin label
// and adding each frame to the gif. The colour scale is fixed over the whole gif so frames compare.
// flipBothAxes: true for data (sky on top, as in the stills), false for sims (already sky-up).
// skyGroundLineY: camera y of the Sky/Ground line; nothing is drawn when it is outside the frame.
void MakeEventGif(const std::vector<std::vector<float>>& frames, int firstTimeBin, const std::string& gifPath,
                  const std::string& fileLabelText, const std::string& eventLabelText,
                  bool flipBothAxes, double skyGroundLineY){
    // brightest pixel over the whole gif sets the top of the colour scale
    float maxPE = 0.0;
    int maxPixel = -1, maxTimeBin = -1;
    for (int frame = 0; frame < (int)frames.size(); frame++){
        for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++){
            if (frames[frame][pixelID] > maxPE){
                maxPE = frames[frame][pixelID];
                maxPixel = pixelID;
                maxTimeBin = firstTimeBin + frame;
            }
        }
    }
    cout << "Brightest pixel in the gif: " << maxPixel << " at " << maxPE << " PE in time bin " << maxTimeBin << endl;
    if (maxPE < 1.0){
        cout << "WARNING: nothing reaches 1 PE in time bins " << firstTimeBin << "-" << firstTimeBin + (int)frames.size() - 1
             << "; this event has no signal in the window" << endl;
        maxPE = 1.0;
    }

    canvas->Clear();
    canvas->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.12); // room for the file name and time bin lines

    // same axes and colour scale floor as the CreateCameraImageFromHist stills
    TH2F camera("camera", ";Pixels along camera x axis;Pixels along camera y axis;PE", 16, -0.5, 15.5, 16, -0.5, 15.5);
    camera.SetDirectory(0);
    camera.SetStats(0);
    camera.GetXaxis()->SetLabelSize(0.03);
    camera.GetYaxis()->SetLabelSize(0.03);
    camera.GetXaxis()->SetTitleOffset(1.1);
    camera.GetYaxis()->SetTitleOffset(1.0);
    camera.GetZaxis()->SetTitleOffset(1.4);
    camera.SetMinimum(-1);
    camera.SetMaximum(maxPE);
    FillCamera(&camera, frames[0], flipBothAxes);
    camera.Draw("COLZ");
    canvas->Update();

    TPaletteAxis* palette = (TPaletteAxis*)camera.GetListOfFunctions()->FindObject("palette");
    if (palette) {
        palette->SetX1NDC(0.86);
        palette->SetX2NDC(0.90);
        palette->SetY1NDC(0.1);
        palette->SetY2NDC(0.88);
    }

    DrawSkyGroundLine(skyGroundLineY);
    plottools->DrawSIABBoundaries();

    TLatex fileLabel(0.5, 0.955, fileLabelText.c_str());
    fileLabel.SetNDC();
    fileLabel.SetTextAlign(22);
    fileLabel.SetTextSize(0.025);
    fileLabel.Draw();

    TLatex timeBinLabel(0.5, 0.915, "");
    timeBinLabel.SetNDC();
    timeBinLabel.SetTextAlign(22);
    timeBinLabel.SetTextSize(0.03);
    timeBinLabel.Draw();

    gSystem->Unlink(gifPath.c_str()); // a rerun would otherwise add its frames to the old gif

    for (int frame = 0; frame < (int)frames.size(); frame++){
        const int timeBin = firstTimeBin + frame;
        FillCamera(&camera, frames[frame], flipBothAxes);
        timeBinLabel.SetTitle(Form("%s   Time bin %d", eventLabelText.c_str(), timeBin));
        canvas->Modified();
        canvas->Update();

        if (PrintFrameDiagnostics){
            int brightestPixel = std::max_element(frames[frame].begin(), frames[frame].end()) - frames[frame].begin();
            float cameraSum = 0.0;
            for (float pe : frames[frame]) cameraSum += pe;
            cout << "  Time bin " << timeBin << ": brightest pixel " << brightestPixel
                 << " at " << Form("%.1f", frames[frame][brightestPixel]) << " PE, camera sum "
                 << Form("%.1f", cameraSum) << " PE" << endl;
        }

        // "gif+N" adds a frame shown for N cs; "gif++N" adds the last frame and makes the gif loop
        const bool lastFrame = (frame == (int)frames.size() - 1);
        canvas->Print(Form("%s%s%d", gifPath.c_str(), lastFrame ? "++" : "+", FrameDelayCentiseconds(frame)));
    }

    canvas->Clear(); // the camera and labels go out of scope here, so take them off the canvas first
    cout << "Wrote " << frames.size() << " frames to " << gifPath << endl;
}

// Every gif one row of the list asks for: opens the file, decides sim or data from the Sim tree,
// sets up the matching pedestal and calibration, and makes one gif per event.
void MakeGifsForRequest(const GifRequest& request){
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

    const bool isSim = FileHasSimTree(f);
    cout << "Data type: " << (isSim ? "simulation (Sim tree found)" : "demonstrator data (no Sim tree)") << endl;

    IEvent *ev = new IEvent();
    testTree->SetBranchAddress("Events", &ev);

    // data: SiPM gains per event and the night's flasher calibration
    const int night = extractYYYYMMDD(request.filePath);
    ISiPM *sipmInfo = nullptr;
    TH2F *flasher = nullptr;
    if (!isSim){
        cout << "Night: " << night << endl;
        sipmInfo = new ISiPM();
        testTree->SetBranchAddress("SiPM", &sipmInfo);
        flasher = LoadFlasherCalibration(night);
        if (!flasher){
            f->Close();
            return;
        }
    }

    // sims: the truth tree gives each event's run name and trigger flag
    TTree *simTree = nullptr;
    SEvent *simEvent = nullptr;
    std::vector<Long64_t> simEntryForTestEntry;
    if (isSim){
        simTree = (TTree*)f->Get("Sim");
        simEvent = new SEvent();
        simTree->SetBranchAddress("GrOptics", &simEvent);
        simEntryForTestEntry = MatchTestToSimEntries(testTree, simTree, simEvent);
    }

    const int firstTimeBin = (request.firstTimeBin >= 0) ? request.firstTimeBin : (isSim ? SimFirstTimeBin : DataFirstTimeBin);
    const int lastTimeBin  = (request.lastTimeBin  >= 0) ? request.lastTimeBin  : (isSim ? SimLastTimeBin  : DataLastTimeBin);
    cout << "Time bins: " << firstTimeBin << "-" << lastTimeBin << " at " << FramesPerSecond << " fps" << endl;

    std::string fileName = request.filePath.substr(request.filePath.find_last_of('/') + 1);

    std::vector<Long64_t> events = EventsToAnimate(request, testTree, simTree, simEvent, simEntryForTestEntry);
    cout << "Events to animate: " << events.size() << endl;

    for (Long64_t testEntry : events){
        cout << "---- Event " << testEntry << " ----" << endl;
        if (testEntry < 0 || testEntry >= testTree->GetEntries()){
            cerr << "Error: event " << testEntry << " is not in the Test tree (" << testTree->GetEntries() << " entries)" << endl;
            continue;
        }
        testTree->GetEntry(testEntry);

        std::vector<std::vector<Int_t>> traces(MaxNofChannels);
        for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++){
            traces[pixelID] = ev->GetSignalValue(pixelID);
        }
        if (firstTimeBin < 0 || lastTimeBin < firstTimeBin || lastTimeBin >= (int)traces[0].size()){
            cerr << "Error: time bins " << firstTimeBin << "-" << lastTimeBin << " do not fit in a "
                 << traces[0].size() << "-sample trace" << endl;
            break;
        }

        std::vector<float> pedestals(MaxNofChannels);
        std::vector<float> calibrationFactors(MaxNofChannels);
        for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++){
            if (isSim){
                pedestals[pixelID] = MeanOfSamples(traces[pixelID], SimPedestalStart, SimPedestalSamples);
                calibrationFactors[pixelID] = 1.0 / SimADCtoPE;
            } else {
                pedestals[pixelID] = MeanOfSamples(traces[pixelID], DataPedestalStart, DataPedestalSamples);
            }
        }
        if (!isSim) calibrationFactors = DataCalibrationFactors(flasher, sipmInfo->GetGain());

        // labels: the file name on top, and the run (sims) or night (data) with the event number below it.
        // Sky/Ground line: data by the night's row ratio, sims by the tilt (the SEvent's own, or the
        // Tilt_ folder for sims converted before PlotCamera stored it, which read back 0).
        std::string eventLabelText = Form("Night %d   Event %lld", night, testEntry);
        double skyGroundLineY = DataSkyGroundLineY(night);
        if (isSim){
            std::string run = "unknown run";
            double tiltDeg = TiltFromPath(request.filePath);
            std::string tiltSource = "Tilt_ folder name";
            if (testEntry < (Long64_t)simEntryForTestEntry.size()){
                simTree->GetEntry(simEntryForTestEntry[testEntry]);
                run = simEvent->GetSumRun();
                cout << "Sim truth: run " << run << ", sim event " << simEvent->GetSimEventNumber()
                     << ", triggered " << simEvent->GetTriggered() << endl;
                if (simEvent->GetTilt() != 0){
                    tiltDeg = simEvent->GetTilt();
                    tiltSource = "SEvent";
                }
            }
            eventLabelText = Form("Run %s   Event %lld", run.c_str(), testEntry);

            skyGroundLineY = HorizonRow(tiltDeg);
            if (tiltDeg == 0.0){
                cout << "WARNING: no tilt in the SEvent or the file path, no horizon drawn" << endl;
            } else {
                cout << "Horizon: tilt " << tiltDeg << " deg (from " << tiltSource << "), camera y "
                     << Form("%.2f", skyGroundLineY) << " of -0.5 to 15.5";
                if (skyGroundLineY <= -0.5)      cout << ", below the frame: whole image is sky, no line";
                else if (skyGroundLineY >= 15.5) cout << ", above the frame: whole image is ground, no line";
                cout << endl;
            }
        } else {
            cout << "Sky/Ground line at camera y " << skyGroundLineY << " (night " << night << ")" << endl;
        }

        std::vector<std::vector<float>> frames;
        for (int timeBin = firstTimeBin; timeBin <= lastTimeBin; timeBin++){
            frames.push_back(CameraAtTimeBin(traces, pedestals, calibrationFactors, timeBin));
        }

        std::string gifPath = Form("%s%s_E%lld.gif", outDir.c_str(), request.output.c_str(), testEntry);
        const bool flipBothAxes = !isSim; // sims are already sky-up
        MakeEventGif(frames, firstTimeBin, gifPath, fileName, eventLabelText, flipBothAxes, skyGroundLineY);
    }

    delete flasher;
    f->Close();
}

int main(int argc, char **argv){
    if(argc < 3){
        cout << "Too few arguments; please include the gif list file and mount flag" << endl;
        cout << "  ./MakeEventGif /DataAnalysis/MakeThingsPretty/EventGiFs/EventGifs.txt /cedar/" << endl;
        return 1;
    }

    std::string ListPath = argv[1];
    std::string mount = argv[2];

    if (mount == "y"){
        std::cout << "using mounted directory path" << std::endl;
        mnt = "/mnt/";
    } else if (mount == "n"){
        std::cout << "using default directory path" << std::endl;
        mnt = "";
    } else {
        std::cout << "using specific directory path: " << mount << std::endl;
        mnt = mount;
    }
    if (mnt != ""){
        outDir = Form("%sDataAnalysis/MakeThingsPretty/EventGiFs/Gifs/", mnt.c_str());
        CalibrationFactorDir = Form("%sDataAnalysis/flasher_calibration/Output/", mnt.c_str());
    }
    gSystem->mkdir(outDir.c_str(), true);
    cout << "Gifs go to: " << outDir << endl;

    util = new IUtilities();
    plottools = new IPlotTools();

    std::vector<std::string> listLines = util->readFileToVectorString(Form("%s%s", mnt.c_str(), ListPath.c_str()));
    if (listLines.size() < 2){
        cout << "The gif list needs a header line and at least one row: " << mnt << ListPath << endl;
        return 1;
    }
    auto rows = splitCSVLines(listLines);

    gErrorIgnoreLevel = kWarning; // ROOT prints an Info line for every gif frame otherwise

    canvas = new TCanvas("Display","Event Display",1000,900);

    for (int i = 1; i < (int)rows.size(); i++){ // row 0 is the header
        if (rows[i].size() < 3){
            cerr << "Error: row needs Output,File,Event: " << listLines[i] << endl;
            continue;
        }

        GifRequest request;
        request.output   = rows[i][0];
        request.filePath = rows[i][1];
        request.event    = rows[i][2];
        if (rows[i].size() >= 5){
            request.firstTimeBin = std::atoi(rows[i][3].c_str());
            request.lastTimeBin  = std::atoi(rows[i][4].c_str());
        }

        cout << "======== " << request.output << " ========" << endl;
        MakeGifsForRequest(request);
    }

    cout << "Completed!" << endl;
}
