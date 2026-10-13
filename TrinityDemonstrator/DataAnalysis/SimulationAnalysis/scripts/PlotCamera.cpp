#include "PlotCamera.h"


int main(int argc, char **argv){
        if(argc < 4){
                cout << "Too few arguments; usage: PlotCamera <YYYYMMDD> <simrun> <mount> [GROPT|CARE] [runFolder]" << endl;
                return 1;
        }
    
        // Get the Arguments
        std::string folString = argv[1];
        std::string simrun = argv[2];

        std::string mount = argv[3];

        // 4th argument selects which simulation stage to read. Defaults to GROPT
        // so the existing calls keep working unchanged.
        if (argc > 4) simType = argv[4];
        std::transform(simType.begin(), simType.end(), simType.begin(), ::toupper);
        if (simType != "GROPT" && simType != "CARE"){
                cout << "Unknown simulation type '" << simType << "'; expected GROPT or CARE" << endl;
                return 1;
        }
        cout << "Simulation type: " << simType << endl;

        // 5th argument restricts the pass to one run folder inside the campaign,
        // e.g. job_000004. This is the unit a condor job gets: a campaign holds
        // tens of run folders and tens of thousands of showers between them, far
        // too much for one node, but each run folder is self contained - its own
        // GROPT/CARE files and its own ISOG file - so they parallelise cleanly.
        // Left off (or "all"), every run folder is processed in one pass, which is
        // the original behaviour.
        std::string runFolder = (argc > 5) ? argv[5] : "";
        if (runFolder == "all") runFolder = "";
        if (!runFolder.empty()) cout << "Restricting to run folder: " << runFolder << endl;

        // Per shower ISOG dump is 22 lines an event, so it is opt in: set
        // ISOGVERBOSE=1 to inspect one run by hand.
        const char* isogVerboseEnv = getenv("ISOGVERBOSE");
        isogVerbose = (isogVerboseEnv && std::string(isogVerboseEnv) != "0");
        if (isogVerbose) cout << "ISOG per shower dump enabled" << endl;

        // Per pixel CARE readout dump, 256 lines an event, so it is opt in the
        // same way: set CAREVERBOSE=1 to check the pedestals and the sample the
        // camera was imaged at, pixel by pixel.
        const char* careVerboseEnv = getenv("CAREVERBOSE");
        careVerbose = (careVerboseEnv && std::string(careVerboseEnv) != "0");
        if (careVerbose) cout << "CARE per pixel readout dump enabled" << endl;

        // ROOT scales the histogram title to fit its box, so the box is widened
        // here; setting the font size alone has no effect on a long title.
        gStyle->SetTitleFontSize(kTitleFontSize);
        gStyle->SetTitleW(0.98);
        gStyle->SetTitleH(0.055);
    
        // The two stages currently live under different simdata folders:
        // GROPT under daviddata, CARE under sofiadata.
        std::string simDataDir = (simType == "CARE") ? "sofiadata" : "sofiadata";

        if (mount == "y"){ // with usingin htcondor you need to have contianers and some use full paths and other use mounts this lets you specify
                std::cout << "using mounted directory path" << std::endl;
                mnt="/mnt/";
                dataDir = Form("/mnt/simdata/%s/",simDataDir.c_str());
                outDir = "/mnt/DataAnalysis/SimulationAnalysis/data/";
        } else if (mount != "n"){
                std::cout << "using specific directory path" << std::endl;
                mnt=mount.c_str();
                dataDir = Form("%ssimdata/%s/",mnt.c_str(),simDataDir.c_str());
                outDir = Form("%sDataAnalysis/SimulationAnalysis/data/",mnt.c_str());
        }

        // GROPT, CARE and ISOG all live side by side under the same run folder.
        // simdata/sofiadata/<date>_<name>/ is a campaign folder now: it holds one
        // sub-folder per run, each with its own GROPT/CARE/ISOG/Tilt_* tree. Every
        // one of them is discovered here and processed in the same pass, so the
        // command line still names the campaign and never the individual runs. A
        // flat run folder, the original layout, resolves to a single unit of work.
        std::string RunDir = Form("%s%s_%s/",dataDir.c_str(),folString.c_str(),simrun.c_str());
        std::vector<RunUnit> runUnits = CollectRunUnits(RunDir, simType);
        if (!runFolder.empty()){
                std::vector<RunUnit> picked;
                for (size_t u = 0; u < runUnits.size(); u++)
                        if (runUnits[u].name == runFolder) picked.push_back(runUnits[u]);
                if (picked.empty()){
                        cout << "Run folder " << runFolder << " not found under " << RunDir << endl;
                        return 1;
                }
                runUnits = picked;
        }
        if (runUnits.empty()){
                cout << "No " << simType << " tilt folder under " << RunDir << ", moving on..." << endl;
                cout << "Completed!"<< endl;
                return 0;
        }
        cout << "Found " << runUnits.size() << " " << simType << " run folder(s) under " << RunDir << endl;
        // GROPT and CARE get their own folder under the night so the two sets of
        // outputs stay separate: <outDir>/<YYYYMMDD>/<GROPT|CARE>/
        mkdir(Form("%s%s",outDir.c_str(),folString.c_str()), 0775); // night folder, may already exist
        std::string outRunDir = Form("%s%s/%s/",outDir.c_str(),folString.c_str(),simType.c_str());
        if (mkdir(Form("%s%s/%s",outDir.c_str(),folString.c_str(),simType.c_str()), 0775) == -1)
        {
                cout << "Error : possibly directory already exists" << endl;

        } else{
                cout << "Directory created" <<  endl;

        }
        cout << "Output Folder: " << outRunDir << endl;

        // Every output file is keyed on this tag. With one condor job per run
        // folder the tag carries the folder name, so the 58 jobs of a campaign
        // write 58 separate sets of files into the same directory instead of
        // overwriting each other. A whole campaign pass keeps the original names.
        std::string outTag = runFolder.empty() ? simrun : simrun + "_" + runFolder;
        cout << "Output tag: " << outTag << endl;

        // Flatten the campaign into one list of files so the event loop below stays
        // a single pass and every run ends up in the same output. The ISOG file is
        // resolved once per run rather than per event. A run whose stage folder came
        // out empty, or that has no ISOG file to describe its showers, is reported
        // and dropped here instead of aborting the whole campaign - half of a night
        // is typically still running when the finished half is analysed.
        std::vector<FileJob> fileJobs;
        int nRunsSkipped = 0;
        for (size_t u = 0; u < runUnits.size(); u++){
                std::vector<std::string> rootFiles = ListFilesWithExt(runUnits[u].stagePath, ".root");
                if (rootFiles.empty()){
                        cout << runUnits[u].name << ": " << simType << " folder empty, skipping" << endl;
                        nRunsSkipped++;
                        continue;
                }
                if (runUnits[u].isogPath.empty()){
                        cout << runUnits[u].name << ": no ISOG file, skipping" << endl;
                        nRunsSkipped++;
                        continue;
                }
                // The shower information is matched to the events by row number, so
                // the file is checked before anything is read from it: the column
                // names have to be the ones this code assigns to the SEvent, and
                // there has to be one row per event. A run that fails either check
                // is dropped rather than written out with shower parameters sitting
                // on the wrong events.
                cout << runUnits[u].name << ": checking ISOG " << runUnits[u].isogPath << endl;
                if (!CheckISOGHeader(runUnits[u].isogPath)){
                        cout << runUnits[u].name << ": ISOG header does not match the expected columns, skipping" << endl;
                        nRunsSkipped++;
                        continue;
                }
                runUnits[u].isogRows = CountISOGRows(runUnits[u].isogPath);
                cout << runUnits[u].name << ": ISOG has " << runUnits[u].isogRows << " shower rows" << endl;
                for (size_t i = 0; i < rootFiles.size(); i++){
                        // A CARE folder also holds <name>_tdem.root, the converted
                        // Sim/Test output of a previous pass rather than a CARE
                        // readout. It has no Events/T0 and would only be opened to be
                        // rejected, so it is left out of the count as well.
                        const std::string tdemSuffix = "_tdem.root";
                        if (rootFiles[i].size() > tdemSuffix.size() &&
                            rootFiles[i].compare(rootFiles[i].size()-tdemSuffix.size(),
                                                 tdemSuffix.size(), tdemSuffix) == 0)
                                continue;
                        FileJob job;
                        job.unit = u;
                        job.file = rootFiles[i];
                        fileJobs.push_back(job);
                }
        }
        cout << "Skipped " << nRunsSkipped << " incomplete run folder(s)" << endl;
        if (fileJobs.empty()){
                cout << "No " << simType << " files to read under " << RunDir << ", moving on..." << endl;
                cout << "Completed!"<< endl;
                return 0;
        }
        cout << "Reading " << fileJobs.size() << " " << simType << " file(s) from "
             << runUnits.size()-nRunsSkipped << " run folder(s)" << endl;

        fileOutputPlots = new TFile(Form("%splots_%s_%s.root",outRunDir.c_str(),folString.c_str(),outTag.c_str()), "RECREATE");
        fileOutputPlots->Close();
        fileOutput = new TFile(Form("%sdatafiles_%s_%s.root",outRunDir.c_str(),folString.c_str(),outTag.c_str()), "RECREATE");

        // The output trees are made once for the whole run. Previously they were
        // recreated inside the file loop, so only the last file survived the Write.
        fileOutput->cd();
        simRun = new SEvent(1);
        treeSims = new TTree("Sim", "Events tree");
        treeSims->Branch("GrOptics","SEvent",&simRun);
        ev = new IEvent();
        treeNew = new TTree("Test", "Events tree");
        treeNew->Branch("Events","IEvent",&ev);

        canvas = new TCanvas("Display","Event Display",2*kPadPixelsW,kPadPixelsH);
        canvas->Divide(2,1);
        for (int ipad = 1; ipad <= 2; ipad++){
                canvas->cd(ipad)->SetRightMargin(kPadRight);
                canvas->cd(ipad)->SetTopMargin(kPadTop);
        }

        // The square PE camera images and the wide two pane views do not sit well
        // in one document, so each goes to its own PDF.
        //   pdf_<date>_<simrun>_PE.pdf      square 16x16 PE camera images, both stages
        //   pdf_<date>_<simrun>_2panel.pdf  wide 1mm / 6mm views, GROPT only
        // The PE canvas width is picked so the plot area comes out square once the
        // margins are removed: W*(1-kPadLeft-kPadRight) = H*(1-kPadBottom-kPadTop),
        // so with the taller header band W = H*0.72/0.75 = 960 for H = 1000.
        // GROPT gets two documents, both built on the wide two pane canvas:
        //   pdf_<date>_<simrun>_1mm6mm.pdf   1mm next to 6mm, in photons
        //   pdf_<date>_<simrun>_6mmPE.pdf    6mm next to the PE pixel map, both in PE
        // CARE has no 6mm view, so it gets the single square PE camera image:
        //   pdf_<date>_<simrun>_PE.pdf
        std::string pdf1mm6mmPath = Form("%spdf_%s_%s_1mm6mm.pdf",outRunDir.c_str(),folString.c_str(),outTag.c_str());
        std::string pdf6mmPEPath  = Form("%spdf_%s_%s_6mmPE.pdf",outRunDir.c_str(),folString.c_str(),outTag.c_str());
        std::string pdfPEPath     = Form("%spdf_%s_%s_PE.pdf",outRunDir.c_str(),folString.c_str(),outTag.c_str());

        // Square canvas for CARE's standalone PE image. Width chosen so the plot
        // area is square once the margins come out:
        // W*(1-kPadLeft-kPadRight) = H*(1-kPadBottom-kPadTop).
        TCanvas* canvasPE = 0;
        if (simType == "CARE"){
                canvasPE = new TCanvas("PE","PE Camera Display",kPadPixelsW,kPadPixelsH);
                canvasPE->SetLeftMargin(kPadLeft);
                canvasPE->SetRightMargin(kPadRight);
                canvasPE->SetBottomMargin(kPadBottom);
                canvasPE->SetTopMargin(kPadTop);
                SetPdfPageToCanvas(canvasPE);
                canvasPE->Print(Form("%s[", pdfPEPath.c_str()));
        } else {
                SetPdfPageToCanvas(canvas);
                canvas->Print(Form("%s[", pdf1mm6mmPath.c_str()));
                canvas->Print(Form("%s[", pdf6mmPEPath.c_str()));
        }

        // CARE also keeps an overlay of every FADC trace. The branch wrote one
        // traces_*.root per input file; here they accumulate into a single run
        // level histogram so a 4000 file night stays manageable.
        TH2F* hTraces = 0;
        TH2F* hTracesCorr = 0;
        if (simType == "CARE"){
                hTraces = new TH2F("hTraces", "All FADC Traces;Timebin;ADC Counts",
                                   200, -0.5, 199.5, 150, -0.5, 1499.5);
                hTraces->SetStats(0);
                hTraces->SetDirectory(0);

                // the same traces after the full correction that create_root_file_care
                // writes into the IEvent: inverted and shifted onto a pedestal+3400
                // baseline, so pulses run downwards the way the demonstrator records
                // them. The range covers the whole band that correction produces.
                hTracesCorr = new TH2F("hTracesCorrected",
                                       "All FADC Traces AFTER correction;Timebin;ADC Counts (inverted, pedestal+3400 baseline)",
                                       200, -0.5, 199.5, 450, -0.5, 4499.5);
                hTracesCorr->SetStats(0);
                hTracesCorr->SetDirectory(0);
        }

        // GROPT stores ray traced photons in T1, CARE stores FADC traces in Events/T0
        std::string treeName = (simType == "CARE") ? "Events/T0" : "T1";

        for(int f = 0; f<static_cast<int>(fileJobs.size()); f++){
                const RunUnit &unit  = runUnits[fileJobs[f].unit];
                std::string FilePath = unit.stagePath + fileJobs[f].file;
                cout << "-------------------------------------------------------" << endl;
                cout << "file " << f+1 << "/" << fileJobs.size() << ", run " << unit.name
                     << " (" << unit.tilt << ")" << endl;
                cout << "file path: " <<FilePath << endl;
                cout << "ISOG file path: " << unit.isogPath << endl;

                if (!util->isBranchPresentInFile(FilePath, treeName)) {
                        continue; // Skip to the next branch if not present

                }

                TFile* file = TFile::Open(FilePath.c_str(), "READ");
                TTree* tree = (TTree*)file->Get(treeName.c_str());
                if (!tree){
                        cout << "No tree " << treeName << " in " << FilePath << endl;
                        file->Close();
                        delete file;
                        continue;
                }

                Int_t iPhotonsInFocalPlane = 0;
                if (simType == "CARE"){
                        // One branch per pixel, each is a vector<int> of 200 timebins
                        for (int pix = 0; pix < kNPixels; pix++) {
                                tree->SetBranchAddress(Form("vFADCTraces%d", pix), &vFADCTraces[pix]);
                        }
                        // CARE's own photon count, the counterpart of the GROPT photon integral
                        if (tree->GetBranch("iPhotonsInFocalPlane"))
                                tree->SetBranchAddress("iPhotonsInFocalPlane", &iPhotonsInFocalPlane);
                } else {
                        tree->SetBranchAddress("photonX", &vPhotonX);
                        tree->SetBranchAddress("photonY", &vPhotonY);
                        tree->SetBranchAddress("photonDcosX", &vPhotonDcosX);
                        tree->SetBranchAddress("photonDcosY", &vPhotonDcosY);
                }

                // titles carry the file name, so trim the extension to leave the
                // title box room to render at a readable size
                std::string fileLabel = fileJobs[f].file;
                size_t dotPos = fileLabel.rfind(".root");
                if (dotPos != std::string::npos) fileLabel = fileLabel.substr(0,dotPos);

                int nEntries = tree->GetEntries();
                if (nEntries==0) {
                        cout << "No events" << endl;
                        file->Close();
                        delete file;
                        continue;
                }
                cout << "Total number of entries: " << nEntries << endl;

                // The ISOG row and the tree entry are the same shower only if the two
                // have the same length. On 20260916 they do (715 rows, 715 entries, in
                // every one of the 58 run folders), but a truncated ISOG would otherwise
                // pair every event with the wrong shower silently, so it is checked here
                // rather than assumed.
                if (unit.isogRows >= 0 && unit.isogRows != nEntries){
                        cout << "ISOG/tree length mismatch in " << unit.name << ": "
                             << unit.isogRows << " ISOG rows against " << nEntries
                             << " tree entries, skipping this file" << endl;
                        file->Close();
                        delete file;
                        continue;
                }

                for(int EventCounter = 0; EventCounter < (nEntries); EventCounter++){
                        cout << "Event number: " << EventCounter << endl;
                        tree->GetEntry(EventCounter);

                        // The PE camera image. Same binning, same axes and same object
                        // name for both pipelines, only the contents are built differently.
                        TH2F* hPhotonDist1xPixels = new TH2F("hPhotonDist1xPixels", Form("Camera View %s, %s, %s;Pixels along camera x axis;Pixels along camera y axis;PE", folString.c_str(),unit.name.c_str(),fileLabel.c_str()), CameraStep, CameraMin, CameraMax, CameraStep, CameraMin, CameraMax);

                        double totalPhotons = 0;
                        double totalPE = 0;   // SIZE: summed PE over the camera
                        bool triggered = false;
                        // Sample the CARE camera image is taken at. Found per
                        // event below; kept out here because the plot labels
                        // further down have to say which sample it was.
                        int careTimeBin = kCareTimeBin;

                        if (simType == "GROPT"){
                                hPhotonX = new TH1F("hist photonX","PhotonX Distribution", 192, -96, 96);
                                hPhotonX->GetXaxis()->SetTitle("Photon X Position (mm)");
                                hPhotonX->GetYaxis()->SetTitle("Number of Photons");
                                hPhotonX->SetStats(0);

                                hPhotonY = new TH1F("hist photonY","PhotonY Distribution", 96, -48, 48);
                                hPhotonY->GetXaxis()->SetTitle("Photon Y Position (mm)");
                                hPhotonY->GetYaxis()->SetTitle("Number of Photons");
                                hPhotonY->SetStats(0);

                                hPhotonDcosX = new TH1F("hist photonDcosX","PhotonDcosX Distribution", 100, -1, 1);
                                hPhotonDcosX->GetXaxis()->SetTitle("PhotonDcosX");
                                hPhotonDcosX->GetYaxis()->SetTitle("Number of Photons");
                                hPhotonDcosX->SetStats(0);

                                hPhotonDcosY = new TH1F("hist photonDcosY","PhotonDcosY Distribution", 100, -1, 1);
                                hPhotonDcosY->GetXaxis()->SetTitle("PhotonDcosY");
                                hPhotonDcosY->GetYaxis()->SetTitle("Number of Photons");
                                hPhotonDcosY->SetStats(0);

                                hPhotonZ = new TH1F("hist photonZ","Photon Incidence Angle Distribution", 180, 0, 90);
                                hPhotonZ->GetXaxis()->SetTitle("Photon Incidence Angle (Deg)");
                                hPhotonZ->GetYaxis()->SetTitle("Number of Photons");
                                hPhotonZ->SetStats(0);

                                hPhotonDist1x = new TH2F("hist1x", Form("Camera View %s, %s, %s, 1mm bins", folString.c_str(),unit.name.c_str(),fileLabel.c_str()), 96, -48.5, 47.5, 96, -48.5, 47.5);
                                hPhotonDist1x->GetXaxis()->SetTitle("X (mm)");
                                hPhotonDist1x->GetYaxis()->SetTitle("Y (mm)");
                                hPhotonDist1x->GetZaxis()->SetTitle("Photons");
                                hPhotonDist1x->GetZaxis()->SetTitleOffset(kZTitleOffset);

                                hPhotonDist1x6mm = new TH2F("hist1x6mm", Form("Camera View %s, %s, %s, 6mm bins", folString.c_str(),unit.name.c_str(),fileLabel.c_str()), 16, -48.5, 47.5, 16, -48.5, 47.5);
                                hPhotonDist1x6mm->GetXaxis()->SetTitle("X (mm)");
                                hPhotonDist1x6mm->GetYaxis()->SetTitle("Y (mm)");
                                hPhotonDist1x6mm->GetZaxis()->SetTitle("Photons");
                                hPhotonDist1x6mm->GetZaxis()->SetTitleOffset(kZTitleOffset);

                                hPhotonDist4x = new TH2F("hist4x", Form("Camera View %s, %s, %s (4x)", folString.c_str(),unit.name.c_str(),fileLabel.c_str()), 32*nResMult, -96.5, 95.5, 16*nResMult, -48.5, 47.5);
                                hPhotonDist4x->GetXaxis()->SetTitle("X (mm)");
                                hPhotonDist4x->GetYaxis()->SetTitle("Y (mm)");
                                hPhotonDist4x->GetZaxis()->SetTitle("Photons");
                                hPhotonDist4x->GetZaxis()->SetTitleOffset(kZTitleOffset);
                                hPhotonDist4x->SetStats(0);

                                for (size_t j = 0; j < vPhotonX->size(); j++)
                                {
                                        hPhotonX->Fill((*vPhotonX)[j]);
                                        hPhotonY->Fill((*vPhotonY)[j]);
                                        hPhotonDcosX->Fill((*vPhotonDcosX)[j]);
                                        hPhotonDcosY->Fill((*vPhotonDcosY)[j]);
                                        hPhotonZ->Fill(TMath::ASin(sqrt(pow((*vPhotonDcosX)[j],2)+pow((*vPhotonDcosY)[j],2)))*RadtoDeg);
                                        hPhotonDist1x->Fill((*vPhotonX)[j], (*vPhotonY)[j]);
                                        hPhotonDist1x6mm->Fill((*vPhotonX)[j], (*vPhotonY)[j]);
                                        hPhotonDist4x->Fill((*vPhotonX)[j], (*vPhotonY)[j]);

                                }

                                totalPhotons = hPhotonDist1x6mm->Integral();
                                // SIZE, summed the same way the PE map is built below:
                                // photons scaled to PE and clipped at the 160 ceiling
                                for (int ix = 1; ix <= hPhotonDist1x6mm->GetNbinsX(); ++ix)
                                        for (int iy = 1; iy <= hPhotonDist1x6mm->GetNbinsY(); ++iy){
                                                double pe = hPhotonDist1x6mm->GetBinContent(ix,iy)*PtoPEconversion;
                                                totalPE += (pe > 160) ? 160 : pe;
                                        }
                                triggered = (hPhotonDist1x6mm->GetMaximum()*PtoPEconversion >= kTriggerPE);

                        } else { // CARE
                                // Ported from the Sofia_CARE_Simanalysis branch: build the
                                // camera image from one FADC sample of every pixel, after
                                // subtracting a pedestal taken from the leading samples.
                                // The pedestals are taken first, for the whole camera, so the
                                // sample to image at can be found before any pixel is filled:
                                // it is where the brightest - the triggering - pixel peaks,
                                // rather than a sample fixed in advance that the shower only
                                // sometimes happens to fall on. See FindCareSnapshot.
                                int pedestals[kNPixels] = {};
                                for (int pix = 0; pix < kNPixels; pix++)
                                        if (vFADCTraces[pix]) pedestals[pix] = GetCarePedestal(*vFADCTraces[pix]);

                                CareSnapshot snapshot = FindCareSnapshot(vFADCTraces, pedestals);
                                careTimeBin = snapshot.bin;
                                if (snapshot.pixel < 0)
                                        cout << "No pixel rises above its pedestal, imaging at the fallback timebin "
                                             << careTimeBin << endl;
                                else
                                        cout << "Brightest pixel " << snapshot.pixel << " peaks at timebin "
                                             << careTimeBin << " with " << snapshot.amplitude << " ADC ("
                                             << snapshot.amplitude/kCareADCtoPE << " PE), imaging the camera there" << endl;

                                for (int pix = 0; pix < kNPixels; pix++) {

                                        if (!vFADCTraces[pix]) continue;

                                        int nSamples = (int)vFADCTraces[pix]->size();
                                        int pedestal = pedestals[pix];
                                        if (careVerbose)
                                                cout << "Pixel " << pix << ": pedestal = " << pedestal
                                                     << " (raw), demonstrator baseline = " << pedestal+3400
                                                     << ", sample " << careTimeBin << " = "
                                                     << ((nSamples > careTimeBin) ? (*vFADCTraces[pix])[careTimeBin] : 0)
                                                     << endl;
                                        float val = (nSamples > careTimeBin)
                                                        ? (*vFADCTraces[pix])[careTimeBin] - pedestal
                                                        : 0.0f;
                                        int nx, ny;
                                        plottools->FindBin(pix, &nx, &ny);
                                        if (val <= 0.0){
                                                val = 0.001;
                                        }
                                        hPhotonDist1xPixels->SetBinContent(nx + 1, ny + 1, val/kCareADCtoPE);

                                        // both trace overlays are filled here, where the
                                        // pedestal for this pixel is already known, and
                                        // before any trigger decision so they stay raw
                                        // readout diagnostics
                                        for (int t = 0; t < nSamples; t++){
                                                hTraces->Fill(t, (*vFADCTraces[pix])[t]);
                                                // exactly what create_root_file_care stores
                                                hTracesCorr->Fill(t, (-1*(*vFADCTraces[pix])[t])+2*pedestal+3400);
                                        }
                                }

                                totalPhotons = iPhotonsInFocalPlane;
                                totalPE = hPhotonDist1xPixels->Integral();
                                triggered = (hPhotonDist1xPixels->GetMaximum() >= kTriggerPE);
                        }

                        // First create the root file for this event
                        // the individual run, not the campaign, so an event can be
                        // traced back to the folder it came out of
                        simRun->SetSimRun(unit.name);
                        simRun->SetDate(std::stoi(folString)); // placeholder for now
                        simRun->SetSimEventNumber(EventCounter);

                        // Flag whether the event has enough signal to trigger. Untriggered
                        // events are still saved to the output, just marked Triggered = false.
                        cout << "Number of photons detected: " << totalPhotons << endl;

                        simRun->SetTriggered(triggered);
                        // tilt angle read straight off the folder name, e.g.
                        // Tilt_90.000000 -> 90.0, so the number of decimals does
                        // not have to be known in advance
                        simRun->SetTilt(atof(unit.tilt.c_str() + kTiltPrefixLen));
                        if (simType == "GROPT"){
                                simRun->SetPosX(*vPhotonX);
                                simRun->SetPosY(*vPhotonY);
                        } else {
                                // CARE has no per photon positions; clear them so the event
                                // does not inherit the previous fill.
                                simRun->SetPosX(std::vector<Float_t>());
                                simRun->SetPosY(std::vector<Float_t>());
                        }

                        // Same ISOG row for both pipelines, keyed on the entry counter
                        LoadISOG(simRun, unit.isogPath, EventCounter);

                        treeSims->Fill();
                        if (!triggered) {
                                cout << "Not enough photons detected in this simulation run, saving as untriggered..." << endl;
                                // Clean up this iteration's histograms before skipping so they
                                // are not left registered in the directory for the next event
                                delete hPhotonDist1xPixels;
                                if (simType == "GROPT"){
                                        delete hPhotonX;
                                        delete hPhotonY;
                                        delete hPhotonDcosX;
                                        delete hPhotonDcosY;
                                        delete hPhotonZ;
                                        delete hPhotonDist1x;
                                        delete hPhotonDist1x6mm;
                                        delete hPhotonDist4x;
                                }
                                continue;
                        }

                        // Emergence angle (azimuth wrapped to [-180,180] from the telescope axis)
                        // and distance to the emergence point (magnitude of the telescope position vector)
                        double emergenceAngle = util->GetEmergenceAngle(simRun->GetAzimuthAngle());
                        double emergenceDistance = util->GetEmergenceDistance(simRun->GetTelescope_Xpos(), simRun->GetTelescope_Ypos(), simRun->GetTelescope_Zpos());
                        cout << "Emergence Angle (from telescope axis): " << emergenceAngle << " deg" << endl;
                        cout << "Distance to Emergence Point: " << emergenceDistance << " (units of telescope coords)" << endl;

                        // Second create the plots for each event and save them to a root file

                        // create a lable with neutrino energy

                        fileOutputPlots = new TFile(Form("%splots_%s_%s.root",outRunDir.c_str(),folString.c_str(),outTag.c_str()), "UPDATE");
                        // two lines: the single line ran off the side of the page
                        TPaveText *ptEnergy = new TPaveText(0.05,kLabelY1,0.95,kLabelY2,"NDC");
                        ptEnergy->SetFillColor(0);
                        ptEnergy->SetTextAlign(12);
                        ptEnergy->SetTextSize(kLabelTextSize);
                        ptEnergy->AddText(Form("Event: %i, Neutrino Energy: %.3g GeV, Emergence Angle: %.3g degrees",
                                               EventCounter, simRun->GetNeutrinoEnergy(), emergenceAngle));
                        ptEnergy->AddText(Form("Shower Distance: %.3g m, Total Photons: %i, SIZE = %.4g PE",
                                               emergenceDistance, int(totalPhotons), totalPE));

                        if (simType == "GROPT"){
                                canvas->Clear();
                                canvas->Divide(2,1);
                                for (int ipad = 1; ipad <= 2; ipad++){
                                        canvas->cd(ipad)->SetRightMargin(kPadRight);
                                        canvas->cd(ipad)->SetTopMargin(kPadTop);
                                }
                                canvas->cd(1);
                                hPhotonDist1x->SetMarkerSize(200.0);
                                hPhotonDist1x->SetStats(0);
                                hPhotonDist1x->Draw("colz");
                                hPhotonDist1x->Write(Form("h_%s_%s_%i_1mm",folString.c_str(),unit.name.c_str(),EventCounter));
                                DrawMUSICBoundaries();
                                ptEnergy->Draw("SAME");
                                canvas->cd(2);
                                hPhotonDist1x6mm->SetMarkerSize(200.0);
                                hPhotonDist1x6mm->SetStats(0);
                                hPhotonDist1x6mm->Draw("colz");
                                hPhotonDist1x6mm->Write(Form("h_%s_%s_%i_6mm",folString.c_str(),unit.name.c_str(),EventCounter));
                                DrawMUSICBoundaries();

                                canvas->Write(Form("C_%s_%s_%i",folString.c_str(),unit.name.c_str(),EventCounter));
                                canvas->Print(pdf1mm6mmPath.c_str()); // 1mm / 6mm page

                                // Scale to PE into the pixel map only. The 6mm histogram is
                                // left in photons so it still reads as 6mm bins when it is
                                // drawn beside the PE map below.
                                // Loop through all bins (excluding underflow/overflow)
                                int nBinsX = hPhotonDist1x6mm->GetNbinsX();
                                int nBinsY = hPhotonDist1x6mm->GetNbinsY();
                                for (int ix = 1; ix <= nBinsX; ++ix) {
                                        for (int iy = 1; iy <= nBinsY; ++iy) {
                                                double content = hPhotonDist1x6mm->GetBinContent(ix, iy)*PtoPEconversion;
                                                if (content > 160){
                                                        content = 160;
                                                }
                                                hPhotonDist1xPixels->SetBinContent(ix,iy,content);
                                        }
                                }
                                // ---- second page: the PE scaled 6mm view beside the PE pixel map,
                                // so the same event can be read off in the same units side by side
                                canvas->Clear();
                                canvas->Divide(2,1);
                                for (int ipad = 1; ipad <= 2; ipad++){
                                        canvas->cd(ipad)->SetRightMargin(kPadRight);
                                        canvas->cd(ipad)->SetTopMargin(kPadTop);
                                }
                                canvas->cd(1);
                                hPhotonDist1x6mm->SetStats(0);
                                hPhotonDist1x6mm->GetXaxis()->SetTitle("X (mm)");
                                hPhotonDist1x6mm->GetYaxis()->SetTitle("Y (mm)");
                                hPhotonDist1x6mm->GetZaxis()->SetTitle("Photons");
                                hPhotonDist1x6mm->GetZaxis()->SetTitleOffset(kZTitleOffset);
                                hPhotonDist1x6mm->Draw("colz");
                                DrawMUSICBoundaries();
                                ptEnergy->Draw("SAME");

                                canvas->cd(2);
                                hPhotonDist1xPixels->SetMarkerSize(200.0);
                                hPhotonDist1xPixels->SetStats(0);
                                hPhotonDist1xPixels->Draw("colz");
                                hPhotonDist1xPixels->Write(Form("h_%s_%s_%iPixelsPE",folString.c_str(),unit.name.c_str(),EventCounter));
                                hPhotonDist1xPixels->GetXaxis()->SetTitle("Pixels along camera x axis");
                                hPhotonDist1xPixels->GetYaxis()->SetTitle("Pixels along camera y axis");
                                hPhotonDist1xPixels->GetZaxis()->SetTitleOffset(kZTitleOffset);
                                hPhotonDist1xPixels->GetZaxis()->SetTitle("PE");
                                plottools->DrawSIABBoundaries();

                                canvas->Write(Form("C_%s_%s_%iPixelsPE",folString.c_str(),unit.name.c_str(),EventCounter));
                                canvas->Print(pdf6mmPEPath.c_str()); // 6mm / PE page

                        } else {
                                // CARE: square canvas, so the 16x16 camera map is not stretched
                                canvasPE->cd();
                                canvasPE->Clear();
                                canvasPE->SetLeftMargin(kPadLeft);
                                canvasPE->SetRightMargin(kPadRight);
                                canvasPE->SetBottomMargin(kPadBottom);
                                canvasPE->SetTopMargin(kPadTop);
                                hPhotonDist1xPixels->SetMarkerSize(200.0);
                                hPhotonDist1xPixels->SetStats(0);
                                hPhotonDist1xPixels->Draw("colz");
                                hPhotonDist1xPixels->Write(Form("h_%s_%s_%iPixelsPE",folString.c_str(),unit.name.c_str(),EventCounter));
                                hPhotonDist1xPixels->GetXaxis()->SetTitle("Pixels along camera x axis");
                                hPhotonDist1xPixels->GetYaxis()->SetTitle("Pixels along camera y axis");
                                hPhotonDist1xPixels->GetZaxis()->SetTitleOffset(kZTitleOffset);
                                hPhotonDist1xPixels->GetZaxis()->SetTitle(Form("PE at timebin %i", careTimeBin));
                                hPhotonDist1xPixels->SetMinimum(0.0);
                                plottools->DrawSIABBoundaries();
                                ptEnergy->Draw("SAME");
                                canvasPE->Write(Form("C_%s_%s_%iPixelsPE",folString.c_str(),unit.name.c_str(),EventCounter));
                                canvasPE->Print(pdfPEPath.c_str()); // PE camera page
                        }
                        gROOT->cd(); // restore root directory before touching the output files

                        // close the plots file for now
                        fileOutputPlots->Close();

                        if (simType == "CARE")
                                create_root_file_care(hPhotonDist1xPixels,treeNew,vFADCTraces);
                        else
                                create_root_file(hPhotonDist1xPixels,treeNew);

                        delete ptEnergy;
                        delete hPhotonDist1xPixels;
                        if (simType == "GROPT"){
                                delete hPhotonX;
                                delete hPhotonY;
                                delete hPhotonDcosX;
                                delete hPhotonDcosY;
                                delete hPhotonZ;
                                delete hPhotonDist1x;
                                delete hPhotonDist1x6mm;
                                delete hPhotonDist4x;
                        }

                }

                file->Close();
                delete file;

        }
        // close the PDFs
        if (simType == "CARE"){
                canvasPE->Print(Form("%s]", pdfPEPath.c_str()));
                util->setFilePermissions(pdfPEPath);
                cout << "PE PDF written: " << pdfPEPath << endl;
        } else {
                canvas->Print(Form("%s]", pdf1mm6mmPath.c_str()));
                util->setFilePermissions(pdf1mm6mmPath);
                cout << "1mm/6mm PDF written: " << pdf1mm6mmPath << endl;
                canvas->Print(Form("%s]", pdf6mmPEPath.c_str()));
                util->setFilePermissions(pdf6mmPEPath);
                cout << "6mm/PE PDF written: " << pdf6mmPEPath << endl;
        }

        if (simType == "CARE"){
                // The trace overlay gets its own file, matching the traces_*.root
                // convention, so it is not buried among the event displays.
                std::string tracesPath = Form("%straces_%s_%s.root",outRunDir.c_str(),folString.c_str(),outTag.c_str());
                TFile* fileTraces = new TFile(tracesPath.c_str(), "RECREATE");
                TCanvas* canvasTraces = new TCanvas("Traces","FADC Traces",1600,800);
                canvasTraces->SetRightMargin(0.15);
                hTraces->Draw("colz");
                hTraces->Write(Form("hTraces_%s_%s",folString.c_str(),outTag.c_str()));
                canvasTraces->Write(Form("CTraces_%s_%s",folString.c_str(),outTag.c_str()));

                // the same overlay after the pedestal correction
                TCanvas* canvasTracesCorr = new TCanvas("TracesCorrected","FADC Traces AFTER correction",1600,800);
                canvasTracesCorr->SetRightMargin(0.15);
                hTracesCorr->Draw("colz");
                hTracesCorr->Write(Form("hTracesCorrected_%s_%s",folString.c_str(),outTag.c_str()));
                canvasTracesCorr->Write(Form("CTracesCorrected_%s_%s",folString.c_str(),outTag.c_str()));
                fileTraces->Close();
                util->setFilePermissions(tracesPath);
                cout << "Trace overlay written: " << tracesPath << endl;
                delete canvasTraces;
                delete canvasTracesCorr;
                delete fileTraces;
                delete hTraces;
                delete hTracesCorr;
        }
        if (canvasPE) delete canvasPE;

        fileOutput->cd();
        treeSims->Write();
        treeNew->Write();
        fileOutput->Close();
        util->setFilePermissions(Form("%sdatafiles_%s_%s.root",outRunDir.c_str(),folString.c_str(),outTag.c_str()));
        util->setFilePermissions(Form("%splots_%s_%s.root",outRunDir.c_str(),folString.c_str(),outTag.c_str()));

        delete simRun;
        delete ev;
        delete fileOutput;
        delete fileOutputPlots;
        delete canvas;
        cout << "Completed!" << endl;   
}


void create_root_file(TH2F* hist, TTree* treeNew) {


    ev->SetCoBoTime(ULong64_t( 174000000000000ULL));
    ev->SetUNIXTime(ULong64_t(0ULL));
    ev->SetTBTime(ULong64_t( 174000000000000ULL));
    ev->SetEventType(ULong64_t(3));
    std::vector<vector<Int_t>> pixels;
    for(int i =0; i < 256; i++){
        std::vector<Int_t> trace;    
        for(int j = 0; j<512; j++) {
            trace.push_back(160.0);
        }
        // get the pixel nx, and ny

        int nx,ny;
        plottools->FindBin(i,&nx,&ny);
        // cout << "pixel:" << i << "bin content: " << hist->GetBinContent(nx+1,ny+1) << endl;
        // trace[239] = 160.0 - hist->GetBinContent(nx+1,ny+1);
        if (hist->GetBinContent(nx+1,ny+1) > 160){
                trace[239] = 160.0 - 160.0;

        } else {
                trace[239] = 160.0 - hist->GetBinContent(nx+1,ny+1);

        }
        // cout << trace[239] << endl;
        pixels.push_back(trace);
    }
    ev->SetSignalValue(pixels);
    int maxBinX, maxBinY, dummy;
    hist->GetMaximumBin(maxBinX, maxBinY, dummy);
    hist->GetBinContent(maxBinX, maxBinY);
//     cout << "Max bin content: " << hist->GetBinContent(maxBinX, maxBinY) << endl;
    int MusicID = plottools->FindMUSIC(maxBinX-1, maxBinY-1);

    ev->SetROIMusicID({MusicID,-1});
    ev->SetROIPixelID({1});
    std::vector<float> temp = {1.0};
    std::string temp_str = "";
    ev->SetParametersFromTimestamp(temp_str,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp);
    
    
    treeNew->Fill();

}


// CARE writes the real FADC waveforms into the IEvent, ported from the
// Sofia_CARE_Simanalysis branch. Each 200 sample trace is inverted and shifted
// so the pulse reads the way the demonstrator records it, then centred inside
// the 512 sample window with the rest held at the pixel baseline.
void create_root_file_care(TH2F* hist, TTree* treeNew, std::vector<int>* traces[]) {

    ev->SetCoBoTime(ULong64_t( 174000000000000ULL));
    ev->SetUNIXTime(ULong64_t(0ULL));
    ev->SetTBTime(ULong64_t( 174000000000000ULL));
    ev->SetEventType(ULong64_t(3));

    const int traceLen = 512;
    const int fadcLen  = 200;
    const int offset   = (traceLen - fadcLen) / 2; // = 156, centers the 200-bin trace in 512

    std::vector<vector<Int_t>> pixels;
    for (int i = 0; i < 256; i++) {

        // pedestal baseline for this pixel's trace, same as the camera image uses
        int nSamples = traces[i] ? (int)traces[i]->size() : 0;
        int nPed = std::min(kCarePedBins, nSamples);
        int pedestal = 0;
        for (int t = 0; t < nPed; t++){
                pedestal += (*traces[i])[t];
        }
        if (nPed > 0) pedestal /= nPed; // average of the first bins

        std::vector<Int_t> trace(traceLen, pedestal+3400); // fill with pedestal baseline

        // Embed the 200-bin FADC trace centered in the 512-bin trace
        for (int t = 0; t < fadcLen && t < nSamples; t++) {
            trace[offset + t] = (-1*(*traces[i])[t])+2*pedestal+3400; // subtract pedestal to center around zero and then shift to match the demonstrator
        }
        pixels.push_back(trace);
    }
    ev->SetSignalValue(pixels);

    int maxBinX, maxBinY, dummy;
    hist->GetMaximumBin(maxBinX, maxBinY, dummy);
    int MusicID = plottools->FindMUSIC(maxBinX-1, maxBinY-1);

    ev->SetROIMusicID({MusicID,-1});
    ev->SetROIPixelID({1});
    std::vector<float> temp = {1.0};
    std::string temp_str = "";
    ev->SetParametersFromTimestamp(temp_str,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp);

    treeNew->Fill();

}


// ROOT always writes an A4 MediaBox for PDF, but gStyle's paper size still sets
// how large the canvas is drawn on it. Left at the default the wide two pane
// canvas is letterboxed into a strip; set larger than the sheet and it runs off
// the edge. So the drawing area is fitted to the canvas aspect inside the usable
// A4 area, which ROOT orients by that same aspect (wide -> landscape).
void SetPdfPageToCanvas(TCanvas *c) {
    double aspect = (double)c->GetWw() / (double)c->GetWh();
    // usable drawing area, mirroring ROOT's own 20x26 default for portrait
    double maxW = (aspect >= 1.0) ? kPdfAreaLongCm  : kPdfAreaShortCm;
    double maxH = (aspect >= 1.0) ? kPdfAreaShortCm : kPdfAreaLongCm;
    double w = std::min(maxW, maxH * aspect);
    gStyle->SetPaperSize(w, w / aspect);
}
