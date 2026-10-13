#include "EventCleaning.h"

// Event cleaning code see Latex for more information 
//https://www.overleaf.com/project/67bf40074acfe901c4446a39
//and the pdf in this folder


int main(int argc, char **argv){
    if(argc < 1){
        cout << "Too few arguments; please include the date (YYYYMMDD) data directory to summarize" << endl;
		return 1;
	}
    
    std::string mount = argv[2];
    std::string filename_argument = "";
    filename_argument = argv[3];

    // Optional trailing arguments, forwarded by run_script.sh: "verbose" for every event,
    // "event <N>" for one event only. Without them the run prints per-file lines and the totals
    // and nothing per event. Everything is processed and written either way.
    ParseDiagnosticArguments(argc, argv);

    // RISE_FRACTION=<percent> picks the rise fraction and TIME_WINDOW=<n> the time-walk window for this
    // run. Read before the output names are built, since they carry both.
    if (!ReadRunSettingsFromEnvironment()) return 1;
    PrintTimeWalkSettings();
    
    if (mount == "y"){ // with usingin htcondor you need to have contianers and some use full paths and other use mounts this lets you specify
        std::cout << "using mounted directory path" << std::endl;
        mnt="/mnt/";
        dataDir = "/mnt/DataAnalysis/MergedData/Output/";
        neighborDir = "/mnt/DataAnalysis/event_cleaning/TimeWalkCleaning/neighbors/";
        CalibrationFactorDir = "/mnt/DataAnalysis/flasher_calibration/Output/";
        outDir = "/mnt/DataAnalysis/event_cleaning/Output/";
    } else if (mount != "n"){
        std::cout << "using specific directory path" << std::endl;
        mnt=mount.c_str();
        dataDir = Form("%sDataAnalysis/MergedData/Output/",mnt.c_str());
        neighborDir = Form("%sDataAnalysis/event_cleaning/TimeWalkCleaning/neighbors/",mnt.c_str());
        CalibrationFactorDir = Form("%sDataAnalysis/flasher_calibration/Output/",mnt.c_str());
        outDir = Form("%sDataAnalysis/event_cleaning/Output/",mnt.c_str());
    }

    

    // Get the Arguments 
    std::string folString = argv[1];

    // Load in all the files
    std::string FolderPath = Form("%s%s/",dataDir.c_str(),folString.c_str());
    cout << "FolderPath: " <<FolderPath << endl;
    // cout << "folString: " << folString.c_str() << endl;
    // Createthe fileNames Vector based on the type of data found in the passes argument
    // Sim - Simulation data (Only Test Branch events)
    // muon - muon data (Forced Branch events)
    // bkg - Background data (Test Branch and HLED Branch events in seperate folder)
    // none - normal data
    std::vector<std::string> fileNamesVec;
    if (folString.find("bkg")==0) {
        whatData = "bkg";
        fileNamesVec=util->readFileToVectorString(Form("%s%s.txt",bkgDir.c_str(),folString.c_str()));
        
    } else if (folString.find("muon")==0) {
        whatData = "muon";
        fileNamesVec=util->readFileToVectorString(Form("%s%s.txt",muonDir.c_str(),folString.c_str()));

    } else if (folString.find("sim")==0) {
        cout << "Simulation data detected" << endl;
        whatData = "sim";
        folString = folString.substr(3, 19); // Extract date from filename
        // cout << "Simulation Information: " << folString << endl;
        // Old flat layout: every CARE file sat directly in <date>/CARE/, next to the
        // plots_/traces_ byproducts that had to be filtered back out.
        // FolderPath = Form("%s%s/CARE/",simDir.c_str(),folString.c_str());
        // cout << "Simulation Folder Path: " << FolderPath << endl;
        // fileNamesVec=util->GetFilesInDirectory(FolderPath,".root");
        // // remove files in the vec that have plots_
        // fileNamesVec.erase(
        //     std::remove_if(
        //         fileNamesVec.begin(),
        //         fileNamesVec.end(),
        //         [](const std::string& name) {
        //             return name.find("plots_") != std::string::npos;
        //         }
        //     ),
        //     fileNamesVec.end()
        // );
        // fileNamesVec.erase(
        //     std::remove_if(
        //         fileNamesVec.begin(),
        //         fileNamesVec.end(),
        //         [](const std::string& name) {
        //             return name.find("traces_") != std::string::npos;
        //         }
        //     ),
        //     fileNamesVec.end()
        // );

        // sofia layout: <date>/<run>/CARE/Tilt_<angle>/<shower>_tdem.root, one run per folder.
        // GetSimTdemFiles walks all the runs and returns the _tdem.root files relative to <date>/.
        FolderPath = Form("%s%s/",simDir.c_str(),folString.c_str());
        cout << "Simulation Folder Path: " << FolderPath << endl;
        fileNamesVec = GetSimTdemFiles(FolderPath);
        cout << "Simulation Files Found: " << fileNamesVec.size() << endl;


    } else {
        if (filename_argument != "n"){ // if the file name not specified then do all files in the directory
        // std::cout << "using specific file name" << std::endl;
        // std::string specificfile = Form("%s%s/%s",dataDir.c_str(),folString.c_str(),filename_argument.c_str());
        // // cout << "Specific File: " << specificfile << endl;
        fileNamesVec.push_back(filename_argument);
    } else {
        // fileNamesVec = read_directory(dirName.c_str());
        fileNamesVec=util->GetFilesInDirectory(FolderPath,".root");
        fileNamesVec.erase(fileNamesVec.begin(), fileNamesVec.begin() + 2);
    }
    }
    
    CreateFileName(folString,whatData);
    
    // Creates a root file that stores all the data    
    file = new TFile(OutputFileRoot.c_str(), "RECREATE");  // "RECREATE" to overwrite if it exists
    file->Close();
    
    PCA::SetVerbose(PrintCutDiagnostics); // the PCA matrix dumps are per-event, so they follow the switch

    plothelp = new PlotHelp(); // initialize the plot help class in scope for whole program
    
    // Create a new canvas with 2x2 sub-panels
    TCanvas* c_cleaned = new TCanvas("CleanedDisplay", "Cleaned CameraPlot", 950, 1000);
    c_cleaned->Divide(2, 2);
    c_cleaned->Print(OutputFilePDFOpen.c_str());


    // new output file -- keep it OPEN so the trees created below are file-resident.
    // Creating a TTree while no file is open makes it memory-resident; large branches
    // like SEvent::posX/posY then fail to flush their baskets mid-Fill ("basket's
    // WriteBuffer failed").
    TFile *fileOutput = new TFile(OutputFileEventCleaningDataRoot.c_str(), "RECREATE");
    TTree *treeSims = new TTree("EventCleaning", "Events tree");
    eventInfo = new EventInfo(1);
    treeSims->Branch("Cleaned", eventInfo); // Store only serializable data members, not TH2F pointers

    // For sim data: copy the input "Sim" (GrOptics / SEvent) tree into the output
    // Data_EventCleaning file exactly as-is, for ALL events (triggered or not, cleaned or not).
    TTree *treeSimOut = nullptr;
    SEvent *simRunOut = new SEvent();
    if (whatData == "sim"){
        treeSimOut = new TTree("Sim", "Events tree");
        treeSimOut->Branch("GrOptics", &simRunOut);
    }

    
    // fileNamesVec.assign(fileNamesVec.begin() + 134, fileNamesVec.begin() + 135);
    
    for(int f = 0; f<10; f++){
    // for(int f = 0; f<static_cast<int>(fileNamesVec.size()); f++){
        cout << "Processing file: " << fileNamesVec[f] << endl;
        if (whatData == "bkg" ||  whatData == "muon") {
            std::string date = fileNamesVec[f].substr(7, 10); // Extract date from filename
            // remove - from the date string
            date.erase(std::remove(date.begin(), date.end(), '-'), date.end());
            folString = date; // Update folString with the extracted date
            FolderPath = Form("%s%s/",dataDir.c_str(),folString.c_str());
        } 

        // checks to make sure the  data  file is  readable
        std::string FilePath = Form("%s%s",FolderPath.c_str(),fileNamesVec[f].c_str());
        // std::string FilePath = Form("%s",fileNamesVec[f].c_str());
        if (!util->isBranchPresentInFile(FilePath, "Test")) {
            continue; // Skip to the next branch if not present
        }

        int TotalEvents;
        int nEntries;
        int simDate = 0; // date pulled from the SEvent class for sim events
        std::string simRun = ""; // sim run number pulled from the SEvent class for sim events
        int simEvent = 0; // sim event number pulled from the SEvent class for sim events
        double fileTilt = 0.0; // sims: tilt from this file's Tilt_ folder, for the horizon overlay; stays 0 for data
        TFile *fSim = nullptr;            // sim file kept open so SEvent can be read per event
        TTree *simTree = nullptr;         // "Sim" tree holding one SEvent entry per THROWN event
        std::vector<Long64_t> simEntryForTestEvent; // "Sim" entry holding the truth for each "Test" event
        SEvent *simRunData = nullptr; // reused across the event loop
        TotalEvents = 0;
        nEntries = 0;
        // loads the events based on the type of data for each file
        if (whatData !="muon" && whatData !="sim"){
            LoadEvents(FilePath, "Test");
            LoadEventsHLED(FilePath, "HLED");
            ev = new IEvent();
            evHLED = new IEvent();
            sipmInfo = new ISiPM();
            SetBranches(ev);
            SetBranchesHLED(evHLED);
            tree->SetBranchAddress("SiPM", &sipmInfo);
            treeHLED->SetBranchAddress("SiPM", &sipmInfo);
            
            
            nEntries = tree->GetEntries();
            int nEntriesHLED = treeHLED->GetEntries();
            if (nEntries==0) {
                continue; // acts like continue in parrellization
            }    
            TotalEvents = nEntries + nEntriesHLED;
            std::cout << "Total Number of Events: " << nEntries << std::endl;
            std::cout << "Total Number of HLED Events: " << nEntriesHLED << std::endl;
            
        } else if (whatData== "muon"){
            
            LoadEvents(FilePath, "Forced");
            ev = new IEvent();
            sipmInfo = new ISiPM();
            SetBranches(ev);
            tree->SetBranchAddress("SiPM", &sipmInfo);
            
            nEntries = tree->GetEntries();
            if (nEntries==0) {
                continue;
            }    
            TotalEvents = nEntries;
            std::cout << "Total Number of Forced Events: " << nEntries << std::endl;
        } else if (whatData == "sim"){
            cout << "Loading Sim File: " << FilePath << endl;
            LoadEvents(FilePath, "Test");
            ev = new IEvent();
            // sipmInfo = new ISiPM();
            SetBranches(ev);
            // tree->SetBranchAddress("SiPM", &sipmInfo);

            // Open the "Sim" tree (GrOptics branch), which holds one SEvent per event.
            // Kept open for the whole event loop so the event number can be read per entry.
            fSim = new TFile(FilePath.c_str(), "READ");
            simTree = (TTree*)fSim->Get("Sim");
            simRunData = new SEvent();
            simTree->SetBranchAddress("GrOptics", &simRunData);
            simTree->GetEntry(0);
            simDate = simRunData->GetDate();
            cout << "Sim Run Date from SEvent: " << simDate << endl;
            fileTilt = TiltFromPath(FilePath);
            cout << "Sim tilt from folder name: " << fileTilt << " deg (horizon row "
                 << Form("%.2f", HorizonRow(fileTilt)) << " of 0-15)" << endl;

            // "Test" entries carry no event identifier, so each one is matched to its "Sim" (truth)
            // entry from the layout of the file. Two layouts exist:
            //   - "Test" holds only the showers that triggered (e.g. 20260827): Test event k is the
            //     k-th "Sim" entry flagged as triggered.
            //   - "Test" holds every thrown shower (e.g. the 20260914 *_tdem files): the two trees are
            //     parallel, so Test event k is "Sim" entry k.
            // Matching the wrong way silently gives many events the same event number and truth.
            std::vector<Long64_t> triggeredSimEntries;
            for (Long64_t s = 0; s < simTree->GetEntries(); s++){
                simTree->GetEntry(s);
                if (simRunData->GetTriggered()) triggeredSimEntries.push_back(s);
            }
            Long64_t nTest = tree->GetEntries();
            cout << "Sim entries: " << simTree->GetEntries()
                 << "  triggered: " << triggeredSimEntries.size()
                 << "  Test entries: " << nTest << endl;

            simEntryForTestEvent.clear();
            if (static_cast<Long64_t>(triggeredSimEntries.size()) == nTest){
                simEntryForTestEvent = triggeredSimEntries;
                cout << "Sim truth matching: Test holds triggered showers only -> Test event k = k-th triggered Sim entry" << endl;
            } else if (simTree->GetEntries() == nTest){
                for (Long64_t s = 0; s < nTest; s++) simEntryForTestEvent.push_back(s);
                cout << "Sim truth matching: Test holds every thrown shower -> Test event k = Sim entry k" << endl;
            } else {
                cout << "WARNING: Test entries (" << nTest << ") match neither the triggered Sim entries ("
                     << triggeredSimEntries.size() << ") nor all Sim entries (" << simTree->GetEntries()
                     << "); these events get no truth info and event number -1." << endl;
            }

            // Several events sharing one (run, event number) make their outputs indistinguishable.
            std::set<std::pair<std::string, int>> seenSimEvents;
            int repeatedSimEvents = 0;
            for (Long64_t s : simEntryForTestEvent){
                simTree->GetEntry(s);
                if (!seenSimEvents.insert({simRunData->GetSumRun(), simRunData->GetSimEventNumber()}).second) repeatedSimEvents++;
            }
            if (repeatedSimEvents > 0){
                cout << "WARNING: " << repeatedSimEvents << " Test events share a (run, event number) with an earlier "
                        "event in this file; this comes from the Sim tree itself." << endl;
            }

            // Copy every entry of the input "Sim" tree into the output Sim tree exactly as-is,
            // so the truth info is saved for ALL events regardless of triggered/cleaning class.
            if (treeSimOut){
                Long64_t nSim = simTree->GetEntries();
                for (Long64_t s = 0; s < nSim; s++){
                    simTree->GetEntry(s);
                    *simRunOut = *simRunData; // copy the SEvent object exactly
                    treeSimOut->Fill();
                }
                cout << "Copied " << nSim << " Sim entries to output" << endl;
            }

            nEntries = tree->GetEntries();
            if (nEntries==0) {
                continue;
            }    
            TotalEvents = nEntries;
            std::cout << "Total Number of Sim Test Events: " << nEntries << std::endl;
        }
        
        
        // loops through each event in a file
        // for(int EventCounter = 0; EventCounter < 1000; EventCounter++){
        for(int EventCounter = 0; EventCounter < TotalEvents; EventCounter++){
            cev = new CEvent(); // creates a new Clean Event object
            eventInfo->Clear();
            cev->SetEventDate(stoi(folString));
            bool simTruthMatched = false; // set when this sim event was matched to its Sim (truth) entry
            if (whatData == "sim"){
                // read the SEvent entry for THIS event so the event number changes per event
                float simEnergy = 0.0;
                if (simTree && EventCounter < static_cast<int>(simEntryForTestEvent.size())){
                    Long64_t simEntry = simEntryForTestEvent[EventCounter];
                    simTree->GetEntry(simEntry);
                    simTruthMatched = true;
                    simEvent = simRunData->GetSimEventNumber();
                    simEnergy = simRunData->GetNeutrinoEnergy(); // switch to GetTauEnergy() for tau energy
                    simDate = simRunData->GetDate();
                    // The run name changes from sub-run to sub-run within one file, so it has to be
                    // re-read for every event rather than cached from the first entry.
                    simRun = simRunData->GetSumRun();
                } else {
                    // no matching truth entry: never reuse the previous event's number and truth
                    simEvent = -1;
                    simRun = "noSimTruth";
                }
                cev->SetEventDate(simDate); // use the date from the SEvent class for sim events
                cev->SetFilename(simRun); // use the run name from the SEvent class for sim events
                cev->SetEventNumber(simEvent);
                plothelp->SetSimEnergy(simEnergy, true); // track this energy for the saved-vs-removed plot
            } else {
                cev->SetFilename(fileNamesVec[f]);
                plothelp->SetSimEnergy(0.0, false); // non-sim: do not track energy

            }

            // Horizon overlay tilt: the SEvent's own, or the Tilt_ folder for sims converted before
            // PlotCamera stored it (those read back 0). Data keeps 0, which draws no horizon.
            HorizonTiltDeg = fileTilt;
            std::string horizonTiltSource = "Tilt_ folder name";
            if (simTruthMatched && simRunData->GetTilt() != 0) {
                HorizonTiltDeg = simRunData->GetTilt();
                horizonTiltSource = "SEvent";
            }

            // Opens this event's diagnostics block. Everything printed below it belongs to this
            // event, and FocusEventNumber is matched against the number printed here.
            PrintEventBanner(EventCounter, (whatData == "sim") ? simEvent : EventCounter, TotalEvents);
            if (whatData == "sim") {
                PrintSimTruth(simTruthMatched ? simRunData : nullptr,
                              simTruthMatched ? simEntryForTestEvent[EventCounter] : -1, EventCounter);
                PrintHorizonTilt(HorizonTiltDeg, horizonTiltSource);
            }

            std::vector<Int_t> TrigMus;
            std::vector<float> AmplitudesTimeBin;  // amplitude per pixel used by all the cuts
            std::vector<int> PeakTimeBin;          // each pixel's own peak bin
            std::vector<float> Pedestals;
            std::vector<float> PeakAmplitudes;     // amplitude at each pixel's own peak bin; used to find the triggered pixel
            std::vector<std::vector<Int_t>> traces(MaxNofChannels); // traces as read from the file, shared by every step below
            Pulse *pulse;

            // Data, HLED, muon and sim traces are all read the same way. Only the expected pulse position
            // differs (239 for the demonstrator, 184 for CARE sims), so the search window is placed around it.
            const int PeakSearchStart = ExpectedPulseBin() - PeakSearchBinsBefore;
            const int PeakSearchEnd   = ExpectedPulseBin() + PeakSearchBinsAfter;

            /**
            --------Pass 1: find every pixel's own peak
            **/
            if(whatData != "muon") {
                if (EventCounter < nEntries) {
                    tree->GetEntry(EventCounter);
                    TrigMus = ev->GetROIMusicID();
                } else {
                    treeHLED->GetEntry(EventCounter-nEntries);
                    TrigMus = evHLED->GetROIMusicID();
                    // EventCounter = EventCounter -nEntries;
                }

                for (int k = 0; k < MaxNofChannels; k++) {
                    if (EventCounter < nEntries) {
                        traces[k] = ev->GetSignalValue(k);
                        eventInfo->SetEventBranch(1);
                    } else {
                        traces[k] = evHLED->GetSignalValue(k);
                        eventInfo->SetEventBranch(9);
                    }
                    // the peak bin is the deepest sample in [PeakSearchStart, PeakSearchEnd); the
                    // amplitude is taken below, against whichever pedestal this data type uses
                    pulse = new Pulse(traces[k], PeakSearchStart, PeakSearchEnd, 512, true);
                    Pedestals.push_back(PedestalForPixel(traces[k], pulse->GetPedestal()));
                    PeakTimeBin.push_back(pulse->GetTimePeak());
                    delete pulse;
                }
                for (int k = 0; k < MaxNofChannels; k++) {
                    PeakAmplitudes.push_back(AmplitudeAtBin(traces[k], Pedestals[k], PeakTimeBin[k]));
                }

            } else if  (whatData == "muon"){
                tree->GetEntry(EventCounter);
                TrigMus = ev->GetROIMusicID();

                for (int k = 0; k < MaxNofChannels; k++) {
                    traces[k] = ev->GetSignalValue(k);
                    pulse = new Pulse(traces[k], PeakSearchStart, PeakSearchEnd, 512, true);

                    Pedestals.push_back(PedestalForPixel(traces[k], pulse->GetPedestal()));
                    PeakTimeBin.push_back(pulse->GetTimePeak());

                    delete pulse;

                }
                for (int k = 0; k < MaxNofChannels; k++) {
                    PeakAmplitudes.push_back(AmplitudeAtBin(traces[k], Pedestals[k], PeakTimeBin[k]));
                }
            }

            cev->SetTriggeredMUSICID(TrigMus[0]);

            // Saturation, measured on the same traces the peaks came from. A railed pulse has a flat
            // bottom, so its peak bin is replaced by the railed sample SaturationPeakMode asks for.
            // Data rails at the fixed ADC floor; a sim pixel rails at its own plateau (see GetRailLevels).
            std::vector<float> RailLevels = GetRailLevels(traces, Pedestals, PeakSearchStart, PeakSearchEnd);
            std::vector<IUtilities::SatSpan> SaturationSpans = FindSaturationSpans(traces, RailLevels, PeakSearchStart, PeakSearchEnd);
            ApplySaturationPeakMode(PeakTimeBin, SaturationSpans);

            // Each pixel's time for the time walk and the timing map is where its rising edge reaches
            // RiseFractionPercent of its max amplitude. CEvent still stores it under the name PeakTimeBin.
            std::vector<int> RiseTimeBin = FindRiseTimeBins(traces, Pedestals, PeakAmplitudes, PeakTimeBin);
            cev->SetPeakTimeBin(RiseTimeBin);

            /**
            --------Find the triggered pixel and the snapshot bin
            **/
            // Showers cross the camera over several ns, so reading every pixel at one fixed bin catches
            // late pixels on their rising edge. The triggered pixel is found from the per-pixel peak
            // amplitudes instead, and its peak bin becomes the moment the whole camera is read at.
            std::vector<float> CalibrationFactors = GetCalibrationFactors();
            std::vector<float> CalibratedPeakAmplitudes = ApplyCalibration(PeakAmplitudes, CalibrationFactors);
            const int TriggeredMUSICID = FindTriggeredMUSIC(CalibratedPeakAmplitudes, TrigMus[0]);
            const int TriggeredPixelID = util->GetMaximumPixelID(CalibratedPeakAmplitudes, TriggeredMUSICID);
            const int SnapshotBin      = PeakTimeBin[TriggeredPixelID];

            /**
            --------Pass 2: amplitude snapshot at the triggered pixel's peak bin
            **/
            // Data, HLED and sims read every pixel at SnapshotBin. Muon amplitudes stay each pixel's own
            // peak amplitude, as they were.
            const bool UsesSnapshot = (whatData != "muon");
            if (UsesSnapshot) {
                AmplitudesTimeBin = TakeAmplitudeSnapshot(traces, Pedestals, SnapshotBin);
            } else {
                AmplitudesTimeBin = PeakAmplitudes;
            }

            // Pixels sitting on the rail at the snapshot bin count as core pixels (see DescribeCorePixels)
            std::vector<int> SaturatedAtSnapshotBin;
            if (UsesSnapshot) {
                for (int pixelID = 0; pixelID < MaxNofChannels; pixelID++) {
                    if (traces[pixelID][SnapshotBin] <= RailLevels[pixelID]) SaturatedAtSnapshotBin.push_back(pixelID);
                }
            }
            StoreCalibratedAmplitudes(AmplitudesTimeBin);
            AmplitudesTimeBin = cev->GetAmplitudeValuesTimeBin(); // reset the  amplitudesTimeBin to be aboslute gain calibrated
            // Where the signal actually is, before any cut has run. For sims this is the first place
            // an empty or mispositioned CARE pulse shows itself.
            PrintTraceSummary(traces, TriggeredPixelID, SnapshotBin);
            PrintBrightestPixels(AmplitudesTimeBin, 10);
            cev->SetAverageAmplitude(util->GetEventAverageAmplitude(cev->GetAmplitudeValuesTimeBin()));
            cev->SetRMS(cev->GetAmplitudeValuesTimeBin());
            plothelp->AddtoRMSratioAVGamp(cev->GetRMSoverAvgAmp());
            // cout <<"Avg amplitude Before: "  << util->GetEventAverageAmplitude(AmplitudesTimeBin) << endl;
            // cout <<"Avg amplitude: "  << cev->GetAverageAmplitude() << endl;
            
            // the amplitude is different when its prior to the october 2024 site visit
            // if (cev->GetEventDate() < 20241001 && whatData != "muon") {
            //     FlasherEventsCutOff = 350; // 350 ADC/8 PE
            // }
            if (whatData != "sim") {
                cev->SetEventNumber(EventCounter);
            }
            if (EventCounter > nEntries && whatData != "muon") {
                cev->SetEventNumber(EventCounter-nEntries);
            }

            if (PrintTimingDiagnostics && DiagnoseThisEvent()) {
                PrintPedestalSummary(Pedestals, traces, TriggeredPixelID);
                PrintRailLevels(RailLevels, SaturationSpans);
                PrintTimingSummary(traces, CalibrationFactors, Pedestals, PeakTimeBin, SaturationSpans,
                                   TrigMus[0], TriggeredMUSICID, TriggeredPixelID, SnapshotBin, UsesSnapshot);
                cout << "Saturated at snapshot bin " << SnapshotBin << " (count as core pixels if they survive):";
                if (SaturatedAtSnapshotBin.empty()) cout << " none";
                for (int pixelID : SaturatedAtSnapshotBin) cout << " " << pixelID;
                cout << endl;
                PrintRiseTimes(Pedestals, PeakAmplitudes, PeakTimeBin, RiseTimeBin, AmplitudesTimeBin, TriggeredPixelID);
            }
            
            // update limit for flasher events cutoff deperciated for now
            // if (cev->GetAverageAmplitude() > (FlasherEventsCutOff/util->GetADCtoPEratio()) && whatData != "muon" && whatData != "sim") {
            //     plothelp->AddHLEDEvent(1);
            //     continue;
            // }
            

            eventInfo->SetDate(cev->GetEventDate());
            eventInfo->SetFile(cev->GetFilename());
            eventInfo->SetEvent(cev->GetEventNumber());
            eventInfo->SetTriggeredMUSICID(cev->GetTriggeredMUSICID());
            eventInfo->SetTriggeredPixelID(cev->GetMaxAmplitudePixelID());
            eventInfo->SetRMS_AvgAmp(cev->GetRMSoverAvgAmp());
            eventInfo->SetRMS(cev->GetRMS());
            eventInfo->SetAvgAmp(cev->GetAverageAmplitude());
            
            
            const bool cameraLitEvenly = (cev->GetRMSoverAvgAmp() < 0.15 && whatData != "muon");
            PrintCutVerdict("RMS / average amplitude", cameraLitEvenly, Form("%.3f", cev->GetRMSoverAvgAmp()),
                            ">= 0.15 (an evenly lit camera is a flasher, not a shower)", 9);
            if (cameraLitEvenly) {
                plothelp->AddEventFlags(9);
                eventInfo->SetEventFlag(9); // HLED  event
                treeSims->Fill();
                delete cev;
                continue;
            }
            
            
            
            // Remove events where the door is closed. This can be done by checking the HV currents of the SiPMs and 
            // if the sky seeing Channels are below 3.8 mA you know that the door is closed.
            // if (ev->Gethvc()[0] < 3.8 && ev->Gethvc()[2] < 3.8){
                //     plothelp->AddCleanedEvent(1);
                //     continue;
                // }
                
                
            // SaturationSpans was measured right after the peak search, on the same traces and in the same
            // window around the expected pulse (HLED events use their own branch there). Every pixel with
            // a railed sample is saturated. This runs before cleaning, so it is the raw list; CompletePanel3 narrows it
            // to the pixels that actually survived. The cut below has to stay on the raw count because
            // there are no surviving pixels yet.
            SaturatedPixelIDsRaw.clear();
            for (std::size_t p = 0; p < SaturationSpans.size(); p++) {
                if (SaturationSpans[p].count > 0) SaturatedPixelIDsRaw.push_back((Int_t)p);
            }
            cev->SetSaturatedPixels((int)SaturatedPixelIDsRaw.size());
            const bool tooManySaturatedPixels = (cev->GetSaturatedPixels() > SaturatedPixelCutoff && whatData != "muon");
            PrintCutVerdict("saturated pixels", tooManySaturatedPixels, Form("%d", cev->GetSaturatedPixels()),
                            Form("<= %d", SaturatedPixelCutoff), 1);
            if (tooManySaturatedPixels) {
                plothelp->AddEventFlags(1);
                eventInfo->SetEventFlag(1);
                treeSims->Fill();
                delete cev;
                continue;
            }
            // cout << "triggered Music pair: " << util->GetMUSICMate(cev->GetTriggeredMUSICID(), true) << endl;
            
            // //check which music has the highest amplitude 
            // cout << "triggered Music amp: " << util->GetMaximum(AmplitudesTimeBin,cev->GetTriggeredMUSICID()) << endl;
            // cout << "triggered Music pair amp: " << util->GetMaximum(AmplitudesTimeBin,util->GetMUSICMate(cev->GetTriggeredMUSICID(), true)) << endl;
            // The MUSIC swap was decided on the per-pixel peak amplitudes (FindTriggeredMUSIC above), so
            // the triggered MUSIC here is the one the snapshot bin came from.
            if (TriggeredMUSICID != cev->GetTriggeredMUSICID()) {
                cev->SetTriggeredMUSICID(TriggeredMUSICID);
                eventInfo->SetTriggeredMUSICID(cev->GetTriggeredMUSICID());
                eventInfo->SetTriggeredPixelID(cev->GetMaxAmplitudePixelID());
                // cout << "Replaced. " << endl;
                // cout << "triggered Music pair: " << util->GetMUSICMate(cev->GetTriggeredMUSICID(), true) << endl;
            }
            
            cev->SetMaxAmplitude(util->GetMaximum(cev->GetAmplitudeValuesTimeBin(),cev->GetTriggeredMUSICID()));
            const float triggeredAmpCutOffPE = TriggeredChannelAmpCutOff/util->GetADCtoPEratio();
            const bool triggeredPixelTooDim = cev->GetMaxAmplitude() < triggeredAmpCutOffPE;
            PrintCutVerdict("triggered pixel amplitude", triggeredPixelTooDim,
                            Form("%.2f PE", cev->GetMaxAmplitude()),
                            Form(">= %.2f PE (in MUSIC %d)", triggeredAmpCutOffPE, cev->GetTriggeredMUSICID()), 1);
            if (triggeredPixelTooDim) {
                plothelp->AddEventFlags(1);
                eventInfo->SetEventFlag(1);
                treeSims->Fill();
                delete cev;
                continue;
            }

            cev->SetPedestalValues(Pedestals);
            // Set explicitly rather than re-derived from the snapshot, so the pixel the time walk starts
            // from is the same pixel that defined SnapshotBin.
            cev->SetMaxAmplitudePixelID(TriggeredPixelID);

            // removed events that have pixel 76 
            // if (cev->GetMaxAmplitudePixelID() == 78 || cev->GetMaxAmplitudePixelID() == 79 || cev->GetMaxAmplitudePixelID() == 77|| cev->GetMaxAmplitudePixelID() == 76 || cev->GetMaxAmplitudePixelID() == 75) {
            //     cout << "Triggered Pixel ID 78 or 79 skipping event." << endl;
            //     plothelp->AddPreCleaned(1);
            //     continue;
            // }

            // removed the hotspot on the worst night 
            const bool onKnownHotSpot = ((cev->GetTriggeredMUSICID() == 18 || cev->GetTriggeredMUSICID() == 19) &&
                                        (folString == "20250305" || folString == "20250331" || folString == "20250328"));
            PrintCutVerdict("hot spot MUSIC 18/19", onKnownHotSpot, Form("MUSIC %d", cev->GetTriggeredMUSICID()),
                            "not MUSIC 18 or 19 on 20250305/0328/0331", 1);
            if (onKnownHotSpot){
                eventInfo->SetEventFlag(1);
                plothelp->AddEventFlags(1);
                treeSims->Fill();
                delete cev;
                continue;
            }
            
       
            cev->SetSumAmplitude(std::accumulate(AmplitudesTimeBin.begin(), AmplitudesTimeBin.end(), 0.0));
            
            
            /** 
            --------Create Panel 1 
            **/
            c_cleaned->cd(1);
            gPad->SetLeftMargin(0.15);  // Increase left margin
            gPad->SetBottomMargin(0.25); // Increase bottom margin
            gPad->SetRightMargin(0.15);
            
            std::string filenameTitle = cev->GetFilename();
            if (whatData != "sim"){
                int pos = filenameTitle.find("T");
                filenameTitle = filenameTitle.substr(pos+1, 5);
            } 
            // create TH2D for Camera
            hcam_panel1 = new TH2F("hcam_panel1", Form("Calibrated Img --N# %i --F# %s-- E# %i ", cev->GetEventDate(), filenameTitle.c_str() ,cev->GetEventNumber()), 16, -0.5, 15.5, 16, -0.5, 15.5);
            hcam_panel1->SetDirectory(0);
            cev->SetPanel1(hcam_panel1);
            DrawSimHorizon(hcam_panel1, kWhite);
            // eventInfo->SetHPanel1(hcam_panel1);
            // cev->SetAverageAmplitudeRMS(cev->GetAmplitudeValuesTimeBin());
            // Add a subtitle using TLatex
            TLatex *subtitle = new TLatex();
            subtitle->SetNDC(); // Set to Normalized Device Coordinates (NDC)
            subtitle->SetTextSize(0.03);
            
            subtitle->DrawLatex(0.1, 0.92, Form("Avg Amp Whole Camera : %.2f - Avg RMS: %.2f", cev->GetAverageAmplitude(), cev->GetRMS()));
            subtitle->DrawLatex(0.25, 0.12, Form("Triggered MUSIC: %i  Triggered Pixel: %i",cev->GetTriggeredMUSICID(),cev->GetMaxAmplitudePixelID()));
            // Saturation on the raw image, before cleaning: the boxes say where, the count says
            // how much. Panel 3 repeats this for the saturated pixels that survived cleaning.
            DrawSaturationMarkers(SaturatedPixelIDsRaw);
            subtitle->DrawLatex(0.25, 0.07, Form("Saturated Pixels: %i", (int)SaturatedPixelIDsRaw.size()));
            subtitle->DrawLatex(0.25, 0.02, ImageTimeBinLabel(UsesSnapshot, SnapshotBin).c_str());
            delete subtitle;
            
            eventInfo->SetHPanel1(hcam_panel1);
            
            /**
            --------Amplitude pre-cut (not displayed)
            **/
            // An event needs PixelSurviveCutOff pixels at or above CorePixelAmpCutOff at the snapshot bin,
            // or it is removed with flag 2. At 0 PE every pixel passes, so this removes nothing; raise
            // CorePixelAmpCutOff to switch it back on. CEvent::SetPanel2 draws its image on the current
            // pad, so it is pointed at pad 2, where the 8 PE walk image replaces it below.
            c_cleaned->cd(2);
            TH2F* hcam_amplitudeCut = (TH2F*)hcam_panel1->Clone("hcam_amplitudeCut");
            hcam_amplitudeCut->SetDirectory(0);
            cev->SetPanel2(hcam_amplitudeCut, AmplitudePreCutPE());
            delete hcam_amplitudeCut;
            const int amplitudeCutPixels = (int)cev->GetSurvivingPixelPanel2().size();
            const bool tooFewAmplitudeCutPixels = amplitudeCutPixels < PixelSurviveCutOff;
            PrintCutVerdict("amplitude pre-cut pixels", tooFewAmplitudeCutPixels, Form("%d of %d", amplitudeCutPixels, MaxNofChannels),
                            Form(">= %d pixels at or above %d PE", PixelSurviveCutOff, AmplitudePreCutPE()), 2);
            if (tooFewAmplitudeCutPixels) {
                eventInfo->SetEventFlag(2);
                plothelp->AddEventFlags(2);
                treeSims->Fill();
                delete hcam_panel1;
                delete cev;
                continue;
            }

            /**
            --------Time walk: panel 2 is the 8 PE seed stage, panel 3 adds the 4 PE extension stage
            **/
            // Both start from the whole snapshot image (panel 1): the walk reads each pixel's signal around
            // its own expected time, so a pixel under the threshold at the snapshot bin can still join.
            // Survivors keep their snapshot amplitude in the images.
            const bool diagnoseTimeWalk = PrintTimingDiagnostics && DiagnoseThisEvent();
            const PixelSignals pixelSignals = {traces, Pedestals, CalibrationFactors, PeakTimeBin};
            RunTimeWalk(pixelSignals, RiseTimeBin, TriggeredPixelID, AmplitudesTimeBin, diagnoseTimeWalk);

            /** 
            --------Create Panel 2: the 8 PE seed stage
            **/
            c_cleaned->cd(2);
            gPad->SetLeftMargin(0.15);  // Increase left margin
            gPad->SetBottomMargin(0.25); // Increase bottom margin
            gPad->SetRightMargin(0.15);
            hcam_panel2 = (TH2F*)hcam_panel1->Clone("hcam_panel2");
            hcam_panel2->SetTitle(Form("%i PE Walk Img --N# %i --F# %s-- E# %i ", SeedThresholdPE,
                                       cev->GetEventDate(), filenameTitle.c_str() ,cev->GetEventNumber()));
            hcam_panel2->SetDirectory(0);
            KeepOnlyPixels(hcam_panel2, SurvivingPixelsPanel2);
            DrawTimeWalkPanel(hcam_panel2);

            /** 
            --------Create Panel 3: the 4 PE extension stage added, the cleaned image
            **/
            // Drawn, with its Hillas ellipse, by CompletePanel3 once the event has passed every cut.
            hcam_panel3 = (TH2F*)hcam_panel1->Clone("hcam_panel3");
            hcam_panel3->SetTitle(Form("%i PE Walk Img --N# %i --F# %s-- E# %i ", ExtensionThresholdPE,
                                       cev->GetEventDate(), filenameTitle.c_str() ,cev->GetEventNumber()));
            hcam_panel3->SetDirectory(0);
            SurvivingTotalAmpPanel3 = KeepOnlyPixels(hcam_panel3, SurvivingPixelsPanel3);
            DescribeCorePixels(hcam_panel3, TriggeredPixelID, SaturatedAtSnapshotBin, diagnoseTimeWalk);


            //find configure events with amp high and above 13
            const bool isConfigureEvent = IsConfigureEvent(SurvivingPixelsPanel3, SurvivingTotalAmpPanel3);
            PrintCutVerdict("configure (single-SIAB) event", isConfigureEvent, isConfigureEvent ? "yes" : "no",
                            "not >=15 bright pixels confined to one SIAB", 3);
            if (isConfigureEvent){
                eventInfo->SetEventFlag(3);
                plothelp->AddEventFlags(3);
                treeSims->Fill();
                delete hcam_panel1;
                delete hcam_panel2;
                delete hcam_panel3;
                delete cev;
                continue;
            }
            const bool tooFewPanel3Pixels = SurvivingPixelsPanel3.size() < static_cast<std::vector<int>::size_type>(PixelSurviveCutOff);
            PrintCutVerdict("panel-3 surviving pixels", tooFewPanel3Pixels, Form("%d", (int)SurvivingPixelsPanel3.size()),
                            Form(">= %d pixels kept by the %i/%i PE time walk", PixelSurviveCutOff, SeedThresholdPE, ExtensionThresholdPE), 3);
            if (tooFewPanel3Pixels) {
                eventInfo->SetEventFlag(3);
                plothelp->AddEventFlags(3);
                treeSims->Fill();
                delete hcam_panel1;
                delete hcam_panel2;
                delete hcam_panel3;
                delete cev;
                continue;
            }
            
            const bool tooSmall = SurvivingTotalAmpPanel3 <= sizeCutoff;
            PrintCutVerdict("panel-3 size", tooSmall, Form("%.2f PE", SurvivingTotalAmpPanel3),
                            Form("> %d PE", sizeCutoff), 3);
            if (tooSmall) {
                eventInfo->SetEventFlag(3);
                plothelp->AddEventFlags(3);
                treeSims->Fill();
                delete hcam_panel1;
                delete hcam_panel2;
                delete hcam_panel3;
                delete cev;
                continue;
            }
            
            /** 
            --------PCA of panel 3
            **/
            //Create the PCA
            PCA pca(2);
            
            // do the Principle componet analysis on the cleaned data
            LoadDataPCA(pca,hcam_panel3,SurvivingTotalAmpPanel3);
            
            TVectorD eigenVals; // remove the use of the pointer because this is already in the PCA class and calulated in LoadDataPCA
            TMatrixD eigenVecs; 
            std::vector<double> EllipicRatio=CreateWLRatio(pca,eigenVals,eigenVecs); // to do make EllipicRatio a member of cev
            //ex. cev->SetEllipticRatio(EllipicRatio);
            const bool tooRound = EllipicRatio[2] > (WLRatioCutOff/100.0);
            PrintCutVerdict("width / length ratio", tooRound, Form("%.3f", EllipicRatio[2]),
                            Form("<= %.2f (length %.2f, width %.2f)", WLRatioCutOff/100.0, EllipicRatio[0], EllipicRatio[1]), 4);
            if (tooRound) {
                eventInfo->SetEventFlag(4);
                plothelp->AddEventFlags(4);
                treeSims->Fill();
                delete hcam_panel1;
                delete hcam_panel2;
                delete hcam_panel3;
                delete cev;
                continue;
            }
            
            
            eventInfo->SetEventFlag(5);
            plothelp->AddEventFlags(5);
            PrintEventOutcome("passed every cut and was kept", 5);
            
            plothelp->AddtoEventDate(cev->GetEventDate());

            c_cleaned->cd(3);
            gPad->SetLeftMargin(0.15);  // Increase left margin
            gPad->SetBottomMargin(0.25); // Increase bottom margin
            gPad->SetRightMargin(0.15);
            CompletePanel3(pca,hcam_panel3,cev,EllipicRatio,eigenVals,eigenVecs,eventInfo);
            treeSims->Fill();

            /** 
            --------Create Panel 4: timing of the panel-3 pixels
            **/
            // Timing heat maps: each pixel's rise bin relative to the triggered pixel's, for the panel-2 and
            // panel-3 survivors. Panel 4 shows the panel-3 map; c_time shows both side by side and goes to
            // the display ROOT file only. Neither is attached to eventInfo (stays out of the data files).
            const int triggeredRiseBin = RiseTimeBin[TriggeredPixelID];
            const int timeMapRange = TimingMapRange(cev, triggeredRiseBin);
            TH2F* hcam_time_panel2 = FillTimingHeatMap(cev, triggeredRiseBin, filenameTitle, SurvivingPixelsPanel2,
                                                       "hcam_time_panel2", Form("%i PE Walk", SeedThresholdPE), timeMapRange);
            TH2F* hcam_time_panel3 = FillTimingHeatMap(cev, triggeredRiseBin, filenameTitle, SurvivingPixelsPanel3,
                                                       "hcam_time_panel3", Form("%i PE Walk", ExtensionThresholdPE), timeMapRange);
            c_cleaned->cd(4);
            gPad->SetBottomMargin(0.25); // as the other panels; DrawTimingMapPad sets the side margins
            DrawTimingMapPad(hcam_time_panel3, TriggeredPixelID);
            
            // For sims, draw the truth info (neutrino energy, distance to and
            // location of the emergence point, distance to the tau decay point)
            // in the empty space under panel 2, clear of the panel 1 subtitle and the panel 3/4 titles.
            // Drawn BEFORE the canvas is written so it is saved in the .root canvas too.
            TPaveText *simInfo = nullptr;
            if (whatData == "sim" && simTruthMatched){
                // Distance to emergence point (magnitude of the telescope position vector)
                double emergenceDistance = util->GetEmergenceDistance(simRunData->GetTelescope_Xpos(),
                                                                      simRunData->GetTelescope_Ypos(),
                                                                      simRunData->GetTelescope_Zpos());

                // Emergence angle (azimuth wrapped to [-180,180] from the telescope axis)
                double emergenceAngle = util->GetEmergenceAngle(simRunData->GetAzimuthAngle());

                // Distance to the tau decay point. The shower start (= decay point) and the
                // telescope position are both stored in the same emergence-point-centred frame,
                // so the telescope-to-decay-point distance is the magnitude of their difference.
                // (|showerStart| on its own is the decay length measured from the emergence point.)
                double decayDx = simRunData->GetTelescope_Xpos() - simRunData->GetShowerStartX();
                double decayDy = simRunData->GetTelescope_Ypos() - simRunData->GetShowerStartY();
                double decayDz = simRunData->GetTelescope_Zpos() - simRunData->GetShowerStartZ();
                double decayDistance = std::sqrt(decayDx*decayDx + decayDy*decayDy + decayDz*decayDz);

                c_cleaned->cd(0);
                simInfo = new TPaveText(0.56, 0.505, 0.96, 0.585, "NDC");
                simInfo->SetFillColorAlpha(0, 0.0);
                simInfo->SetFillStyle(0);
                simInfo->SetBorderSize(0);
                simInfo->SetTextSize(0.018);
                simInfo->AddText(Form("Neutrino Energy: %.3g GeV", simRunData->GetNeutrinoEnergy()));
                simInfo->AddText(Form("Emergence Angle: %.3g deg", emergenceAngle));
                simInfo->AddText(Form("Distance to Emergence Point: %.3g m", emergenceDistance));
                simInfo->AddText(Form("Distance to Decay Point: %.3g m", decayDistance));
                simInfo->Draw();
                c_cleaned->Update();
            }

            DrawTimingHeatMaps(hcam_time_panel2, hcam_time_panel3, TriggeredPixelID);

            // Peak-amplitude display: the same four panels and surviving pixels as c_cleaned, but every
            // pixel shows its amplitude at its own peak time bin instead of at the triggered pixel's.
            std::vector<TH2F*> peakAmpMaps = DrawPeakAmplitudeDisplay(pixelSignals, filenameTitle, hcam_time_panel3, eventInfo);
            if (simInfo) {
                c_peakamp->cd(0);
                simInfo->Draw(); // the same truth box as on c_cleaned, since this is the display in the PDF
            }
            c_cleaned->cd(); // DrawTimingHeatMaps and DrawPeakAmplitudeDisplay switched the current pad

            // Everything this event puts in the display ROOT file, all named after the 2x2 display N...:
            //   Pre...TH2F  panel 1 (calibrated image)        Post...TH2F  panel 3 (cleaned image)
            //   TimePanel2...TH2F / TimePanel3...TH2F  timing maps of panels 2 and 3 (panel 3's is panel 4)
            //   Time...  canvas with both timing maps   PeakAmp...  peak-amplitude display (+ Pre/Post TH2F)
            const std::string eventName = Form("N%i-F%s-E%i", cev->GetEventDate(), filenameTitle.c_str(), cev->GetEventNumber());
            file = new TFile(OutputFileRoot.c_str(), "UPDATE");
            c_cleaned->Write(eventName.c_str());
            hcam_panel1->Write(("Pre" + eventName + "TH2F").c_str());
            hcam_panel3->Write(("Post" + eventName + "TH2F").c_str());
            hcam_time_panel2->Write(("TimePanel2" + eventName + "TH2F").c_str());
            hcam_time_panel3->Write(("TimePanel3" + eventName + "TH2F").c_str());
            c_time->Write(("Time" + eventName).c_str());
            c_peakamp->Write(("PeakAmp" + eventName).c_str());
            peakAmpMaps[0]->Write(("PeakAmpPre" + eventName + "TH2F").c_str());
            peakAmpMaps[2]->Write(("PeakAmpPost" + eventName + "TH2F").c_str());
            file->Close();

            // The PDF gets the peak-amplitude display; the snapshot display stays in the ROOT file only
            if (SnapshotDisplayInPDF) c_cleaned->Print(OutputFilePDF.c_str());
            c_peakamp->Print(OutputFilePDF.c_str());
            if (simInfo){
                c_cleaned->GetListOfPrimitives()->Remove(simInfo); // keep it off later pages
                c_peakamp->GetListOfPrimitives()->Remove(simInfo);
                delete simInfo;
            }


            delete hcam_panel1;
            delete hcam_panel2;
            delete hcam_panel3;
            c_time->Clear("D"); // drop the maps from c_time's pads before they are deleted
            delete hcam_time_panel2;
            delete hcam_time_panel3;
            c_peakamp->Clear("D"); // same for c_peakamp
            for (TH2F* peakAmpMap : peakAmpMaps) delete peakAmpMap;
            delete cev;

        }
     
        if (whatData != "muon"){
            delete treeHLED;
            delete evHLED;
        }
        if (whatData == "sim"){
            delete simRunData; // close the Sim tree opened for this file
            delete fSim;       // also deletes simTree (owned by the file)
            simRunData = nullptr;
            simTree = nullptr;
            fSim = nullptr;
        }
        delete tree;
        delete ev;
    }
    cout << "we here" << endl;
    fileOutput->cd();
    treeSims->Write();
    if (treeSimOut) treeSimOut->Write(); // sims: input Sim tree copied as-is for all events
    fileOutput->Close();
    file = new TFile(OutputFileRoot.c_str(), "UPDATE");
    plothelp->PlotEventFlags(c_cleaned, OutputFilePDF);
    plothelp->PlothSimEnergySavedVsNotSaved(c_cleaned, OutputFilePDF); // sims only: saved vs removed energy
    plothelp->PlothWL(c_cleaned,OutputFilePDF,outDir,to_string(cev->GetEventDate()));
    plothelp->PlothdistLandW(c_cleaned,OutputFilePDF);
    plothelp->PlothSize(c_cleaned,OutputFilePDF);    
    plothelp->PlothSizeConc(c_cleaned,OutputFilePDF);    
    plothelp->PlothSizeSPC(c_cleaned,OutputFilePDF);    
    plothelp->PlothCoreRatio(c_cleaned,OutputFilePDF);    
    plothelp->PlothNumberofCores(c_cleaned,OutputFilePDF);    
    plothelp->PlothCoreRatioandNumberofCores(c_cleaned,OutputFilePDF);    
    plothelp->PlothWLConc(c_cleaned,OutputFilePDF);    
    plothelp->PlothWLSPC(c_cleaned,OutputFilePDF); 
    plothelp->PlothSIZEWL(c_cleaned,OutputFilePDF);   
    plothelp->PlothLvW(c_cleaned,OutputFilePDF);
    plothelp->PlothSPCConc(c_cleaned,OutputFilePDF);    
    plothelp->PlothCOG(c_cleaned,OutputFilePDF);    
    plothelp->PlothTPID(c_cleaned,OutputFilePDF);
    plothelp->PlothRMSratioAvg(c_cleaned, OutputFilePDF);
    plothelp->PlothEventDateTriggeredPixel(c_cleaned, OutputFilePDF);
    plothelp->PlothncoreSIZE(c_cleaned, OutputFilePDF);
    plothelp->PlothnCoreSPC(c_cleaned, OutputFilePDF);
    plothelp->PlothCRSPC(c_cleaned, OutputFilePDF);
    plothelp->PlothSIZEWandL(c_cleaned, OutputFilePDF);
    plothelp->PlotPixelsDistanceToMajorAxis(c_cleaned, OutputFilePDF);
    plothelp->PlotPixelsRatioDistanceToMajorAxis(c_cleaned, OutputFilePDF);
    plothelp->PlotdistRMSandWeightedRMS(c_cleaned, OutputFilePDF);
    plothelp->PlotPixelsOnMajorAxis(c_cleaned, OutputFilePDF);
    plothelp->PlotPixelsOffMajorAxis(c_cleaned, OutputFilePDF);
    plothelp->PlothRMSvsRatioDistance(c_cleaned, OutputFilePDF);
    plothelp->PlothWRMSvsRatioDistance(c_cleaned, OutputFilePDF);
    plothelp->PlothWLvsRatioPixelsMajorAxis(c_cleaned, OutputFilePDF);
    plothelp->PlothDistancevsRatioPixelsMajorAxis(c_cleaned, OutputFilePDF);
    plothelp->PlothCoreOverSPCvsRatioPixelsMajorAxis(c_cleaned, OutputFilePDF);
    plothelp->PlothSPCvsRatioPixelsMajorAxis(c_cleaned, OutputFilePDF);
    plothelp->PlothSurvivingPixelsX(c_cleaned, OutputFilePDF);
    plothelp->PlothSurvivingPixelsY(c_cleaned, OutputFilePDF);
    plothelp->PlothCOGx(c_cleaned, OutputFilePDF);
    plothelp->PlothCOGy(c_cleaned, OutputFilePDF);
    plothelp->PlothAngle(c_cleaned, OutputFilePDF);
    plothelp->PlothUpDownSym(c_cleaned, OutputFilePDF);
    plothelp->PlothSaturatedPixels(c_cleaned, OutputFilePDF);
    plothelp->PlothSaturatedOverSurviving(c_cleaned, OutputFilePDF);
    plothelp->PlothSaturatedPixelsX(c_cleaned, OutputFilePDF);
    plothelp->PlothSaturatedPixelsY(c_cleaned, OutputFilePDF);
    //plothelp->PlothOnOffMajorAxisvsratio(c_cleaned, OutputFilePDF);

    file->Close();
    c_cleaned->Print(OutputFilePDFClose.c_str());
    util->setFilePermissions(OutputFileRoot.c_str());
    util->setFilePermissions(OutputFileEventCleaningDataRoot.c_str());
    cout << "Completed!" << endl;   
    // cout << "Total Events: " << plothelp->GetTotalEvent() << endl;
    delete plothelp;


    
}






