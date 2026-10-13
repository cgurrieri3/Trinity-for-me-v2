#include "CEvent.h"
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <queue>
#include <set>
#include <string>

// Constructor
CEvent::CEvent() : TriggeredMUSICID(-1) {
    // Initialize TriggeredMUSICID to -1 (default value)
}

void CEvent::LoadFlasherCalibration(std::string CalibrationFactorDir, std::string folString, double biasVoltage) {
    // Start by creating a histogram for panel 1 then use that to copy and create the 2nd panel
    // create the first panel, apply the Gain, apply the ADC to PE 
    TFile *flasher = new TFile(Form("%s%s_FlasherCalibration_Factor.root", CalibrationFactorDir.c_str(), folString.c_str()));
    if (!flasher || flasher->IsZombie()) {
        std::cerr << "Error: Could not open file " << Form("%s%s_FlasherCalibration_Factor.root", CalibrationFactorDir.c_str(), folString.c_str()) << std::endl;
    }

    std::string biasString = Form("CamFlasher%s", std::to_string(static_cast<int>(biasVoltage)).c_str());
    hflasher = (TH2F*)flasher->Get(biasString.c_str());
    if (!hflasher || hflasher->IsZombie()) {
        std::cerr << "Error: Could not find '" << biasString << "' in the file or it is invalid." << std::endl;
    }
    hflasher->SetDirectory(0);
    flasher->Close();
    delete flasher;
}

// this is for when you set all your gain and such when outside the script
void CEvent::SetAmplitudeValues(std::vector<float> amp) {
    // sets the peak amplitude  and assumes that gain calibration as already been set to the CEvent object
    AmplitudeValues = amp;
}

// setting the abs gain and the flasher calibration within the CEvent object for when you pull a time bin
void CEvent::SetAmplitudeValuesTimeBin(std::vector<float> amp, std::string CalibrationFactorDir, std::string folString, std::vector<double> absgain, double biasVoltage) {
    // sets the peak amplitude  and assumes that gain calibration as already been set to the CEvent object
    // Load the flasher calibration histogram
    LoadFlasherCalibration(CalibrationFactorDir, folString, biasVoltage);
    TH2F* hf = GetFlasherCalibration();
    
    // set the absolute gain from the ISiPM calss
    SetAbsoluteGain(absgain);

    AmplitudeValuesTimeBin = amp;
    for (int i = 0; i < MaxNofChannels; i++) {
        int nx, ny;
        IPlotTools::FindBin(i, &nx, &ny);
        double flashercalibratedamp = amp[i]*(1.0/hf->GetBinContent(nx + 1, ny + 1));
        // std::cout << "flasher calib: " << flashercalibratedamp << std::endl;
        AmplitudeValuesTimeBin[i] = flashercalibratedamp*(1.0/AbsoluteGain[i]);
        // std::cout << "AmplitudeValuesTimeBin[" << i << "]: " << AmplitudeValuesTimeBin[i] << std::endl;
    }
    delete hf;
}

void CEvent::SetRMS(std::vector<float> amp){
    float sumsq = 0;
    for (std::size_t i = 0; i < amp.size(); i++){
        sumsq += amp[i]*amp[i];
        // cout << "sumsq: " << sumsq << endl;
        RMS = sqrt(sumsq)/MaxNofChannels;
        // cout << "RMS: " << RMS << endl;
    }
}

void CEvent::SetPanel1(TH2F* hcam_panel1) {
    hcam_panel1->SetStats(0);
    hcam_panel1->GetXaxis()->SetLabelSize(0.03);
    hcam_panel1->GetYaxis()->SetLabelSize(0.03);
    hcam_panel1->SetXTitle("Pixel Bin along x axis of Camera");
    hcam_panel1->SetYTitle("Pixel Bin along y axis of Camera");
    hcam_panel1->GetXaxis()->SetTitleOffset(1.1); // Adjust X-axis title offset
    hcam_panel1->GetYaxis()->SetTitleOffset(1.0); // Adjust Y-axis title offset
    hcam_panel1->SetMinimum(-1);
    
    for(int j = 0; j<MaxNofChannels; j++){
        int nx, ny;
        IPlotTools::FindBin(j, &nx, &ny);
        hcam_panel1->SetBinContent(nx + 1, ny + 1, AmplitudeValuesTimeBin[j]);
    }
    
    hcam_panel1->Draw("colz");
    IPlotTools::DrawMUSICBoundaries();
    
}

void CEvent::SetPanel2(TH2F* hcam_panel2, int CorePixelAmpCutOff) {

    
    for(int j = 0; j<MaxNofChannels; j++){
        
        int nx, ny;
        IPlotTools::FindBin(j, &nx, &ny);
        float binContent= hcam_panel2->GetBinContent(nx+1,ny+1);
        if (binContent >= CorePixelAmpCutOff) { 
            hcam_panel2->SetBinContent(nx + 1, ny + 1, binContent);
            // std::cout << "Surviving Panel 2 Pixel: " << j << " Amplitude: " << AmplitudeValuesTimeBin[j] << std::endl;
            SurvivingPixelPanel2.push_back(j);
        } else {
            hcam_panel2->SetBinContent(nx + 1, ny + 1, 0);
        }
    }

    hcam_panel2->Draw("colz");
    hcam_panel2->SetMinimum(7.95);
    IPlotTools::DrawMUSICBoundaries();
}


// identifies if a neighbor to the triggered pixel is a core and then only checks the surrounding pixels of the core to be saved

void CEvent::SetPanel3(TH2F* hcam_panel3, std::string neighborDir, int NumberOfCores, float CoreCutOff , float Crosstalkallowance, int maxTimeStep, int maxTimeDrift, bool verbose) {
    // std::vector<int> PixelPeakTimes = PeakTimeBin;
    // int TriggeredPixelPeakTime = PixelPeakTimes[MaxAmplitdePixelID];
    // std::cout << "Triggered Pixel Peak Time: " << TriggeredPixelPeakTime << std::endl;
    // Get neighboring pixels to triggered pixel
    
    // get the 2nd core pixel 
    
    // Check to see if there is a pixel that is connected to the triggered pixel that is at least PE_second > PE_Trig - 1.5sqrt(PE_Trig)
    isSecondCorePixel(hcam_panel3,MaxAmplitdePixelID, neighborDir, CoreCutOff ,NumberOfCores, verbose);
    if (SecondCoreID == -1){
        SurvivingPixelPanel3 = {MaxAmplitdePixelID};
        if (verbose) {
            std::cout << "Time walk skipped: no second core pixel, only the triggered pixel is kept" << std::endl;
        }
        return;
    }

    if (verbose) std::cout << "Trig Id: " << MaxAmplitdePixelID << " Core ID; " << SecondCoreID << std::endl;

    // is you are using core avg
    // double coreAverageAmp = (MaxEventAmplitde +SecondCoreAmp)/2.0;
    // std::cout << "Core Average Amplitude: " << coreAverageAmp << std::endl;
    // set  them to the surving 
    // std::vector<int> neighborsurvivingPanel3Pixels = {MaxAmplitdePixelID, SecondCoreID};

    // max between the two core pixels
    double maxcoreAmp = std::max(MaxEventAmplitde, SecondCoreAmp);

    
    
    // SetCoreIDs({-300});
    // std::cout<< "Try getting Cores:" <<updateCore[0] << std::endl;

    // go around triggered pixel and find those that are 30% of core avg
    // go  around second pixel and find those that are   30% ofthe second core pixel 
    // CheckCloseNeigborsOnly(hcam_panel3, MaxAmplitdePixelID, SecondCoreID, 0.3 * maxcoreAmp, neighborDir);
    CheckCloseNeigborsOnly(hcam_panel3, MaxAmplitdePixelID, GetCoreIDs(), Crosstalkallowance * maxcoreAmp, neighborDir, verbose);
    
    

    std::vector<int> updateCore =  SecondCoreIDs;
    if (verbose) {
        std::cout << "Core pixels:";
        for (int coreID : updateCore) std::cout << " " << coreID;
        std::cout << std::endl;
    }

    // Timing + connectivity in one pass. Showers cross the camera over several ns, so pixels far
    // from the triggered pixel peak in later FADC bins. Instead of one fixed time window for every
    // pixel, the cluster grows outward from the triggered pixel (all 8 neighbours of a core, broadside
    // neighbours of any other pixel) and each new pixel is compared with the pixel it was reached
    // from. Pixels not reached this way (including a second core that is not connected in time) are
    // cleaned in the loop below.
    if (verbose) {
        int pixelsSharingBin = 0;
        int connectedSharingBin = CountSameTimeBinSurvivors(hcam_panel3, neighborDir, pixelsSharingBin);
        std::cout << "No time steps: " << pixelsSharingBin << " pixels sit in the triggered pixel's own time bin ("
                  << PeakTimeBin[MaxAmplitdePixelID] << "), " << connectedSharingBin
                  << " of them connect to it" << std::endl;
    }
    std::vector<int> neighborsurvivingPanel3Pixels;
    TimePropagatedNeighborhoodChecker(hcam_panel3, MaxAmplitdePixelID, neighborsurvivingPanel3Pixels, neighborDir, maxTimeStep, maxTimeDrift, verbose);
    SurvivingPixelPanel3={};
    SurvivingPixelTotalAmpPanel3 = 0.0;
    totalCoreAmp = 0;
    for(int j = 0; j < MaxNofChannels; j++){
        int nx, ny;
        IPlotTools::FindBin(j, &nx, &ny);
        float binContent = hcam_panel3->GetBinContent(nx+1,ny+1);
        
        if ( binContent > 0.0 && std::find(neighborsurvivingPanel3Pixels.begin(), neighborsurvivingPanel3Pixels.end(), j) != neighborsurvivingPanel3Pixels.end()) {
            // std::cout << "Surviving Panel 3 Pixel: " << j << " Amplitude: " << AmplitudeValuesTimeBin[j] << std::endl;
            SurvivingPixelPanel3.push_back(j);
            SurvivingPixelTotalAmpPanel3 += binContent;
            hcam_panel3->SetBinContent(nx + 1, ny + 1, binContent);
            if (std::find(updateCore.begin(), updateCore.end(), j) != updateCore.end()){
                totalCoreAmp += binContent;
            }
        } else {
            if (j == MaxAmplitdePixelID || j == SecondCoreID) {
                if (verbose) {
                    std::cout << "Core pixel " << j << " is not in the panel-3 cluster (amplitude " << binContent
                              << ") -> MaxAmplitude and SecondCoreAmp set to 0, event will be removed" << std::endl;
                }
                MaxEventAmplitde = 0.0;
                SecondCoreAmp = 0.0;
            }
            hcam_panel3->SetBinContent(nx + 1, ny + 1, 0);
            // TO DO:fix the core output int to reset the coreIDS to only include those plotted
        }
    }
    if (verbose) std::cout << "Total Core amp: " << totalCoreAmp << std::endl;
    // Check to see if there is a pixel that is connected to the triggered pixel that is at least PE_second > PE_Trig - 1.5sqrt(PE_Trig)
    // isSecondCorePixel(SurvivingPixelPanel3, hcam_panel3,MaxAmplitdePixelID, neighborDir);
    
    hcam_panel3->Draw("colz");
    IPlotTools::DrawMUSICBoundaries();

}


void CEvent::SetPanel3AllNeighboringPixels(TH2F* hcam_panel3, std::string neighborDir) {
    // std::vector<int> PixelPeakTimes = PeakTimeBin;
    // int TriggeredPixelPeakTime = PixelPeakTimes[MaxAmplitdePixelID];
    // std::cout << "Triggered Pixel Peak Time: " << TriggeredPixelPeakTime << std::endl;
    // Get neighboring pixels to triggered pixel
    std::vector<int> neighborsurvivingPanel3Pixels;
    std::vector<int> tempSurvivingPixelPanel3;

    // std::cout << "Neightsbors panel2:  ";
    // for (auto i: SurvivingPixelPanel2)
    //     std::cout << i;
    // std::cout << std::endl;
    
    TriggeredPixelNeighborhoodChecker(SurvivingPixelPanel2,MaxAmplitdePixelID,neighborsurvivingPanel3Pixels, neighborDir);
    
    // std::cout << "Neightsbors before before: ";
    // for (auto i: neighborsurvivingPanel3Pixels)
    //     std::cout << i;
    // std::cout << std::endl;
    
    
    for(int j = 0; j < MaxNofChannels; j++){
        int nx, ny;
        IPlotTools::FindBin(j, &nx, &ny);
        if ((PeakTimeBin[j] == 239 || PeakTimeBin[j] == 238 || PeakTimeBin[j] == 237) && std::find(std::begin(neighborsurvivingPanel3Pixels), std::end(neighborsurvivingPanel3Pixels), j) != std::end(neighborsurvivingPanel3Pixels)) {
            // std::cout << "Surviving Panel 3 Pixel: " << j << " Amplitude: " << AmplitudeValuesTimeBin[j] << std::endl;
            tempSurvivingPixelPanel3.push_back(j);
            //hcam_panel3->SetBinContent(nx + 1, ny + 1, AmplitudeValuesTimeBin[j]);
        } //else {
            // hcam_panel3->SetBinContent(nx + 1, ny + 1, 0);
            //}
    }
    
    // std::cout << "Neightsbors before 1: ";
    // for (auto i: tempSurvivingPixelPanel3)
    //     std::cout << i;
    // std::cout << std::endl;
    neighborsurvivingPanel3Pixels={};
    TriggeredPixelNeighborhoodChecker(tempSurvivingPixelPanel3,MaxAmplitdePixelID,neighborsurvivingPanel3Pixels, neighborDir);
    tempSurvivingPixelPanel3 = {};
    // Remove elements from neighbors that are not in surviving2
    for (int p = 0; p < neighborsurvivingPanel3Pixels.size(); p++) {
        // std::cout << "Pixel " << neighborsurvivingPanel3Pixels[p] << std::endl;
        for (int s=0; s<SurvivingPixelPanel2.size(); s++) {
            if (neighborsurvivingPanel3Pixels[p] == SurvivingPixelPanel2[s]) {
                // std::cout << "Found surviving pixel: " << neighborsurvivingPanel3Pixels[p] << std::endl;
                tempSurvivingPixelPanel3.push_back(neighborsurvivingPanel3Pixels[p]);
            }
            
        }
        
    }
    neighborsurvivingPanel3Pixels = tempSurvivingPixelPanel3;
    // std::cout << "Neightsbors after 1: ";
    // for (auto i: neighborsurvivingPanel3Pixels)
    //     std::cout << i;
    // std::cout << std::endl;
    for(int j = 0; j < MaxNofChannels; j++){
        int nx, ny;
        IPlotTools::FindBin(j, &nx, &ny);
        float binContent = hcam_panel3->GetBinContent(nx+1,ny+1);
        if ((PeakTimeBin[j] ==237 || PeakTimeBin[j] ==238 || PeakTimeBin[j] == 239) && std::find(std::begin(neighborsurvivingPanel3Pixels), std::end(neighborsurvivingPanel3Pixels), j) != std::end(neighborsurvivingPanel3Pixels)) {
            std::cout << "Surviving Panel 3 Pixel: " << j << " Amplitude: " << AmplitudeValuesTimeBin[j] << std::endl;
            SurvivingPixelPanel3.push_back(j);
            SurvivingPixelTotalAmpPanel3 += binContent;
            hcam_panel3->SetBinContent(nx + 1, ny + 1, binContent);
        } else {
            hcam_panel3->SetBinContent(nx + 1, ny + 1, 0);
        }
    }
    
    // Check to see if there is a pixel that is connected to the triggered pixel that is at least PE_second > PE_Trig - 1.5sqrt(PE_Trig)
    // isSecondCorePixel(hcam_panel3,MaxAmplitdePixelID, neighborDir, NumberOfCores);
    
    hcam_panel3->Draw("colz");
    IPlotTools::DrawMUSICBoundaries();

}

void CEvent::NeighborhoodCheckerHelper(int ID, std::vector<int>& IDSelected,std::vector<int>& cleaned_pixels, std::vector<int>& visited_ids, std::string neighborDir) {
    bool already_searched = std::find(visited_ids.begin(), visited_ids.end(), ID) != visited_ids.end();
    visited_ids.push_back(ID);
    if (already_searched == false) {
        for(std::vector<int>::size_type i=0; i < IDSelected.size(); i++){
            //cout << "ID: " << ID << " IDSelected: " << IDSelected[i] << endl;

            if (ID == IDSelected[i] ){
                cleaned_pixels.push_back(IDSelected[i]);

                // cout << "ID selected: " << ID << endl;
                // check recursilely the board side of these pixels
                // Removes the ID that has been seleceted in the array so that the loops end
                

                for (int arrayID: GetNeighborArray(ID,Form("%sneighborbroadside.csv",neighborDir.c_str()))) {
                    // cout << "ArrayID: NCH: " << arrayID << endl;
                    NeighborhoodCheckerHelper(arrayID,IDSelected, cleaned_pixels,visited_ids, neighborDir);
                }
            }

        }
    }
    //cout << ID <<", ";
}


void CEvent::TriggeredPixelNeighborhoodChecker(std::vector<int>& IDSelected, int triggeredpixel,std::vector<int>& cleaned_pixels, std::string neighborDir){
    // we already know that the triggered pixel has a high enough amplitude 
    // We have a list of amplitudes that survied the first part of the the cleaning the straight amplitude cut across the whole camera
    // that will have 256 pixels and if it is to low the pixel will be 0 otherwise it will have the correct amplitude value if it meets the amplitude requirement
    //IDSelected.erase(std::remove(IDSelected.begin(), IDSelected.end(), triggeredpixel), IDSelected.end());
    
    std::vector<int> visited_ids;
    // get the Trigger pixels neighbors
    cleaned_pixels.push_back(triggeredpixel);
    for (int arrayID : GetNeighborArray(triggeredpixel,Form("%sneighborbroadside.csv",neighborDir.c_str()))) {
        // cout << "arrayID: Triggered() " << arrayID << endl;
        NeighborhoodCheckerHelper(arrayID,IDSelected,cleaned_pixels,visited_ids, neighborDir);
        // cout << "completed ID: "<< endl;
    }



}

// Builds the panel-3 cluster by walking outward from the triggered pixel, carrying the time
// reference along with the walk. A core pixel (the triggered pixel or any of GetCoreIDs()) looks at
// all 8 of its neighbours, the same neighbor.csv the crosstalk step uses to pick the boundary pixels
// around each core, so a boundary pixel or second core sitting diagonally off a core can join. Every
// other pixel steps only to its 4 broadside neighbours (neighborbroadside.csv).
//
// A neighbour joins the cluster when all of these hold:
//   - it is not already in the cluster
//   - it still has amplitude in hcam_panel3 (it survived the amplitude and crosstalk cuts)
//   - its peak bin is within maxStep bins of the pixel it was reached from
//   - its peak bin is within maxDrift bins of the triggered pixel's peak bin (maxDrift < 0 = no cap)
//
// A pixel is only marked as visited once it joins. A pixel that fails the step test from one
// neighbour can still join through another neighbour, so the final cluster is the same whatever
// order the neighbours are searched in.
void CEvent::TimePropagatedNeighborhoodChecker(TH2F* hcam_panel3, int triggeredpixel, std::vector<int>& cleaned_pixels, std::string neighborDir, int maxStep, int maxDrift, bool verbose) {
    const int triggeredPeakBin = PeakTimeBin[triggeredpixel];
    const std::string broadsideNeighborFile = Form("%sneighborbroadside.csv", neighborDir.c_str());
    const std::string allNeighborFile       = Form("%sneighbor.csv", neighborDir.c_str());

    std::vector<bool> isCore(MaxNofChannels, false);
    isCore[triggeredpixel] = true;
    for (int coreID : GetCoreIDs()) {
        if (coreID >= 0 && coreID < MaxNofChannels) isCore[coreID] = true;
    }

    std::vector<bool> inCluster(MaxNofChannels, false);
    std::vector<bool> rejectedOnTime(MaxNofChannels, false);
    std::queue<int> pixelsToSearchFrom;

    cleaned_pixels.push_back(triggeredpixel);
    inCluster[triggeredpixel] = true;
    pixelsToSearchFrom.push(triggeredpixel);

    if (verbose) {
        std::cout << "Time walk (step +/-" << maxStep << ", max drift "
                  << (maxDrift < 0 ? std::string("off") : std::to_string(maxDrift))
                  << ") starting at pixel " << triggeredpixel << ", bin " << triggeredPeakBin << ":" << std::endl;
    }

    while (!pixelsToSearchFrom.empty()) {
        int parentPixelID = pixelsToSearchFrom.front();
        pixelsToSearchFrom.pop();
        const std::string& neighborFile = isCore[parentPixelID] ? allNeighborFile : broadsideNeighborFile;

        for (int neighborID : GetNeighborArray(parentPixelID, neighborFile)) {
            if (neighborID < 0 || neighborID >= MaxNofChannels) continue;

            int nx, ny;
            IPlotTools::FindBin(neighborID, &nx, &ny);
            int stepFromParent   = PeakTimeBin[neighborID] - PeakTimeBin[parentPixelID];
            int driftFromTrigger = PeakTimeBin[neighborID] - triggeredPeakBin;

            bool alreadyVisited      = inCluster[neighborID];
            bool passedAmplitudeCuts = hcam_panel3->GetBinContent(nx + 1, ny + 1) > 0.0;
            bool withinStep          = std::abs(stepFromParent) <= maxStep;
            bool withinDrift         = (maxDrift < 0) || (std::abs(driftFromTrigger) <= maxDrift);

            // Pixels already removed by the amplitude/crosstalk cuts are not a timing decision,
            // so they are skipped silently to keep the printout about timing only.
            if (alreadyVisited || !passedAmplitudeCuts) continue;

            if (withinStep && withinDrift) {
                cleaned_pixels.push_back(neighborID);
                inCluster[neighborID] = true;
                pixelsToSearchFrom.push(neighborID);
            } else {
                rejectedOnTime[neighborID] = true;
            }

            if (verbose) {
                std::cout << "  " << (withinStep && withinDrift ? "admitted" : "rejected")
                          << "  pixel " << neighborID
                          << "  bin " << PeakTimeBin[neighborID]
                          << "  delta " << std::showpos << driftFromTrigger
                          << "  step " << stepFromParent << std::noshowpos
                          << "  from " << (isCore[parentPixelID] ? "core " : "parent ") << parentPixelID;
                if (!withinStep)  std::cout << "  (step > " << maxStep << ")";
                if (!withinDrift) std::cout << "  (drift > " << maxDrift << ")";
                std::cout << std::endl;
            }
        }
    }

    if (verbose) {
        // a pixel rejected from one parent but admitted through another is not counted as rejected
        int rejectedCount = 0;
        for (int j = 0; j < MaxNofChannels; j++) {
            if (rejectedOnTime[j] && !inCluster[j]) rejectedCount++;
        }
        std::cout << "Time walk result: " << cleaned_pixels.size() << " pixels admitted (incl. triggered pixel), "
                  << rejectedCount << " rejected on time" << std::endl;
    }
}

// Diagnostic only, affects no cut: the cluster that would survive if no time steps were allowed at
// all, i.e. only pixels whose time bin is exactly the triggered pixel's. PeakTimeBin holds whatever
// clock the time walk runs on (the 50% rise bin as EventCleaning fills it), so this says how much of
// what the walk keeps it owes to the +/-maxStep allowance rather than to pixels genuinely in step.
int CEvent::CountSameTimeBinSurvivors(TH2F* hcam_panel3, std::string neighborDir, int& pixelsSharingBin) {
    const int triggeredTimeBin = PeakTimeBin[MaxAmplitdePixelID];

    std::vector<int> sameTimeBin;
    for (int j = 0; j < MaxNofChannels; j++) {
        int nx, ny;
        IPlotTools::FindBin(j, &nx, &ny);
        bool hasAmplitude = hcam_panel3->GetBinContent(nx + 1, ny + 1) > 0.0;
        if (hasAmplitude && PeakTimeBin[j] == triggeredTimeBin) sameTimeBin.push_back(j);
    }
    pixelsSharingBin = (int)sameTimeBin.size();

    std::vector<int> connected;
    TriggeredPixelNeighborhoodChecker(sameTimeBin, MaxAmplitdePixelID, connected, neighborDir);

    // the checker always adds the triggered pixel (and can add it twice), so count unique pixels
    // that share the bin themselves
    std::set<int> kept;
    for (int pixelID : connected) {
        if (std::find(sameTimeBin.begin(), sameTimeBin.end(), pixelID) != sameTimeBin.end()) {
            kept.insert(pixelID);
        }
    }
    return (int)kept.size();
}

void CEvent::isSecondCorePixel(TH2F* hcam_panel3, int triggeredpixel,std::string neighborDir,float CoreCutOff, int NumberOfCores, bool verbose) {
    // gather all pixels that match the core requirement 
    std::vector<int> cores = {};
    int nx, ny;
    IPlotTools::FindBin(triggeredpixel, &nx, &ny);
    float TriggererContent = hcam_panel3->GetBinContent(nx+1,ny+1);
    if (verbose) std::cout << "Triggered pixel: " << triggeredpixel << " amplitude: " << TriggererContent << " PE" << std::endl;
    // A pixel is a core when its amplitude is within CoreCutOff*sqrt(trigger) of the triggered pixel, or
    // when it is railed at the snapshot bin. A railed pixel's amplitude is just its rail value after
    // calibration, which differs from pixel to pixel, so it cannot be compared by amplitude.
    const float coreAmplitudeCutoff = TriggererContent - CoreCutOff*sqrt(TriggererContent);
    if (verbose) std::cout << "Core amplitude cutoff: " << coreAmplitudeCutoff << " PE (trigger - "
                           << CoreCutOff << "*sqrt(trigger))" << std::endl;
    for(int j = 0; j< MaxNofChannels; j++) {
        IPlotTools::FindBin(j, &nx, &ny);
        float binContent = hcam_panel3->GetBinContent(nx+1,ny+1);
        bool aboveCoreCutoff     = binContent >= coreAmplitudeCutoff;
        bool saturatedAtSnapshot = std::find(SaturatedAtSnapshot.begin(), SaturatedAtSnapshot.end(), j) != SaturatedAtSnapshot.end();
        if (aboveCoreCutoff || saturatedAtSnapshot) {
            if (verbose) std::cout << " core: " << j << " Amplitude: " << binContent
                      << (aboveCoreCutoff ? "" : "  (below cutoff, core because saturated at snapshot bin)") << std::endl;
            cores.push_back(j);
        }
    }
    
    std::vector<int> neighborsurvivingPanel3Pixels={};
    // then ensure they are all attached to the triggered pixel via broadside and adjecent neighborhold to eachother as fasr out at core pixels needed 
    TriggeredPixelNeighborhoodChecker(cores,MaxAmplitdePixelID,neighborsurvivingPanel3Pixels, neighborDir);
    
    SetCoreIDs(cores);
    SetNumberofCorePixels(GetCoreIDs().size());
    sort(cores.begin(), cores.end());


    float secondcoreamp = -1;
    for (int j = 0; j < GetCoreIDs().size(); j++){
        if (GetCoreIDs()[j] == triggeredpixel){
            continue;
        }
        IPlotTools::FindBin(GetCoreIDs()[j], &nx, &ny);
        float binContent = hcam_panel3->GetBinContent(nx+1,ny+1);
        if (verbose) std::cout << " candidate second core: " << GetCoreIDs()[j] << " amplitude: " << binContent << " PE" << std::endl;
        // will set the first pixel that has a bin content as the highest core pixel if there are multiple at 160
        if (secondcoreamp < binContent) {
                secondcoreamp = binContent;
                SetSecondCoreID(GetCoreIDs()[j]);
                SetSecondCoreAmp(binContent);
        }
            
    }


    SetCoreRatio((1.0*TriggererContent)/(1.0*GetSecondCoreAmp()));
    if (verbose) {
        std::cout << "Second core pixel: " << GetSecondCoreID() << " amplitude: " << GetSecondCoreAmp()
                  << " PE, core ratio " << GetCoreRatio() << ", cores found " << GetCoreIDs().size()
                  << " (need " << NumberOfCores << ")" << std::endl;
    }
    if (GetCoreIDs().size() < NumberOfCores || TriggererContent == 0) {
        if (verbose) std::cout << "No second core pixel: only the triggered pixel is kept, so the event"
                                  " cannot pass the surviving-pixel cut" << std::endl;
        SurvivingPixelPanel3={triggeredpixel};
        SetSecondCoreID(-1);
        SetSecondCoreAmp(-1);
    }
    
    // Get  setsecond core through second max pixel

    
    // // check all the neighbors of the triggered pixel to see if there is another pixel that is within the range of the triggered pxiel to be a second core pixel
    // int nx, ny;
    // IPlotTools::FindBin(triggeredpixel, &nx, &ny);
    // float TriggererContent = hcam_panel3->GetBinContent(nx+1,ny+1);
    // // std::cout << "triggered pXiel: " << triggeredpixel << " Amp:  " << TriggererContent << std::endl;
    // int scount = 0;
    // float secondcoreamp = 0;
    
    // for (int arrayID : GetNeighborArray(triggeredpixel,Form("%sneighbor.csv",neighborDir.c_str()))) {
    //     // std::cout << "arrayID: Second() " << arrayID << std::endl;
    //     int nx, ny;
    //     IPlotTools::FindBin(arrayID, &nx, &ny);
    //     float binContent = hcam_panel3->GetBinContent(nx+1,ny+1);
    //     // std::cout << "binContent: " << binContent <<std::endl;
    //     if (binContent >= (TriggererContent - CoreCutOff*sqrt(TriggererContent))) {
    //         std::cout << "Second core: " << arrayID << " Amplitude: " << binContent << std::endl;
    //         scount++;
    //         if (secondcoreamp < binContent) {
    //             secondcoreamp = binContent;
    //             SetSecondCoreID(arrayID);
    //             SetSecondCoreAmp(binContent);
    //         }
    //     }

    // }
    // SetNumberofCorePixels(scount+1);
    // if (NumberOfCores == 0 && scount == 0) {
    //     SetSecondCoreID(MaxAmplitdePixelID);
    //     SetSecondCoreAmp(TriggererContent);       
    // }

    // std::cout << "scount: " << scount << std::endl;
    // if (scount < NumberOfCores || TriggererContent == 0) {
    //     std::cout << "Second core pixel was not found: "<< std::endl;
    //     SurvivingPixelPanel3={triggeredpixel};
    //     SetSecondCoreID(-1);
    //     SetSecondCoreAmp(-1);
    // }

    // SetCoreRatio(TriggererContent/secondcoreamp);
}

// Crosstalk step before the time walk. A pixel is kept when its amplitude is at least cutoff AND it
// is the triggered pixel, one of the core pixels, or a neighbour (neighbor.csv, incl. diagonals) of
// the triggered pixel or of a core pixel. Every other pixel is set to 0 in hcam_panel3.
//
// The core pixels are kept explicitly because a pixel is never in its own neighbour list: a second
// core that is not next to the triggered pixel used to be zeroed here, which then zeroed both core
// amplitudes in SetPanel3 and removed the event.
void CEvent::CheckCloseNeigborsOnly(TH2F* hcam_panel3, int pixel,std::vector<int> cores, double cutoff,std::string neighborDir, bool verbose) {
    const std::string neighborFile = Form("%sneighbor.csv", neighborDir.c_str());

    // pixels allowed to keep their amplitude: triggered pixel, cores, and all their neighbours
    std::set<int> allowedPixels = {pixel};
    for (int neighborID : GetNeighborArray(pixel, neighborFile)) allowedPixels.insert(neighborID);
    for (int coreID : cores) {
        allowedPixels.insert(coreID);
        for (int neighborID : GetNeighborArray(coreID, neighborFile)) allowedPixels.insert(neighborID);
    }

    std::vector<int> removedNotNextToCore;
    for(int j = 0; j < MaxNofChannels; j++){
        int nx, ny;
        IPlotTools::FindBin(j, &nx, &ny);
        float binContent = hcam_panel3->GetBinContent(nx+1,ny+1);

        bool aboveCutoff       = binContent >= cutoff;
        bool nearTriggerOrCore = allowedPixels.count(j) > 0;

        if (aboveCutoff && nearTriggerOrCore) {
            hcam_panel3->SetBinContent(nx + 1, ny + 1, binContent);
        } else {
            if (aboveCutoff) removedNotNextToCore.push_back(j);
            hcam_panel3->SetBinContent(nx + 1, ny + 1, 0);
        }
    }

    if (verbose) {
        std::cout << "Crosstalk step (cutoff " << cutoff << " PE), core pixels kept:";
        for (int coreID : cores) std::cout << " " << coreID;
        std::cout << std::endl;
        std::cout << "Crosstalk step removed above cutoff but not next to triggered/core pixel:";
        if (removedNotNextToCore.empty()) std::cout << " none";
        for (int pixelID : removedNotNextToCore) std::cout << " " << pixelID;
        std::cout << std::endl;
    }
}


// Function to open a lookup table for the neighbors of each pixel ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int> CEvent::GetNeighborArray(int id, std::string filename) {

	// give name of lookup table file
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open the neighbor.csv file for reading." << std::endl;
        return std::vector<int>();
    }

    std::string line;
	// load each line from file
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        char delimiter = ',';

        int neighborID;
        std::string neighborArrayStr;

        if (iss >> neighborID >> delimiter >> neighborArrayStr) {
            if (neighborID == id) {
                std::vector<int> neighborArray;

                // Handle the first value separately
                size_t pos = 0;
                int value;

                try {
                    value = std::stoi(neighborArrayStr, &pos);
                    neighborArray.push_back(value);
                } catch (const std::invalid_argument& e) {
                 //   std::cerr << "Invalid argument: " << e.what() << std::endl;
                } catch (const std::out_of_range& e) {
                    std::cerr << "Out of range: " << e.what() << std::endl;
                }

                // Continue with the rest of the values
                while (pos < neighborArrayStr.length()) {
                    size_t nextPos;
                    try {
                        value = std::stoi(neighborArrayStr.substr(pos + 1), &nextPos);
                        neighborArray.push_back(value);
                        pos += nextPos + 1;
                    } catch (const std::invalid_argument& e) {
                     //   std::cerr << "Invalid argument: " << e.what() << std::endl;
                        break;
                    } catch (const std::out_of_range& e) {
                        std::cerr << "Out of range: " << e.what() << std::endl;
                        break;
                    }
                }

                file.close();
                return neighborArray;
            }
        }
    }

    file.close();
    return std::vector<int>(); // returns array of neighbors for given pixel
}

int CEvent::isConfigureEvent(){
    // std::cout << "Pixels: " << SurvivingPixelPanel3.size() <<" And size: " << SurvivingPixelTotalAmpPanel3 << std::endl;
    if (SurvivingPixelPanel3.size() < 15 || SurvivingPixelTotalAmpPanel3 < 2080){
        // std::cout << "Exitting" << std::endl;
        return 0;
    }
    std::vector<int> survivingSIABS; 
    for(unsigned int i=0; i < SurvivingPixelPanel3.size(); i++){
        int nx, ny;
        IPlotTools::FindBin(SurvivingPixelPanel3[i],&nx,&ny);
        survivingSIABS.push_back(IPlotTools::FindSIAB(nx,ny));
    }

    std::sort(survivingSIABS.begin(), survivingSIABS.end());
    auto unique=std::unique(survivingSIABS.begin(), survivingSIABS.end());
    survivingSIABS.erase(unique, survivingSIABS.end());

    if (survivingSIABS.size() > 1){
        return 0;
    }
    return 1;

}



