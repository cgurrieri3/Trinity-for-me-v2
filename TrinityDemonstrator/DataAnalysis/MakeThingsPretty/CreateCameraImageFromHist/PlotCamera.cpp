#include "PlotCamera.h"

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

TH2F* makeSkyTopHistogram(const TH2F* source){
    TH2F *flipped = (TH2F*)source->Clone(Form("%s_sky_top", source->GetName()));
    flipped->SetDirectory(0);
    flipped->Reset("ICES");

    // Flip the camera image in both axes so the orientation matches the desired
    // camera view while keeping the original histogram unchanged.
    //
    // Each bin is copied to its mirror position in both x and y:
    // old left column -> new right column and old bottom row -> new top row.
    // The same mirrored copy is also applied to bin errors so plotted
    // uncertainties stay with their corresponding pixels.
    const int nBinsX = source->GetNbinsX();
    const int nBinsY = source->GetNbinsY();
    for (int xBin = 0; xBin <= nBinsX + 1; xBin++){
        for (int yBin = 0; yBin <= nBinsY + 1; yBin++){
            const int flippedXBin = nBinsX + 1 - xBin;
            const int flippedYBin = nBinsY + 1 - yBin;
            flipped->SetBinContent(flippedXBin, flippedYBin, source->GetBinContent(xBin, yBin));
            flipped->SetBinError(flippedXBin, flippedYBin, source->GetBinError(xBin, yBin));
        }
    }

    flipped->SetEntries(source->GetEntries());
    return flipped;
}

int main(int argc, char **argv){
    if(argc < 3){
        cout << "Too few arguments; please include the date data directory and mount flag" << endl;
        return 1;
    }

    std::string FilePath = argv[1];
    std::string mount = argv[2];

    if (mount == "y"){
        std::cout << "using mounted directory path" << std::endl;
        mnt = "/mnt/";
        outDir = "/mnt/DataAnalysis/MakeThingsPretty/Images/";
    } else if (mount == "n"){
        std::cout << "using default directory path" << std::endl;
        mnt = "";
    } else {
        std::cout << "using specific directory path: " << mount << std::endl;
        mnt = mount.c_str();
        outDir = Form("%sDataAnalysis/MakeThingsPretty/Images2/",mnt.c_str());
    }
    
    std::vector<string> tempfiles = util->readFileToVectorString(Form("%s%s", mnt.c_str(), FilePath.c_str()));
    cout << tempfiles[0] << endl;

    auto files = splitCSVLines(tempfiles);

    for (const auto& row : files) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }

    canvas = new TCanvas("Display","Event Display",1000,900);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);

    // Int_t ci = TColor::GetFreeColorIndex();
    // TColor *color = new TColor(ci, 238.f/255.f, 243.f/255.f, 246.f/255.f);
    // Color_t thecolor = color->GetNumber();

    // canvas->SetFillColor(thecolor);

    for(int i = 1; i < (int)files.size(); i++){
        FilePath = Form("%s",(files[i][3]).c_str());
        cout << "file path: " << FilePath << endl;

        canvas->cd();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);

        TFile *f = new TFile(FilePath.c_str());
        if (!f || f->IsZombie()){
            cerr << "Error opening file: " << FilePath << endl;
            continue;
        }

        std::string HistTitle = files[i][2];
        std::replace(HistTitle.begin(), HistTitle.end(), '%', ' ');
        cout << HistTitle << endl;

        TH2F *sourceHist = (TH2F*)f->Get(HistTitle.c_str());
        if (!sourceHist){
            cerr << "Error: histogram not found: " << HistTitle << endl;
            f->Close();
            continue;
        }

        TH2F *h2 = makeSkyTopHistogram(sourceHist);

        if(files[i][4] == "PE"){
            h2->SetTitle("Temp;Pixels along camera x axis;Pixels along camera y axis;PE");
        } else {
            h2->SetTitle("Temp;Pixels along camera x axis;Pixels along camera y axis;Photons");
        }

        std::string Title = files[i][1];
        std::replace(Title.begin(), Title.end(), '_', ' ');
        h2->SetTitle("");
        h2->GetZaxis()->SetTitleOffset(1.4);

        TFile *outFile = new TFile(Form("%s%s.root", outDir.c_str(), (files[i][0]).c_str()), "RECREATE");
        canvas->Print(Form("%s%s.pdf[", outDir.c_str(), (files[i][0]).c_str()));

        h2->Draw("COLZ");
        canvas->Update();

        TPaletteAxis* palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
        if (palette) {
            palette->SetX1NDC(0.86);
            palette->SetX2NDC(0.90);
            palette->SetY1NDC(0.1);
            palette->SetY2NDC(0.9);
            canvas->Modified();
            canvas->Update();
        }

        const double yMin = h2->GetYaxis()->GetXmin();
        // Sofia's commissioning-paper note gives the horizon by camera-row
        // ratios, not by an exact pointing-change day:
        // - Oct. 2024 through Sept. 2025: 12 ground rows and 4 sky rows.
        // - Oct. 2025 through Dec. 2025: 9 ground rows and 7 sky rows.
        //
        // With 16 y bins spanning -0.5 to 15.5, the row boundary is:
        // - 12:4 ground:sky -> line at y=11.5, between rows 11 and 12.
        // -  9:7 ground:sky -> line at y=8.5, between rows 8 and 9.
        //
        // If Sofia later gives a different first day for the 9:7 period, only
        // update this constant.
        const int newHorizonRatioDate = 20251001;
        const int eventDate = extractYYYYMMDD(HistTitle + " " + FilePath);
        const double oldGroundSkyBoundary = yMin + 12.0;
        const double newGroundSkyBoundary = yMin + 9.0;
        const double skyGroundLineY = (eventDate >= newHorizonRatioDate)
            ? newGroundSkyBoundary
            : oldGroundSkyBoundary;

        TLine *l = new TLine(-0.5, skyGroundLineY, 15.5, skyGroundLineY);
        l->SetLineWidth(6);
        l->SetLineColor(kGray);

        TText *t1 = new TText(13.1, skyGroundLineY + 0.75, "Sky");
        t1->SetTextAlign(22);
        t1->SetTextSize(0.05);
        t1->SetTextColor(kGray);

        TText *t2 = new TText(13.1, skyGroundLineY - 0.6, "Ground");
        t2->SetTextAlign(22);
        t2->SetTextSize(0.05);
        t2->SetTextColor(kGray);

        if (files[i][0].find("Cleaned") != std::string::npos || files[i][0].find("Sim") != std::string::npos) {
            l->SetLineColor(kBlack);
            t1->SetTextColor(kBlack);
            t2->SetTextColor(kBlack);
        }

        if (files[i][0].find("Sim") == std::string::npos){
            l->Draw();
            t1->Draw();
            t2->Draw();
        }

        plottools->DrawSIABBoundaries();
        canvas->Update();

        outFile->cd();
        h2->Write("hist");
        canvas->Write(h2->GetTitle());

        canvas->Print(Form("%s%s.pdf",  outDir.c_str(), (files[i][0]).c_str()));
        canvas->SaveAs(Form("%s%s.png", outDir.c_str(), (files[i][0]).c_str()));
        canvas->Print(Form("%s%s.pdf]", outDir.c_str(), (files[i][0]).c_str()));

        f->Close();
        outFile->Close();
        delete h2;
    }

    cout << "Completed!" << endl;
}
