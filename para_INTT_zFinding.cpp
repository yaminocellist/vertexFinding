#include "Fit/BinData.h"
#include <TStopwatch.h>
#include "TApplication.h"
#include <thread>
#include <chrono>
#include <mutex>

#include "headers/histogramPlotting.h"
#include "headers/CommonParaFuncs.h"
#include "headers/zFinding.h"

std::mutex m_mutex;
TH1D *h;
using namespace rootBranchEssentials;

void para_histoAccumulator (
    const int &lower_bound, 
    const int &upper_bound, 
    const std::vector<myTrackletMemberLite> &t0, 
    const std::vector<myTrackletMemberLite> &t1
) {
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        // std::cout << "Processed range: " << lower_bound << " to " << upper_bound << "," << event << "," << NClus << std::endl;
        int sizeMax = t1.size();
        double dPhi, foundZ_revised, nearest_d_revisedSQUARED;
        std::pair<double, double> nearestZResults_revised;
        for (int i = lower_bound; i < upper_bound; i++) {
            for (int j = 0; j < sizeMax; j++) {
                dPhi = t0[i].phi - t1[j].phi;
                if (dPhi >= -0.1 && dPhi <= +0.1) {
                    myPoint3D P1 = {t0[i].x, t0[i].y, t0[i].z};
                    myPoint3D P2 = {t1[j].x, t1[j].y, t1[j].z};
                    nearestZResults_revised         = nearestZ_revised(P1, P2);
                    foundZ_revised           = nearestZResults_revised.first;
                    nearest_d_revisedSQUARED = nearestZResults_revised.second;
                    if (foundZ_revised >= zmin && foundZ_revised <= zmax && nearest_d_revisedSQUARED <= DCA_cutSQUARED) {
                        h -> Fill(foundZ_revised);
                    }
                }
            }
        }
    }
}

double singleFinding(const int &target) {
    int sizeOfLayer0;   double foundZ;
    std::vector<int> threadRanges;
    // b_event->GetEntry(target);
    // b_MBD_centrality->GetEntry(target);
    b_NClus->GetEntry(target);
    b_ClusLayer->GetEntry(target);
    b_ClusX->GetEntry(target);
    b_ClusY->GetEntry(target);
    b_ClusZ->GetEntry(target);
    // b_ClusR->GetEntry(target);
    b_ClusPhi->GetEntry(target);
    // b_ClusEta->GetEntry(target);
    std::vector<myTrackletMemberLite> tracklet_layer_0, tracklet_layer_1;
        for (int i = 0; i < ClusX->size(); i++) {
            if (ClusLayer->at(i) == 3 || ClusLayer->at(i) == 4) {
                // std::cout << std::atan2(ClusY->at(i), ClusX->at(i)) - ClusPhi->at(i) << std::endl;
                for (int k = 0; k <= 32; k++) {
                    tracklet_layer_0.push_back({ClusX->at(i), ClusY->at(i), ClusZ->at(i) + k*0.05,
                                                ClusPhi->at(i)});
                }
            }
            else {
                for (int k = 0; k <= 32; k++) {
                    tracklet_layer_1.push_back({ClusX->at(i), ClusY->at(i), ClusZ->at(i) + k*0.05,
                                                ClusPhi->at(i)});
                }
            }
        }
        sizeOfLayer0 = tracklet_layer_0.size();
        threadRanges.push_back(0);
        threadRanges.push_back(sizeOfLayer0/8);
        threadRanges.push_back(sizeOfLayer0/4);
        threadRanges.push_back(sizeOfLayer0/8*3);
        threadRanges.push_back(sizeOfLayer0/2);
        threadRanges.push_back(sizeOfLayer0/8*5);
        threadRanges.push_back(sizeOfLayer0/8*6);
        threadRanges.push_back(sizeOfLayer0/8*7);
        threadRanges.push_back(sizeOfLayer0);
        // Use 'sysctl -n hw.logicalcpu' to determine the max number of threads:
        h = new TH1D("", "", bins, zmin - scanstep/2, zmax + scanstep/2);
        h->Sumw2(kFALSE);
        std::thread thsafe[8];
        for(int i = 0; i < 8; ++i)     thsafe[i] = std::thread(para_histoAccumulator,threadRanges[i],threadRanges[i+1],tracklet_layer_0,tracklet_layer_1);
        for(int i = 0; i < 8; ++i)     thsafe[i].join();

        TSpectrum *s = new TSpectrum();
        TH1 *bg = s -> Background(h, 20 , "nosmoothing"); 
        bg -> SetLineColor(kRed);
        h->Add(bg, -1.);
        Int_t nfound = s->Search(h, 7., "", 0.001); // total number of xpeaks
        Double_t *xpeaks = s->GetPositionX();
        foundZ = xpeaks[0];
        // foundZ = (xpeaks[0] > xpeaks[1]) ? xpeaks[0] : xpeaks[1];
        // printRed(xpeaks[0], xpeaks[1]);

    tracklet_layer_0.clear();   tracklet_layer_1.clear();
    h -> Reset("ICESM");    delete h;   h = 0;
    // TSpectrumOutcomes(h, bg, target, NClus, foundZ, MBD_z_vtx);
    // printRed(target, NClus, foundZ, MBD_z_vtx);
    return foundZ;
}

int main(int argc, char* argv[]) {
    // Start the stopwatch:
    TStopwatch timer;   timer.Start();
    // Get and print the current PC time
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    std::cout << "Run started at: " << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S") << std::endl;
    TApplication theApp("App", &argc, argv);

    std::vector<std::string> method;
    // <METHOD> <TARGET> <CEN_LOW> <CEN_HIGH> <Z_LOW> <Z_HIGH>
    int target       = 100;
    double phi_low   = 0.0;
    double phi_high  = 0.7;
    double z_low     = 25.;
    double z_high    = 5.; 

    if (argc > 1) {
        if (argc > 2) target   = std::stoi(argv[2]);
        if (argc > 3) phi_low  = std::stod(argv[3]);
        if (argc > 4) phi_high = std::stod(argv[4]);
        if (argc > 5) z_low    = std::stod(argv[5]);
        if (argc > 6) z_high   = std::stod(argv[6]);
        for (int i = 1; i < argc; i++) {
            std::string substring = argv[i];
            method.push_back(substring);
        }
    }
    std::vector<std::string> ops = splitString(method[0], '_');

    TFile *file = TFile::Open("../../External/Data_CombinedNtuple_Run54280_20241113.root");
    objectExistenceCheck(file);
    TTree *tree = (TTree*)file->Get("EventTree");
    objectExistenceCheck(tree);

    TObjArray* branches = tree->GetListOfBranches();
    std::string branchName;
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch* branch = (TBranch*)branches->At(i);
        branchName = branch->GetName();
        if (branchName=="event")        idx_event = i;
        if (branchName=="MBD_z_vtx")    idx_MBD_z_vtx = i;
        if (branchName=="MBD_centrality")    idx_MBD_centrality = i;
        if (branchName=="MBD_charge_sum")    idx_MBD_charge_sum = i;
        if (branchName=="NClus")        idx_NClus = i;
        if (branchName=="ClusLayer")    idx_ClusLayer = i;
        if (branchName=="ClusX")        idx_ClusX = i;
        if (branchName=="ClusY")        idx_ClusY = i;
        if (branchName=="ClusZ")        idx_ClusZ = i;
        if (branchName=="ClusR")        idx_ClusR = i;
        if (branchName=="ClusPhi")      idx_ClusPhi = i;
        if (branchName=="ClusEta")      idx_ClusEta = i;
    }
    b_event     = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_event));
    b_MBD_z_vtx = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_MBD_z_vtx));
    b_MBD_centrality = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_MBD_centrality));
    b_MBD_charge_sum = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_MBD_charge_sum));
    b_NClus     = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_NClus));
    b_ClusLayer = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_ClusLayer));
    b_ClusX     = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_ClusX));
    b_ClusY     = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_ClusY));
    b_ClusZ     = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_ClusZ));
    b_ClusR     = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_ClusR));
    b_ClusPhi   = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_ClusPhi));
    b_ClusEta   = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_ClusEta));
    b_event    ->SetAddress(&event);
    b_MBD_z_vtx->SetAddress(&MBD_z_vtx);
    b_MBD_centrality->SetAddress(&MBD_centrality);
    b_MBD_charge_sum->SetAddress(&MBD_charge_sum);
    b_NClus    ->SetAddress(&NClus);
    b_ClusLayer->SetAddress(&ClusLayer);
    b_ClusX    ->SetAddress(&ClusX);
    b_ClusY    ->SetAddress(&ClusY);
    b_ClusZ    ->SetAddress(&ClusZ);
    b_ClusR    ->SetAddress(&ClusR);
    b_ClusPhi  ->SetAddress(&ClusPhi);
    b_ClusEta  ->SetAddress(&ClusEta);

    ROOT::EnableThreadSafety();
    double results;
    if (method[0] == "single") {
        b_MBD_z_vtx->GetEntry(target);
        results = singleFinding(std::stoi(method[1]));
    }
    if (method[0] == "all") {
        gROOT->SetBatch(kTRUE); gErrorIgnoreLevel = kWarning;   // Supressing TCanvas auto creation;
        // std::string filePath = "../zFindingResults/DCAfit_0_16_-001_001_step2e-1_z_-40_30.txt";
        std::string filePath = "../zFindingResults/DummyZ.txt";
        std::ofstream outputFile(filePath, std::ios_base::app);
        if (!outputFile.is_open()) {
            std::cout << "Unable to open the file to be written." << std::endl;
            system("read -n 1 -s -p \"Press any key to continue...\" echo");
            exit(1);
        }
        std::ifstream checkFile(filePath, std::ios::ate); // Check if the text file is empty;
        if (!checkFile.is_open()) {
            std::cerr << "Failed to open file for reading.\n";
            system("read -n 1 -s -p \"Press any key to continue...\" echo");
            exit(1);
        }
        bool isEmpty = (checkFile.tellg() == 0);
        if(isEmpty) outputFile << "index,event,NClus,foundZ,trueZ,centrality" << std::endl;

        Long64_t nEntries = tree -> GetEntries();
        for (Int_t i = 0; i < nEntries; i++) {
            b_MBD_z_vtx->GetEntry(i);
            b_event    ->GetEntry(i);
            if (MBD_z_vtx > -40 && MBD_z_vtx < 30) {
                results = singleFinding(i);
                outputFile << i << "," << event << "," << NClus << "," << results << "," << MBD_z_vtx << std::endl;
            }
        }
        printBlue("Finished!");
    }
    if (method[0] == "centrality") {
        // ↓↓ foundZ data access;
        std::string filePath = "../zFindingResults/DCAfit_0_16_-001_001_step2e-1_z_-40_30.txt";
        std::ifstream myfile(filePath);
        fileExistenceCheck(myfile);
        std::string line, value;
        std::vector<int> index, evt, NHits;
        std::vector<double> foundZ, MBD_true_z;

        getline(myfile, line);
        while (getline(myfile, line)) {
            std::stringstream data(line);
            getline(data, value, ',');  int i     = std::stoi(value);
            getline(data, value, ',');  int e     = std::stoi(value);
            getline(data, value, ',');  int N     = std::stoi(value);
            getline(data, value, ',');  double f  = std::stod(value);
            getline(data, value, ',');  double Mz = std::stod(value);

            index.push_back(i);     evt.push_back(e);         NHits.push_back(N);
            foundZ.push_back(f);    MBD_true_z.push_back(Mz);   
        }
        // ↑↑ foundZ data access;

        TH1D *MBDChargeSum = new TH1D("", "", 501, 0, 2000);
        Long64_t nEntries = tree -> GetEntries();
        for (Int_t i = 0; i < nEntries; i++) {
            b_MBD_z_vtx->GetEntry(i);
            if (MBD_z_vtx > -40 && MBD_z_vtx < 30) {
                b_MBD_charge_sum->GetEntry(i);
                MBDChargeSum->Fill(MBD_charge_sum);
            }
        }
        // ↓↓ Write to a new file;
        std::string outFile = "../zFindingResults/DCAfit_0_16_-001_001_step2e-1_z_-40_30_centrality_calculated.txt";
        std::ofstream outputFile(outFile);
        if (!outputFile.is_open()) {
            std::cout << "Unable to open the file to be written." << std::endl;
            system("read -n 1 -s -p \"Press any key to continue...\" echo");
            exit(1);
        }
        std::ifstream checkFile(outFile, std::ios::ate); // Check if the text file is empty;
        if (!checkFile.is_open()) {
            std::cerr << "Failed to open file for reading.\n";
            system("read -n 1 -s -p \"Press any key to continue...\" echo");
            exit(1);
        }
        bool isEmpty = (checkFile.tellg() == 0);
        if(isEmpty) outputFile << "index,event,NClus,foundZ,MBDZ,calculated_centrality,MBD_charge_sum" << std::endl;
        // ↑↑ Write to a new file;
        TCanvas *c = new TCanvas("c", "c", 0,50,2000,1000); // It's a must!
        MBDChargeSum->Draw();
        gPad -> SetLogy();
        gPad -> Update();
        int minBin = MBDChargeSum->FindFixBin(0);
        int maxBin = MBDChargeSum->FindFixBin(1900);
        double sum = MBDChargeSum->Integral(minBin, maxBin, "");
        printBlue(sum);
        for (int i = 0; i < index.size(); i++) {
            b_MBD_charge_sum->GetEntry(index[i]);
            double MBD_calculated_centrality = (MBDChargeSum->Integral(MBDChargeSum->FindFixBin(MBD_charge_sum),maxBin,""))/sum;
            // printRed(index[i], evt[i], NHits[i], foundZ[i], MBD_true_z[i], MBD_charge_sum, MBD_calculated_centrality);
            outputFile << index[i] << "," << evt[i] << "," << NHits[i] << "," << foundZ[i] << ","
                       << MBD_true_z[i] << "," << MBD_calculated_centrality << "," << MBD_charge_sum << std::endl;
        }
        
        // for (int bin = minBin; bin <= maxBin; bin++) {
        //     double boundary0 = MBDChargeSum->Integral(bin, maxBin, "");
        //     if (boundary0/sum <= 0.05) {
        //         printRed(bin, MBDChargeSum->GetBinCenter(bin));
        //         break;
        //     }
        // }
    }
    if (method[0] == "analyze") {
        // ↓↓ foundZ data access;
        std::string filePath = "../zFindingResults/DCAfit_0_16_-001_001_step2e-1_z_-40_30_centrality_calculated.txt";
        std::ifstream myfile(filePath);
        fileExistenceCheck(myfile);
        std::string line, value;
        std::vector<int> index, evt, NHits;
        std::vector<double> foundZ, MBD_true_z, calculated_cen;

        getline(myfile, line);
        while (getline(myfile, line)) {
            std::stringstream data(line);
            getline(data, value, ',');  int i     = std::stoi(value);
            getline(data, value, ',');  int e     = std::stoi(value);
            getline(data, value, ',');  int N     = std::stoi(value);
            getline(data, value, ',');  double f  = std::stod(value);
            getline(data, value, ',');  double Mz = std::stod(value);
            getline(data, value, ',');  double cc = std::stod(value);

            index.push_back(i);     evt.push_back(e);           NHits.push_back(N);
            foundZ.push_back(f);    MBD_true_z.push_back(Mz);   calculated_cen.push_back(cc);
        }
        // ↑↑ foundZ data access;

        TCanvas *c = new TCanvas("c", "c", 0,50,2000,1000); // It's a must!
        TH1D *h_resolution = new TH1D("found z resolution", Form("INTT found z resolution of %lu events;dZ [cm];# of counts", index.size()), 1001, -1e2, 1e2);
        TH1D *h_foundZ = new TH1D("found z", ";Z vtx [cm];# of counts", 401, -1e2, 1e2);
        TH1D *h_MBDZ = new TH1D("MBD z", Form("INTT MBD z of %lu events;Z vtx [cm];# of counts", index.size()), 401, -1e2, 1e2);
        h_resolution -> GetXaxis() -> CenterTitle(true); h_resolution -> GetYaxis() -> CenterTitle(true);
        h_resolution -> SetFillStyle(3003); h_resolution -> SetFillColor(kBlue - 7);
        h_resolution -> SetLineWidth(2);    h_resolution->GetXaxis()->SetRangeUser(-20,20);
        h_foundZ -> GetXaxis() -> CenterTitle(true); h_foundZ -> GetYaxis() -> CenterTitle(true);
        h_foundZ -> SetLineWidth(3);    h_foundZ->GetXaxis()->SetRangeUser(-60,60); h_MBDZ -> SetLineWidth(3);
        TLine *l = new TLine(0, 0, 0, h_resolution->GetBinContent(h_resolution->GetMaximumBin()));
        l -> Draw("same"); l -> SetLineColor(kRed); l -> SetLineWidth(2);
        TH2D *h = new TH2D("", ";MBD Z;INTT found Z", 100, -250, -50, 100, -250, -50);
        TGraph *g_NFunction    = new TGraph();
        g_NFunction->SetTitle(Form("INTT found z resolution of %lu events;dZ [cm];# of counts", index.size()));
        g_NFunction->GetYaxis()->SetTitle("dZ [cm]"); g_NFunction->GetXaxis()->SetTitle("NClus");
        g_NFunction->GetXaxis()->CenterTitle(true); g_NFunction->GetYaxis()->CenterTitle(true);
        TGraph *g_ZFunction    = new TGraph();
        TGraph *g_CFunction    = new TGraph();
        TGraph *g_ZCorrelation = new TGraph();
        TGraph *g_Calibration  = new TGraph();
        double z_resolution;
        for (int i = 0; i < index.size(); i++) {
            h -> Fill(MBD_true_z[i], foundZ[i]);
            // if (NHits[i] > 100) {
                h_foundZ       -> Fill(foundZ[i]);
                h_MBDZ         -> Fill(MBD_true_z[i]);
                g_ZCorrelation -> SetPoint(g_ZCorrelation->GetN(), foundZ[i], MBD_true_z[i]);
            // }
            
            z_resolution = foundZ[i] - MBD_true_z[i];
            h_resolution   -> Fill(z_resolution);
            g_NFunction    -> SetPoint(g_NFunction->GetN(), NHits[i], z_resolution);
            g_ZFunction    -> SetPoint(g_ZFunction->GetN(), MBD_true_z[i], z_resolution);
            g_CFunction    -> SetPoint(g_CFunction->GetN(), calculated_cen[i], z_resolution);
            g_Calibration  -> SetPoint(g_Calibration->GetN(), MBD_true_z[i], MBD_true_z[i]);
        }
        // h_resolution -> Draw("");
        h_foundZ->SetTitle(Form("INTT found z/MBD z correlation of %d events",static_cast<int>(h_foundZ->GetEntries())));
        h_foundZ -> Draw(); h_foundZ->SetLineColor(kRed);
        h_MBDZ -> Draw("same"); 
        // h_foundZ->GetYaxis()->SetRangeUser(0,1400);
        static TLegend *lg = new TLegend(0.12, 0.8, 0.38, 0.9);
        lg -> AddEntry(h_foundZ, "found Z", "l");
        lg -> AddEntry(h_MBDZ, "MBD_z_vtx", "l");
        // lg -> AddEntry(h_MBDZ, "NClus > 100", "p");
        lg -> SetTextSize(.032);
        lg->Draw("same");

        g_ZCorrelation->SetTitle(Form("INTT found z/MBD z of %d events, NClus > 100", g_ZCorrelation->GetN()));
        g_ZCorrelation->GetXaxis()->SetTitle("found Z [cm]"); g_ZCorrelation->GetYaxis()->SetTitle("MBD_z_vtx [cm]");
        g_ZCorrelation->GetXaxis()->CenterTitle(true); g_ZCorrelation->GetYaxis()->CenterTitle(true);
        // g_ZCorrelation->Draw("AP");
        g_ZCorrelation->GetXaxis()->SetLimits(-40, 30);
        g_ZCorrelation->SetMinimum(-40); // Setting y range;
        g_ZCorrelation->SetMaximum(30);  // Setting y range;
        g_ZCorrelation->SetMarkerStyle(29); g_ZCorrelation -> SetMarkerColor(12); g_ZCorrelation -> SetMarkerSize(1.6);
        // TLine *lg = new TLine(0, 0, 9000, 0);
        // lg -> Draw("same"); lg -> SetLineColor(kRed); lg -> SetLineWidth(2);
        TH1D *h_projX = h->ProjectionX("projX", h->GetXaxis()->FindBin(-250), h->GetXaxis()->FindBin(-50));
        gPad->SetGrid(1, 1); gPad -> Update();
    }
    
    // Stop the stopwatch and print the runtime:
    timer.Stop();   timer.Print();

    theApp.Run();
    return 0;
}