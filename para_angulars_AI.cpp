#include "Fit/BinData.h"
#include <TStopwatch.h>
#include "TApplication.h"
#include <thread>
#include <chrono>
#include <mutex>
#include <omp.h>
#include <TSpectrum.h>

#include "headers/histogramPlotting.h"
#include "headers/CommonParaFuncs.h"
#include "headers/zFinding.h"
#include "ROOT/TThreadedObject.hxx"

using namespace rootBranchEssentials;

std::mutex m_mutex;
int N = 1000;
double range_min = -M_PI;
double range_max = M_PI;
double bin_width = (range_max - range_min) / N;
TH1D *h_dPhi_nomix;
TH1D *h_ZonOne[20];
TH1D *h_CenonOne[14];
std::vector<double> Phi0, Phi1;
TH1D *h_Background_dPhi = new TH1D("", "", N, range_min, range_max);

void para_dPhiAccumulator (const int &lower, const int &upper) {
    double dPhi;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int k = lower; k < upper; k++) {
        for (int l = 0; l < Phi1.size(); l++) {
            dPhi = Phi0[k] - Phi1[l];
            if (dPhi > M_PI)    dPhi -= TWO_PI;
            if (dPhi < -M_PI)   dPhi += TWO_PI;
            h_dPhi_nomix->Fill(dPhi); 
        }
    }
}

void dPhi_in_bins_of_Centrality (
    const int &id,
    const int &num,
    const std::vector<int>    &index,
    const std::vector<double> &MBD_true_z,
    const std::vector<double> &cen
) {
    h_CenonOne[id] = new TH1D(Form("dPhi of %.2f to %.2f", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), Form("dPhi of %.2f to %.2f;dPhi;# of counts", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), N, range_min, range_max);
    double phi, dPhi;
    int ClusPhiSize, Phi0Size, Phi1Size, layer;
    // std::vector<double> Phi0, Phi1;
    std::unique_lock<std::mutex> lock(m_mutex);
    double cen_lower_range  = static_cast<double>(id)*0.05;
    double cen_higher_range = static_cast<double>(id+1)*0.05;
    for (int i = 0; i < num; i++) {
        if (cen[i] >= cen_lower_range && cen[i] <= cen_higher_range) {
            b_ClusPhi->GetEntry(index[i]);   // ClusPhi;
            b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
            ClusPhiSize = ClusPhi->size();
            for (int j = 0; j < ClusPhiSize; j++) {
                phi = ClusPhi->at(j);
                layer = ClusLayer->at(j);
                if (layer == 3 || layer == 4) {
                    Phi0.push_back(phi);
                }
                else {
                    Phi1.push_back(phi);
                }
            }
            Phi0Size = Phi0.size();
            Phi1Size = Phi1.size();
            for (int k = 0; k < Phi0Size; k++) {
                for (int l = 0; l < Phi1Size; l++) {
                    dPhi = Phi0[k] - Phi1[l];
                    if (dPhi > M_PI)    dPhi -= TWO_PI;
                    if (dPhi < -M_PI)   dPhi += TWO_PI;
                    h_CenonOne[id] -> Fill(dPhi);
                }
            }

            Phi0.clear();   Phi1.clear();
        }
    }
    lock.unlock();
}

void dPhi_in_bins_of_Z_vtx (
    const int &id,
    const int &num,
    const std::vector<int> &index,
    const std::vector<double> &MBD_true_z,
    const std::vector<double> &cen
) {
    double phi, dPhi;
    int layer, Phi0Size, Phi1Size, ClusPhiSize;
    // std::vector<double> Phi0, Phi1;
    std::unique_lock<std::mutex> lock(m_mutex);
    double z_lower_range  = -6 - static_cast<double>(id);
    double z_higher_range = -5 - static_cast<double>(id);
    h_ZonOne[id] = new TH1D(Form("dPhi of %1.0f to %1.0f", z_lower_range, z_higher_range), Form("dPhi of %1.0f to %1.0f;dPhi;# of counts", z_lower_range, z_higher_range), N, range_min, range_max);
    for (int i = 0; i < num; i++) {
        if (MBD_true_z[i] >= z_lower_range && MBD_true_z[i] <= z_higher_range) {
            b_ClusPhi->GetEntry(index[i]);   // ClusPhi;
            b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
            ClusPhiSize = ClusPhi->size();
            for (int j = 0; j < ClusPhiSize; j++) {
                phi = ClusPhi->at(j);
                layer = ClusLayer->at(j);
                if (layer == 3 || layer == 4) {
                    Phi0.push_back(phi);
                }
                else {
                    Phi1.push_back(phi);
                }
            }

            Phi0Size = Phi0.size();
            Phi1Size = Phi1.size();
            for (int k = 0; k < Phi0Size; k++) {
                for (int l = 0; l < Phi1Size; l++) {
                    dPhi = Phi0[k] - Phi1[l];
                    if (dPhi > M_PI)    dPhi -= TWO_PI;
                    if (dPhi < -M_PI)   dPhi += TWO_PI;
                    h_ZonOne[id] -> Fill(dPhi);
                }
            }

            Phi0.clear();   Phi1.clear();
        }
    }
    lock.unlock();
}

int main(int argc, char* argv[]) {
    // Start the stopwatch:
    TStopwatch timer;   timer.Start();
    // Get and print the current PC time
    printSeparation();
    auto now = std::chrono::system_clock::now(); auto now_c = std::chrono::system_clock::to_time_t(now);
    printBlueSingle("Run started at: "); printBlueSingle(std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S\n"));
    TApplication theApp("App", &argc, argv);

    std::vector<std::string> method;
    // <METHOD> <num> <CEN_LOW> <CEN_HIGH> <Z_LOW> <Z_HIGH>
    int num         = 100;
    double cen_low  = 0.0;
    double cen_high = 0.7;
    double z_low    = 25.;
    double z_high   = 5.; 
    int target      = 4;

    if (argc > 1) {
        if (argc > 2) num      = std::stoi(argv[2]);
        if (argc > 3) cen_low  = std::stod(argv[3]);
        if (argc > 4) cen_high = std::stod(argv[4]);
        if (argc > 5) z_low    = std::stod(argv[5]);
        if (argc > 6) z_high   = std::stod(argv[6]);
        if (argc > 7) target   = std::stod(argv[7]);
        for (int i = 1; i < argc; i++) {
            std::string substring = argv[i];
            method.push_back(substring);
        }
    }
    std::vector<std::string> opts = splitString(method[0], '_');

    // TFile *file = TFile::Open("../../External/Data_CombinedNtuple_Run54280_20241113.root");
    TFile *file = TFile::Open("../../External/Sim_Ntuple_HIJING_ana443_20241102.root");
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
    b_event     ->SetAddress(&event);
    b_MBD_z_vtx ->SetAddress(&MBD_z_vtx);
    b_MBD_centrality->SetAddress(&MBD_centrality);
    b_MBD_charge_sum->SetAddress(&MBD_charge_sum);
    b_NClus     ->SetAddress(&NClus);
    b_ClusLayer ->SetAddress(&ClusLayer);
    b_ClusX     ->SetAddress(&ClusX);
    b_ClusY     ->SetAddress(&ClusY);
    b_ClusZ     ->SetAddress(&ClusZ);
    b_ClusR     ->SetAddress(&ClusR);
    b_ClusPhi   ->SetAddress(&ClusPhi);
    b_ClusEta   ->SetAddress(&ClusEta);

    // ↓↓ foundZ data access;
    // std::string filePath = "../zFindingResults/DCAfit_0_16_-001_001_step2e-1_z_-40_30_centrality_calculated.txt";
    std::string filePath = "../zFindingResults/dummyForFit_simulated_data.txt";
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

        if (opts[0] == "mix") {
            if (Mz >= z_low && Mz <= z_high && cc >= cen_low && cc <= cen_high && N > 100) {
                index.push_back(i);     evt.push_back(e);           NHits.push_back(N);
                foundZ.push_back(f);    MBD_true_z.push_back(Mz);   calculated_cen.push_back(cc);    
            }
        }
        else if (N > 100) {
            index.push_back(i);     evt.push_back(e);           NHits.push_back(N);
            foundZ.push_back(f);    MBD_true_z.push_back(Mz);   calculated_cen.push_back(cc);
        }
    }
    // ↑↑ foundZ data access;
    num = num > index.size() ? index.size() : num;
    ROOT::EnableThreadSafety();
    if (method[0] == "nomix") {
        h_dPhi_nomix = new TH1D("dPhi values", ";dPhi;# of counts", N, range_min, range_max);
        double phi, dPhi;
        int layer, layerSize, Phi0Size, Phi1Size;
        for (int i = 0; i < num; i++) {
            if (NHits[i] > 100) {
                b_ClusLayer->GetEntry(index[i]);
                b_ClusPhi  ->GetEntry(index[i]);
                for (int j = 0; j < NHits[i]; j++) {
                    phi = ClusPhi->at(j);
                    layer = ClusLayer->at(j);
                    if (layer == 3 || layer == 4) {
                        Phi0.push_back(phi);
                    } else {
                        Phi1.push_back(phi);
                    }
                }
                Phi0Size = Phi0.size();
                Phi1Size = Phi1.size();
                for (int k = 0; k < Phi0Size; k++) {
                    for (int l = 0; l < Phi1Size; l++) {
                        dPhi = Phi0[k] - Phi1[l];
                        if (dPhi > M_PI)    dPhi -= TWO_PI;
                        if (dPhi < -M_PI)   dPhi += TWO_PI;
                        h_dPhi_nomix->Fill(dPhi); 
                    }
                }
            }

            Phi0.clear(); Phi1.clear(); 
        }
        angularPlot1D(h_dPhi_nomix, method, Form("SIMULATION DATA dPhi of unmixed for %d events", num));
    }
    else if (method[0] == "para1") {
        h_dPhi_nomix = new TH1D("dPhi values", ";dPhi;# of counts", N, range_min, range_max);
        double phi, dPhi;
        int layer, layerSize;
        std::vector<int> boundaries;
        std::thread thsafe[8];
        for (int i = 0; i < num; i++) {
            if (NHits[i] > 100) {
                b_ClusLayer->GetEntry(index[i]);
                b_ClusPhi  ->GetEntry(index[i]);

                #pragma omp parallel
                {
                    std::vector<double> local_Phi0, local_Phi1;

                    #pragma omp for nowait
                    for (int j = 0; j < NHits[i]; j++) { // Adjust NHits index as necessary
                        double phi = ClusPhi->at(j);
                        int layer = ClusLayer->at(j);
                        if (layer == 3 || layer == 4) {
                            local_Phi0.push_back(phi);
                        } else {
                            local_Phi1.push_back(phi);
                        }
                    }

                    // Merge local vectors into the shared ones
                    #pragma omp critical
                    {
                        Phi0.insert(Phi0.end(), local_Phi0.begin(), local_Phi0.end());
                        Phi1.insert(Phi1.end(), local_Phi1.begin(), local_Phi1.end());
                    }
                }

                for (int k = 0; k < Phi0.size(); k++) {
                    for (int l = 0; l < Phi1.size(); l++) {
                        dPhi = Phi0[k] - Phi1[l];
                        if (dPhi > M_PI)    dPhi -= TWO_PI;
                        if (dPhi < -M_PI)   dPhi += TWO_PI;
                        h_dPhi_nomix->Fill(dPhi); 
                    }
                }
            }

            Phi0.clear(); Phi1.clear(); 
        }
        angularPlot1D(h_dPhi_nomix, method, "dPhi of unmixed");
    }
    else if (method[0] == "para2") {
        ROOT::TThreadedObject<TH1D> h_dPhi_para("dPhi values", ";dPhi;# of counts", N, range_min, range_max);
        double phi, dPhi;
        int layer, layerSize;
        std::vector<int> boundaries;
        std::thread thsafe[8];
        for (int i = 0; i < num; i++) {
            if (NHits[i] > 100) {
                b_ClusLayer->GetEntry(index[i]);
                b_ClusPhi  ->GetEntry(index[i]);
                for (int j = 0; j < NHits[i]; j++) {
                    phi = ClusPhi->at(j);
                    layer = ClusLayer->at(j);
                    if (layer == 3 || layer == 4) {
                        Phi0.push_back(phi);
                    } else {
                        Phi1.push_back(phi);
                    }
                }
            
                auto paraFill = [&](UInt_t id) {
                    auto h = h_dPhi_para.Get();
                    for (int k = 0; k < Phi0.size(); k++) {
                        for (int l = 0; l < Phi1.size(); l++) {
                            dPhi = Phi0[k] - Phi1[l];
                            if (dPhi > M_PI)    dPhi -= TWO_PI;
                            if (dPhi < -M_PI)   dPhi += TWO_PI;
                            h->Fill(dPhi); 
                        }
                    }
                };
                std::vector<std::thread> threads;
                for (auto threadID : ROOT::TSeqI(nThreads))
                    threads.emplace_back(paraFill, threadID);
                for (auto &&th : threads)
                    th.join();

                Phi0.clear(); Phi1.clear(); threads.clear();   
            }
        }
        auto h_final = h_dPhi_para.Merge();
        TCanvas *can1 = new TCanvas("c1d","c1d",0,50,1920,1056);
        h_final -> DrawCopy();
        can1 -> Update();
        // angularPlot1D(h_dPhi_nomix, method, "dPhi of unmixed");
    }
    else if (method[0] == "cen") {
        std::thread thsafe[14];
        std::cout << "dPhi of different centralities" << std::endl;
        for (int i = 0; i < 14; i++)
            thsafe[i] = std::thread(dPhi_in_bins_of_Centrality,i,num,std::cref(index),std::cref(MBD_true_z),std::cref(calculated_cen));

        for (int i = 0; i < 14; i++)
            thsafe[i].join();

        std::vector<TH1D*> h(h_CenonOne, h_CenonOne + 14);
        ArrayPlot1D_Rescale(h, method, Form("dPhi per centralities rescale with %d events", num));
        // ArrayPlot1D_Rescale(h, method, Form("dPhi_per_centralities_rescale_with_dEta_cut_%d_events",num));
    }
    else if (method[0] == "z") {
        std::thread thsafe[20];
        std::cout << "dPhi of different Z vertices" << std::endl;
        for (int i = 0; i < 20; i++)
            thsafe[i] = std::thread(dPhi_in_bins_of_Z_vtx,i,num,std::cref(index),std::cref(MBD_true_z),std::cref(calculated_cen));

        for (int i = 0; i < 20; i++)
            thsafe[i].join();

        std::vector<TH1D*> h(h_ZonOne, h_ZonOne + 20);
        ArrayPlot1D_Rescale_ver2(h, method, Form("dPhi_per_Z_vtx_rescale_%d_events", num));
    }
    else if (opts[0] == "mix") {
        if (opts[1] == "wo") {
            printRed("Mixing events without any cuts, ", index.size(), " events selected, ", num, " events involved.");
            double phi, dPhi; 
            int clus_layer, ClusPhiSize;
            std::vector<double> Phi0, Phi1;
            std::vector<std::vector <double>> event_Phi0, event_Phi1;
            // Signal dPhi is unmixed events' dPhi:
            TH1D *h_Signal_dPhi = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %d events;dPhi value;# of counts", num), N, range_min, range_max);
            for (int i = 0; i < num; i++) {  // Loop over events;
                b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
                b_ClusPhi->GetEntry(index[i]);     // ClusPhi;
                event_Phi0.push_back(std::vector <double>());   
                event_Phi1.push_back(std::vector <double>());
                ClusPhiSize = ClusPhi->size();
                for (int j = 0; j < ClusPhiSize; j++) { // Loop inside one event, over all hits;
                    phi        = ClusPhi->at(j);
                    clus_layer = ClusLayer->at(j);
                    if (clus_layer == 3 || clus_layer == 4) {
                        Phi0.push_back(phi);
                        event_Phi0[i].push_back(phi);
                    }
                    else {
                        Phi1.push_back(phi);
                        event_Phi1[i].push_back(phi);
                    }
                }
                for (int k = 0; k < Phi0.size(); k++) {
                    for (int l = 0; l < Phi1.size(); l++) {
                        dPhi = Phi0[k] - Phi1[l];
                        if (dPhi > M_PI)    dPhi -= TWO_PI;
                        if (dPhi < -M_PI)   dPhi += TWO_PI;
                        h_Signal_dPhi -> Fill(dPhi);
                    }
                }
                Phi0.clear();   Phi1.clear();
            }

            for (int m = 0; m < event_Phi0.size(); m++) {
                for (int n = m; n < event_Phi1.size(); n++) {
                    // Process pairs from event_Phi0[i] and event_Phi1[j]
                    for (double phi0 : event_Phi0[m]) {
                        for (double phi1 : event_Phi1[n]) {
                            dPhi = phi0 - phi1;
                            if (dPhi > M_PI)  dPhi -= TWO_PI;
                            if (dPhi < -M_PI) dPhi += TWO_PI;
                            h_Background_dPhi->Fill(dPhi);
                        }
                    }
                    // Process pairs from event_Phi0[j] and event_Phi1[i] to ensure symmetry
                    for (double phi0 : event_Phi0[n]) {
                        for (double phi1 : event_Phi1[m]) {
                            dPhi = phi0 - phi1;
                            if (dPhi > M_PI)  dPhi -= TWO_PI;
                            if (dPhi < -M_PI) dPhi += TWO_PI;
                            h_Background_dPhi->Fill(dPhi);
                        }
                    }
                }
            }

            backgroundCancelling_dPhi(h_Background_dPhi, h_Signal_dPhi, method, num);
        }
        else if (opts[1] == "single") {
            printBlue("Mixing events with |Eta| < ", Eta_range, ", ", index.size(), " events selected, ", num, " events involved, ", target, "th event's under processing");
            double phi, dPhi, found_z, dZ, R, halfTheta, eta, abs_eta;
            int clus_layer, ClusPhiSize;
            std::vector<double> Phi0, Phi1;
            std::vector<std::vector <double>> event_Phi0, event_Phi1;
            // Signal dPhi is unmixed events' dPhi:
            TH1D *h_Signal_dPhi = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %d events;dPhi value;# of counts", num), N, range_min, range_max);
            for (int i = 0; i < num; i++) {  // Loop over events;
                b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
                b_ClusPhi  ->GetEntry(index[i]);   // ClusPhi;
                b_ClusR    ->GetEntry(index[i]);   // ClusR;
                b_ClusZ    ->GetEntry(index[i]);   // ClusPhi;
                event_Phi0.push_back(std::vector <double>());   
                event_Phi1.push_back(std::vector <double>());
                ClusPhiSize = ClusPhi->size();
                // found_z = foundZ[i];
                found_z = MBD_true_z[i];
                for (int j = 0; j < ClusPhiSize; j++) { // Loop inside one event, over all hits;
                    phi        = ClusPhi->at(j);
                    clus_layer = ClusLayer->at(j);
                    dZ         = ClusZ->at(j) - found_z;
                    R          = ClusR->at(j);
                    halfTheta  = (std::atan2(R, dZ))/2.;
                    if (dZ >= 0)    eta = -std::log(std::tan(halfTheta));
                    if (dZ <  0)    eta = std::log(std::tan(halfPI - halfTheta));
                    abs_eta = std::abs(eta);
                    if (abs_eta < Eta_range) {
                        if (clus_layer == 3 || clus_layer == 4) {
                            if (index[i] == target)   Phi0.push_back(phi);
                            event_Phi0[i].push_back(phi);
                        }
                        else {
                            if (index[i] == target)   Phi1.push_back(phi);
                            event_Phi1[i].push_back(phi);
                        }
                    }
                }
                for (int k = 0; k < Phi0.size(); k++) {
                    for (int l = 0; l < Phi1.size(); l++) {
                        dPhi = Phi0[k] - Phi1[l];
                        if (dPhi > M_PI)    dPhi -= TWO_PI;
                        if (dPhi < -M_PI)   dPhi += TWO_PI;
                        h_Signal_dPhi -> Fill(dPhi);
                    }
                }
                Phi0.clear();   Phi1.clear();
            }
            printBlue(event_Phi0.size(), ", ", event_Phi1.size());
            for (int m = 0; m < event_Phi0.size(); m++) {
                for (int n = m; n < event_Phi1.size(); n++) {
                    // Process pairs from event_Phi0[i] and event_Phi1[j]
                    for (double phi0 : event_Phi0[m]) {
                        for (double phi1 : event_Phi1[n]) {
                            dPhi = phi0 - phi1;
                            if (dPhi > M_PI)  dPhi -= TWO_PI;
                            if (dPhi < -M_PI) dPhi += TWO_PI;
                            h_Background_dPhi->Fill(dPhi);
                        }
                    }
                    // Process pairs from event_Phi0[j] and event_Phi1[i] to ensure symmetry
                    for (double phi0 : event_Phi0[n]) {
                        for (double phi1 : event_Phi1[m]) {
                            dPhi = phi0 - phi1;
                            if (dPhi > M_PI)  dPhi -= TWO_PI;
                            if (dPhi < -M_PI) dPhi += TWO_PI;
                            h_Background_dPhi->Fill(dPhi);
                        }
                    }
                }
            }
            backgroundCancelling_dPhi(h_Background_dPhi, h_Signal_dPhi, method, num);
            // TCanvas *c = new TCanvas("c", "c", 0,50,2000,1000); // It's a must!
            // h_Signal_dPhi->Draw();
            // gPad -> SetGrid(1,1); gPad -> Update();                 // It's a must;
        }
        else if (opts[1] == "E") {
            printRed("Mixing events with |Eta| < ", Eta_range, ", ", index.size(), " events selected, ", num, " events involved.");
            double phi, dPhi, found_z, dZ, R, halfTheta, eta, abs_eta;
            int clus_layer, ClusPhiSize;
            std::vector<double> Phi0, Phi1;
            std::vector<std::vector <double>> event_Phi0, event_Phi1;
            // Signal dPhi is unmixed events' dPhi:
            TH1D *h_Signal_dPhi = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %d events;dPhi value;# of counts", num), N, range_min, range_max);
            for (int i = 0; i < num; i++) {  // Loop over events;
                b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
                b_ClusPhi  ->GetEntry(index[i]);   // ClusPhi;
                b_ClusR    ->GetEntry(index[i]);   // ClusR;
                b_ClusZ    ->GetEntry(index[i]);   // ClusPhi;
                event_Phi0.push_back(std::vector <double>());   
                event_Phi1.push_back(std::vector <double>());
                ClusPhiSize = ClusPhi->size();
                // found_z = foundZ[i];
                found_z = MBD_true_z[i];
                for (int j = 0; j < ClusPhiSize; j++) { // Loop inside one event, over all hits;
                    phi        = ClusPhi->at(j);
                    clus_layer = ClusLayer->at(j);
                    dZ         = ClusZ->at(j) - found_z;
                    R          = ClusR->at(j);
                    halfTheta  = (std::atan2(R, dZ))/2.;
                    if (dZ >= 0)    eta = -std::log(std::tan(halfTheta));
                    if (dZ <  0)    eta = std::log(std::tan(halfPI - halfTheta));
                    abs_eta = std::abs(eta);
                    if (abs_eta < Eta_range) {
                        if (clus_layer == 3 || clus_layer == 4) {
                            Phi0.push_back(phi);
                            event_Phi0[i].push_back(phi);
                        }
                        else {
                            Phi1.push_back(phi);
                            event_Phi1[i].push_back(phi);
                        }
                    }
                }
                for (int k = 0; k < Phi0.size(); k++) {
                    for (int l = 0; l < Phi1.size(); l++) {
                        dPhi = Phi0[k] - Phi1[l];
                        if (dPhi > M_PI)    dPhi -= TWO_PI;
                        if (dPhi < -M_PI)   dPhi += TWO_PI;
                        h_Signal_dPhi -> Fill(dPhi);
                    }
                }
                Phi0.clear();   Phi1.clear();
            }

            for (int m = 0; m < event_Phi0.size(); m++) {
                for (int n = m; n < event_Phi1.size(); n++) {
                    // Process pairs from event_Phi0[i] and event_Phi1[j]
                    for (double phi0 : event_Phi0[m]) {
                        for (double phi1 : event_Phi1[n]) {
                            dPhi = phi0 - phi1;
                            if (dPhi > M_PI)  dPhi -= TWO_PI;
                            if (dPhi < -M_PI) dPhi += TWO_PI;
                            h_Background_dPhi->Fill(dPhi);
                        }
                    }
                    // Process pairs from event_Phi0[j] and event_Phi1[i] to ensure symmetry
                    for (double phi0 : event_Phi0[n]) {
                        for (double phi1 : event_Phi1[m]) {
                            dPhi = phi0 - phi1;
                            if (dPhi > M_PI)  dPhi -= TWO_PI;
                            if (dPhi < -M_PI) dPhi += TWO_PI;
                            h_Background_dPhi->Fill(dPhi);
                        }
                    }
                }
            }
            backgroundCancelling_dPhi(h_Background_dPhi, h_Signal_dPhi, method, num);
        }
    }
    // Stop the stopwatch and print the runtime:
    timer.Stop();   timer.Print();

    theApp.Run();
    return 0;
}