#include "../headers/globalDefinitions.h"

using namespace rootBranchEssentials;
// void printTreeInfo(const char* filename, const char* treeName) {
void rootFileChecker_AI (std::string opt = "") {
    // TFile *file = TFile::Open("../External/Data_CombinedNtuple_Run54280_20241113.root");
    TFile *file = TFile::Open("../External/Sim_Ntuple_HIJING_ana443_20241102.root");
    objectExistenceCheck(file);
    TTree* tree = (TTree*)file->Get("EventTree");
    objectExistenceCheck(tree);

    TObjArray* branches = tree->GetListOfBranches();
    std::string branchName;
    int idx_event, idx_MBD_z_vtx, idx_MBD_centrality, idx_NClus, idx_ClusLayer,
        idx_ClusX, idx_ClusY, idx_ClusZ, idx_ClusR, idx_ClusPhi, idx_ClusEta;
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

    if (opt == "print") {
        for (int i = 0; i < branches->GetEntries(); ++i) {
            TBranch* branch = (TBranch*)branches->At(i);
            branchName = branch->GetName();
            TLeaf* leaf = branch->GetLeaf(branchName.c_str());
            std::string leafType = (leaf) ? leaf->GetTypeName() : "Unknown";
            std::cout << "Branch " << i << ": Name = " << branch->GetName()
                    << ", Type = " << leafType << std::endl;
        }
    }
    if (opt == "size") {
        TBranch *b_event     = (TBranch*)tree->GetListOfBranches()->At(idx_event);
        TBranch *b_MBD_z_vtx = (TBranch*)tree->GetListOfBranches()->At(idx_MBD_z_vtx);
        TBranch *b_MBD_centrality = (TBranch*)tree->GetListOfBranches()->At(idx_MBD_centrality);
        TBranch *b_NClus     = (TBranch*)tree->GetListOfBranches()->At(idx_NClus);
        TBranch *b_ClusLayer = (TBranch*)tree->GetListOfBranches()->At(idx_ClusLayer);
        TBranch *b_ClusX     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusX);
        TBranch *b_ClusY     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusY);
        TBranch *b_ClusZ     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusZ);
        TBranch *b_ClusR     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusR);
        TBranch *b_ClusPhi   = (TBranch*)tree->GetListOfBranches()->At(idx_ClusPhi);
        TBranch *b_ClusEta   = (TBranch*)tree->GetListOfBranches()->At(idx_ClusEta);
        b_event    ->SetAddress(&event);
        b_MBD_z_vtx->SetAddress(&MBD_z_vtx);
        b_MBD_centrality->SetAddress(&MBD_centrality);
        b_NClus    ->SetAddress(&NClus);
        b_ClusLayer->SetAddress(&ClusLayer);
        b_ClusX    ->SetAddress(&ClusX);
        b_ClusY    ->SetAddress(&ClusY);
        b_ClusZ    ->SetAddress(&ClusZ);
        b_ClusR    ->SetAddress(&ClusR);
        b_ClusPhi  ->SetAddress(&ClusPhi);
        b_ClusEta  ->SetAddress(&ClusEta);
        // ***********************************************
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            b_event->GetEntry(i);
            b_MBD_z_vtx->GetEntry(i);
            b_MBD_centrality->GetEntry(i);
            b_NClus->GetEntry(i);
            b_ClusLayer->GetEntry(i);
            b_ClusX->GetEntry(i);
            b_ClusY->GetEntry(i);
            b_ClusZ->GetEntry(i);
            b_ClusR->GetEntry(i);
            b_ClusPhi->GetEntry(i);
            b_ClusEta->GetEntry(i);
            if (NClus!=ClusX->size() || ClusX->size()!=ClusLayer->size() || ClusX->size()!=ClusY->size()
               || ClusX->size()!=ClusZ->size() || NClus!=ClusR->size()
               || NClus!=ClusPhi->size() || NClus!=ClusEta->size()) {
                printRed("Error!\n");
                exit(1);
            }
            // std::cout <<i<<","<<event<< ","<< MBD_centrality << std::endl;
        }
    }
    if (opt == "phi") {
        TBranch *b_event     = (TBranch*)tree->GetListOfBranches()->At(idx_event);
        TBranch *b_MBD_z_vtx = (TBranch*)tree->GetListOfBranches()->At(idx_MBD_z_vtx);
        TBranch *b_NClus     = (TBranch*)tree->GetListOfBranches()->At(idx_NClus);
        TBranch *b_ClusLayer = (TBranch*)tree->GetListOfBranches()->At(idx_ClusLayer);
        TBranch *b_ClusX     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusX);
        TBranch *b_ClusY     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusY);
        TBranch *b_ClusZ     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusZ);
        TBranch *b_ClusR     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusR);
        TBranch *b_ClusPhi   = (TBranch*)tree->GetListOfBranches()->At(idx_ClusPhi);
        TBranch *b_ClusEta   = (TBranch*)tree->GetListOfBranches()->At(idx_ClusEta);
        b_event    ->SetAddress(&event);
        b_MBD_z_vtx->SetAddress(&MBD_z_vtx);
        b_NClus    ->SetAddress(&NClus);
        b_ClusLayer->SetAddress(&ClusLayer);
        b_ClusX    ->SetAddress(&ClusX);
        b_ClusY    ->SetAddress(&ClusY);
        b_ClusZ    ->SetAddress(&ClusZ);
        b_ClusR    ->SetAddress(&ClusR);
        b_ClusPhi  ->SetAddress(&ClusPhi);
        b_ClusEta  ->SetAddress(&ClusEta);
        Long64_t nEntries = tree->GetEntries();
        TH1D *h = new TH1D("phi sanity checker",";calculated_phi - ClusPhi;# of counts", 100, -M_PI, M_PI);
        double phi;
        for (Long64_t i = 0; i < nEntries; ++i) {
            b_event->GetEntry(i);
            b_MBD_z_vtx->GetEntry(i);
            b_NClus->GetEntry(i);
            b_ClusLayer->GetEntry(i);
            b_ClusX->GetEntry(i);
            b_ClusY->GetEntry(i);
            b_ClusZ->GetEntry(i);
            b_ClusR->GetEntry(i);
            b_ClusPhi->GetEntry(i);
            b_ClusEta->GetEntry(i);
            for (int j = 0; j < ClusX->size(); j++) {
                    phi = std::atan2(ClusY->at(j),ClusX->at(j));
                if (phi != ClusPhi->at(j)) {
                    // std::cout << i << "," << event << "," << phi <<","<< ClusPhi->at(j) <<","<<phi - ClusPhi->at(j) << std::endl;
                    h->Fill(phi - ClusPhi->at(j));
                }
            }
        }
        h->Draw();
    }
    if (opt == "eta") {
        TBranch *b_event     = (TBranch*)tree->GetListOfBranches()->At(idx_event);
        TBranch *b_MBD_z_vtx = (TBranch*)tree->GetListOfBranches()->At(idx_MBD_z_vtx);
        TBranch *b_NClus     = (TBranch*)tree->GetListOfBranches()->At(idx_NClus);
        TBranch *b_ClusLayer = (TBranch*)tree->GetListOfBranches()->At(idx_ClusLayer);
        TBranch *b_ClusX     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusX);
        TBranch *b_ClusY     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusY);
        TBranch *b_ClusZ     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusZ);
        TBranch *b_ClusR     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusR);
        TBranch *b_ClusPhi   = (TBranch*)tree->GetListOfBranches()->At(idx_ClusPhi);
        TBranch *b_ClusEta   = (TBranch*)tree->GetListOfBranches()->At(idx_ClusEta);
        b_event    ->SetAddress(&event);
        b_MBD_z_vtx->SetAddress(&MBD_z_vtx);
        b_NClus    ->SetAddress(&NClus);
        b_ClusLayer->SetAddress(&ClusLayer);
        b_ClusX    ->SetAddress(&ClusX);
        b_ClusY    ->SetAddress(&ClusY);
        b_ClusZ    ->SetAddress(&ClusZ);
        b_ClusR    ->SetAddress(&ClusR);
        b_ClusPhi  ->SetAddress(&ClusPhi);
        b_ClusEta  ->SetAddress(&ClusEta);
        Long64_t nEntries = tree->GetEntries();
        TH1D *h = new TH1D("eta sanity check",";calculated_eta - ClusEta;#", 100, -2*M_PI, 2*M_PI);
        double R, z_vtx, dZ, theta_half, eta;
        for (Long64_t i = 0; i < nEntries; ++i) {
            b_event->GetEntry(i);
            b_MBD_z_vtx->GetEntry(i);
            b_NClus->GetEntry(i);
            b_ClusLayer->GetEntry(i);
            b_ClusX->GetEntry(i);
            b_ClusY->GetEntry(i);
            b_ClusZ->GetEntry(i);
            b_ClusR->GetEntry(i);
            b_ClusPhi->GetEntry(i);
            b_ClusEta->GetEntry(i);
            for (int j = 0; j < ClusX->size(); j++) {
                R = ClusR->at(j);
                dZ = ClusZ->at(j);
                theta_half = (std::atan2(R, dZ))/2;
                if (dZ >= 0)    eta = -std::log(std::tan(theta_half));
                if (dZ <  0)    eta = std::log(std::tan(halfPI - theta_half));
                if (eta != ClusEta->at(j)) {
                    // std::cout << i << "," << event << "," << phi <<","<< ClusPhi->at(j) <<","<<phi - ClusPhi->at(j) << std::endl;
                    h->Fill(eta - ClusEta->at(j));
                }
            }
        }
        h->Draw();
    }
    if (opt == "list") {
        std::string filePath = "./zFindingResults/dummyForFit_simulated_data.txt";
        std::ofstream outputFile(filePath);
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
        if(isEmpty) outputFile << "index,event,NClus,foundZ,trueZ,centrality,chargeSum" << std::endl;

        TBranch *b_event     = (TBranch*)tree->GetListOfBranches()->At(idx_event);
        TBranch *b_MBD_z_vtx = (TBranch*)tree->GetListOfBranches()->At(idx_MBD_z_vtx);
        TBranch *b_MBD_centrality = (TBranch*)tree->GetListOfBranches()->At(idx_MBD_centrality);
        b_MBD_charge_sum = static_cast<TBranch*>(tree->GetListOfBranches()->At(idx_MBD_charge_sum));
        TBranch *b_NClus     = (TBranch*)tree->GetListOfBranches()->At(idx_NClus);
        TBranch *b_ClusLayer = (TBranch*)tree->GetListOfBranches()->At(idx_ClusLayer);
        TBranch *b_ClusX     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusX);
        TBranch *b_ClusY     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusY);
        TBranch *b_ClusZ     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusZ);
        TBranch *b_ClusR     = (TBranch*)tree->GetListOfBranches()->At(idx_ClusR);
        TBranch *b_ClusPhi   = (TBranch*)tree->GetListOfBranches()->At(idx_ClusPhi);
        TBranch *b_ClusEta   = (TBranch*)tree->GetListOfBranches()->At(idx_ClusEta);
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
        Long64_t nEntries = tree->GetEntries();
        double R, z_vtx, dZ, theta_half, eta;
        for (Long64_t i = 0; i < nEntries; ++i) {
            b_event->GetEntry(i);
            b_MBD_z_vtx->GetEntry(i);
            b_MBD_centrality->GetEntry(i);
            b_MBD_charge_sum->GetEntry(i);
            b_NClus->GetEntry(i);
            b_ClusLayer->GetEntry(i);
            b_ClusX->GetEntry(i);
            b_ClusY->GetEntry(i);
            b_ClusZ->GetEntry(i);
            b_ClusR->GetEntry(i);
            b_ClusPhi->GetEntry(i);
            b_ClusEta->GetEntry(i);
            // if (std::abs(MBD_z_vtx)<2 && NClus >= 100) 
            if (MBD_z_vtx > -40 && MBD_z_vtx < 30)
                outputFile << i << "," << event << "," << NClus << "," << std::nan("") << "," << MBD_z_vtx << "," << MBD_centrality << "," << MBD_charge_sum << std::endl;
        }
    }

    // Clean up
    // file->Close();
}
