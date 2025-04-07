#include "Fit/BinData.h"
#include "HFitInterface.h"
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/WrappedParamFunction.h"
#include "TROOT.h"
#include "TVirtualFitter.h"
#include <TStopwatch.h>
#include "headers/histogramPlotting.h"
#include "headers/CommonParaFuncs.h"
#include "TCanvas.h"
#include "TApplication.h"
#include <thread>
#include <mutex>
#include <chrono>

using namespace rootBranchEssentials;

TH1 *h[100];
TH1D *h_ZonOne[20];
TH1D *h_CenonOne[14];
std::mutex m_mutex;
int batch_number = 1250;   // 28147*8 = 225176
int N = 1000;
double range_min = -M_PI;
double range_max = M_PI;
double bin_width = (range_max - range_min) / N;
TH1D *h_dPhi_nomix = new TH1D("dPhi values", ";dPhi;# of counts", N, range_min, range_max);
TH1D *h_Background_dPhi = new TH1D("", "", N, range_min, range_max);

void dPhiAccumulator (
    const int &id,
    const int &batch_size,
    const std::vector<int> &index,
    const std::vector<double> &MBD_cen,
    TBranch *b_ClusLayer,
    TBranch* b_ClusPhi,
    const std::vector<float>* const ClusPhi,
    const std::vector<int>* const ClusLayer
) {
    double phi, dPhi;
    std::vector<double> Phi0, Phi1;
    
    int local_index;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < batch_size; i++) {
        local_index = i + batch_size*id;
        b_ClusPhi->GetEntry(index[local_index]);   b_ClusLayer->GetEntry(index[local_index]);
        // std::cout << MBD_true_z[local_index] << std::endl;
        // std::cout << id << ", " << local_index << std::endl;
        // std::cout << id << ", " << local_index << ", " << index[local_index] << ", " << MBD_true_z[local_index] << ", " << MBD_cen[local_index] << ", " << ClusPhi->size() << ", " << ClusLayer->size() << std::endl;
        Phi0.reserve(ClusPhi->size());
        Phi1.reserve(ClusPhi->size());
        for (int j = 0; j < ClusPhi->size(); j++) {
            phi = ClusPhi->at(j);
            if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                Phi0.push_back(phi);
            }
            else {
                Phi1.push_back(phi);
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
        Phi0.clear();   Phi1.clear();
    }
    lock.unlock();
}

void dPhiNoMixWithEta (
    const int &id,
    const int &batch_size,
    const std::vector<int> &index,
    const std::vector<double> &MBD_true_z,
    TBranch *b_ClusLayer,
    TBranch *b_ClusZ,
    TBranch *b_ClusR,
    TBranch* b_ClusPhi,
    const std::vector<float>* const ClusZ,
    const std::vector<float>* const ClusR,
    const std::vector<float>* const ClusPhi,
    const std::vector<int>*   const ClusLayer
) {
    double phi, dPhi, eta, dEta;
    double z_vtx, dZ, R, halfTheta;
    std::vector<EtaWithPhi> Phi0, Phi1;
    
    int local_index;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < batch_size; i++) {
        local_index = i + batch_size*id;
        b_ClusLayer->GetEntry(index[local_index]);   // ClusLayer;
        b_ClusZ->GetEntry(index[local_index]);   // ClusZ;
        b_ClusR->GetEntry(index[local_index]);   // ClusR;
        b_ClusPhi->GetEntry(index[local_index]);   // ClusPhi;
        z_vtx = MBD_true_z[local_index];
        for (int j = 0; j < ClusPhi->size(); j++) {
            dZ       = ClusZ->at(j) - z_vtx;
            R        = ClusR->at(j);
            halfTheta    = (std::atan2(R, dZ))/2.;
            phi      = ClusPhi->at(j);
            if (dZ >= 0)    eta = -std::log(std::tan(halfTheta));
            if (dZ <  0)    eta = std::log(std::tan(halfPI - halfTheta));
            if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                Phi0.emplace_back(eta, phi);
            }
            else {
                Phi1.emplace_back(eta, phi);
            }
        }
        for (int k = 0; k < Phi0.size(); k++) {
            for (int l = 0; l < Phi1.size(); l++) {
                dEta = Phi0[k].eta_value - Phi1[l].eta_value;
                dPhi = Phi0[k].phi_value - Phi1[l].phi_value;
                if (dPhi > M_PI)    dPhi -= TWO_PI;
                if (dPhi < -M_PI)   dPhi += TWO_PI;
                if (std::abs(dEta) < dEta_cut)
                    h_dPhi_nomix->Fill(dPhi);
            }
        }
        Phi0.clear();   Phi1.clear();
    }
    lock.unlock();
}

void dPhi_in_bins_of_Centrality (
    const int &id,
    const int &target,
    const std::vector<int>    &index,
    const std::vector<double> &MBD_true_z,
    const std::vector<double> &MBD_cen,
    TBranch *b_ClusLayer,
    TBranch* b_ClusPhi,
    const std::vector<float>* const ClusPhi,
    const std::vector<int>*   const ClusLayer
) {
    h_CenonOne[id] = new TH1D(Form("dPhi of %.2f to %.2f", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), Form("dPhi of %.2f to %.2f;dPhi;# of counts", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), N, range_min, range_max);
    double phi, dPhi;
    std::vector<double> Phi0, Phi1;
    std::unique_lock<std::mutex> lock(m_mutex);
    double cen_lower_range  = static_cast<double>(id)*0.05;
    double cen_higher_range = static_cast<double>(id+1)*0.05;
    for (int i = 0; i < target; i++) {
        if (MBD_cen[i] >= cen_lower_range && MBD_cen[i] <= cen_higher_range) {
            b_ClusPhi->GetEntry(index[i]);   // ClusPhi;
            b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
            Phi0.reserve(ClusPhi->size());  Phi1.reserve(ClusPhi->size());
            for (int j = 0; j < ClusPhi->size(); j++) {
                double phi = ClusPhi->at(j);
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Phi0.push_back(phi);
                }
                else {
                    Phi1.push_back(phi);
                }
            }

            for (int k = 0; k < Phi0.size(); k++) {
                for (int l = 0; l < Phi1.size(); l++) {
                    dPhi = Phi0[k] - Phi1[l];
                    if (dPhi > M_PI)    dPhi = dPhi - TWO_PI;
                    if (dPhi < -M_PI)   dPhi = dPhi + TWO_PI;
                    h_CenonOne[id] -> Fill(dPhi);
                }
            }

            Phi0.clear();   Phi1.clear();
        }
    }
    lock.unlock();
}

void dPhi_in_bins_of_Centrality_with_dEta_cut (
    const int &id,
    const int &target,
    const std::vector<int>    &index,
    const std::vector<double> &MBD_true_z,
    const std::vector<double> &MBD_cen,
    TBranch *b_ClusLayer,
    TBranch *b_ClusZ,
    TBranch *b_ClusR,
    TBranch *b_ClusPhi,
    TBranch *b_MBD_z_vtx,
    float &MBD_z_vtx,
    const std::vector<int>*   const ClusLayer,
    const std::vector<float>* const ClusZ,
    const std::vector<float>* const ClusR,
    const std::vector<float>* const ClusPhi
) {
    h_CenonOne[id] = new TH1D(Form("dPhi of %.2f to %.2f", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), Form("dPhi of %.2f to %.2f;dPhi;# of counts", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), N, range_min, range_max);
    double phi, dPhi, eta, dEta;
    double z_vtx, dZ, R, halfTheta;
    std::vector<EtaWithPhi> Phi0, Phi1;
    std::unique_lock<std::mutex> lock(m_mutex);
    double cen_lower_range  = static_cast<double>(id)*0.05;
    double cen_higher_range = static_cast<double>(id+1)*0.05;
    for (int i = 0; i < target; i++) {
        if (MBD_cen[i] >= cen_lower_range && MBD_cen[i] <= cen_higher_range) {
            b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
            b_ClusZ->GetEntry(index[i]);   // ClusZ;
            b_ClusR->GetEntry(index[i]);   // ClusR;
            b_ClusPhi->GetEntry(index[i]);   // ClusPhi;
            b_MBD_z_vtx->GetEntry(index[i]);   // MBD_z_vtx;
            z_vtx = MBD_true_z[i];
            // if (std::abs(z_vtx - MBD_z_vtx) > 1e-4) {
            //     printRed("Error!!");
            //     exit(1);
            // }
            for (int j = 0; j < ClusPhi->size(); j++) {
                dZ       = ClusZ->at(j) - z_vtx;
                R        = ClusR->at(j);
                halfTheta    = (std::atan2(R, dZ))/2.;
                phi      = ClusPhi->at(j);
                if (dZ >= 0)    eta = -std::log(std::tan(halfTheta));
                if (dZ <  0)    eta = std::log(std::tan(halfPI - halfTheta));
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Phi0.emplace_back(eta, phi);
                }
                else {
                    Phi1.emplace_back(eta, phi);
                }
            }

            for (int k = 0; k < Phi0.size(); k++) {
                for (int l = 0; l < Phi1.size(); l++) {
                    dEta = Phi0[k].eta_value - Phi1[l].eta_value;
                    dPhi = Phi0[k].phi_value - Phi1[l].phi_value;
                    if (dPhi > M_PI)    dPhi -= TWO_PI;
                    if (dPhi < -M_PI)   dPhi += TWO_PI;
                    if (std::abs(dEta) < dEta_cut)
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
    const int &target,
    const std::vector<int> &index,
    const std::vector<double> &MBD_true_z,
    const std::vector<double> &MBD_cen,
    TBranch *b_ClusLayer,
    TBranch* b_ClusPhi,
    const std::vector<float>* const ClusPhi,
    const std::vector<int>* const ClusLayer
) {
    // h_ZonOne[id] = new TH1D(Form("dPhi of %.2f to %.2f", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), Form("dPhi of %.2f to %.2f;dPhi;# of counts", static_cast<double>(id)*0.05, static_cast<double>(id + 1)*0.05), N, range_min, range_max);
    double phi, dPhi;
    std::vector<double> Phi0, Phi1;
    std::unique_lock<std::mutex> lock(m_mutex);
    double z_lower_range  = -6 - static_cast<double>(id);
    double z_higher_range = -5 - static_cast<double>(id);
    h_ZonOne[id] = new TH1D(Form("dPhi of %1.0f to %1.0f", z_lower_range, z_higher_range), Form("dPhi of %1.0f to %1.0f;dPhi;# of counts", z_lower_range, z_higher_range), N, range_min, range_max);
    for (int i = 0; i < target; i++) {
        if (MBD_true_z[i] >= z_lower_range && MBD_true_z[i] <= z_higher_range) {
            b_ClusPhi->GetEntry(index[i]);   // ClusPhi;
            b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
            Phi0.reserve(ClusPhi->size());  Phi1.reserve(ClusPhi->size());
            for (int j = 0; j < ClusPhi->size(); j++) {
                double phi = ClusPhi->at(j);
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Phi0.push_back(phi);
                }
                else {
                    Phi1.push_back(phi);
                }
            }

            for (int k = 0; k < Phi0.size(); k++) {
                for (int l = 0; l < Phi1.size(); l++) {
                    dPhi = Phi0[k] - Phi1[l];
                    if (dPhi > M_PI)    dPhi = dPhi - TWO_PI;
                    if (dPhi < -M_PI)   dPhi = dPhi + TWO_PI;
                    h_ZonOne[id] -> Fill(dPhi);
                }
            }

            Phi0.clear();   Phi1.clear();
        }
    }
    lock.unlock();
}

void dPhi_in_bins_of_Z_vtx_with_dEta_cut (
    const int &id,
    const int &target,
    const std::vector<int> &index,
    const std::vector<double> &MBD_true_z,
    const std::vector<double> &MBD_cen,
    TBranch *b_ClusLayer,
    TBranch *b_ClusZ,
    TBranch *b_ClusR,
    TBranch* b_ClusPhi,
    const std::vector<int>*   const ClusLayer,
    const std::vector<float>* const ClusZ,
    const std::vector<float>* const ClusR,
    const std::vector<float>* const ClusPhi
) {
    double phi, dPhi, eta, dEta;
    double z_vtx, dZ, R, halfTheta;
    std::vector<EtaWithPhi> Phi0, Phi1;
    std::unique_lock<std::mutex> lock(m_mutex);
    double z_lower_range  = -30 + static_cast<double>(id)*3;
    double z_higher_range = -27 + static_cast<double>(id)*3;
    h_ZonOne[id] = new TH1D(Form("dPhi of %1.0f to %1.0f", z_lower_range, z_higher_range), Form("dPhi of %1.0f to %1.0f;dPhi;# of counts", z_lower_range, z_higher_range), N, range_min, range_max);
    for (int i = 0; i < target; i++) {
        z_vtx = MBD_true_z[i];
        if (z_vtx >= z_lower_range && z_vtx <= z_higher_range) {
            b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
            b_ClusZ->GetEntry(index[i]);   // ClusZ;
            b_ClusR->GetEntry(index[i]);   // ClusR;
            b_ClusPhi->GetEntry(index[i]);   // ClusPhi;
            for (int j = 0; j < ClusPhi->size(); j++) {
                dZ       = ClusZ->at(j) - z_vtx;
                R        = ClusR->at(j);
                halfTheta    = (std::atan2(R, dZ))/2.;
                phi      = ClusPhi->at(j);
                if (dZ >= 0)    eta = -std::log(std::tan(halfTheta));
                if (dZ <  0)    eta = std::log(std::tan(halfPI - halfTheta));
                if (ClusLayer->at(j) == 3 || ClusLayer->at(j) == 4) {
                    Phi0.emplace_back(eta, phi);
                }
                else {
                    Phi1.emplace_back(eta, phi);
                }
            }

            for (int k = 0; k < Phi0.size(); k++) {
                for (int l = 0; l < Phi1.size(); l++) {
                    dEta = Phi0[k].eta_value - Phi1[l].eta_value;
                    dPhi = Phi0[k].phi_value - Phi1[l].phi_value;
                    if (dPhi > M_PI)    dPhi -= TWO_PI;
                    if (dPhi < -M_PI)   dPhi += TWO_PI;
                    if (std::abs(dEta) < dEta_cut) {
                        h_ZonOne[id] -> Fill(dPhi);
                    }
                }
            }

            Phi0.clear();   Phi1.clear();
        }
    }
    lock.unlock();
}

void dPhi_mixing (
    const int &id,
    const int &chunck_size,
    const std::vector<std::vector<double>> &event_Phi0,
    const std::vector<std::vector<double>> &event_Phi1
) {
    int local_index;
    double dPhi;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < chunck_size; i++) {
        local_index = i + chunck_size*id;
        for (int j = local_index; j < event_Phi1.size(); j++) {
            // Process pairs from event_Phi0[i] and event_Phi1[j]
            for (double phi0 : event_Phi0[local_index]) {
                for (double phi1 : event_Phi1[j]) {
                    dPhi = phi0 - phi1;
                    if (dPhi > M_PI)  dPhi -= TWO_PI;
                    if (dPhi < -M_PI) dPhi += TWO_PI;
                    h_Background_dPhi->Fill(dPhi);
                }
            }
            // Process pairs from event_Phi0[j] and event_Phi1[i] to ensure symmetry
            for (double phi0 : event_Phi0[j]) {
                for (double phi1 : event_Phi1[local_index]) {
                    dPhi = phi0 - phi1;
                    if (dPhi > M_PI)  dPhi -= TWO_PI;
                    if (dPhi < -M_PI) dPhi += TWO_PI;
                    h_Background_dPhi->Fill(dPhi);
                }
            }
        }
    }
    lock.unlock();
}

// Using an averaged load for each thread:
void dPhi_mixing_with_dEta_cut (
    const int &id,
    const int &chunck_size,
    const std::vector<std::vector<EtaWithPhi>> &event_Phi0,
    const std::vector<std::vector<EtaWithPhi>> &event_Phi1
) {
    int local_index;
    double dPhi, dEta;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < chunck_size; i++) {
        local_index = i + chunck_size*id;
        for (int j = local_index; j < event_Phi1.size(); j++) {
            // Process pairs from event_Phi0[i] and event_Phi1[j]
            for (const EtaWithPhi& phi0 : event_Phi0[local_index]) {
                for (const EtaWithPhi& phi1 : event_Phi1[j]) {
                    dPhi = phi0.phi_value - phi1.phi_value;
                    dEta = phi0.eta_value - phi1.eta_value;
                    if (dPhi > M_PI)  dPhi -= TWO_PI;
                    if (dPhi < -M_PI) dPhi += TWO_PI;
                    if (std::abs(dEta) < dEta_cut) h_Background_dPhi->Fill(dPhi);
                }
            }
            // Process pairs from event_Phi0[j] and event_Phi1[i] to ensure symmetry
            for (const EtaWithPhi& phi0 : event_Phi0[j]) {
                for (const EtaWithPhi& phi1 : event_Phi1[local_index]) {
                    dPhi = phi0.phi_value - phi1.phi_value;
                    dEta = phi0.eta_value - phi1.eta_value;
                    if (dPhi > M_PI)  dPhi -= TWO_PI;
                    if (dPhi < -M_PI) dPhi += TWO_PI;
                    if (std::abs(dEta) < dEta_cut) h_Background_dPhi->Fill(dPhi);
                }
            }
        }
    }
    lock.unlock();
}

// Using a calculated, balanced (thought to be) load for each thread:
void dPhi_mixing_with_dEta_cut_ver2 (
    const int &starter,
    const int &ender,
    const std::vector<std::vector<EtaWithPhi>> &event_Phi0,
    const std::vector<std::vector<EtaWithPhi>> &event_Phi1
) {
    double dPhi, dEta;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = starter; i <= ender; i++) {
        for (int j = i; j < event_Phi1.size(); j++) {
            // Process pairs from event_Phi0[i] and event_Phi1[j]
            for (const EtaWithPhi& phi0 : event_Phi0[i]) {
                for (const EtaWithPhi& phi1 : event_Phi1[j]) {
                    dPhi = phi0.phi_value - phi1.phi_value;
                    dEta = phi0.eta_value - phi1.eta_value;
                    if (dPhi > M_PI)  dPhi -= TWO_PI;
                    if (dPhi < -M_PI) dPhi += TWO_PI;
                    if (std::abs(dEta) < dEta_cut) h_Background_dPhi->Fill(dPhi);
                }
            }
            // Process pairs from event_Phi0[j] and event_Phi1[i] to ensure symmetry
            for (const EtaWithPhi& phi0 : event_Phi0[j]) {
                for (const EtaWithPhi& phi1 : event_Phi1[i]) {
                    dPhi = phi0.phi_value - phi1.phi_value;
                    dEta = phi0.eta_value - phi1.eta_value;
                    if (dPhi > M_PI)  dPhi -= TWO_PI;
                    if (dPhi < -M_PI) dPhi += TWO_PI;
                    if (std::abs(dEta) < dEta_cut) h_Background_dPhi->Fill(dPhi);
                }
            }
        }
    }
    lock.unlock();
}

void dPhi_mixing_with_Eta_range (
    const int &id,
    const int &chunck_size,
    const std::vector<std::vector<EtaWithPhi>> &event_Phi0,
    const std::vector<std::vector<EtaWithPhi>> &event_Phi1
) {
    int local_index;
    double dPhi, dEta, eta1, eta2;
    std::unique_lock<std::mutex> lock(m_mutex);
    for (int i = 0; i < chunck_size; i++) {
        local_index = i + chunck_size*id;
        for (int j = local_index; j < event_Phi1.size(); j++) {
            // Process pairs from event_Phi0[i] and event_Phi1[j]
            for (const EtaWithPhi& phi0 : event_Phi0[local_index]) {
                for (const EtaWithPhi& phi1 : event_Phi1[j]) {
                    dPhi = phi0.phi_value - phi1.phi_value;
                    if (dPhi > M_PI)  dPhi -= TWO_PI;
                    if (dPhi < -M_PI) dPhi += TWO_PI;
                    h_Background_dPhi->Fill(dPhi);
                }
            }
            // Process pairs from event_Phi0[j] and event_Phi1[i] to ensure symmetry
            for (const EtaWithPhi& phi0 : event_Phi0[j]) {
                for (const EtaWithPhi& phi1 : event_Phi1[local_index]) {
                    dPhi = phi0.phi_value - phi1.phi_value;
                    if (dPhi > M_PI)  dPhi -= TWO_PI;
                    if (dPhi < -M_PI) dPhi += TWO_PI;
                    h_Background_dPhi->Fill(dPhi);
                }
            }
        }
    }
    lock.unlock();
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
    double cen_low   = 0.0;
    double cen_high  = 0.7;
    double z_low     = 25.;
    double z_high    = 5.; 

    if (argc > 1) {
        if (argc > 2) target   = std::stoi(argv[2]);
        if (argc > 3) cen_low  = std::stod(argv[3]);
        if (argc > 4) cen_high = std::stod(argv[4]);
        if (argc > 5) z_low    = std::stod(argv[5]);
        if (argc > 6) z_high   = std::stod(argv[6]);
        for (int i = 1; i < argc; i++) {
            std::string substring = argv[i];
            method.push_back(substring);
        }
    }
    std::vector<std::string> cmds = splitString(method[0], '_');

    // TFile *file = TFile::Open("../../External/Data_CombinedNtuple_Run20869_HotDead_BCO_ADC_Survey.root");
    TFile *file = TFile::Open("../../External/Sim_Ntuple_HIJING_ana443_20241102.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        exit(1);
    }
    TTree *EventTree = (TTree*)file->Get("EventTree");
    if (!EventTree) {
        std::cerr << "Error getting TTree!" << std::endl;
        exit(1);
    }

    // std::string filePath = "../zFindingResults/new_foundZ_DCAfit_0_16_-001_001_step2e-1_z_25_5.txt";
    std::string filePath = "../zFindingResults/dummy.txt";
    std::ifstream myfile(filePath);
    if (!myfile.is_open()){
		std::cout << "Unable to open linelabel" << std::endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
 	}

    std::string line, value;
    std::vector<int> index, evt, NHits;
    std::vector<double> foundZ, MBD_true_z, MBD_cen;

    getline(myfile, line);
    while (getline(myfile, line)) {
        std::stringstream data(line);
        getline(data, value, ',');  int i     = std::stoi(value);
        getline(data, value, ',');  int e     = std::stoi(value);
        getline(data, value, ',');  int N     = std::stoi(value);
        getline(data, value, ',');  double f  = std::stod(value);
        getline(data, value, ',');  double Mz = std::stod(value);
        getline(data, value, ',');  double Mc = std::stod(value);

        if (cmds[0] == "mix") {
            if (Mz >= z_low && Mz <= z_high && Mc >= cen_low && Mc <= cen_high) {
                index.push_back(i);     evt.push_back(e);         NHits.push_back(N);
                foundZ.push_back(f);    MBD_true_z.push_back(Mz);   MBD_cen.push_back(Mc);
            }
        }
        else {
            // if (Mz >= -25. && Mz <= -5. && Mc <= 0.7) {
            // if (Mz >= -5. && Mz <= 5. && Mc <= 0.7) {
                index.push_back(i);     evt.push_back(e);         NHits.push_back(N);
                foundZ.push_back(f);    MBD_true_z.push_back(Mz);   MBD_cen.push_back(Mc);
            // }
        }
    }
    Long64_t nEntries = EventTree -> GetEntries();
    // int event25, NClus;
    // float MBD_centrality, MBD_z_vtx;
    // std::vector<float> *ClusX     = nullptr;
    // std::vector<float> *ClusY     = nullptr;
    // std::vector<float> *ClusZ     = nullptr;
    // std::vector<int>   *ClusLayer = nullptr;
    // std::vector<float> *ClusR     = nullptr;    // FYI only;
    // std::vector<float> *ClusPhi   = nullptr;    // FYI only;
    // std::vector<float> *ClusEta   = nullptr;    // FYI only;
    TObjArray* branches = EventTree->GetListOfBranches();
    std::string branchName;
    int idx_event, idx_MBD_z_vtx, idx_MBD_centrality, idx_NClus, idx_ClusLayer,
        idx_ClusX, idx_ClusY, idx_ClusZ, idx_ClusR, idx_ClusPhi, idx_ClusEta;
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch* branch = (TBranch*)branches->At(i);
        branchName = branch->GetName();
        if (branchName=="event")        idx_event = i;
        if (branchName=="MBD_z_vtx")    idx_MBD_z_vtx = i;
        if (branchName=="MBD_centrality")    idx_MBD_centrality = i;
        if (branchName=="NClus")        idx_NClus = i;
        if (branchName=="ClusLayer")    idx_ClusLayer = i;
        if (branchName=="ClusX")        idx_ClusX = i;
        if (branchName=="ClusY")        idx_ClusY = i;
        if (branchName=="ClusZ")        idx_ClusZ = i;
        if (branchName=="ClusR")        idx_ClusR = i;
        if (branchName=="ClusPhi")      idx_ClusPhi = i;
        if (branchName=="ClusEta")      idx_ClusEta = i;
    }
    TBranch *b_event     = (TBranch*)EventTree->GetListOfBranches()->At(idx_event);
    TBranch *b_MBD_z_vtx = (TBranch*)EventTree->GetListOfBranches()->At(idx_MBD_z_vtx);
    TBranch *b_MBD_centrality = (TBranch*)EventTree->GetListOfBranches()->At(idx_MBD_centrality);
    TBranch *b_NClus     = (TBranch*)EventTree->GetListOfBranches()->At(idx_NClus);
        TBranch *b_ClusLayer = (TBranch*)EventTree->GetListOfBranches()->At(idx_ClusLayer);
        TBranch *b_ClusX     = (TBranch*)EventTree->GetListOfBranches()->At(idx_ClusX);
        TBranch *b_ClusY     = (TBranch*)EventTree->GetListOfBranches()->At(idx_ClusY);
        TBranch *b_ClusZ     = (TBranch*)EventTree->GetListOfBranches()->At(idx_ClusZ);
        TBranch *b_ClusR     = (TBranch*)EventTree->GetListOfBranches()->At(idx_ClusR);
        TBranch *b_ClusPhi   = (TBranch*)EventTree->GetListOfBranches()->At(idx_ClusPhi);
        TBranch *b_ClusEta   = (TBranch*)EventTree->GetListOfBranches()->At(idx_ClusEta);
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

    // TBranch *branch25 = (TBranch*)EventTree->GetListOfBranches()->At(25);    // event25;
    // TBranch *branch10 = (TBranch*)EventTree->GetListOfBranches()->At(10);    // NClus;
    // TBranch *b_ClusLayer = (TBranch*)EventTree->GetListOfBranches()->At(11);    // ClusLayer;
    // TBranch *branch12 = (TBranch*)EventTree->GetListOfBranches()->At(12);    // ClusX;
    // TBranch *branch13 = (TBranch*)EventTree->GetListOfBranches()->At(13);    // ClusY;
    // TBranch *b_ClusZ = (TBranch*)EventTree->GetListOfBranches()->At(14);    // ClusZ;
    // TBranch *branch30 = (TBranch*)EventTree->GetListOfBranches()->At(30);    // MBD_centrality;
    // TBranch *b_MBD_z_vtx = (TBranch*)EventTree->GetListOfBranches()->At(31);    // MBD_z_vtx;
    // TBranch *b_ClusR = (TBranch*)EventTree->GetListOfBranches()->At(15);    // ClusR;
    // TBranch *b_ClusPhi = (TBranch*)EventTree->GetListOfBranches()->At(16);    // ClusPhi;
    // TBranch *branch17 = (TBranch*)EventTree->GetListOfBranches()->At(17);    // ClusEta;
    
    // branch25->SetAddress(&event25);
    // branch10->SetAddress(&NClus);
    // b_ClusLayer->SetAddress(&ClusLayer);
    // branch12->SetAddress(&ClusX);
    // branch13->SetAddress(&ClusY);
    // b_ClusZ->SetAddress(&ClusZ);
    // b_ClusR->SetAddress(&ClusR);
    // b_ClusPhi->SetAddress(&ClusPhi);
    // branch17->SetAddress(&ClusEta);
    // branch30->SetAddress(&MBD_centrality);
    // b_MBD_z_vtx->SetAddress(&MBD_z_vtx);

    std::cout << "Total size: " << evt.size() << std::endl;
    target = target > evt.size() ? evt.size() : target;
    std::cout << "Selected size: " << target << std::endl;

    std::cout << target << "," << cen_low << "," << cen_high << "," << z_low << "," << z_high << std::endl;
    
    ROOT::EnableThreadSafety();
    // Use 'sysctl -n hw.logicalcpu' to determine max number of threads:
    if (method[0] == "nomix") {
        std::thread thsafe[8];
        std::cout<<"multi-thready safe:"<<std::endl;
        // for(int i = 0; i < 8; ++i)     thsafe[i]= std::thread(dPhiAccumulator,i,target/8,std::cref(index),std::cref(MBD_cen),b_ClusLayer,b_ClusPhi,ClusPhi,ClusLayer);
        for(int i = 0; i < 8; ++i)     thsafe[i]= std::thread(dPhiNoMixWithEta,i,target/8,std::cref(index),std::cref(MBD_true_z),b_ClusLayer,b_ClusZ,b_ClusR,b_ClusPhi,ClusZ,ClusR,ClusPhi,ClusLayer);
        for(int i = 0; i < 8; ++i)     thsafe[i].join();

        h_dPhi_nomix -> SetTitle(Form("dPhi of %d un-mixed events, with dEta cut %0.2f", target, dEta_cut));
        angularPlot1D(h_dPhi_nomix, method, "dPhi of unmixed with dEta cut");
        // h_dPhi_nomix -> SetTitle(Form("dPhi of %d un-mixed events", target));
        // angularPlot1D(h_dPhi_nomix, method, "dPhi of unmixed");
        TFile* file = TFile::Open("dPhi.root", "UPDATE");
        if (!file || file -> IsZombie()) {
            std::cerr << "You messed up!" << std::endl;
            exit(1);
        }
        h_dPhi_nomix -> Write("nomix", TObject::kOverwrite);
    }
    else if (method[0] == "perCen") {
        std::thread thsafe[14];
        std::cout << "safe dPhi of different centralities" << std::endl;
        for (int i = 0; i < 14; i++)
            thsafe[i] = std::thread(dPhi_in_bins_of_Centrality_with_dEta_cut,i,target,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),b_ClusLayer,b_ClusZ,b_ClusR,b_ClusPhi,b_MBD_z_vtx,std::ref(MBD_z_vtx),ClusLayer,ClusZ,ClusR,ClusPhi);
            // thsafe[i] = std::thread(dPhi_in_bins_of_Centrality,i,target,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),b_ClusLayer,b_ClusPhi,ClusPhi,ClusLayer);

        for (int i = 0; i < 14; i++)
            thsafe[i].join();

        std::vector<TH1D*> h(h_CenonOne, h_CenonOne + 14);
        // ArrayPlot1D_Rescale(h, method, "dPhi_per_centralities_rescale");
        ArrayPlot1D_Rescale(h, method, Form("dPhi_per_centralities_rescale_with_dEta_cut_%d_events",target));
    }
    else if (method[0] == "perZ") {
        std::thread thsafe[20];
        std::cout << "safe dPhi of different Z vertices" << std::endl;
        for (int i = 0; i < 20; i++)
            thsafe[i] = std::thread(dPhi_in_bins_of_Z_vtx_with_dEta_cut,i,target,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),b_ClusLayer,b_ClusZ,b_ClusR,b_ClusPhi,ClusLayer,ClusZ,ClusR,ClusPhi);
            // thsafe[i] = std::thread(dPhi_in_bins_of_Z_vtx,i,target,std::cref(index),std::cref(MBD_true_z),std::cref(MBD_cen),b_ClusLayer,b_ClusPhi,ClusPhi,ClusLayer);

        for (int i = 0; i < 20; i++)
            thsafe[i].join();

        std::vector<TH1D*> h(h_ZonOne, h_ZonOne + 20);
        ArrayPlot1D_Rescale_ver2(h, method, Form("dPhi_per_Z_vtx_rescale_%d_events", target));
        // ArrayPlot1D_Rescale_ver2(h, method, "dPhi_per_Z_vtx_rescale_with_dEta_cut");
    }
    // else if (method[0] == "mix") {
    //     double phi, dPhi, phi_0, phi_1; int clus_layer;
    //     std::vector<double> Phi0, Phi1;
    //     std::vector<std::vector <double>> event_Phi0, event_Phi1;
    //     // Signal dPhi is unmixed events' dPhi:
    //     TH1D *h_Signal_dPhi = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %d events;dPhi value;# of counts", target), N, range_min, range_max);
    //     // Fill Signals:
    //     for (int i = 0; i < target; i++) {
    //         b_ClusPhi->GetEntry(index[i]);   b_ClusLayer->GetEntry(index[i]);
    //         event_Phi0.push_back(std::vector <double>());   event_Phi1.push_back(std::vector <double>());
    //         for (int j = 0; j < ClusPhi->size(); j++) {
    //             phi        = ClusPhi->at(j);
    //             clus_layer = ClusLayer->at(j);
    //             if (clus_layer == 3 || clus_layer == 4) {
    //                 Phi0.push_back(phi);
    //                 event_Phi0[i].push_back(phi);
    //             }
    //             else {
    //                 Phi1.push_back(phi);
    //                 event_Phi1[i].push_back(phi);
    //             }
    //         }
    //         for (int k = 0; k < Phi0.size(); k++) {
    //             for (int l = 0; l < Phi1.size(); l++) {
    //                 dPhi = Phi0[k] - Phi1[l];
    //                 if (dPhi > M_PI)    dPhi = dPhi - TWO_PI;
    //                 if (dPhi < -M_PI)   dPhi = dPhi + TWO_PI;
    //                 h_Signal_dPhi -> Fill(dPhi);
    //             }
    //         }
    //         Phi0.clear();   Phi1.clear();
    //     }

    //     std::thread thsafe[8];
    //     std::cout << "Mixing events: \n" << std::endl;
    //     for (int i = 0; i < 8; i++)
    //         thsafe[i] = std::thread(dPhi_mixing,i,target/8,event_Phi0,event_Phi1);
    //     for (int i = 0; i < 8; i++)
    //         thsafe[i].join();


    //     TCanvas *c1 = new TCanvas("c1", "dPhi Histogram", 1920, 1056);
    //     double phi_range_low = -2.4, phi_range_high = -1.8;
    //     int bin_range_low = h_Background_dPhi->FindBin(phi_range_low), bin_range_high = h_Background_dPhi->FindBin(phi_range_high);
    //     double max_unmixed = -1, max_mixed = -1, current_binContent;
    //     for (int bin = bin_range_low; bin <= bin_range_high; bin++) {
    //         current_binContent = h_Signal_dPhi->GetBinContent(bin);
    //         if (max_unmixed < current_binContent)  max_unmixed = current_binContent;
    //         current_binContent = h_Background_dPhi->GetBinContent(bin);
    //         if (max_mixed < current_binContent)  max_mixed = current_binContent;
    //     }
    //     h_Signal_dPhi -> Add(h_Signal_dPhi, max_mixed/max_unmixed - 1);
    //     h_Signal_dPhi -> Draw("SAME");

    //     max_mixed   = h_Background_dPhi->GetBinContent(h_Background_dPhi->GetMaximumBin());
    //     max_unmixed = h_Signal_dPhi->GetBinContent(h_Signal_dPhi->GetMaximumBin());
    //     std::cout << h_Signal_dPhi->GetBinContent(h_Signal_dPhi->GetMaximumBin()) << ", " << h_Background_dPhi->GetBinContent(h_Background_dPhi->GetMaximumBin())  << std::endl;
    //     if (max_unmixed > max_mixed) {
    //         h_Signal_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_unmixed*1.2);
    //         h_Background_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_unmixed*1.2);
    //     }   
    //     else {
    //         h_Signal_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_mixed*1.2);
    //         h_Background_dPhi -> GetYaxis() -> SetRangeUser(1e7, max_mixed*1.2);
    //     }                      
    //     h_Background_dPhi -> Draw("SAME");
    //     h_Background_dPhi -> SetLineColor(2);

    //     double pi = TMath::Pi();
    //     int bin_min  = 1;  // The first bin
    //     int bin_max  = h_Signal_dPhi->GetNbinsX();  // The last bin
    //     // Calculate bin positions for each label
    //     int bin_pi   = bin_max;
    //     int bin_0    = bin_min + (bin_max - bin_min)/2;
    //     int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    //     int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    //     // Set the labels at the calculated positions
    //     h_Background_dPhi->GetXaxis()->SetBinLabel(bin_0, "0");
    //     h_Background_dPhi->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    //     h_Background_dPhi->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    //     h_Background_dPhi->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    //     h_Background_dPhi->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    //     // Ensure the custom labels are displayed by setting the number of divisions
    //     h_Background_dPhi->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    //     h_Background_dPhi->GetXaxis()->SetLabelSize(0.04);
    //     // Update histogram to refresh the axis
    //     // h[0]->Draw("HIST");
    //     h_Background_dPhi->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    //     h_Signal_dPhi->GetXaxis()->SetBinLabel(bin_0, "0");
    //     h_Signal_dPhi->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    //     h_Signal_dPhi->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    //     h_Signal_dPhi->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    //     h_Signal_dPhi->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    //     // Ensure the custom labels are displayed by setting the number of divisions
    //     h_Signal_dPhi->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    //     h_Signal_dPhi->GetXaxis()->SetLabelSize(0.04);
    //     // Update histogram to refresh the axis
    //     // h[0]->Draw("HIST");
    //     h_Signal_dPhi->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    //     h_Signal_dPhi -> GetXaxis() -> CenterTitle(true);
    //     h_Signal_dPhi -> GetYaxis() -> CenterTitle(true);
    //     h_Background_dPhi -> GetXaxis() -> CenterTitle(true);
    //     h_Background_dPhi -> GetYaxis() -> CenterTitle(true);

    //     // TLegend *lg = new TLegend(0.12, 0.8, 0.33, 0.9);
    //     // lg -> AddEntry(h_Background_dPhi, Form("z_vtx between %2.2f and %2.2f", z_lower_range, z_upper_range), "l");
    //     // gStyle -> SetLegendTextSize(.023);
    //     // lg->Draw("same");
    
    //     c1 -> SaveAs("../../External/zFindingPlots/dPhi_mixed.png");
    //     backgroundCancelling_dPhi(h_Background_dPhi, h_Signal_dPhi, method, target);
    // }
    
    if (cmds[0] == "mix") {
        if (cmds[1] == "wo") {
            double phi, dPhi; 
            int clus_layer;
            std::vector<double> Phi0, Phi1;
            std::vector<std::vector <double>> event_Phi0, event_Phi1;
            // Signal dPhi is unmixed events' dPhi:
            TH1D *h_Signal_dPhi = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %d events;dPhi value;# of counts", target), N, range_min, range_max);
            for (int i = 0; i < target; i++) {  // Loop over events;
                b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
                b_ClusPhi->GetEntry(index[i]);   // ClusPhi;
                event_Phi0.push_back(std::vector <double>());   
                event_Phi1.push_back(std::vector <double>());
                for (int j = 0; j < ClusPhi->size(); j++) { // Loop inside one event, over all hits;
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
            std::thread thsafe[8];
            std::cout << "Mixing events without dEta cut: \n" << std::endl;
            for (int i = 0; i < 8; i++)
                thsafe[i] = std::thread(dPhi_mixing,i,target/8,event_Phi0,event_Phi1);
            for (int i = 0; i < 8; i++)
                thsafe[i].join();

            backgroundCancelling_dPhi(h_Background_dPhi, h_Signal_dPhi, method, target);
        } else if (cmds[1] == "wdE") {
            double phi, dPhi, eta, dEta;
            double z_vtx, dZ, R, theta;
            int clus_layer;
            std::vector<EtaWithPhi> Phi0, Phi1;
            std::vector<std::vector <EtaWithPhi>> event_Phi0, event_Phi1;
            // Signal dPhi is unmixed events' dPhi:
            TH1D *h_Signal_dPhi = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %d events;dPhi value;# of counts", target), N, range_min, range_max);
            for (int i = 0; i < target; i++) {  // Loop over events;
                b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
                b_ClusZ->GetEntry(index[i]);   // ClusZ;
                b_ClusR->GetEntry(index[i]);   // ClusR;
                b_ClusPhi->GetEntry(index[i]);   // ClusPhi;
                event_Phi0.emplace_back(std::vector <EtaWithPhi>());   
                event_Phi1.emplace_back(std::vector <EtaWithPhi>());
                z_vtx = MBD_true_z[i];
                for (int j = 0; j < ClusPhi->size(); j++) { // Loop inside one event, over all hits;
                    dZ         = ClusZ->at(j) - z_vtx;
                    R          = ClusR->at(j);
                    theta      = std::atan2(R, dZ);
                    phi        = ClusPhi->at(j);
                    clus_layer = ClusLayer->at(j);
                    if (dZ >= 0)    eta = -std::log(std::tan(theta/2));
                    if (dZ <  0)    eta = std::log(std::tan((M_PI - theta)/2));
                    if (clus_layer == 3 || clus_layer == 4) {
                        Phi0.emplace_back(eta, phi);
                        event_Phi0[i].emplace_back(eta, phi);
                    }
                    else {
                        Phi1.emplace_back(eta, phi);
                        event_Phi1[i].emplace_back(eta, phi);
                    }
                }

                for (int k = 0; k < Phi0.size(); k++) {
                    for (int l = 0; l < Phi1.size(); l++) {
                        dEta = Phi0[k].eta_value - Phi1[l].eta_value;
                        dPhi = Phi0[k].phi_value - Phi1[l].phi_value;
                        if (dPhi > M_PI)    dPhi -= TWO_PI;
                        if (dPhi < -M_PI)   dPhi += TWO_PI;
                        if (std::abs(dEta) < dEta_cut)  h_Signal_dPhi -> Fill(dPhi);
                    }
                }
                Phi0.clear();   Phi1.clear();
            }

            // std::vector<int> boundaries = readCsvToVector("../../codeGarage/boundaries.csv");
            // std::vector<int> boundaries = readCsvToVector("/Users/yaminocellist/codeGarage/boundaries.csv");

            std::thread thsafe[8];
            std::cout << "Mixing events with dEta cut: \n" << std::endl;
            for (int i = 0; i < 8; i++)
                thsafe[i] = std::thread(dPhi_mixing_with_dEta_cut,i,target/8,event_Phi0,event_Phi1);
                // thsafe[i] = std::thread(dPhi_mixing_with_dEta_cut_ver2,boundaries[2*i],boundaries[2*i+1],event_Phi0,event_Phi1);
            for (int i = 0; i < 8; i++)
                thsafe[i].join();

            backgroundCancelling_dPhi(h_Background_dPhi, h_Signal_dPhi, method, target);
        } else if (cmds[1] == "wE") {
            double phi, dPhi, eta, abs_eta, dEta, eta1, eta2;
            double z_vtx, dZ, R, halfTheta;
            int clus_layer;
            std::vector<EtaWithPhi> Phi0, Phi1;
            std::vector<std::vector <EtaWithPhi>> event_Phi0, event_Phi1;
            // Signal dPhi is unmixed events' dPhi:
            TH1D *h_Signal_dPhi = new TH1D("dPhi of unmixed", Form("dPhi of unmixed %d events;dPhi value;# of counts", target), N, range_min, range_max);
            for (int i = 0; i < target; i++) {  // Loop over events;
                b_ClusLayer->GetEntry(index[i]);   // ClusLayer;
                b_ClusZ->GetEntry(index[i]);   // ClusZ;
                b_ClusR->GetEntry(index[i]);   // ClusR;
                b_ClusPhi->GetEntry(index[i]);   // ClusPhi;
                event_Phi0.emplace_back(std::vector <EtaWithPhi>());   
                event_Phi1.emplace_back(std::vector <EtaWithPhi>());
                z_vtx = MBD_true_z[i];
                for (int j = 0; j < ClusPhi->size(); j++) { // Loop inside one event, over all hits;
                    dZ         = ClusZ->at(j) - z_vtx;
                    R          = ClusR->at(j);
                    halfTheta      = (std::atan2(R, dZ))/2.;
                    phi        = ClusPhi->at(j);
                    clus_layer = ClusLayer->at(j);
                    if (dZ >= 0)    eta = -std::log(std::tan(halfTheta));
                    if (dZ <  0)    eta = std::log(std::tan(halfPI - halfTheta));
                    abs_eta = std::abs(eta);
                    if (abs_eta < Eta_range) {
                        if (clus_layer == 3 || clus_layer == 4) {
                            Phi0.emplace_back(eta, phi);
                            event_Phi0[i].emplace_back(eta, phi);
                        }
                        else {
                            Phi1.emplace_back(eta, phi);
                            event_Phi1[i].emplace_back(eta, phi);
                        }
                    }
                }

                for (int k = 0; k < Phi0.size(); k++) {
                    for (int l = 0; l < Phi1.size(); l++) {
                        dPhi = Phi0[k].phi_value - Phi1[l].phi_value;
                        if (dPhi > M_PI)    dPhi -= TWO_PI;
                        if (dPhi < -M_PI)   dPhi += TWO_PI;
                        h_Signal_dPhi -> Fill(dPhi);
                    }
                }
                Phi0.clear();   Phi1.clear();
            }

            std::thread thsafe[8];
            std::cout << "Mixing events with Eta range: " << Eta_range << std::endl;
            for (int i = 0; i < 8; i++)
                thsafe[i] = std::thread(dPhi_mixing_with_Eta_range,i,target/8,event_Phi0,event_Phi1);
            for (int i = 0; i < 8; i++)
                thsafe[i].join();

            backgroundCancelling_dPhi(h_Background_dPhi, h_Signal_dPhi, method, target);
        }
    }

    // Stop the stopwatch and print the runtime:
    timer.Stop();   timer.Print();

    theApp.Run();
    return 0;
}
