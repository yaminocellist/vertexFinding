#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H
#include <vector>
#include <cmath>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <functional>
#include <chrono>
#include <filesystem>

#include <TH1.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TPolyMarker.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TFile.h>
#include <TTree.h>
#include <TMultiGraph.h>
#include <TProfile.h>
#include <TParameter.h>
#include "ROOT/TSeq.hxx"

/*********************************************************************
 *                      GLOBAL VARIABLES;
 * ******************************************************************/
double zmin = -45;
double zmax = 35;        // From sPHENIX paper, the stave's length is aroung 27.12 cm;
double scanstep = 0.2;  // unit: cm; 
int bins = (zmax - zmin)/scanstep + 1;
const double dPhi_cut  = 0.01;
const double dEta_cut  = 1;
const double Eta_range = 1;
const double TWO_PI    = 2*M_PI;
const double halfPI    = M_PI/2;
double DCA_cut        = 0.2;    // unit: cm;
double DCA_cutSQUARED = 0.04;   // unit: cm;
double MBD_lower = 0., MBD_upper = 10.;
double abs_fit_range = M_PI/60;
const UInt_t nThreads = 8U;
// double abs_fit_range = 0.05;

 // ANSI escape code for red text
const std::string COLOR_GREEN="\033[0;32m";
const std::string COLOR_RED="\033[0;31m";
const std::string COLOR_YELLOW="\033[0;33m";
const std::string COLOR_OCHRE="\033[38;5;95m";
const std::string COLOR_BLUE="\033[0;34m";
const std::string COLOR_WHITE="\033[0;37m";
const std::string COLOR_RESET="\033[0m";

/*  Genearting unequal bin ranges for a historgram  */
// Double_t binEdges[481];
// double dx_1 = 0.01;
// binEdges[0] = -2;
// for (int i = 1; i <= 190; i++) {
//     binEdges[i] = binEdges[0] + i*dx_1;
// }
// double dx_2 = 0.002;
// for (int j = 1; j <= 100; j++) {
//     binEdges[j+190] = binEdges[190] + j*dx_2;
// }
// for (int l = 1; l <= 190; l++) {
//     binEdges[l + 290] = binEdges[290] + l*dx_1;
// }

//////////////////////////////////////////////////////////////////////

struct myPoint3D {
    double x, y, z;
};

struct myTrackletMemberLite {
    Double_t x, y, z;
    Double_t phi;
};

struct myTrackletMember {
    Double_t x, y, z, r;
    Double_t eta, phi;
    Int_t layer;
};

struct myTrackletMemberExtended : myTrackletMember {
    Int_t trackID;
};

struct EtaWithPhi {
    double eta_value;
    double phi_value;

    // You can add constructors, methods, etc., to enhance functionality
    EtaWithPhi(double e, double p) : eta_value(e), phi_value(p) {}
};

namespace rootBranchEssentials {
    inline int event, NClus;
    inline float MBD_z_vtx, MBD_centrality, MBD_charge_sum;
    inline std::vector<int>   *ClusLayer = nullptr;
    inline std::vector<float> *ClusX    = nullptr;
    inline std::vector<float> *ClusY    = nullptr;
    inline std::vector<float> *ClusZ    = nullptr;
    inline std::vector<float> *ClusR    = nullptr;
    inline std::vector<float> *ClusPhi  = nullptr;
    inline std::vector<float> *ClusEta  = nullptr;

    inline int idx_event, idx_MBD_z_vtx, idx_MBD_centrality, idx_MBD_charge_sum, idx_NClus, idx_ClusLayer,
               idx_ClusX, idx_ClusY, idx_ClusZ, idx_ClusR, idx_ClusPhi, idx_ClusEta;

    inline TBranch *b_event, *b_MBD_z_vtx, *b_MBD_centrality, *b_MBD_charge_sum, *b_NClus, *b_ClusLayer, *b_ClusX, *b_ClusY, *b_ClusZ, *b_ClusR, *b_ClusPhi, *b_ClusEta;
}

bool isInteger(const std::string& s) {
    std::istringstream iss(s);
    int x;
    return (iss >> x) && (iss.eof());
}

std::pair<double, double> nearestZ (const myPoint3D &p1, const myPoint3D &p2) {
    double numeratorZ = (p1.y - p2.y)*(p1.y*p2.z - p2.y*p1.z) - (p1.x - p2.x)*(p2.x*p1.z - p1.x*p2.z);
    double denominatorZ = (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
    double nearest_z = numeratorZ/denominatorZ;

    double part123 = (p1.x*p1.x + p2.x*p2.x + p1.y*p1.y + p2.y*p2.y - 2*p1.x*p2.x - 2*p1.y*p2.y)*nearest_z*nearest_z +
                    2*((p1.x*p2.x + p1.y*p2.y - p2.x*p2.x - p2.y*p2.y)*p1.z + (p1.x*p2.x + p1.y*p2.y - p1.x*p1.x - p1.y*p1.y)*p2.z)*nearest_z +
                    p1.y*p1.y*p2.z*p2.z + p2.y*p2.y*p1.z*p1.z + p1.x*p1.x*p2.z*p2.z + p2.x*p2.x*p1.z*p1.z + p1.x*p1.x*p2.y*p2.y + p2.x*p2.x*p1.y*p1.y -
                    2*p1.x*p2.x*p1.y*p2.y - 2*p1.y*p2.y*p1.z*p2.z - 2*p1.x*p2.x*p1.z*p2.z;

    double denominator = sqrt(p1.x*p1.x + p2.x*p2.x + p1.y*p1.y + p2.y*p2.y + p1.z*p1.z + p2.z*p2.z - 2*p1.x*p2.x - 2*p1.y*p2.y - 2*p1.z*p2.z);

    return std::make_pair(nearest_z, sqrt(part123)/denominator);
}

std::pair<double, double> nearestZ_revised (const myPoint3D &p1, const myPoint3D &p2) {
    double x1 = p1.x, y1 = p1.y, z1 = p1.z;
    double dx = p2.x - x1, dy = p2.y - y1, dz = p2.z - z1;
    double bPlusde            = 2*(x1*dx + y1*dy);
    double fourAMinusdSQUARED = 4*(dx*dx + dy*dy); 

    double z_m = z1 + (-2)*dz*bPlusde/fourAMinusdSQUARED;
    // double dist_m = std::sqrt(x1*x1 + y1*y1 - bPlusde*bPlusde/fourAMinusdSQUARED);
    double dist_mSQUARED = x1*x1 + y1*y1 - bPlusde*bPlusde/fourAMinusdSQUARED;

    return std::make_pair(z_m, dist_mSQUARED);
}

std::vector<std::string> splitString(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);

    while (std::getline(tokenStream, token, delim)) {
        tokens.push_back(token);
    }

    return tokens;
}

/**
 * @brief 
 * 
 * @return ** void 
 */

void current_PC_time () {
    auto now   = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    std::cout << "Run started at: " << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S") << std::endl;
}

// template <typename T>
// T readCsvToVector(const std::string& filename) {
//     std::vector<int> result;
//     std::ifstream file(filename);
//     if (!file.is_open()) {
//         std::cerr << "Failed to open file: " << filename << std::endl;
//         return result;
//     }
    
//     std::string line;
//     if (std::getline(file, line)) {
//         std::stringstream ss(line);
//         std::string value;
        
//         while (std::getline(ss, value, ',')) {
//             result.push_back(std::stoi(value));
//         }
//     }
    
//     file.close();
//     return result;
// }

// Function to read a CSV file and return a vector of vectors containing the data
template <typename T>
T readCsvToVector(const std::string& filename) {
    T result;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return result;
    }
    
    std::string line;
    bool isHeader = true; // Skip the first row if it contains headers
    while (std::getline(file, line)) {
        if (isHeader) {
            isHeader = false;
            // If headers need to be processed, you can handle them here
            continue;
        }
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row;
        
        while (std::getline(ss, value, ',')) {
            try {
                row.push_back(std::stod(value)); // Convert value to double
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid value: " << value << " in file " << filename << std::endl;
                row.push_back(0.0); // Optional: Handle invalid values gracefully
            }
        }
        result.push_back(row);
    }
    
    file.close();
    return result;
}

// Helper function to print the read contents of .csv
void printCsvData(const std::vector<std::vector<double>>& data) {
    for (const auto& row : data) {
        for (const auto& value : row) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
}

std::string findRootFile (const std::string &directory, const std::string &prefix) {
    for (const auto &entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            std::string fileName = entry.path().filename().string();
            if (fileName.find(prefix) == 0 && fileName.find(".root") != std::string::npos) {
                return entry.path().string();
            }
        }
    }
    return "";
}

std::vector<std::string> findRootFiles (const std::string &directory, const std::string &prefix) {
    std::vector<std::string> matchingFiles;
    for (const auto &entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            std::string fileName = entry.path().filename().string();
            if (fileName.find(prefix) == 0 && fileName.find(".root") != std::string::npos) {
                matchingFiles.push_back(entry.path().string());
            }
        }
    }
    std::sort(matchingFiles.begin(), matchingFiles.end());
    return matchingFiles;
}

template <typename T>
void printRedSingle(const T &content) {
    std::cout << COLOR_RED << content << COLOR_RESET;
}
template <typename First, typename... Rest>
void printRed(const First& first, const Rest&... rest){
    printRedSingle(first);
    if constexpr (sizeof...(rest) > 0) {
        printRedSingle("");
        printRed(rest...);
    }
    else
        std::cout << std::endl;
}

template <typename T>
void printWhite(const T &content) {
    std::cout << COLOR_WHITE << content << COLOR_RESET << std::endl;
}

template <typename T>
void printBlueSingle(const T &content) {
    std::cout << COLOR_BLUE << content << COLOR_RESET;
}
template <typename First, typename... Rest>
void printBlue(const First& first, const Rest&... rest){
    printBlueSingle(first);
    if constexpr (sizeof...(rest) > 0) {
        printBlueSingle("");
        printBlue(rest...);
    }
    else
        std::cout << std::endl;
}

void printSeparation () {
    printWhite("=================================================================================================================");
}

template <typename T>
void fileExistenceCheck (const T &file) {
    if (!file.is_open()){
		std::cout << "Unable to open linelabel" << std::endl;
		system("read -n 1 -s -p \"Press any key to continue...\" echo");
		exit(1);
 	}
}

template <typename T>
void objectExistenceCheck (const T &file) {
    if (!file) {
        std::cerr << "Error: Could not open/find the object " << std::endl;
        exit(1);
    }
}

#endif