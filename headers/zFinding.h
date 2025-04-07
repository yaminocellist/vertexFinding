#include <iomanip>
#include "globalDefinitions.h"
#include <TSpectrum.h>
#include <TLatex.h>

// #define __PEAKS_C_FIT_AREAS__ 1 /* fit peaks' areas */

Int_t npeaks = 10;
Double_t fpeaks(Double_t *x, Double_t *par) {
   Double_t result = par[0] + par[1]*x[0];
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm  = par[3*p+2]; // "height" or "area"
      Double_t mean  = par[3*p+3];
      Double_t sigma = par[3*p+4];
      
      #if defined(__PEAKS_C_FIT_AREAS__)
        norm /= sigma * (TMath::Sqrt(TMath::TwoPi())); // "area"
      #endif /* defined(__PEAKS_C_FIT_AREAS__) */
      
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}

// Original DCA (Distance of the Closest Approach) method:
// The histogram is processed directly after Fill, without filtering its background;
double nearest_z_method (const int &evt, const std::vector<myTrackletMemberLite> &t0, const std::vector<myTrackletMemberLite> &t1,
                         const double &eta_cut_low, const double &eta_cut_high,
                         const double &phi_cut_low, const double &phi_cut_high, const double &trueZ) {
    // TH1D *h = new TH1D("", "", bins, zmin - scanstep/2, zmax + scanstep/2);
    auto h = std::make_unique<TH1D>("", "", bins, zmin - scanstep/2, zmax + scanstep/2);
    h -> Sumw2(kFALSE);     // Disable the storage of the sum of squares of weights;
                            // reduce both memory and run time;
    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t1.size(); j++) {
            double dPhi = t0[i].phi - t1[j].phi;
            if (dPhi >= phi_cut_low && dPhi <= phi_cut_high) {
                myPoint3D P1 = {t0[i].x, t0[i].y, t0[i].z};
                myPoint3D P2 = {t1[j].x, t1[j].y, t1[j].z};
                double foundZ = nearestZ(P1, P2).first;
                double nearest_d = nearestZ(P1, P2).second;
                if (foundZ >= zmin && foundZ <= zmax && nearest_d <= DCA_cut) {
                    h -> Fill(foundZ);
                }
            }
        }
    }

    // TFitResultPtr fitResult = h->Fit("gaus", "SQ");  // Get the fit result
    // if (fitResult && fitResult->Status() == 0) {
    //     double mean = fitResult->Parameter(1);  // Get the mean
    //     double sigma = fitResult->Parameter(2);  // Get the standard deviation
    // }
    
    int binmax      = h -> GetMaximumBin(); 
    double xlabel0  = h -> GetXaxis() -> GetBinCenter(binmax);      double entry0  = h -> GetBinContent(binmax);
    double xlabel_1 = h -> GetXaxis() -> GetBinCenter(binmax - 1);  double entry_1 = h -> GetBinContent(binmax - 1);
    double xlabel1  = h -> GetXaxis() -> GetBinCenter(binmax + 1);  double entry1  = h -> GetBinContent(binmax + 1);
    double ctz      = (xlabel0*entry0 + xlabel_1*entry_1 + xlabel1*entry1) / (entry0 + entry1 + entry_1);
    // double ctz = mean;
// /*
    TCanvas *can2 = new TCanvas("c2","c2",0,50,1800,550);
    h -> Draw();
    h -> SetFillColor(kYellow - 7);
    h -> SetLineWidth(1);
    h -> SetFillStyle(1001);
    h -> GetXaxis() -> SetTitle("Z coordinate [cm]");
    h -> GetXaxis() -> SetTitleSize(.05);
    h -> GetXaxis() -> SetLabelSize(.04);
    h -> GetXaxis() -> CenterTitle(true);
    h -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h -> GetXaxis() -> SetTitleOffset(.8);
    h -> GetYaxis() -> SetTitle("# of Counts");
    h -> GetYaxis() -> SetTitleSize(.05);
    h -> GetYaxis() -> SetLabelSize(.04);
    h -> GetYaxis() -> SetTitleOffset(.62);
    h -> GetYaxis() -> CenterTitle(true);
    h -> SetTitle(Form("Found Z of Event %d, #bf{Nearest Z Method}", evt));
    h -> SetMinimum(0);
    int max_entry = h -> GetBinContent(h -> FindBin(ctz));
    // TLine *l1 = new TLine(ctz, 0, ctz, max_entry);
    static TLine *l1 = new TLine();
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(2);
    l1 -> SetLineWidth(4);
    // l1 -> Draw("same");
    // TLine *l2 = new TLine(trueZ, 0, trueZ, max_entry);
    static TLine *l2 = new TLine();
    l2 -> SetLineColor(kBlue);
    l2 -> SetLineStyle(1);
    l2 -> SetLineWidth(4);
    // l2 -> Draw("same");
    // TLegend *lg = new TLegend(0.12, 0.8, 0.42, 0.9);
    static TLegend *lg = new TLegend(0.12, 0.8, 0.42, 0.9);
    lg -> AddEntry(h.get(), Form("%lu hits, found z = %fcm, true z = %fcm", t0.size() + t1.size(), ctz, trueZ), "f");
    lg -> SetTextSize(.03);
    // Update line positions and redraw for the current event
    l1->SetX1(ctz); l1->SetX2(ctz);
    l1->SetY1(0);   l1->SetY2(max_entry);
    l1->Draw("same");

    l2->SetX1(trueZ); l2->SetX2(trueZ);
    l2->SetY1(0);     l2->SetY2(max_entry);
    l2->Draw("same");

    // Clear and update the legend for the current event
    lg->Clear();
    lg->AddEntry(h.get(), Form("%lu hits, found z = %fcm, true z = %fcm", t0.size() + t1.size(), ctz, trueZ), "f");
    lg->Draw("same");

    // lg -> Draw("same");
    gPad -> SetGrid(1,1); gPad -> Update();
    // gPad -> SetLogy();
    can2 -> SaveAs(Form("../External/zFindingPlots/nearest_foundz_doublebins_%d.png", evt));
    delete can2;
    // delete h;   h = 0; 
    // delete l1; delete l2; delete lg;
// */
    return ctz;
}

//  DCA (Distance of Closest Approach), integrating npeaks.C tutorial:
double DCA_npeaks_fitLite (const int &evt, const std::vector<myTrackletMemberLite> &t0, const std::vector<myTrackletMemberLite> &t1,
                         const double &eta_cut_low, const double &eta_cut_high,
                         const double &phi_cut_low, const double &phi_cut_high, const double &trueZ) {
    // auto h = std::make_unique<TH1D>("", "", bins, zmin - scanstep/2, zmax + scanstep/2);
    TH1D *h = new TH1D("", Form("Found Z of Event %d, NHits = %lu, #bf{DCA with fit};z position [cm];# of counts", evt, (t0.size()+t1.size())/33), bins, zmin - scanstep/2, zmax + scanstep/2);
    // Disable the storage of the sum of squares of weights; 
    // Save both memory and run time;
    h->Sumw2(kFALSE);
    std::pair<double, double> nearestZResults_revised;
    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t1.size(); j++) {
            double dPhi = t0[i].phi - t1[j].phi;
            if (dPhi >= phi_cut_low && dPhi <= phi_cut_high) {
                myPoint3D P1 = {t0[i].x, t0[i].y, t0[i].z};
                myPoint3D P2 = {t1[j].x, t1[j].y, t1[j].z};
                nearestZResults_revised         = nearestZ_revised(P1, P2);
                double foundZ_revised           = nearestZResults_revised.first;
                double nearest_d_revisedSQUARED = nearestZResults_revised.second;
                if (foundZ_revised >= zmin && foundZ_revised <= zmax && nearest_d_revisedSQUARED <= DCA_cutSQUARED) {
                    h -> Fill(foundZ_revised);
                }
            }
        }
    }

    TSpectrum *s = new TSpectrum();
    TH1 *bg = s -> Background(h, 40, "nosmoothing"); 
    bg -> SetLineColor(kRed);
    TCanvas *c = new TCanvas("c", "c", 0,50,2000,1000);
    h->Add(bg, -1.);
    Int_t nfound = s->Search(h, 7., "", 0.001);
    Double_t *xpeaks;
    xpeaks = s->GetPositionX();
    double ctz = xpeaks[0];
    
    delete c;   delete h;   delete bg;
    h = 0;  bg = 0;
    return ctz;
}

double DCA_npeaks_fit (const int &evt, const std::vector<myTrackletMemberLite> &t0, const std::vector<myTrackletMemberLite> &t1,
                         const double &eta_cut_low, const double &eta_cut_high,
                         const double &phi_cut_low, const double &phi_cut_high, const double &trueZ) {
    // auto h = std::make_unique<TH1D>("", "", bins, zmin - scanstep/2, zmax + scanstep/2);
    TH1D *h = new TH1D("", Form("Found Z of Event %d, NHits = %lu, #bf{DCA with fit};z position [cm];# of counts", evt, (t0.size()+t1.size())/33), bins, zmin - scanstep/2, zmax + scanstep/2);
    // Disable the storage of the sum of squares of weights; 
    // Save both memory and run time;
    h->Sumw2(kFALSE);
    std::pair<double, double> nearestZResults, nearestZResults_revised;
    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t1.size(); j++) {
            double dPhi = t0[i].phi - t1[j].phi;
            if (dPhi >= phi_cut_low && dPhi <= phi_cut_high) {
                myPoint3D P1 = {t0[i].x, t0[i].y, t0[i].z};
                myPoint3D P2 = {t1[j].x, t1[j].y, t1[j].z};
                // double foundZ = nearestZ(P1, P2).first;
                // double nearest_d = nearestZ(P1, P2).second;
                nearestZResults  = nearestZ(P1, P2);
                double foundZ    = nearestZResults.first;
                double nearest_d = nearestZResults.second;
                nearestZResults_revised  = nearestZ_revised(P1, P2);
                double foundZ_revised    = nearestZResults_revised.first;
                double nearest_d_revisedSQUARED = nearestZResults_revised.second;
                if (std::abs(foundZ - foundZ_revised) > 1e-12 || std::abs(nearest_d*nearest_d - nearest_d_revisedSQUARED) > 1e-07)   std::cout << foundZ - foundZ_revised << ", " << nearest_d*nearest_d - nearest_d_revisedSQUARED << std::endl;
                if (foundZ_revised >= zmin && foundZ_revised <= zmax && nearest_d_revisedSQUARED <= DCA_cutSQUARED) {
                    h -> Fill(foundZ_revised);
                }
            }
        }
    }

    TSpectrum *s = new TSpectrum();
    TH1 *bg = s -> Background(h, 40, "nosmoothing"); 
    bg -> SetLineColor(kRed);
    TCanvas *c = new TCanvas("c", "c", 0,50,2000,1000);
    h->Add(bg, -1.);
    Int_t nfound = s->Search(h, 7., "", 0.001);
    Double_t *xpeaks;
    xpeaks = s->GetPositionX();
    double ctz = xpeaks[0];
    
    if (gPad) gPad->SetGrid(1, 1);
    // h -> DrawCopy("");
    bg -> Draw("SAME");
    c -> Update();
    h -> GetXaxis() -> CenterTitle(true);   h -> GetYaxis() -> CenterTitle(true);   
    h -> Draw("same");
    h -> SetFillColor(kYellow - 7);
    h -> SetLineWidth(1);
    h -> SetFillStyle(1001);

    static TLine *l1 = new TLine();
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(2);
    l1 -> SetLineWidth(4);
    static TLine *l2 = new TLine();
    l2 -> SetLineColor(kBlue);
    l2 -> SetLineStyle(1);
    l2 -> SetLineWidth(4);
    static TLegend *lg = new TLegend(0.12, 0.85, 0.46, 0.9);
    lg -> AddEntry(h, Form("%lu tracklets, found z = %fcm, true z = %fcm", t0.size() + t1.size(), ctz, trueZ), "f");
    lg -> SetTextSize(.02);
    int max_entry = h -> GetBinContent(h -> FindBin(ctz));
    // Update line positions and redraw
    l1->SetX1(ctz); l1->SetX2(ctz);
    l1->SetY1(0);   l1->SetY2(max_entry);
    l1->Draw("same");
    l2->SetX1(trueZ); l2->SetX2(trueZ);
    l2->SetY1(0);     l2->SetY2(max_entry);
    l2->Draw("same");

    // Clear and update the legend for the current event
    lg->Clear();
    lg->AddEntry(h, Form("%lu hits, found z = %fcm, true z = %fcm", t0.size() + t1.size(), ctz, trueZ), "f");
    lg->Draw("same");
    gPad -> SetGrid(1,1); gPad -> Update();
    c -> SaveAs(Form("zFindingPlots/DCA_npeaks_fit_%d.png", evt));
    // delete c;   delete h;   delete bg;
    // h = 0;  bg = 0;


    /*
    TSpectrum *s = new TSpectrum();
    TH1 *bg = s -> Background(h, 40, "nosmoothing");   // estimate the background;
    bg -> SetLineColor(kRed);
    TCanvas *c = new TCanvas("c", "c", 0,50,2000,1000);
    c -> Divide(1, 2);
    c -> cd(1);
    h -> DrawCopy("");
    bg -> Draw("SAME");
    if (gPad) gPad->SetGrid(1, 1);

    c->cd(2);
    h->Add(bg, -1.);
    TH1F *h2 = (TH1F*)h -> Clone("h2");
    Int_t nfound = s->Search(h2, 7., "", 0.001);
    Double_t par[302];
    TF1 *fline = new TF1("fline", "pol1", zmin*2, zmax+25.);
    h -> Fit("fline", "qn"); // 'q': quiet; 'n': no drawing
    // Loop on all found peaks. Eliminate peaks at the background level:
    par[0] = fline->GetParameter(0);
    par[1] = fline->GetParameter(1);
    npeaks = 0;
    Double_t *xpeaks;
    xpeaks = s->GetPositionX();
   
    for (Int_t p = 0; p < nfound; p++) {
        Double_t xp = xpeaks[p];
        std::cout << xp << std::endl;
        Int_t bin = h->GetXaxis()->FindBin(xp);
        Double_t yp = h->GetBinContent(bin);
        if (yp - TMath::Sqrt(yp) < fline->Eval(xp)) continue;
        par[3*npeaks+2] = yp;   // height
        par[3*npeaks+3] = xp;   // mean
        par[3*npeaks+4] = 3;    // sigma
        #if defined(__PEAK_C_FIT_AREAS__)
           par[3*npeaks+2] *= par[3*npeaks+4]*(TMath::Sqrt(TMath::TwoPi()));
        #endif
        npeaks++;
    }
    TF1 *fit = new TF1("fit", fpeaks, zmin*2, zmax+25., 2+3*npeaks);
    // We may have more than the default 25 parameters:
    TVirtualFitter::Fitter(h2, 20+3*npeaks);
    fit -> SetParameters(par);
    fit -> SetNpx(1000);
    h2 -> Fit("fit");
    double ctz = xpeaks[0];
    int max_entry = h2 -> GetBinContent(h2 -> FindBin(ctz));
    h2 -> SetFillColor(kYellow - 7);
    h2 -> SetLineWidth(1);
    h2 -> SetFillStyle(1001);
    h2 -> SetTitle(Form("Found Z of Event %d, #bf{DCA with multiple peaks fit}", evt));
    static TLine *l1 = new TLine();
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(2);
    l1 -> SetLineWidth(4);
    static TLine *l2 = new TLine();
    l2 -> SetLineColor(kBlue);
    l2 -> SetLineStyle(1);
    l2 -> SetLineWidth(4);
    static TLegend *lg = new TLegend(0.12, 0.8, 0.42, 0.9);
    lg -> AddEntry(h2, Form("%lu hits, found z = %fcm, true z = %fcm", t0.size() + t1.size(), ctz, trueZ), "f");
    lg -> SetTextSize(.04);
    // Update line positions and redraw for the current event
    l1->SetX1(ctz); l1->SetX2(ctz);
    l1->SetY1(0);   l1->SetY2(max_entry);
    l1->Draw("same");

    l2->SetX1(trueZ); l2->SetX2(trueZ);
    l2->SetY1(0);     l2->SetY2(max_entry);
    l2->Draw("same");

    // Clear and update the legend for the current event
    lg->Clear();
    lg->AddEntry(h, Form("%lu hits, found z = %fcm, true z = %fcm", t0.size() + t1.size(), ctz, trueZ), "f");
    lg->Draw("same");

    gPad -> SetGrid(1,1); gPad -> Update();
    // gPad -> SetLogy();
    c -> SaveAs(Form("../External/zFindingPlots/DCA_npeaks_%d.png", evt));
    // delete c;
// */
    return ctz;
}

double nearest_z_with_error (const int &evt, std::vector<myTrackletMemberExtended> t0, std::vector<myTrackletMemberExtended> t1,
                             const double &eta_cut_low, const double &eta_cut_high,
                             const double &phi_cut_low, const double &phi_cut_high, const double &trueZ) {
    TH1D *h_weighted = new TH1D("", "", bins, zmin - scanstep/2, zmax + scanstep/2);
    double sum_of_weight = 0;
    double sum_of_weighted_mean = 0;
    for (int i = 0; i < t0.size(); i++) {
        for (int j = 0; j < t1.size(); j++) {
            double dPhi = t0[i].phi - t1[j].phi;
            if (dPhi < phi_cut_high && dPhi > phi_cut_low) {
                myPoint3D P1 = {t0[i].x, t0[i].y, t0[i].z};
                myPoint3D P2 = {t1[j].x, t1[j].y, t1[j].z};
                // P1_l is the possible position of P1 with maximum uncertainty,
                // if it was on the left boundary of the silicon strip;
                myPoint3D P1_l = {t0[i].x, t0[i].y, t0[i].z - 0.8};
                myPoint3D P1_r = {t0[i].x, t0[i].y, t0[i].z + 0.8};
                myPoint3D P2_l = {t1[j].x, t1[j].y, t1[j].z - 0.8};
                myPoint3D P2_r = {t1[j].x, t1[j].y, t1[j].z + 0.8};

                double foundZ = nearestZ(P1, P2).first;
                double foundZ_m1 = nearestZ(P1_l, P2_r).first;
                double foundZ_m2 = nearestZ(P1_r, P2_l).first;
                
                if (foundZ_m1 >= zmin && foundZ_m1 <= zmax && foundZ_m2 >= zmin && foundZ_m2 <= zmax) {
                    double zError = std::abs(foundZ_m1 - foundZ_m2);
                    h_weighted -> Fill(foundZ, 1/zError);
                    sum_of_weight += 1/zError;
                    sum_of_weighted_mean += foundZ/zError;
                }   
            }
        }
    }
    double weightedSum = 0;
    double totalWeight = 0;

    int nBins = h_weighted->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
        double binCenter = h_weighted->GetBinCenter(i);
        double weight = h_weighted->GetBinContent(i);
        if (weight >= 1/40) {
            weightedSum += binCenter * weight;
            totalWeight += weight;
        }
    }

    double weightedAverage = 0;
    if (totalWeight != 0) {
        weightedAverage = weightedSum / totalWeight;
    } else {
        std::cerr << "Total weight is zero, can't compute weighted average" << std::endl;
    }

    double ctz = weightedAverage;
    // double ctz = sum_of_weighted_mean/sum_of_weight;

    TCanvas *can = new TCanvas("c2","c2",0,50,1800,550);
    h_weighted -> Draw();
    h_weighted -> SetFillColor(kYellow - 7);
    h_weighted -> SetLineWidth(1);
    h_weighted -> SetFillStyle(1001);
    h_weighted -> GetXaxis() -> SetTitle("Z coordinate [cm]");
    h_weighted -> GetXaxis() -> SetTitleSize(.05);
    h_weighted -> GetXaxis() -> SetLabelSize(.04);
    h_weighted -> GetXaxis() -> CenterTitle(true);
    h_weighted -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_weighted -> GetXaxis() -> SetTitleOffset(.8);
    h_weighted -> GetYaxis() -> SetTitle("# of Counts");
    h_weighted -> GetYaxis() -> SetTitleSize(.05);
    h_weighted -> GetYaxis() -> SetLabelSize(.04);
    h_weighted -> GetYaxis() -> SetTitleOffset(.62);
    h_weighted -> GetYaxis() -> CenterTitle(true);
    h_weighted -> SetTitle(Form("Found Z of Event %d, #bf{Nearest Z with Error}", evt));
    h_weighted -> SetMinimum(0);
    int max_entry = h_weighted -> GetBinContent(h_weighted->GetMaximumBin());
    TLine *l1 = new TLine(ctz, 0, ctz, max_entry*1.5);
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(2);
    l1 -> SetLineWidth(4);
    l1 -> Draw("same");
    TLine *l2 = new TLine(trueZ, 0, trueZ, max_entry*1.5);
    l2 -> SetLineColor(kBlue);
    l2 -> SetLineStyle(1);
    l2 -> SetLineWidth(4);
    l2 -> Draw("same");
    
    TLatex latex1;
    latex1.SetNDC(); // Use normalized coordinates for positioning
    latex1.SetTextSize(0.04); // Set text size
    latex1.SetTextAlign(13); // Set text alignment
    TString text1 = Form("Weighted average: %.2f", ctz);
    // The coordinates (0.2, 0.8) are in NDC (Normalized Device Coordinates)
    // (0,0) is bottom left, (1,1) is top right
    double textX1 = (ctz + 40.)/40.; // someOffset depends on how far from the line you want the text
    double textY1 = 0.5; // Choose an appropriate y-value
    latex1.DrawLatex(textX1, textY1, text1);
    TLatex latex2;
    latex2.SetNDC(); // Use normalized coordinates for positioning
    latex2.SetTextSize(0.04); // Set text size
    latex2.SetTextAlign(13); // Set text alignment
    TString text2 = Form("True Z: %.2f", trueZ);
    // The coordinates (0.2, 0.8) are in NDC (Normalized Device Coordinates)
    // (0,0) is bottom left, (1,1) is top right
    double textX2 = (trueZ + 40.)/40.; // someOffset depends on how far from the line you want the text
    double textY2 = 0.5; // Choose an appropriate y-value
    latex2.DrawLatex(textX2, textY2, text2);

    TLegend *lg = new TLegend(0.12, 0.8, 0.42, 0.9);
    lg -> AddEntry(h_weighted, Form("%lu hits, found z = %fcm, true z = %fcm", t0.size()+t1.size(), ctz, trueZ), "f");
    // if (globalPhi.size() > ls_cut) {
    //     lg -> AddEntry(h, Form("%2.4fcm, Scan_step=%1.0f#mum, threshold_eta |#Delta#eta|<%1.5f, |#Delta#phi|<%1.5f", ctz, scanstep*1e4, threshold_eta_l, threshold_phi_l), "f");
    // }else {
    //     lg -> AddEntry(h, Form("%2.4fcm, Scan_step=%1.0f#mum, threshold_eta |#Delta#eta|<%1.5f, |#Delta#phi|<%1.5f", ctz, scanstep*1e4, threshold_eta_s, threshold_phi_s), "f");
    // }
    gStyle -> SetLegendTextSize(.03);
    lg -> Draw("same");
    gPad -> SetGrid(1,1); gPad -> Update();
    // gPad -> SetLogy();
    can -> SaveAs(Form("../External/zFindingPlot/nearest_z_with_error_%d.png", evt));
    // delete can;
    // delete h_weighted;   h_weighted = 0; 
    return ctz;

}

double zScan (const int &evt, std::vector<myTrackletMemberExtended> t0, std::vector<myTrackletMemberExtended> t1,
                        const double &eta_cut_low, const double &eta_cut_high,
                        const double &phi_cut_low, const double &phi_cut_high, const double &trueZ) {
    
    TH1D *h = new TH1D("", "", bins, zmin - scanstep/2, zmax + scanstep/2);
    for (int step = 0; step < bins; step++) {
        double z_candidate = zmin + step*scanstep;
        for (int i = 0; i < t0.size(); i++) {
            double z = t0[i].z - z_candidate;
            double theta = std::atan2(t0[i].r, z);
            double eta;
            if (z >= 0) {
                eta = -std::log(std::tan(theta/2));
            } else {
                eta = std::log(std::tan((M_PI - theta)/2));
            }
            t0[i].eta = eta;

            for (int j = 0; j < t1.size(); j++) {
                double z = t1[j].z - z_candidate;
                double theta = std::atan2(t1[j].r, z);
                double eta;
                if (z >= 0) {
                    eta = -std::log(std::tan(theta/2));
                } else {
                    eta = std::log(std::tan((M_PI - theta)/2));
                }
                t1[j].eta = eta;

                double dEta = t0[i].eta - t1[j].eta;
                double dPhi = t0[i].phi - t1[j].phi;
                if (dEta >= eta_cut_low && dEta <= eta_cut_high && dPhi >= phi_cut_low && dPhi <= phi_cut_high) {
                    h -> Fill(z_candidate);
                }
            }
        }
    }
    int binmax = h -> GetMaximumBin(); 
    float xlabel0 = h -> GetXaxis() -> GetBinCenter(binmax);
    int entry0 = h -> GetBinContent(binmax);
    float xlabel_1 = h -> GetXaxis() -> GetBinCenter(binmax - 1); double entry_1 = h -> GetBinContent(binmax - 1);
    float xlabel1 = h -> GetXaxis() -> GetBinCenter(binmax + 1); double entry1 = h -> GetBinContent(binmax + 1);
    float ctz = (xlabel0*entry0 + xlabel_1*entry_1 + xlabel1*entry1) / (entry0 + entry1 + entry_1);

    TCanvas *can = new TCanvas("c2","c2",0,50,1800,550);
    h -> Draw();
    h -> SetFillColor(kYellow - 7);
    h -> SetLineWidth(1);
    h -> SetFillStyle(1001);
    h -> GetXaxis() -> SetTitle("Z coordinate [cm]");
    h -> GetXaxis() -> SetTitleSize(.05);
    h -> GetXaxis() -> SetLabelSize(.04);
    h -> GetXaxis() -> CenterTitle(true);
    h -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h -> GetXaxis() -> SetTitleOffset(.8);
    h -> GetYaxis() -> SetTitle("# of Counts");
    h -> GetYaxis() -> SetTitleSize(.05);
    h -> GetYaxis() -> SetLabelSize(.04);
    h -> GetYaxis() -> SetTitleOffset(.62);
    h -> GetYaxis() -> CenterTitle(true);
    // h -> SetTitle("Displacement at the farthest end is 0.7 mm");
    h -> SetTitle(Form("Z Scan method with z_extension for Event %d, %f < dPhi < %f", evt, phi_cut_low, phi_cut_high));
    h -> SetMinimum(0);
    int max_entry = h -> GetBinContent(h -> FindBin(ctz));
    TLine *l1 = new TLine(ctz, 0, ctz, max_entry);
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(2);
    l1 -> SetLineWidth(4);
    l1 -> Draw("same");
    TLine *l2 = new TLine(trueZ, 0, trueZ, max_entry);
    l2 -> SetLineColor(kBlue);
    l2 -> SetLineStyle(1);
    l2 -> SetLineWidth(4);
    l2 -> Draw("same");
    TLegend *lg = new TLegend(0.12, 0.8, 0.42, 0.9);
    lg -> AddEntry(h, Form("%lu hits, found z = %fcm, true z = %fcm", t0.size() + t1.size(), ctz, trueZ), "f");
    gStyle -> SetLegendTextSize(.03);
    lg -> Draw("same");
    gPad -> SetGrid(1,1); gPad -> Update();
    // gPad -> SetLogy();
    can -> SaveAs(Form("../External/zFindingPlot/zscan_foundz_%d.png", evt));
    // delete can;
    // delete h;   h = 0; 
    return ctz;
}

double zScan_with_error (const int &evt, std::vector<myTrackletMemberExtended> t0, std::vector<myTrackletMemberExtended> t1,
                        const double &eta_cut_low, const double &eta_cut_high,
                        const double &phi_cut_low, const double &phi_cut_high, const double &trueZ) {
    double a = 0;
    double b = 0;
    TH1D *h_weighted = new TH1D("", "", bins, zmin - scanstep/2, zmax + scanstep/2);
    for (int step = 0; step < bins; step++) {
        double z_candidate = zmin + step*scanstep;
        for (int i = 0; i < t0.size(); i++) {
            double z_i = t0[i].z - z_candidate;
            double theta_i = std::atan2(t0[i].r, z_i);
            double eta_i;
            if (z_i >= 0) {
                eta_i = -std::log(std::tan(theta_i/2));
            } else {
                eta_i = std::log(std::tan((M_PI - theta_i)/2));
            }
            t0[i].eta = eta_i;

            for (int j = 0; j < t1.size(); j++) {
                double z_j = t1[j].z - z_candidate;
                double theta_j = std::atan2(t1[j].r, z_j);
                double eta_j;
                if (z_j >= 0) {
                    eta_j = -std::log(std::tan(theta_j/2));
                } else {
                    eta_j = std::log(std::tan((M_PI - theta_j)/2));
                }
                t1[j].eta = eta_j;

                double dEta = t0[i].eta - t1[j].eta;
                double dPhi = t0[i].phi - t1[j].phi;
                if (dEta >= eta_cut_low && dEta <= eta_cut_high && dPhi >= phi_cut_low && dPhi <= phi_cut_high) {
                    double weight = std::abs((t0[i].r - t1[j].r)/(1.6*(t0[i].r + t1[j].r)));
                    h_weighted -> Fill(z_candidate, weight);
                    a += z_candidate*weight;
                    b += weight;
                }
            }
        }
    }
    double weightedSum = 0;
    double totalWeight = 0;
    int nBins = h_weighted->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
        double binCenter = h_weighted->GetBinCenter(i);
        double weight = h_weighted->GetBinContent(i);
            weightedSum += binCenter * weight;
            totalWeight += weight;
        
    }

    double weightedAverage = 0;
    if (totalWeight != 0) {
        weightedAverage = weightedSum / totalWeight;
    } else {
        std::cerr << "Total weight is zero, can't compute weighted average" << std::endl;
    }

    double ctz = a/b;

    TCanvas *can = new TCanvas("c2","c2",0,50,1800,550);
    h_weighted -> Draw();
    h_weighted -> SetFillColor(kYellow - 7);
    h_weighted -> SetLineWidth(1);
    h_weighted -> SetFillStyle(1001);
    h_weighted -> GetXaxis() -> SetTitle("Z coordinate [cm]");
    h_weighted -> GetXaxis() -> SetTitleSize(.05);
    h_weighted -> GetXaxis() -> SetLabelSize(.04);
    h_weighted -> GetXaxis() -> CenterTitle(true);
    h_weighted -> GetXaxis() -> SetNdivisions(31, 5, 0);
    h_weighted -> GetXaxis() -> SetTitleOffset(.8);
    h_weighted -> GetYaxis() -> SetTitle("# of Counts");
    h_weighted -> GetYaxis() -> SetTitleSize(.05);
    h_weighted -> GetYaxis() -> SetLabelSize(.04);
    h_weighted -> GetYaxis() -> SetTitleOffset(.62);
    h_weighted -> GetYaxis() -> CenterTitle(true);
    h_weighted -> SetTitle(Form("Found Z of Event %d, #bf{Z Scan with Errors}", evt));
    h_weighted -> SetMinimum(0);
    int max_entry = h_weighted -> GetBinContent(h_weighted->GetMaximumBin());
    TLine *l1 = new TLine(ctz, 0, ctz, max_entry*1.5);
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(2);
    l1 -> SetLineWidth(4);
    l1 -> Draw("same");
    TLine *l2 = new TLine(trueZ, 0, trueZ, max_entry*1.5);
    l2 -> SetLineColor(kBlue);
    l2 -> SetLineStyle(1);
    l2 -> SetLineWidth(4);
    l2 -> Draw("same");
    
    TLatex latex1;
    latex1.SetNDC(); // Use normalized coordinates for positioning
    latex1.SetTextSize(0.04); // Set text size
    latex1.SetTextAlign(13); // Set text alignment
    TString text1 = Form("Weighted average: %.2f", ctz);
    // The coordinates (0.2, 0.8) are in NDC (Normalized Device Coordinates)
    // (0,0) is bottom left, (1,1) is top right
    double textX1 = (ctz + 40.)/40.; // someOffset depends on how far from the line you want the text
    double textY1 = 0.5; // Choose an appropriate y-value
    latex1.DrawLatex(textX1, textY1, text1);
    TLatex latex2;
    latex2.SetNDC(); // Use normalized coordinates for positioning
    latex2.SetTextSize(0.04); // Set text size
    latex2.SetTextAlign(13); // Set text alignment
    TString text2 = Form("True Z: %.2f", trueZ);
    // The coordinates (0.2, 0.8) are in NDC (Normalized Device Coordinates)
    // (0,0) is bottom left, (1,1) is top right
    double textX2 = (trueZ + 40.)/40.; // someOffset depends on how far from the line you want the text
    double textY2 = 0.5; // Choose an appropriate y-value
    latex2.DrawLatex(textX2, textY2, text2);

    TLegend *lg = new TLegend(0.12, 0.8, 0.42, 0.9);
    lg -> AddEntry(h_weighted, Form("%lu hits, found z = %fcm, true z = %fcm", t0.size() + t1.size(), ctz, trueZ), "f");
    // if (globalPhi.size() > ls_cut) {
    //     lg -> AddEntry(h, Form("%2.4fcm, Scan_step=%1.0f#mum, threshold_eta |#Delta#eta|<%1.5f, |#Delta#phi|<%1.5f", ctz, scanstep*1e4, threshold_eta_l, threshold_phi_l), "f");
    // }else {
    //     lg -> AddEntry(h, Form("%2.4fcm, Scan_step=%1.0f#mum, threshold_eta |#Delta#eta|<%1.5f, |#Delta#phi|<%1.5f", ctz, scanstep*1e4, threshold_eta_s, threshold_phi_s), "f");
    // }
    gStyle -> SetLegendTextSize(.03);
    lg -> Draw("same");
    gPad -> SetGrid(1,1); gPad -> Update();
    // gPad -> SetLogy();
    can -> SaveAs(Form("../External/zFindingPlot/zScan_with_error_%d.png", evt));
    // delete can;
    // delete h_weighted;   h_weighted = 0; 
    return ctz;
}

void dEtaPhiAnalysis(const int &evt, const std::vector<double> &xL, const std::vector<double> &yL, const std::vector<double> &zL, const std::vector<double> &rL, const std::vector<int> &hitL, const std::vector<double> &phi0, const std::vector<double> &phi1, const double &eta_cut, const double &phi_cut) {
    TH2D *h_eta_phi = new TH2D("", "", 81, -4 - .05, 4 + .05, 71, -3.5 - .05, 3.5 + .05);

    for (int step = 0; step < bins; step++) {
        double z_candidate = zmin + step*scanstep;
        // TH2D *h_eta_phi = new TH2D("", "", 80, -4, 4, 70, -3.5, 3.5);
        std::vector<double> eta_0, eta_1;
        for (int i = 0; i < xL.size(); i++) {
            // double x = hitsForOneEvent[i].x_local;
            // double y = hitsForOneEvent[i].y_local;
            double z = zL[i] - z_candidate;
            // double phi = std::atan2(y, x);
            double theta = std::atan2(rL[i], z);
            double eta;
            if (z >= 0) {
                eta = -std::log(std::tan(theta/2));
            } else {
                eta = std::log(std::tan((M_PI - theta)/2));
            }
            if (hitL[i] == 0) {
                eta_0.push_back(eta);
            } else {
                eta_1.push_back(eta);
            }

        }

        for (int i = 0; i < eta_0.size(); i++) {
            for (int j = 0; j < eta_1.size(); j++) {
                double dEta = std::abs(eta_0[i] - eta_1[j]);
                double dPhi = std::abs(phi0[i] - phi1[j]);
                h_eta_phi -> Fill(eta_0[i] - eta_1[j], phi0[i] - phi1[j]);
            }
        }
        // TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,550);
        // h_eta_phi -> Draw("COLZ");
        // can1 -> SaveAs(Form("histoZScan/dEtaPhi_%d.png", step));
        // delete can1;
        // delete h_eta_phi;   h_eta_phi = 0; 
        
        eta_0.clear();  eta_1.clear();
    }
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,550);
    h_eta_phi -> Draw("COLZ");
    TLine *l = new TLine(-4, 0, 4, 0);
    l -> SetLineColor(kRed - 7);
    l -> SetLineWidth(2);
    l -> Draw("same");
    can1 -> SaveAs(Form("histoZScan/dEtaPhiAll_%d.png", evt));
}



void dEtaPhiCheck (const int &evt, const std::vector<double> &xL, const std::vector<double> &yL, const std::vector<double> &zL, 
                     const std::vector<double> &rL, const std::vector<int> &hitL, const std::vector<double> &phi, 
                     const double &phi_cut_lower, const double &phi_cut_upper) {
    TH1D *h_d_phi = new TH1D("", "", 70, -3.5 - .05, 3.5 + .05);

    for (int i = 0; i < xL.size(); i++) {
        for (int j = i + 1; j < xL.size(); j++) {
            h_d_phi -> Fill(phi[i] - phi[j]);
        }
    }

    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,550);
    h_d_phi -> Draw();
    gPad -> SetLogy();
    can1 -> SaveAs(Form("histoZScan/dPhiAll_%d.png", evt));
}

void dEtaPhiCheck_with_TrueZ (const int &evt, const std::vector<double> &xL, const std::vector<double> &yL, const std::vector<double> &zL, 
                               const std::vector<double> &rL, const int &N, const double &TrueX, const double &TrueY, const double &TrueZ) {
    TH2D *h_eta_phi = new TH2D("", "", 80, -4 - .05, 4 + .05, 70, -3.5 - .05, 3.5 + .05);
    TH1D *h_phi = new TH1D("", "", 280, -3.5 - .05, 3.5 + .05);

    std::vector<double> Eta, Phi;

    for (int i = 0; i < xL.size(); i++) {
        double x = xL[i] - TrueX;
        double y = yL[i] - TrueY;
        double phi = std::atan2(y, x);
        double z = zL[i] - TrueZ;
        double theta = std::atan2(rL[i], z);
        double eta;
        if (z >= 0) {
            eta = -std::log(std::tan(theta/2));
        } else {
            eta = std::log(std::tan((M_PI - theta)/2));
        }
        Eta.push_back(eta); Phi.push_back(phi);
    }
    for (int i = 0; i < xL.size(); i++) {
        for (int j = i + 1; j < xL.size(); j++) {
            h_eta_phi -> Fill(Eta[i] - Eta[j], Phi[i] - Phi[j]);
            h_phi -> Fill(Phi[i] - Phi[j]);
        }
    }
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1800,550);
    h_phi -> Draw();
    // h_eta_phi -> Draw("COLZ");
    // h_eta_phi -> GetXaxis() -> SetTitle("Delta Eta");
    // h_eta_phi -> GetYaxis() -> SetTitle("Delta Phi");
    // h_eta_phi -> SetTitle(Form("dEtaPhi checking for Event %d", evt));
    // TLine *l = new TLine(-4, 0, 4, 0);
    // l -> SetLineColor(kRed - 7);
    // l -> SetLineWidth(2);
    // l -> Draw("same");
    // can1 -> SaveAs(Form("histoZScan/dEtaPhiAll_%d.png", evt));
}

double zScan_norm (const int &evt, std::vector<myTrackletMemberExtended> t0, std::vector<myTrackletMemberExtended> t1,
                        const double &eta_cut_low, const double &eta_cut_high,
                        const double &phi_cut_low, const double &phi_cut_high, const double &trueZ) {
    
    TH1D *h_sig = new TH1D("signal", "signal", bins, zmin - scanstep/2, zmax + scanstep/2);
    TH1D *h_bg = new TH1D("background", "background", bins, zmin - scanstep/2, zmax + scanstep/2);
    for (int step = 0; step < bins; step++) {
        double z_candidate = zmin + step*scanstep;
        for (int i = 0; i < t0.size(); i++) {
            double z = t0[i].z - z_candidate;
            double theta = std::atan2(t0[i].r, z);
            double eta;
            if (z >= 0) {
                eta = -std::log(std::tan(theta/2));
            } else {
                eta = std::log(std::tan((M_PI - theta)/2));
            }
            t0[i].eta = eta;

            for (int j = 0; j < t1.size(); j++) {
                double z = t1[j].z - z_candidate;
                double theta = std::atan2(t1[j].r, z);
                double eta;
                if (z >= 0) {
                    eta = -std::log(std::tan(theta/2));
                } else {
                    eta = std::log(std::tan((M_PI - theta)/2));
                }
                t1[j].eta = eta;

                double dEta = t0[i].eta - t1[j].eta;
                double dPhi = t0[i].phi - t1[j].phi;
                if (dEta >= eta_cut_low && dEta <= eta_cut_high && dPhi >= phi_cut_low && dPhi <= phi_cut_high) {
                    h_sig -> Fill(z_candidate);
                } else if (dEta > eta_cut_high && dEta <= eta_cut_high*5 && dPhi > phi_cut_high && dPhi <= 0.2) {
                    h_bg -> Fill(z_candidate);
                }
            }
        }
    }

    double N1 = h_bg -> Integral(h_bg->FindFixBin(-2), h_bg->FindFixBin(-0.1), "") + 
                h_bg -> Integral(h_bg->FindFixBin(0.1), h_bg->FindFixBin(2), "");
    double N2 = h_sig -> Integral(h_sig->FindFixBin(-2), h_sig->FindFixBin(-0.1), "") + 
                h_sig -> Integral(h_sig->FindFixBin(0.1), h_sig->FindFixBin(2), "");
    double N  = N2/N1;
    TH1D* HNormalized = (TH1D*) h_bg->Clone("Normalized HBackground");
    HNormalized -> Add(h_bg, N-1);
    TH1D* hDiff = (TH1D*) h_sig->Clone("Subtracted Signal");
    hDiff -> Add(HNormalized, -1);
    int binmax = hDiff -> GetMaximumBin(); 
    float xlabel0 = hDiff -> GetXaxis() -> GetBinCenter(binmax);
    int entry0 = hDiff -> GetBinContent(binmax);
    float xlabel_1 = hDiff -> GetXaxis() -> GetBinCenter(binmax - 1); double entry_1 = hDiff -> GetBinContent(binmax - 1);
    float xlabel1 = hDiff -> GetXaxis() -> GetBinCenter(binmax + 1); double entry1 = hDiff -> GetBinContent(binmax + 1);
    float ctz = (xlabel0*entry0 + xlabel_1*entry_1 + xlabel1*entry1) / (entry0 + entry1 + entry_1);

    TCanvas *can = new TCanvas("c2","c2",0,50,1800,550);
    hDiff -> Draw();
    hDiff -> SetFillColor(kYellow - 7);
    hDiff -> SetLineWidth(1);
    hDiff -> SetFillStyle(1001);
    hDiff -> GetXaxis() -> SetTitle("Z coordinate [cm]");
    hDiff -> GetXaxis() -> SetTitleSize(.05);
    hDiff -> GetXaxis() -> SetLabelSize(.04);
    hDiff -> GetXaxis() -> CenterTitle(true);
    hDiff -> GetXaxis() -> SetNdivisions(31, 5, 0);
    hDiff -> GetXaxis() -> SetTitleOffset(.8);
    hDiff -> GetYaxis() -> SetTitle("# of Counts");
    hDiff -> GetYaxis() -> SetTitleSize(.05);
    hDiff -> GetYaxis() -> SetLabelSize(.04);
    hDiff -> GetYaxis() -> SetTitleOffset(.62);
    hDiff -> GetYaxis() -> CenterTitle(true);
    // hDiff -> SetTitle("Displacement at the farthest end is 0.7 mm");
    hDiff -> SetTitle(Form("Z Scan Normalized for Event %d", evt));
    hDiff -> SetMinimum(0);
    int max_entry = hDiff -> GetBinContent(hDiff -> FindBin(ctz));
    TLine *l1 = new TLine(ctz, 0, ctz, max_entry);
    l1 -> SetLineColor(kRed);
    l1 -> SetLineStyle(2);
    l1 -> SetLineWidth(4);
    l1 -> Draw("same");
    TLine *l2 = new TLine(trueZ, 0, trueZ, max_entry);
    l2 -> SetLineColor(kBlue);
    l2 -> SetLineStyle(1);
    l2 -> SetLineWidth(4);
    l2 -> Draw("same");
    TLegend *lg = new TLegend(0.12, 0.8, 0.42, 0.9);
    lg -> AddEntry(hDiff, Form("%lu hits, found z = %fcm, true z = %fcm", t0.size() + t1.size(), ctz, trueZ), "f");
    gStyle -> SetLegendTextSize(.03);
    lg -> Draw("same");
    gPad -> SetGrid(1,1); gPad -> Update();
    // gPad -> SetLogy();
    can -> SaveAs(Form("../External/zFindingPlot/zscan_normalized_%d.png", evt));
    // delete can;
    // delete hDiff;   hDiff = 0; 
    delete h_bg;    delete h_sig;   delete HNormalized;
    h_bg = 0;   h_sig = 0;  HNormalized = 0;
    return ctz;
}

