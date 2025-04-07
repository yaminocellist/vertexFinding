#include "globalDefinitions.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"

void angularPlot1D (TH1D* const histo, std::vector<std::string> method, const std::string &fileTitle) {
    histo->SetTitle(fileTitle.c_str());
    int maxBin = histo->GetMaximumBin();
    double maxBinCenter = histo->GetBinCenter(maxBin);
    int maxEntry = histo -> GetBinContent(maxBin);
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,1920,1056);
    histo -> Draw();
    histo -> SetFillColor(kYellow - 7);
    histo -> SetLineWidth(1);
    histo -> SetFillStyle(1001);
    histo -> GetXaxis() -> SetTitleSize(.05);
    histo -> GetXaxis() -> SetLabelSize(.03);
    histo -> GetXaxis() -> CenterTitle(true);
    histo -> GetXaxis() -> SetNdivisions(31, 5, 0);
    histo -> GetXaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> SetTitleSize(.05);
    histo -> GetYaxis() -> SetLabelSize(.03);
    histo -> GetYaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> CenterTitle(true);

    if (histo->GetXaxis()->GetXmin() < -M_PI) {
        double two_pi = TMath::Pi()*2;
        int bin_min = 1;                   // The first bin
        int bin_max = histo->GetNbinsX();  // The last bin

        // Calculate bin positions for each label
        int bin_pi   = bin_max;
        int bin_0    = bin_min + (bin_max - bin_min)/2;
        int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
        int bin_pi_2 = bin_min + (bin_max - bin_min)/4;

        // Set the labels at the calculated positions
        histo->GetXaxis()->SetBinLabel(bin_0, "0");
        histo->GetXaxis()->SetBinLabel(bin_pi_2, "-#pi");
        histo->GetXaxis()->SetBinLabel(binPi_2, "#pi");
        histo->GetXaxis()->SetBinLabel(bin_pi, "2#pi");
        histo->GetXaxis()->SetBinLabel(bin_min, "-2#pi");
        histo->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

        // Ensure the custom labels are displayed by setting the number of divisions
        histo->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
        histo->GetXaxis()->SetLabelSize(0.04);
    }
    else {
        double pi = TMath::Pi();
        int bin_min = 1;                   // The first bin
        int bin_max = histo->GetNbinsX();  // The last bin

        // Calculate bin positions for each label
        int bin_pi   = bin_max;
        int bin_0    = bin_min + (bin_max - bin_min)/2;
        int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
        int bin_pi_2 = bin_min + (bin_max - bin_min)/4;

        // Set the labels at the calculated positions
        histo->GetXaxis()->SetBinLabel(bin_0, "0");
        histo->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
        histo->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
        // histo->GetXaxis()->SetBinLabel(bin_3pi_4, "#frac{3#pi}{4}");
        histo->GetXaxis()->SetBinLabel(bin_pi, "#pi");
        histo->GetXaxis()->SetBinLabel(bin_min, "-#pi");
        histo->GetXaxis()->LabelsOption("h"); // Draw the labels vertically
        // histo->GetXaxis()->SetBinLabel(bin_max - 3*(bin_max - bin_min)/4, "-#frac{3#pi}{4}");
        // histo->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/2, "-#frac{#pi}{2}");
        // histo->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/4, "-#frac{#pi}{4}");

        // Ensure the custom labels are displayed by setting the number of divisions
        histo->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
        histo->GetXaxis()->SetLabelSize(0.04);
    }
    
    TLine *l = new TLine(0, 0, 0, maxEntry);
	l -> Draw("same"); 
    l -> SetLineColor(kRed);
    TLegend *lg = new TLegend(0.12, 0.82, 0.46, 0.9);
    lg -> AddEntry(histo, Form("%.2f <= z_vtx <= %.2f cm, centrality <= 0.70", std::stod(method[4]), std::stod(method[5])), "f");
    gStyle -> SetLegendTextSize(.025);
    lg->Draw("same");
    // gPad -> SetLogy();
    can1 -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void angularPlot2D(TH2D *const h2D, std::vector<std::string> method, const std::string &fileTitle) {
    TCanvas *can = new TCanvas("c2d","c2d",0,50,1920,1056);
    // h2D -> SetTitle(Form("All dPhi values (no event mixed) with %d bins of Z vtx", h2D->GetNbinsY()));
    // h2D -> Draw("lego2");

    if (h2D->GetXaxis()->GetXmin() < -M_PI) {
        double two_pi = TMath::Pi()*2;
        int bin_min = 1;                   // The first bin
        int bin_max = h2D->GetNbinsX();  // The last bin

        // Calculate bin positions for each label
        int bin_pi   = bin_max;
        int bin_0    = bin_min + (bin_max - bin_min)/2;
        int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
        int bin_pi_2 = bin_min + (bin_max - bin_min)/4;

        // Set the labels at the calculated positions
        h2D->GetXaxis()->SetBinLabel(bin_0, "0");
        h2D->GetXaxis()->SetBinLabel(bin_pi_2, "-#pi");
        h2D->GetXaxis()->SetBinLabel(binPi_2, "#pi");
        h2D->GetXaxis()->SetBinLabel(bin_pi, "2#pi");
        h2D->GetXaxis()->SetBinLabel(bin_min, "-2#pi");
        h2D->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

        // Ensure the custom labels are displayed by setting the number of divisions
        h2D->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
        h2D->GetXaxis()->SetLabelSize(0.04);
    }
    else {
        double pi = TMath::Pi();
        int bin_min = 1;                   // The first bin
        int bin_max = h2D->GetNbinsX();  // The last bin

        // Calculate bin positions for each label
        int bin_pi   = bin_max;
        int bin_0    = bin_min + (bin_max - bin_min)/2;
        int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
        int bin_pi_2 = bin_min + (bin_max - bin_min)/4;

        // Set the labels at the calculated positions
        h2D->GetXaxis()->SetBinLabel(bin_0, "0");
        h2D->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
        h2D->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
        // h2D->GetXaxis()->SetBinLabel(bin_3pi_4, "#frac{3#pi}{4}");
        h2D->GetXaxis()->SetBinLabel(bin_pi, "#pi");
        h2D->GetXaxis()->SetBinLabel(bin_min, "-#pi");
        h2D->GetXaxis()->LabelsOption("h"); // Draw the labels vertically
        // h2D->GetXaxis()->SetBinLabel(bin_max - 3*(bin_max - bin_min)/4, "-#frac{3#pi}{4}");
        // h2D->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/2, "-#frac{#pi}{2}");
        // h2D->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/4, "-#frac{#pi}{4}");

        // Ensure the custom labels are displayed by setting the number of divisions
        h2D->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
        h2D->GetXaxis()->SetLabelSize(0.04);
    }

    h2D -> Draw("colz");
    h2D -> GetXaxis() -> CenterTitle(true);
    h2D -> GetYaxis() -> CenterTitle(true);
    h2D -> GetXaxis() -> SetTitleOffset(1.4);   h2D -> GetYaxis() -> SetTitleOffset(1.4);
    can -> Modified();
    can -> Update();
    can -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void angularPlot3D(TH2D * h_dPhi_Z, std::vector<std::string> method, const std::string &fileTitle) {
    TCanvas *can1 = new TCanvas("c3d","c3d",0,50,1920,1056);
    h_dPhi_Z -> SetTitle(Form("All dPhi values (no event mixed) with %d bins of Z vtx", h_dPhi_Z->GetNbinsY()));
    h_dPhi_Z -> Draw("lego2");
    h_dPhi_Z -> GetXaxis() -> CenterTitle(true);
    h_dPhi_Z -> GetYaxis() -> CenterTitle(true);
    h_dPhi_Z -> GetXaxis() -> SetTitleOffset(1.8);   h_dPhi_Z -> GetYaxis() -> SetTitleOffset(1.8);
    can1 -> Modified();
    can1 -> Update();

    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void ZResolutionSinglePlot (TH1D* const histo, std::vector<std::string> method, const std::string &fileTitle) {
    if (isInteger(method[1]))   int lowerRange = stoi(method[1]);
    if (isInteger(method[2]))   int upperRange = stoi(method[2]);
    int maxBin = histo->GetMaximumBin();
    double maxBinCenter = histo->GetBinCenter(maxBin);
    int maxEntry = histo -> GetBinContent(maxBin);
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,1920,1056);
    histo -> Draw();
    histo -> SetFillColor(kYellow - 7);
    histo -> SetLineWidth(1);
    histo -> SetFillStyle(1001);
    histo -> GetXaxis() -> SetTitleSize(.05);
    histo -> GetXaxis() -> SetLabelSize(.03);
    histo -> GetXaxis() -> CenterTitle(true);
    histo -> GetXaxis() -> SetNdivisions(31, 5, 0);
    histo -> GetXaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> SetTitleSize(.05);
    histo -> GetYaxis() -> SetLabelSize(.03);
    histo -> GetYaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> CenterTitle(true);

    TLine *l = new TLine(0, 0, 0, maxEntry);
	l -> Draw("same"); 
    l -> SetLineColor(kRed);
    TLegend *lg = new TLegend(0.12, 0.8, 0.42, 0.9);
    lg -> AddEntry(histo, "Fiducial cut as -250 < MBD_z_vtx < -50 mm", "f");
    gStyle -> SetLegendTextSize(.02);
    lg->Draw("same");
    // gPad -> SetLogy();
    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void TGraphSinglePlot (TGraph* g0, const char *title, const char *Xtitle, const std::string &fileTitle) {
    TCanvas *can1 = new TCanvas("cg","cg",0,50,2100,1200);
    g0 -> SetMarkerStyle(29);
    g0 -> SetMarkerSize(1.1);
    g0 -> SetMarkerColor(kBlue - 7);
    g0 -> SetLineWidth(3);
    g0 -> SetLineColor(kWhite);
    gStyle -> SetTitleW(0.7);  //per cent of the pad width
    gStyle -> SetTitleH(0.08); //per cent of the pad height
    g0 -> SetTitle(title);
    // g0 -> GetXaxis() -> SetTitle("# of hits");
    g0 -> GetXaxis() -> SetTitle(Xtitle);
    g0 -> GetXaxis() -> SetTitleSize(0.05);
    g0 -> GetXaxis() -> SetLabelSize(0.04);
    g0 -> GetXaxis() -> CenterTitle(true);
    g0 -> GetYaxis() -> SetTitle("found z resolution [mm]");
    g0 -> GetYaxis() -> SetTitleSize(0.05);
    g0 -> GetYaxis() -> SetLabelSize(0.025);
    g0 -> GetYaxis() -> CenterTitle(true);
    // g0 -> SetMinimum(-200); // Setting y range;
    // g0 -> SetMaximum(200);  // Setting y range;
    g0 -> GetYaxis() -> SetTitleOffset(0.8); 
    g0 -> GetXaxis() -> SetTitleOffset(0.8); 
    // g0 -> GetXaxis() -> SetLimits(0, 6000); // Setting x range;
    g0 -> Draw("AP SAME");
    gPad->SetGrid(5, 2); gPad->Update();

    TLine *l = new TLine(-250, 0, 8000, 0);
	l -> Draw("same"); 
    l -> SetLineColor(kRed);
    TLegend *lg = new TLegend(0.55, 0.85, 0.9, 0.9);
    lg -> AddEntry(g0, "Fiducial cut as -250 < MBD_z_vtx < -50 mm, centrality <= 0.7", "f");
    gStyle -> SetLegendTextSize(.02);
    lg->Draw("same");
    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void TGraphSinglePlot_Squared (TGraph* g0, const char *title, const char *Xtitle, const char *Ytitle, const std::string &fileTitle) {
    TCanvas *can1 = new TCanvas("cgs","cgs",0,50,1200,1200);
    g0 -> SetMarkerStyle(29);
    g0 -> SetMarkerSize(1.1);
    g0 -> SetMarkerColor(kBlue - 7);
    g0 -> SetLineWidth(3);
    g0 -> SetLineColor(kWhite);
    gStyle -> SetTitleW(0.7);  //per cent of the pad width
    gStyle -> SetTitleH(0.08); //per cent of the pad height
    g0 -> SetTitle(title);
    g0 -> GetXaxis() -> SetTitle(Xtitle);   g0 -> GetYaxis() -> SetTitle(Ytitle);
    g0 -> GetXaxis() -> SetTitleSize(0.05);
    g0 -> GetXaxis() -> SetLabelSize(0.04);
    g0 -> GetXaxis() -> CenterTitle(true);
    g0 -> GetYaxis() -> SetTitleSize(0.05);
    g0 -> GetYaxis() -> SetLabelSize(0.025);
    g0 -> GetYaxis() -> CenterTitle(true);
    // g0 -> SetMinimum(-200); // Setting y range;
    // g0 -> SetMaximum(200);  // Setting y range;
    g0 -> GetYaxis() -> SetTitleOffset(0.8); 
    g0 -> GetXaxis() -> SetTitleOffset(0.8); 
    // g0 -> GetXaxis() -> SetLimits(0, 6000); // Setting x range;
    g0 -> Draw("AP SAME");
    gPad->SetGrid(5, 2); gPad->Update();

    TLine *l = new TLine(-250, -250, -50, -50);
	l -> Draw("same"); 
    l -> SetLineColor(kRed);
    TLegend *lg = new TLegend(0.36, 0.85, 0.9, 0.9);
    lg -> AddEntry(g0, "Fiducial cut as -250 < MBD_z_vtx < -50 mm, centrality <= 0.7", "p");
    gStyle -> SetLegendTextSize(.02);
    lg->Draw("same");
    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void TGraphMultiPlot (TMultiGraph *mg, const char *title, const char *Xtitle, const char *Ytitle, const std::string &fileTitle) {
    TList *graphs = mg->GetListOfGraphs();
    int nGraphs = graphs->GetSize();    // number of graphs;
    TIter next(graphs);
    TGraph *graph;
    for (int i = 0; i < nGraphs; i++) {
        graph = (TGraph *)next();
        graph -> SetMarkerStyle(29);
        graph -> SetMarkerSize(1.1);
        graph -> SetMarkerColor(kBlue - 7);
        if (i == 1)     graph -> SetMarkerColor(2);
        graph -> SetLineWidth(3);
        graph -> SetLineColor(kWhite);
        graph -> GetXaxis() -> SetTitleSize(0.05);
        graph -> GetXaxis() -> SetLabelSize(0.04);
        graph -> GetXaxis() -> CenterTitle(true);
        graph -> GetYaxis() -> SetTitleSize(0.05);
        graph -> GetYaxis() -> SetLabelSize(0.025);
        graph -> GetYaxis() -> CenterTitle(true);
        // graph -> SetMinimum(-200); // Setting y range;
        // graph -> SetMaximum(200);  // Setting y range;
        graph -> GetYaxis() -> SetTitleOffset(0.8); 
        graph -> GetXaxis() -> SetTitleOffset(0.8); 
    }

    // Use "xdpyinfo | grep dimensions" in terminal to get your display's max size;
    TCanvas *can1 = new TCanvas("cgm","cgm",0,50,1056,1056);
    mg->SetTitle(title);
    mg->GetXaxis()->SetTitle(Xtitle);
    mg->GetYaxis()->SetTitle(Ytitle);
    mg -> Draw("AP SAME");
    gStyle -> SetTitleW(0.7);  //per cent of the pad width
    gStyle -> SetTitleH(0.08); //per cent of the pad height
    gPad->SetGrid(5, 2); gPad->Update();

    // TLine *l = new TLine(-250, -250, -50, -50);
	// l -> Draw("same"); 
    // l -> SetLineColor(kRed);
        // TLegend *lg = new TLegend(0.36, 0.85, 0.9, 0.9);
        // lg -> AddEntry(graphs[0], "Fiducial cut as -250 < MBD_z_vtx < -50 mm, centrality <= 0.7", "p");
        // gStyle -> SetLegendTextSize(.02);
        // lg->Draw("same");
    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void EtaPhiSinglePlot (TH1D* const histo, std::vector<std::string> method, Int_t const & target) {
    if (isInteger(method[1]))   int lowerRange = stoi(method[1]);
    if (isInteger(method[2]))   int upperRange = stoi(method[2]);
    int maxBin = histo->GetMaximumBin();
    double maxBinCenter = histo->GetBinCenter(maxBin);
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1920,1056);
    histo -> Draw();
    histo -> SetFillColor(kYellow - 7);
    histo -> SetLineWidth(1);
    histo -> SetFillStyle(1001);
    histo -> GetXaxis() -> SetTitleSize(.05);
    histo -> GetXaxis() -> SetLabelSize(.03);
    histo -> GetXaxis() -> CenterTitle(true);
    histo -> GetXaxis() -> SetNdivisions(31, 5, 0);
    histo -> GetXaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> SetTitleSize(.05);
    histo -> GetYaxis() -> SetLabelSize(.03);
    histo -> GetYaxis() -> SetTitleOffset(.8);
    histo -> GetYaxis() -> CenterTitle(true);

    double pi = TMath::Pi();
    int bin_min = 1;  // The first bin
    int bin_max = histo->GetNbinsX();  // The last bin

    // Calculate bin positions for each label
    int bin_pi = bin_max;
    int bin_0 = bin_min + (bin_max - bin_min)/2;
    int binPi_2 = bin_min + 3*(bin_max - bin_min)/4;
    // int bin_pi_2 = bin_min + (bin_max - bin_min)/2;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;

    // Set the labels at the calculated positions
    histo->GetXaxis()->SetBinLabel(bin_0, "0");
    histo->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    histo->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    // histo->GetXaxis()->SetBinLabel(bin_3pi_4, "#frac{3#pi}{4}");
    histo->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    histo->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // histo->GetXaxis()->SetBinLabel(bin_max - 3*(bin_max - bin_min)/4, "-#frac{3#pi}{4}");
    // histo->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/2, "-#frac{#pi}{2}");
    // histo->GetXaxis()->SetBinLabel(bin_max - (bin_max - bin_min)/4, "-#frac{#pi}{4}");

    // Ensure the custom labels are displayed by setting the number of divisions
    histo->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    histo->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    histo->Draw("HIST");
    histo->GetXaxis()->LabelsOption("h"); // Draw the labels vertically
    // histo -> GetYaxis() -> SetRangeUser(5300e3, 6130e3);
    // histo -> SetTitle(Form("dPhi data of all events whose found z vtx is ~ [-%d, -%d] mm, centered at %0.4f", lowerRange, upperRange, maxBinCenter));
    // gPad -> SetLogy();
    // can1 -> SaveAs(Form("../External/xyFindingPlots/dPhi_all_%d_%d.png", lowerRange, upperRange));
    // histo-> SetTitle(Form("dPhi of %d events mixed up in range of [-209.375, -207.5] mm", target));
    // can1 -> SaveAs(Form("../External/xyFindingPlots/dPhi_mixed_%d.png", target));
}

void doublePlot (TH1D* const hBackground, TH1D* const hSignal, std::vector<std::string> method, Int_t const & target) {
    int upperRange = stoi(method[2]);
    int lowerRange = stoi(method[1]);
    int maxBin = hBackground->GetMaximumBin();
    double maxBinCenter = hBackground->GetBinCenter(maxBin);
    TCanvas *can1 = new TCanvas("c1","c1",0,50,1920,1056);

    // Draw the first histogram
    hBackground -> Draw();
    hBackground -> SetFillColor(kYellow - 7);
    hBackground -> SetLineWidth(1);
    hBackground -> SetFillStyle(1001);
    hBackground -> GetXaxis() -> SetTitle("Phi value");
    hBackground -> GetXaxis() -> SetTitleSize(.05);
    hBackground -> GetXaxis() -> SetLabelSize(.03);
    hBackground -> GetXaxis() -> CenterTitle(true);
    hBackground -> GetXaxis() -> SetNdivisions(31, 5, 0);
    hBackground -> GetXaxis() -> SetTitleOffset(.8);
    hBackground -> GetYaxis() -> SetTitle("# of Counts");
    hBackground -> GetYaxis() -> SetTitleSize(.05);
    hBackground -> GetYaxis() -> SetLabelSize(.03);
    hBackground -> GetYaxis() -> SetTitleOffset(.8);
    hBackground -> GetYaxis() -> CenterTitle(true);
    hBackground -> GetYaxis() -> SetRangeUser(4000e3, 6200e3); // Setting x range;

    // Draw the second histogram on the same canvas
    hSignal -> Draw("SAME");
    hSignal -> SetFillColor(kBlue - 7);
    hSignal -> SetLineWidth(1);
    hSignal -> SetFillStyle(1001);

    // hBackground -> SetTitle(Form("dPhi of %d events mix-up/non-mix-up in range of [-209.375, -207.5] mm", target));
    // hBackground-> SetTitle(Form("dPhi of %d events mixed up in range of [-209.375, -207.5] mm", target));
    // can1 -> SaveAs(Form("../External/xyFindingPlots/dPhi_mixed_%d.png", target));
}

void backgroundCancelling (TH1D* const hBackground, TH1D* const hSignal, std::vector<std::string> method, Int_t const & target) {
    TCanvas *can1 = new TCanvas("csub","csub",0,50,1920,1056);
    can1 -> Divide(1, 2);
    can1 -> cd(1);
    hBackground -> Draw();
    hBackground -> SetFillColor(kYellow - 7);
    hBackground -> SetLineWidth(1);
    hBackground -> SetFillStyle(1001);
    hBackground -> GetXaxis() -> SetTitle("dPhi");
    hBackground -> GetXaxis() -> SetTitleSize(.05);
    hBackground -> GetXaxis() -> SetLabelSize(.03);
    hBackground -> GetXaxis() -> CenterTitle(true);
    hBackground -> GetXaxis() -> SetNdivisions(31, 5, 0);
    hBackground -> GetXaxis() -> SetTitleOffset(.8);
    hBackground -> GetYaxis() -> SetTitle("# of Counts");
    hBackground -> GetYaxis() -> SetTitleSize(.05);
    hBackground -> GetYaxis() -> SetLabelSize(.03);
    hBackground -> GetYaxis() -> SetTitleOffset(.8);
    hBackground -> GetYaxis() -> CenterTitle(true);
    hBackground -> GetYaxis() -> SetRangeUser(4300e3, 6500e3); // Setting x range;

    // Draw the second histogram on the same canvas
    hSignal -> Draw("SAME");
    hSignal -> SetFillColor(kBlue - 7);
    hSignal -> SetLineWidth(1);
    hSignal -> SetFillStyle(1001);

    hBackground -> SetTitle(Form("dPhi of %d events mix-up/non-mix-up in range of [-209.375, -207.5] mm", target));

    can1 -> cd(2);
    double centralPeak = 0.05;
    double N1 = hBackground -> Integral(hBackground->FindFixBin(-M_PI), hBackground->FindFixBin(-centralPeak), "") + 
                hBackground -> Integral(hBackground->FindFixBin(centralPeak), hBackground->FindFixBin(M_PI), "");
    double N2 = hSignal -> Integral(hSignal->FindFixBin(-M_PI), hSignal->FindFixBin(-centralPeak), "") + 
                hSignal -> Integral(hSignal->FindFixBin(centralPeak), hSignal->FindFixBin(M_PI), "");
    double N  = N2/N1;

    TH1D* hNormalized = (TH1D*) hBackground->Clone("Normalized hBackground");
    hNormalized -> Add(hBackground, N-1);
    TH1D* hDiff = (TH1D*) hSignal->Clone("Background Subtracted Signal");
    hDiff -> Add(hNormalized, -1);
    hDiff -> Sumw2();
    hDiff -> Draw("HIST SAME");
    hDiff -> Draw("e1psame");
    hDiff -> GetXaxis() -> SetTitle("dPhi");
    hDiff -> GetXaxis() -> CenterTitle(true);
    hDiff -> GetXaxis() -> SetTitleSize(.05);
    hDiff -> SetFillColor(kYellow - 7);
    hDiff -> SetFillStyle(1001);
    //hDiff -> GetXaxis() -> SetRange(hDiff -> FindFixBin(-0.5),hDiff -> FindFixBin(0.5));
    double peak  = hDiff -> Integral(hDiff->FindFixBin(-centralPeak), hDiff->FindFixBin(centralPeak), "");
    // double Ratio = peak/events;
    
    // hDiff -> SetTitle(Form("Subtracted Signal for %5.0f Events", events));
    hDiff -> SetTitle("Background Subtracted dPhi");
    gStyle->SetEndErrorSize(6);
    gStyle->SetErrorX(0.5);
    
    gStyle -> SetTitleFont(100,"t");
    gStyle -> SetTitleSize(0.065,"t");

    //gPad -> SetLogy(1);
    TLine *l = new TLine(-M_PI,0,M_PI,0);
	l -> Draw("same"); 
    l -> SetLineColor(kGreen);

    can1 -> SaveAs("../../External/zFindingPlots/dPhi_subtraction.png");
}

void backgroundCancelling_dPhi (TH1D* const hBackground, TH1D* const hSignal, std::vector<std::string> method, Int_t const & target) {
    double pi = TMath::Pi();
    double signalRatio;
    int bin_min  = 1;                     // The first bin;
    int bin_max  = hSignal->GetNbinsX();  // The last bin;
    // Calculate bin positions for each label
    int bin_pi   = bin_max;
    int bin_0    = bin_min + (bin_max - bin_min)/2;
    int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    // Set the labels at the calculated positions
    hBackground->GetXaxis()->SetBinLabel(bin_0, "0");
    hBackground->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    hBackground->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    hBackground->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    hBackground->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions:
    hBackground->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    // Draw the labels vertically:
    hBackground->GetXaxis()->LabelsOption("h"); 

    hSignal->GetXaxis()->SetBinLabel(bin_0, "0");
    hSignal->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    hSignal->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    hSignal->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    hSignal->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    hSignal->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    hSignal->GetXaxis()->LabelsOption("h");
    
    std::vector<std::string> options = splitString(method[0], '_');
    TCanvas *can1 = new TCanvas("csub","csub",0,50,2560,1440);
    can1 -> Divide(1, 2);
    can1 -> cd(1);
    double phi_range_low = -2.4, phi_range_high = -1.8;
    double phi_range_low2 = 1.01, phi_range_high2 = 1.55;
    int bin_range_low = hBackground->FindBin(phi_range_low), bin_range_high = hBackground->FindBin(phi_range_high);
    double max_unmixed = -1, max_mixed = -1, min_unmixed = std::numeric_limits<double>::max(), min_mixed = std::numeric_limits<double>::max(), current_binContent;
    for (int bin = bin_range_low; bin <= bin_range_high; bin++) {
        current_binContent = hSignal->GetBinContent(bin);
        if (max_unmixed < current_binContent)   max_unmixed = current_binContent;
        current_binContent = hBackground->GetBinContent(bin);
        if (max_mixed < current_binContent)     max_mixed = current_binContent;
    }
    std::cout << (hBackground->Integral(hBackground->FindFixBin(phi_range_low),hBackground->FindFixBin(phi_range_high),""))/(hSignal->Integral(hSignal->FindFixBin(phi_range_low),hSignal->FindFixBin(phi_range_high),"")) << std::endl;
    std::cout << (hBackground->Integral(hBackground->FindFixBin(phi_range_low),hBackground->FindFixBin(phi_range_high),"")+hBackground->Integral(hBackground->FindFixBin(phi_range_low2),hBackground->FindFixBin(phi_range_high2),""))/(hSignal->Integral(hSignal->FindFixBin(phi_range_low),hSignal->FindFixBin(phi_range_high),"")+hSignal->Integral(hSignal->FindFixBin(phi_range_low2),hSignal->FindFixBin(phi_range_high2),"")) << std::endl;
    std::cout << (hBackground->Integral(hBackground->FindFixBin(-M_PI),hBackground->FindFixBin(phi_range_high),"")+hBackground->Integral(hBackground->FindFixBin(phi_range_low2),hBackground->FindFixBin(M_PI),""))/(hSignal->Integral(hSignal->FindFixBin(-M_PI),hSignal->FindFixBin(phi_range_high),"")+hSignal->Integral(hSignal->FindFixBin(phi_range_low2),hSignal->FindFixBin(M_PI),"")) << std::endl;
    signalRatio = (hBackground->Integral(hBackground->FindFixBin(-M_PI),hBackground->FindFixBin(phi_range_high),"")+hBackground->Integral(hBackground->FindFixBin(phi_range_low2),hBackground->FindFixBin(M_PI),""))/(hSignal->Integral(hSignal->FindFixBin(-M_PI),hSignal->FindFixBin(phi_range_high),"")+hSignal->Integral(hSignal->FindFixBin(phi_range_low2),hSignal->FindFixBin(M_PI),""));
    hSignal -> Add(hSignal, signalRatio - 1);

    // Loop again for y-axis plotting range:
    phi_range_low = -M_PI;
    phi_range_high = M_PI;
    bin_range_low = hBackground->FindBin(phi_range_low);
    bin_range_high = hBackground->FindBin(phi_range_high);
    for (int bin = bin_range_low; bin <= bin_range_high; bin++) {
        current_binContent = hSignal->GetBinContent(bin);
        if (max_unmixed < current_binContent)   max_unmixed = current_binContent;
        if (min_unmixed > current_binContent && current_binContent > 0)   min_unmixed = current_binContent;
        current_binContent = hBackground->GetBinContent(bin);
        if (max_mixed < current_binContent)     max_mixed = current_binContent;
        if (min_mixed > current_binContent && current_binContent > 0)     min_mixed = current_binContent;
    }

    double max_y_range = max_unmixed, min_y_range = min_mixed;
    if (max_unmixed < max_mixed)    max_y_range = max_mixed;
    if (min_unmixed < min_mixed)    min_y_range = min_unmixed;

    hSignal -> GetYaxis() -> SetRangeUser(min_y_range*0.9, max_y_range*1.1);
    hBackground -> GetYaxis() -> SetRangeUser(min_y_range*0.9, max_y_range*1.1);

    hSignal -> Draw("SAME");
    hSignal -> GetXaxis()->SetLabelSize(0.06);
    hSignal -> GetXaxis()->SetTitleSize(0.05);
    hSignal -> GetXaxis() -> CenterTitle(true);
    hSignal -> GetXaxis() -> SetTitleOffset(.9);
    hSignal -> GetYaxis() -> CenterTitle(true);
    hSignal -> GetYaxis() -> SetTitleSize(0.05);
    hBackground -> Draw("SAME");
    hBackground -> SetLineColor(2);
    hBackground -> GetXaxis() -> CenterTitle(true);
    hBackground -> GetYaxis() -> CenterTitle(true);
    if (options[1] == "dE") {
        hSignal -> SetTitle(Form("%d events with |dEta| < %.2f, %2.2fcm < z vtx < %2.2fcm, %1.2f < centrality < %1.2f", target, dEta_cut, std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])));
        hBackground -> SetTitle(Form("%d events with |dEta| < %.2f, %2.2fcm < z vtx < %2.2fcm, %1.2f < centrality < %1.2f", target, dEta_cut, std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])));
    }
    else if (options[1] == "E") {
        hSignal -> SetTitle(Form("%d events with |Eta| < %.2f,  %2.2fcm < z vtx < %2.2fcm, %1.2f < centrality < %1.2f", target, Eta_range, std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])));
        hBackground -> SetTitle(Form("%d events with |Eta| < %.2f,  %2.2fcm < z vtx < %2.2fcm, %1.2f < centrality < %1.2f", target, Eta_range, std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])));
    }

    TLegend *lg = new TLegend(0.12, 0.8, 0.35, 0.9);
    lg -> AddEntry(hSignal, "Unmixed Events' dPhi", "l");
    lg -> AddEntry(hBackground, Form("Mixed Events' dPhi, multiplied by %.2f", signalRatio), "l");
    gStyle -> SetLegendTextSize(.043);
    lg->Draw("same");
    if (options[3] == "s") {
        // Assuming you have a canvas already created and the plots are drawn
        TLatex text2;
        text2.SetTextSize(0.1); // Adjust the size as needed
        text2.SetTextColor(kRed); // Set color to red
        text2.SetTextAlign(22);   // Center the text
        text2.DrawLatex(halfPI, max_y_range, "Simulated Data"); // Adjust the coordinates (x, y) to place the text appropriately
    }

    // TLine *l2 = new TLine(phi_range_low2, 0, phi_range_low2, max_y_range);
	// l2 -> Draw("same"); 
    // l2 -> SetLineColor(kRed);
    // TLine *l3 = new TLine(phi_range_high2, 0, phi_range_high2, max_y_range);
	// l3 -> Draw("same"); 
    // l3 -> SetLineColor(kRed);
    can1 -> cd(2);
    double centralPeak = 0.05;
    double N1 = hBackground -> Integral(hBackground->FindFixBin(-M_PI), hBackground->FindFixBin(-centralPeak), "") + 
                hBackground -> Integral(hBackground->FindFixBin(centralPeak), hBackground->FindFixBin(M_PI), "");
    double N2 = hSignal -> Integral(hSignal->FindFixBin(-M_PI), hSignal->FindFixBin(-centralPeak), "") + 
                hSignal -> Integral(hSignal->FindFixBin(centralPeak), hSignal->FindFixBin(M_PI), "");
    double N  = N2/N1;

    TH1D* hNormalized = (TH1D*) hBackground->Clone("Normalized hBackground");
    hNormalized -> Add(hBackground, N-1);
    TH1D* hDiff = (TH1D*) hSignal->Clone("Background Subtracted Signal");
    hDiff -> Add(hNormalized, -1);
    hDiff -> Add(hDiff, 1./signalRatio - 1);
    hDiff -> Sumw2();
    hDiff -> Draw("HIST SAME");
    hDiff -> Draw("e1psame");
    hDiff -> GetXaxis() -> SetTitle("dPhi value");
    hDiff -> GetXaxis() -> CenterTitle(true);
    hDiff -> GetXaxis() -> SetTitleSize(.05);
    hDiff -> GetXaxis() -> SetLabelSize(0.06);
    hDiff -> GetXaxis() -> SetTitleOffset(.9);
    hDiff -> GetYaxis() -> SetTitleSize(.05);
    hDiff -> SetFillColor(kYellow - 7);
    hDiff -> SetFillStyle(1001);
    //hDiff -> GetXaxis() -> SetRange(hDiff -> FindFixBin(-0.5),hDiff -> FindFixBin(0.5));
    double peak  = hDiff -> Integral(hDiff->FindFixBin(-centralPeak), hDiff->FindFixBin(centralPeak), "");
    // double Ratio = peak/events;
    
    // if (options[1] == "wdE")
    //     hDiff -> SetTitle(Form("Subtracted Signal for %d Events, with |dEta| < %.2f", target, dEta_cut));
    // else if (options[1] == "wE")
    //     hDiff -> SetTitle(Form("Subtracted Signal for %d Events, with |Eta| < %.2f", target, Eta_range));
    // else
    //     hDiff -> SetTitle(Form("Subtracted Signal for %d Events", target));
    gStyle->SetEndErrorSize(6);
    gStyle->SetErrorX(0.5);
    
    gStyle -> SetTitleFont(100,"t");
    gStyle -> SetTitleSize(0.065,"t");

    //gPad -> SetLogy(1);
    TLine *l = new TLine(-M_PI,0,M_PI,0);
	l -> Draw("same"); 
    l -> SetLineColor(kGreen);

    int eventCount = std::stoi(method[1]);
    double method4 = std::stod(method[4]);
    double method5 = std::stod(method[5]);
    double method2 = std::stod(method[2]);
    double method3 = std::stod(method[3]);

    std::string filePrefix;
    if (options[1] == "wo") {
        hDiff -> SetTitle(Form("Subtracted Signal for %d Events, %.2f < centrality < %.2f", target, method2, method3));
        filePrefix = Form("%devents_%2.2f_%2.2f_%1.2f_%1.2f", eventCount, method4, method5, method2, method3);
    } else if (options[1] == "dE") {
        hDiff -> SetTitle(Form("Subtracted Signal for %d Events, with |dEta| < %.2f, %.2f < centrality < %.2f", target, dEta_cut, method2, method3));
        filePrefix = Form("with_dEta_cut_%.2f_%devents_%2.2f_%2.2f_%1.2f_%1.2f", dEta_cut, eventCount, method4, method5, method2, method3);
    } else {
        hDiff -> SetTitle(Form("Subtracted Signal for %d Events, with |Eta| < %.2f, %.2f < centrality < %.2f", target, Eta_range, method2, method3));
        filePrefix = Form("with_Eta_range_%.2f_%devents_%2.2f_%2.2f_%1.2f_%1.2f", Eta_range, eventCount, method4, method5, method2, method3);
    }
    if (options[3] == "s") {
        // Loop again for y-axis plotting range:
        bin_range_low = hDiff->FindBin(-M_PI);
        bin_range_high = hDiff->FindBin(M_PI);
        double max_content = -1;
        for (int bin = bin_range_low; bin <= bin_range_high; bin++) {
            current_binContent = hDiff->GetBinContent(bin);
            if (max_content < current_binContent)   max_content = current_binContent;
        }
        // Assuming you have a canvas already created and the plots are drawn
        TLatex text;
        text.SetTextSize(0.1); // Adjust the size as needed
        text.SetTextColor(kRed); // Set color to red
        text.SetTextAlign(22);   // Center the text
        text.DrawLatex(halfPI, max_content/2, "Simulated Data"); // Adjust the coordinates (x, y) to place the text appropriately
        can1->SaveAs(("../../External/zFindingPlots/Simulation_data_dPhi_mixedsubtract_" + filePrefix + ".png").c_str());
    }
    else {
        can1->SaveAs(("../../External/zFindingPlots/dPhi_mixedsubtract_" + filePrefix + ".png").c_str());
    }

    if (options[2] == "f") {
        std::string rootFileName;
        TParameter<double>* ratioParam = new TParameter<double>("signal_multiplied_ratio", signalRatio);
        TParameter<double>* numParam   = new TParameter<double>("number_of_events", target);
        TParameter<double>* lCenParam  = new TParameter<double>("lower_centrality", method2);
        TParameter<double>* hCenParam  = new TParameter<double>("higher_centrality", method3);
        TParameter<double>* EtaParam   = new TParameter<double>("eta_range", Eta_range);
        if (options[3] == "s")
            rootFileName = "../../External/forFit/Simulation_data_hDiff_" + filePrefix + ".root";
        else
            rootFileName = "../../External/forFit/hDiff_" + filePrefix + ".root";
        TFile *outputFile = new TFile(rootFileName.c_str(), "RECREATE");
        hDiff     ->Write();
        ratioParam->Write();
        numParam  ->Write();
        lCenParam ->Write();
        hCenParam ->Write();
        EtaParam  ->Write();
        outputFile->Close();
    }

    // if (options[1] == "wo") {
    //     can1 -> SaveAs(Form("../../External/zFindingPlots/dPhi_mixedsubtract_%devents_%2.2f_%2.2f_%1.2f_%1.2f.png", std::stoi(method[1]), std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])));

    //     if (options[2] == "f") {
    //         TFile *outputFile = new TFile(Form("../../External/zFindingPlots/hDiff_%devents_%2.2f_%2.2f_%1.2f_%1.2f.root", std::stoi(method[1]), std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])), "RECREATE");
    //         hDiff->Write();
    //         outputFile->Close();
    //     }
    // } else if (options[1] == "wdE") {
    //     can1 -> SaveAs(Form("../../External/zFindingPlots/dPhi_mixedsubtract_with_dEta_cut_%.2f_%devents_%2.2f_%2.2f_%1.2f_%1.2f.png", dEta_cut, std::stoi(method[1]), std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])));

    //     if (options[2] == "f") {
    //         TFile *outputFile = new TFile(Form("../../External/zFindingPlots/hDiff_with_dEta_cut_%.2f_%devents_%2.2f_%2.2f_%1.2f_%1.2f.root", dEta_cut, std::stoi(method[1]), std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])), "RECREATE");
    //         hDiff->Write();
    //         outputFile->Close();
    //     }
    // } else {
    //     can1 -> SaveAs(Form("../../External/zFindingPlots/dPhi_mixedsubtract_with_Eta_range_%.2f_%devents_%2.2f_%2.2f_%1.2f_%1.2f.png", Eta_range, std::stoi(method[1]), std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])));

    //     if (options[2] == "f") {
    //         TFile *outputFile = new TFile(Form("../../External/zFindingPlots/hDiff_with_Eta_range_%.2f_%devents_%2.2f_%2.2f_%1.2f_%1.2f.root", Eta_range, std::stoi(method[1]), std::stod(method[4]), std::stod(method[5]), std::stod(method[2]), std::stod(method[3])), "RECREATE");
    //         hDiff->Write();
    //         outputFile->Close();
    //     }
    // }        
}

void ArrayPlot1D_Logy (const std::vector<TH1D*>& h, std::vector<std::string> method, const std::string &fileTitle) {
    double entries = 0;
    for (int i = 0; i < h.size(); i++) {
        entries += h[i] -> GetEntries();
    }
    entries /= 14;
    // h[0] -> Scale(1.0 / entries);
    // h[0]->SetFillStyle(3025);

    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2100,1200);
    h[0] -> SetLineWidth(3);
    h[0] -> SetLineColor(2);
    h[0] -> Draw();
    for (int d = 1; d < h.size(); d++) {
        // entries = h[d] -> GetEntries();
        // h[d] -> Scale(1.0 / entries);
        // h[d]->SetLineColor(d + 1); // ROOT colors start from 1 (0 is white)
        // h[d]->SetFillColor(d + 1);
        // h[d]->SetFillStyle(3001 + d); 
        h[d] -> SetLineWidth(3);
        h[d] -> Draw("SAME");
    }
    h[0] -> GetYaxis() -> SetRangeUser(5e4, h[1]->GetEntries()/2e2); // Setting x range;
    h[0] -> GetXaxis() -> CenterTitle(true);    h[0] -> GetYaxis() -> CenterTitle(true);

    double pi = TMath::Pi();
    int bin_min  = 1;  // The first bin
    int bin_max  = h[0]->GetNbinsX();  // The last bin
    // Calculate bin positions for each label
    int bin_pi   = bin_max;
    int bin_0    = bin_min + (bin_max - bin_min)/2;
    int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    // Set the labels at the calculated positions
    h[0]->GetXaxis()->SetBinLabel(bin_0, "0");
    h[0]->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h[0]->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h[0]->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h[0]->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h[0]->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h[0]->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    // h[0]->Draw("HIST");
    h[0]->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    gPad -> SetLogy(1);
    can1 -> SaveAs(Form("../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void ArrayPlot1D_Logy_ver2 (const std::vector<TH1D*>& h, std::vector<std::string> method, const std::string &fileTitle) {
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2100,1200);
    // h[0] -> Draw();
    // h[1] -> Draw("SAME");
    // h[2] -> Draw("SAME");
    // h[3] -> Draw("SAME");
    // h[4] -> Draw("SAME");
    // h[5] -> Draw("SAME");
    // h[6] -> Draw("SAME");
    // h[7] -> Draw("SAME");
    // h[8] -> Draw("SAME");
    // h[9] -> Draw("SAME");
    // h[10] -> Draw("SAME");
    // h[11] -> Draw("SAME");
    // h[12] -> Draw("SAME");
    // h[13] -> Draw("SAME");
    // h[14] -> Draw("SAME");
    // h[15] -> Draw("SAME");
    // h[16] -> Draw("SAME");
    // h[17] -> Draw("SAME");
    // h[18] -> Draw("SAME");
    // h[19] -> Draw("SAME");

    double pi = TMath::Pi();
    int bin_min  = 1;  // The first bin
    int bin_max  = h[16]->GetNbinsX();  // The last bin
    // Calculate bin positions for each label
    int bin_pi   = bin_max;
    int bin_0    = bin_min + (bin_max - bin_min)/2;
    int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    // Set the labels at the calculated positions
    h[16]->GetXaxis()->SetBinLabel(bin_0, "0");
    h[16]->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h[16]->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h[16]->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h[16]->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h[16]->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h[16]->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    h[16]->Draw("HIST");
    h[16]->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    for (int i = 0; i < h.size(); i++) {
        // std::cout << h[i] -> GetEntries() << std::endl;
        h[i] -> Draw("SAME");
        h[i] -> SetLineWidth(3);
        h[i] -> GetYaxis() -> SetRangeUser(1e5, h[16]->GetEntries()/5e2);
        h[i] -> GetXaxis() -> CenterTitle(true);    h[i] -> GetYaxis() -> CenterTitle(true);
    }
    gPad -> SetLogy(1);
}

void ArrayPlot1D_Rescale (const std::vector<TH1D*>& h, std::vector<std::string> method, const std::string &fileTitle) {
    std::vector<double> max_entries;
    double min_y = std::numeric_limits<double>::max();
    for (int i = 0; i < h.size(); i++) {
        max_entries.push_back(h[i]->GetBinContent(h[i]->GetMaximumBin()));
        h[i] -> Add(h[i], 1/(max_entries[i]+1) - 1);

        // Find minimum bin content across all histograms
        for (int bin = 1; bin <= h[i]->GetNbinsX(); ++bin) {
            double bin_content = h[i]->GetBinContent(bin);
            if (bin_content > 0 && bin_content < min_y) {
                min_y = bin_content;
            }
        }
        // h[i] -> GetYaxis() -> SetRangeUser(0.6, 1.1);
    }
    TLegend *lg = new TLegend(0.12, 0.7, 0.27, 0.90);
    gStyle -> SetLegendTextSize(.018);
    for (int i = 0; i < h.size(); i++) {
        h[i] -> GetYaxis() -> SetRangeUser(min_y*0.9, 1.1);
        h[i] -> SetLineColor(32 + i);
        lg -> AddEntry(h[i], Form("Centrality ~ [%.2f, %.2f]", static_cast<double>(i)*0.05, static_cast<double>(i + 1)*0.05), "l");
        h[i] -> SetTitle(fileTitle.c_str());
    }
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2560,1600);
    h[0] -> SetLineWidth(2);
    h[0] -> SetLineColor(2);
    h[0] -> Draw();
    for (int d = 1; d < h.size(); d++) {
        // entries = h[d] -> GetEntries();
        // h[d] -> Scale(1.0 / entries);
        // h[d]->SetLineColor(d + 1); // ROOT colors start from 1 (0 is white)
        // h[d]->SetFillColor(d + 1);
        // h[d]->SetFillStyle(3001 + d); 
        h[d] -> SetLineWidth(2);
        h[d] -> Draw("SAME");
    }
    // h[0] -> GetYaxis() -> SetRangeUser(5e4, h[1]->GetEntries()/2e2); // Setting x range;
    h[0] -> GetXaxis() -> CenterTitle(true);    h[0] -> GetYaxis() -> CenterTitle(true);

    double pi = TMath::Pi();
    int bin_min  = 1;  // The first bin
    int bin_max  = h[0]->GetNbinsX();  // The last bin
    // Calculate bin positions for each label
    int bin_pi   = bin_max;
    int bin_0    = bin_min + (bin_max - bin_min)/2;
    int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    // Set the labels at the calculated positions
    h[0]->GetXaxis()->SetBinLabel(bin_0, "0");
    h[0]->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h[0]->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h[0]->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h[0]->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h[0]->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h[0]->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    // h[0]->Draw("HIST");
    h[0]->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    lg->Draw("same");
    can1 -> Update();

    can1 -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void ArrayPlot1D_Rescale_ver2 (const std::vector<TH1D*>& h, std::vector<std::string> method, const std::string &fileTitle) {
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2100,1200);
    std::vector<double> max_entries;
    double min_y = std::numeric_limits<double>::max();
    for (int i = 0; i < h.size(); i++) {
        max_entries.push_back(h[i]->GetBinContent(h[i]->GetMaximumBin()));
        h[i] -> Add(h[i], 1/(max_entries[i]+1) - 1);
        for (int bin = 1; bin <= h[i]->GetNbinsX(); ++bin) {
            double bin_content = h[i]->GetBinContent(bin);
            if (bin_content > 0 && bin_content < min_y) {
                min_y = bin_content;
            }
        }
    }
    h[16] -> Draw();
    TLegend *lg = new TLegend(0.12, 0.73, 0.28, 0.9);
    TLegend *lg2 = new TLegend(0.28 , 0.73, 0.44, 0.9);
    gStyle -> SetLegendTextSize(.015);
    for (int i = 0; i < h.size(); i++) {
        if (i < 10) {
            lg -> AddEntry(h[i], Form("Z VTX between %.2f and %.2f cm", -30 + static_cast<double>(i)*3, -27 + static_cast<double>(i)*3), "l");
        } else {
            lg2 -> AddEntry(h[i], Form("Z VTX between %.2f and %.2f cm", -30 + static_cast<double>(i)*3, -27 + static_cast<double>(i)*3), "l");
        }
        h[i] -> GetYaxis() -> SetRangeUser(min_y*0.95, 1.05);
        h[i] -> GetXaxis() -> CenterTitle(true);    h[i] -> GetYaxis() -> CenterTitle(true);
        h[i] -> SetLineWidth(2);
        h[i] -> SetLineColor(29 + i);
        h[i] -> Draw("SAME");
        h[i] -> SetTitle(fileTitle.c_str());
    }

    h[16] -> SetLineColor(2);
    double pi = TMath::Pi();
    int bin_min  = 1;  // The first bin
    int bin_max  = h[0]->GetNbinsX();  // The last bin
    // Calculate bin positions for each label
    int bin_pi   = bin_max;
    int bin_0    = bin_min + (bin_max - bin_min)/2;
    int binPi_2  = bin_min + 3*(bin_max - bin_min)/4;
    int bin_pi_2 = bin_min + (bin_max - bin_min)/4;
    // Set the labels at the calculated positions
    h[16]->GetXaxis()->SetBinLabel(bin_0, "0");
    h[16]->GetXaxis()->SetBinLabel(bin_pi_2, "#frac{-#pi}{2}");
    h[16]->GetXaxis()->SetBinLabel(binPi_2, "#frac{#pi}{2}");
    h[16]->GetXaxis()->SetBinLabel(bin_pi, "#pi");
    h[16]->GetXaxis()->SetBinLabel(bin_min, "-#pi");
    // Ensure the custom labels are displayed by setting the number of divisions
    h[16]->GetXaxis()->SetNdivisions(9, 0, 0, kFALSE);
    h[16]->GetXaxis()->SetLabelSize(0.04);
    // Update histogram to refresh the axis
    // h[0]->Draw("HIST");
    h[16]->GetXaxis()->LabelsOption("h"); // Draw the labels vertically

    lg->Draw("same");   lg2->Draw("same");
    can1 -> Update();

    can1 -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void ArrayPlot1D_Rescale_dEta (const std::vector<TH1D*>& h, std::vector<std::string> method, const std::string &fileTitle) {
    double minRange = h[0]->GetXaxis()->GetXmin();
    double maxRange = h[0]->GetXaxis()->GetXmax();
    int index = 0;
    double max = 0;
    TCanvas *can1 = new TCanvas("c1d","c1d",0,50,2100,1200);
    for (int i = 0; i < h.size(); i++) {
        double N = h[i] -> Integral(h[i]->FindFixBin(minRange), h[i]->FindFixBin(maxRange), "");
        h[i] -> Add(h[i], 1/(N/1000+1) - 1);
        h[i] -> GetXaxis() -> SetRangeUser(minRange*5/8, maxRange*5/8);
        h[i] -> GetXaxis() -> CenterTitle(true);    h[i] -> GetYaxis() -> CenterTitle(true);
        int maxEntry = h[i] -> GetBinContent(h[i]->GetMaximumBin());
        if (maxEntry > max) {
            max = maxEntry;
            index = i;
        } 
        // h[i] -> SetLineColor(30+i);
        h[i] -> SetLineWidth(3);
        // h[i] -> Draw("SAME");
    }
    std::cout << index << ", " << max << std::endl;
    h[index] -> Draw("SAME");
    h[index] -> GetYaxis() -> SetRangeUser(0, max*1.05);
    for (int i = 0; i < h.size(); i++) {
        if (i != index)     h[i] -> Draw("SAME");
    }
    h[0] -> SetLineColor(2);
    h[1] -> SetLineColor(4);
    h[2] -> SetLineColor(6);
    h[3] -> SetLineColor(8);
    h[4] -> SetLineColor(9);
    h[5] -> SetLineColor(12);
    h[6] -> SetLineColor(28);
    h[7] -> SetLineColor(30);
    h[8] -> SetLineColor(31);
    h[9] -> SetLineColor(32);
    h[10] -> SetLineColor(38);
    h[11] -> SetLineColor(40);
    h[12] -> SetLineColor(42);
    h[13] -> SetLineColor(46);

    can1 -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileTitle.c_str()));
}

void ArrayPlot_dEta_1D_Logy (
    const std::vector<TH1D*> &h,
    const std::vector<std::string> &method,
    const std::string &fileName
) {
    TCanvas *c1 = new TCanvas("c1dEta", "dPhi Histogram", 1920, 1056);
    // h[0]->Draw("same");
    // h[1]->Draw("same");   h[1]->SetLineColor(kGreen+2);
    // h[2]->Draw("same");
    // h[3]->Draw("same");
    // h[4]->Draw("same");   h[4]->SetLineColor(2);
    // h[5]->Draw("same");
    // h[6]->Draw("same");
    // h[7]->Draw("same");
    // h[8]->Draw("same");
    // h[9]->Draw("same");
    // h[10]->Draw("same");  h[10]->SetLineColor(2);
    // h[11]->Draw("same");
    // h[12]->Draw("same");  h[12]->SetLineColor(kMagenta);
    // h[13]->Draw("same");

    int maxEntry = 0;
    for (int i = 0; i < h.size(); i++) {
            h[i]->Draw("same");
            h[i]->GetXaxis()->CenterTitle(true);    h[i]->GetYaxis()->CenterTitle(true);
            h[i]->SetLineWidth(3);
            h[i]->SetLineColor(20+i);
            int temp = h[i]->GetBinContent(h[i]->GetMaximumBin());
            if (temp >= maxEntry)   maxEntry = temp;
    }
    for (int i = 0; i < h.size(); i++)
        h[i]->GetYaxis()->SetRangeUser(1, maxEntry*1.2);
    TLine *l = new TLine(0, 0, 0, maxEntry);
	l -> Draw("same"); 
    l -> SetLineColor(kRed);
    l -> SetLineWidth(3);
    gPad->SetLogy();
    c1->Update();    c1->Modified();

    c1 -> SaveAs(Form("../../External/zFindingPlots/%s.png", fileName.c_str()));
}

/**
 * @brief 
 * 
 */

// void plotMultiCen () {
//     TCanvas *c1 = new TCanvas("c1","c1",0,50,1920,1056);
//     TProfile *Prof_Cen_Multi = new TProfile("Centrality v.s. Multiplicity", "Centrality v.s. Multiplicity, Fit with #bf{Least Square};Centrality;Multiplicity Density", 7, 0.0, 0.7, "S");
//     Prof_Cen_Multi->Fill(.0, 1129837.10);
//     Prof_Cen_Multi->Fill(.1, 890113.85);
//     Prof_Cen_Multi->Fill(.2, 554845.81);
//     Prof_Cen_Multi->Fill(.3, 352526.54);
//     Prof_Cen_Multi->Fill(.4, 200179.80);
//     Prof_Cen_Multi->Fill(.5, 111203.69);
//     Prof_Cen_Multi->Fill(.6, 52319.40);
//     // TProfile *Prof_Cen_Multi = new TProfile("Centrality v.s. Multiplicity", "Centrality v.s. Multiplicity, Fit with #bf{Total Least Square};Centrality;Multiplicity Density", 7, 0.0, 0.7, "S");
//     // Prof_Cen_Multi->Fill(.0, 960965.19);
//     // Prof_Cen_Multi->Fill(.1, 901032.76);
//     // Prof_Cen_Multi->Fill(.2, 546477.03);
//     // Prof_Cen_Multi->Fill(.3, 377208.87);
//     // Prof_Cen_Multi->Fill(.4, 156769.72);
//     // Prof_Cen_Multi->Fill(.5, 108235.72);
//     // Prof_Cen_Multi->Fill(.6, 54046.14);

//     Prof_Cen_Multi->SetLineColor(kBlue);
//     Prof_Cen_Multi->SetLineWidth(2);
//     Prof_Cen_Multi->SetMarkerSize(kBlue);
//     Prof_Cen_Multi->SetMarkerStyle(20);
//     Prof_Cen_Multi->SetMarkerSize(1.2);
//     Prof_Cen_Multi->GetXaxis()->CenterTitle(true);
//     Prof_Cen_Multi->GetYaxis()->CenterTitle(true);

//     Prof_Cen_Multi->Draw();

//     double bincontent, bincenter;
// 	auto graph_CenMulti = new TGraph();
// 	for (int b=1; b<=Prof_Cen_Multi->GetNbinsX(); b++) {
// 	   bincontent = Prof_Cen_Multi->GetBinContent(b);
// 	   if (bincontent) {
// 	      bincenter = Prof_Cen_Multi->GetBinCenter(b);
// 	      graph_CenMulti->AddPoint(bincenter,bincontent);
// 	   }
// 	}

//     graph_CenMulti->SetLineColor(kRed);
//     graph_CenMulti->SetLineWidth(2);
//     graph_CenMulti->SetMarkerColor(kRed);
//     graph_CenMulti->SetMarkerStyle(21);
//     graph_CenMulti->SetMarkerSize(1.2);
// 	graph_CenMulti->Draw("LP");

//     auto legend = new TLegend(0.12, 0.2, 0.35, 0.35); // Adjust coordinates as needed
//     legend->AddEntry(Prof_Cen_Multi, "INTT Centrality vs. Multiplicity", "lp");
//     legend->AddEntry(graph_CenMulti, "Least Square Fit", "lp");
//     // legend->AddEntry(graph_CenMulti, "Total Least Square Fit", "lp");
//     legend->Draw("same");
//     legend->SetTextSize(0.03);
//     legend->SetTextColor(kBlack);
//     legend->SetLineColor(kWhite);
//     legend->SetFillColor(kWhite);

//     c1->SaveAs("../External/zFindingPlots/Cen_Multi_LS.png");
//     // c1->SaveAs("../External/zFindingPlots/Cen_Multi_TLS.png");
// }

void plotMultiCen() {
    TCanvas *c1 = new TCanvas("c1", "Centrality vs. Multiplicity Comparison", 0, 50, 1920, 1056);
    
    // First dataset: Least Square Fit
    TProfile *Prof_Cen_Multi_LS = new TProfile("Centrality v.s. Multiplicity LS", "Centrality v.s. Multiplicity, different methods comparison;Centrality;Multiplicity Density", 7, 0.0, 0.7, "S");
    Prof_Cen_Multi_LS->Fill(.0, 1129837.10);
    Prof_Cen_Multi_LS->Fill(.1, 890113.85);
    Prof_Cen_Multi_LS->Fill(.2, 554845.81);
    Prof_Cen_Multi_LS->Fill(.3, 352526.54);
    Prof_Cen_Multi_LS->Fill(.4, 200179.80);
    Prof_Cen_Multi_LS->Fill(.5, 111203.69);
    Prof_Cen_Multi_LS->Fill(.6, 52319.40);
    Prof_Cen_Multi_LS->Draw();
    Prof_Cen_Multi_LS->GetXaxis()->CenterTitle(true);
    Prof_Cen_Multi_LS->GetYaxis()->CenterTitle(true);

    // Create a TGraph from the TProfile data for Least Square Fit
    auto graph_LS = new TGraph();
    for (int b = 1; b <= Prof_Cen_Multi_LS->GetNbinsX(); b++) {
        double bincontent = Prof_Cen_Multi_LS->GetBinContent(b);
        double bincenter = Prof_Cen_Multi_LS->GetBinCenter(b);
        graph_LS->AddPoint(bincenter, bincontent);
    }
    graph_LS->SetLineColor(6);
    graph_LS->SetLineWidth(2);
    graph_LS->SetMarkerColor(6);
    graph_LS->SetMarkerStyle(20);
    graph_LS->SetMarkerSize(1.2);

    // Second dataset: Total Least Square Fit
    TProfile *Prof_Cen_Multi_TLS = new TProfile("Centrality v.s. Multiplicity TLS", "Centrality v.s. Multiplicity, Fit with #bf{Total Least Square};Centrality;Multiplicity Density", 7, 0.0, 0.7, "S");
    Prof_Cen_Multi_TLS->Fill(.0, 960965.19);
    Prof_Cen_Multi_TLS->Fill(.1, 901032.76);
    Prof_Cen_Multi_TLS->Fill(.2, 546477.03);
    Prof_Cen_Multi_TLS->Fill(.3, 377208.87);
    Prof_Cen_Multi_TLS->Fill(.4, 156769.72);
    Prof_Cen_Multi_TLS->Fill(.5, 108235.72);
    Prof_Cen_Multi_TLS->Fill(.6, 54046.14);
    Prof_Cen_Multi_TLS->Draw("same");

    // Create a TGraph from the TProfile data for Total Least Square Fit
    auto graph_TLS = new TGraph();
    for (int b = 1; b <= Prof_Cen_Multi_TLS->GetNbinsX(); b++) {
        double bincontent = Prof_Cen_Multi_TLS->GetBinContent(b);
        double bincenter = Prof_Cen_Multi_TLS->GetBinCenter(b);
        graph_TLS->AddPoint(bincenter, bincontent);
    }
    graph_TLS->SetLineColor(7);
    graph_TLS->SetLineWidth(2);
    graph_TLS->SetMarkerColor(7);
    graph_TLS->SetMarkerStyle(21);
    graph_TLS->SetMarkerSize(1.2);

    // Draw the first graph (Least Square) and set axis titles
    graph_LS->Draw("LP SAME");
    graph_LS->GetXaxis()->SetTitle("Centrality");
    graph_LS->GetYaxis()->SetTitle("Multiplicity Density");

    // Draw the second graph (Total Least Square) on the same canvas
    graph_TLS->Draw("LP SAME");

    // Create and customize the legend
    auto legend = new TLegend(0.12, 0.2, 0.35, 0.35); // Adjust coordinates as needed
    legend->AddEntry(graph_LS, "Least Square Fit", "lp");
    legend->AddEntry(graph_TLS, "Total Least Square Fit", "lp");
    legend->SetTextSize(0.03);
    legend->SetTextColor(kBlack);
    legend->SetLineColor(kWhite);
    legend->SetFillColor(kWhite);
    legend->Draw();

    // Save the canvas with both profiles
    c1->SaveAs("../External/zFindingPlots/Cen_Multi_LS_TLS_Comparison.png");
}

void Npart () {
    TCanvas *c1 = new TCanvas("c1", "Npart vs. dN/d#eta/(<Npart>/2)", 0, 50, 1920, 1056);

    std::vector<std::vector<double>> data = readCsvToVector<std::vector<std::vector<double>>>("../External/dNchdEtaNpart.csv");
    auto graph_paper = new TGraph();
    graph_paper->AddPoint(353, 3.45);   // 0-3
    graph_paper->AddPoint(329, 3.34);   // 3-6
    graph_paper->AddPoint(291, 3.25);   // 6-10
    graph_paper->AddPoint(252, 3.16);   // 10-15
    graph_paper->AddPoint(215, 3.12);   // 15-20
    graph_paper->AddPoint(180, 3.08);   // 20-25
    graph_paper->AddPoint(149, 3.03);   // 25-30
    graph_paper->AddPoint(122, 3.00);   // 30-35
    graph_paper->AddPoint(102, 2.91);   // 35-40
    graph_paper->AddPoint(83, 2.87);    // 40-45

    graph_paper->SetLineColor(8);
    graph_paper->SetLineWidth(2);
    graph_paper->SetMarkerColor(8);
    graph_paper->SetMarkerStyle(20);
    graph_paper->SetMarkerSize(1.2);
    graph_paper->GetXaxis()->SetTitle("N_{part}");
    graph_paper->GetXaxis()->CenterTitle(true);
    graph_paper->GetYaxis()->SetTitle("dN_{ch}/d#eta/(N_{part}/2)");
    graph_paper->GetYaxis()->CenterTitle(true);

    auto graph_0052 = new TGraph();
    graph_0052->AddPoint(353, 941.965/353);   // 0-3
    graph_0052->AddPoint(329, 963.663/329);   // 3-6
    graph_0052->AddPoint(291, 899.042/291);   // 6-10
    graph_0052->AddPoint(252, 788.956/252);   // 10-15
    graph_0052->AddPoint(215, 653.359/215);   // 15-20
    graph_0052->AddPoint(180, 552.067/180);   // 20-25
    graph_0052->AddPoint(149, 445.922/149);   // 25-30
    graph_0052->AddPoint(122, 363.475/122);   // 30-35
    graph_0052->AddPoint(102, 293.811/102);   // 35-40
    graph_0052->AddPoint(83, 227.804/83);     // 40-45

    graph_0052->SetLineColor(2);
    graph_0052->SetLineWidth(2);
    graph_0052->SetMarkerColor(2);
    graph_0052->SetMarkerStyle(20);
    graph_0052->SetMarkerSize(1.2);
    graph_0052->GetXaxis()->SetTitle("N_{part}");
    graph_0052->GetXaxis()->CenterTitle(true);
    graph_0052->GetYaxis()->SetTitle("dN_{ch}/d#eta/(N_{part}/2)");  // Set the title of the x-axis
    graph_0052->GetYaxis()->CenterTitle(true);     // Center the x-axis title

    auto graph_020 = new TGraph();
    graph_020->AddPoint(353, 1385.11/353);   // 0-3
    graph_020->AddPoint(329, 1255.29/329);   // 3-6
    graph_020->AddPoint(291, 1078.82/291);   // 6-10
    graph_020->AddPoint(252, 927.49/252);   // 10-15
    graph_020->AddPoint(215, 805.64/215);   // 15-20
    graph_020->AddPoint(180, 654.26/180);   // 20-25
    graph_020->AddPoint(149, 557.81/149);   // 25-30
    graph_020->AddPoint(122, 452.34/122);   // 30-35
    graph_020->AddPoint(102, 367.35/102);   // 35-40
    graph_020->AddPoint(83, 278.72/83);     // 40-45
    // printCsvData(data);

    // graph_020->SetTitle("Hahahaha");
    graph_020->SetLineColor(4);
    graph_020->SetLineWidth(2);
    graph_020->SetMarkerColor(4);
    graph_020->SetMarkerStyle(20);
    graph_020->SetMarkerSize(1.2);
    graph_020->GetXaxis()->SetTitle("N_{part}");
    graph_020->GetXaxis()->CenterTitle(true);
    graph_020->GetYaxis()->SetTitle("dN_{ch}/d#eta/(N_{part}/2)");
    graph_020->GetYaxis()->CenterTitle(true);
    graph_020->GetXaxis()->SetRangeUser(80, 360);
    graph_020->GetYaxis()->SetRangeUser(2.6, 4.2);

    graph_paper->SetTitle("Pseudorapidity Density as a function of Npart");
    graph_020->SetTitle("Pseudorapidity Density as a function of Npart");
    graph_0052->SetTitle("Pseudorapidity Density as a function of Npart");
    graph_020->Draw("");
    graph_paper->Draw("LP same");
    // graph_0052->Draw("LP same");

    // Create and customize the legend
    auto legend = new TLegend(0.72, 0.75, 0.88, 0.89); // Adjust coordinates as needed
    legend->AddEntry(graph_paper, "Data from Paper", "lp");
    legend->AddEntry(graph_020, "Fit range: dPhi ~ [-0.2, 0.2]", "lp");
    legend->AddEntry(graph_0052, "Fit range: dPhi ~ [-#pi/6, #pi/6]", "lp");
    legend->SetTextSize(0.025);
    legend->SetTextColor(kBlack);
    legend->SetLineColor(kWhite);
    legend->SetFillColor(kWhite);
    legend->Draw();

    // Save the canvas with both profiles
    c1->SaveAs("../External/zFindingPlots/Pseudo_Npart.png");
}

// For showing DCA zFinding results:
void TSpectrumOutcomes (
    TH1D *const h,
    TH1  *const bg,
    const int &idx,
    const int &NClus,
    const double &foundZ,
    const double &MBD_z_vtx
) {
    TCanvas *c = new TCanvas("c", "c", 0,50,2000,1000); // It's a must!
    if (gPad) gPad->SetGrid(1, 1);
    // h -> DrawCopy("");
    c -> Update();
    h -> GetXaxis() -> CenterTitle(true);   h -> GetYaxis() -> CenterTitle(true);
    h -> SetTitle(Form("Found Z of Event #%d, NHits = %d, #bf{DCA with fit};z position [cm];# of counts", idx, NClus));
    h -> Draw("same");
    bg -> Draw("SAME");
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
    int max_entry = h -> GetBinContent(h -> FindBin(foundZ));
    // Update line positions and redraw
    l1->SetX1(foundZ); l1->SetX2(foundZ);
    l1->SetY1(0);   l1->SetY2(max_entry);
    l1->Draw("same");
    l2->SetX1(MBD_z_vtx); l2->SetX2(MBD_z_vtx);
    l2->SetY1(0);     l2->SetY2(max_entry);
    l2->Draw("same");
    static TLegend *lg = new TLegend(0.12, 0.85, 0.46, 0.9);
    lg -> AddEntry(h, Form("found z = %.4fcm, true z = %.4fcm", foundZ, MBD_z_vtx), "f");
    lg -> SetTextSize(.028);
    lg->Draw("same");
    gPad -> SetGrid(1,1); gPad -> Update();                 // It's a must;
}