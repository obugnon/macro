/*
 *  DrawFiguresForAcceptanceEfficiency.C
 *
 *  Created by Ophelie Bugnon on 04/06/20.
 *
*/
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "TMath.h"

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsAcceptanceEfficiency.C>


Double_t pT_shape(Double_t *x, Double_t *par)
{
    // power-law
    Double_t A = par[0];
    Double_t B = par[1];
    Double_t n1 = par[2];
    Double_t n2 = par[3];
    // pt
    Double_t pt = x[0];
    return (A * pt) / pow(1 + pow(pt / B, n1), n2);
}

Double_t y_shape(Double_t *x, Double_t *par)
{
    Double_t A = par[0];
    Double_t B = par[1];
    //y
    Double_t y = x[0];
    return A*exp(-0.5*pow((y/B),2));
}

Int_t centrality[6]={0, 10, 30, 50, 70, 90};
Int_t iteration[5]={2, 2, 2, 3, 2};
int lcolor[5]={12, 9, 32, 42, 46};

std::vector<std::vector<Double_t>> PtParameters = {
    {1, 3.02547, 2.35018, 3.03864},
    {1, 3.20589, 2.20675, 3.22158},
    {1, 3.71266, 1.99116, 3.7466},
    {1, 4.11143, 1.8616, 4.14776},
    {1, 8.63433, 1.23962, 10.5712},
};
std::vector<std::vector<Double_t>> YParameters = {
    {1, 2.14206},
    {1, 2.19626},
    {1, 2.38394},
    {1, 2.48649},
    {1, 2.65303},
}; 

void DrawGenerationFonctions()
{
    //Graph vs Pt
    TCanvas *cPtFunction = new TCanvas("cPtFunction", "", 900, 800);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);    
    cPtFunction->SetLogy();
    TH1F* hPtFunction=new TH1F("hPtFunction", " ;#it{p}_{T} (GeV/#it{c});dN_{J/#psi}/d#it{p}_{T}", 150, 0, 15);
    hPtFunction->Draw("SAME");

    TLatex * textPt = new TLatex (10,0.1,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    textPt->Draw("SAME");
    textPt->SetTextSizePixels(24);
    TLatex * text1Pt = new TLatex (10,0.05,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    text1Pt->SetTextSizePixels(22);
    text1Pt->Draw("SAME");
    TLatex * text2Pt = new TLatex (10,0.01,"A#times#frac{#it{p}_{T}}{ (1+ (#it{p}_{T}/B)^{n_{1}})^{n_{2}}  }");
    text2Pt->SetTextSizePixels(22);
    text2Pt->Draw("SAME");

    //Graph vs Y
    TCanvas *cYFunction = new TCanvas("cYFunction", "", 900, 800);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TH1F* hYFunction=new TH1F("hYFunction", " ;#it{y};dN_{J/#psi}/d#it{y}", 15, -4, -2.5);
    hYFunction->Draw();

    TLatex * textY = new TLatex (-3.9,0.92,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    textY->Draw("SAME");
    textY->SetTextSizePixels(22);
    TLatex * text1Y = new TLatex (-3.9,0.88,"J/#psi #rightarrow #mu^{+}#mu^{-}, 0.3 < #it{p}_{T} < 15 GeV/#it{c}");
    text1Y->SetTextSizePixels(20);
    text1Y->Draw("SAME");
    TLatex * text2Y = new TLatex (-3.9,0.8,"A#times e^{-0.5 #times #left(#frac{y}{B} #right)^{2}}");
    text2Y->SetTextSizePixels(22);
    text2Y->Draw("SAME");

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);

    
    for(int i=0; i<5; i++)
    {
        TF1 *fPtDistribution = new TF1("fPtDistribution", pT_shape, 0, 15, 4);
        TF1 *fRapDistribution = new TF1("fRapDistribution", y_shape, -4, -2.5, 2);

        for(int j=0; j<4; j++)
        {
        fPtDistribution->SetParameter(j, PtParameters[i][j]);
        cout << "Parameter " << j << " = " <<  PtParameters[i][j] << endl;
        }
        for (int k=0; k<2; k++)
        {
        fRapDistribution->SetParameter(k, YParameters[i][k]);
        }
        
        fPtDistribution->SetLineColor(lcolor[i]);
        fRapDistribution->SetLineColor(lcolor[i]);
        legend->AddEntry(fPtDistribution,Form("%i-%i %% : iteration %i",centrality[i], centrality[i+1], iteration[i]),"l");

        cPtFunction->cd();
        fPtDistribution->Draw("SAME");

        cYFunction->cd();
        fRapDistribution->Draw("same");        
        
    }
    cPtFunction->cd();
    legend->Draw();

    cYFunction->cd();
    legend->Draw();
}


// std::vector<std::vector<Double_t>> PtRange;
// xRange.push_back({0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5});
// xRange.push_back({0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11, 13.5});
// xRange.push_back({0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11});

// std::vector<std::vector<Double_t>> dpTRange;
// dxRange.push_back({0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5});
// dxRange.push_back({0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1.5});
// dxRange.push_back({0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1});

// Double_t YRange[6] = {-3.875, -3.625, -3.375, -3.125, -2.875, -2.625};
// Double_t dYRange[6] = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125};
void DrawAcceptanceEfficiency()
{
    //Graph vs Pt
    TCanvas *cPtFunction = new TCanvas("cPtFunction", "", 900, 800);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);    
    TH1F* hPtFunction=new TH1F("hPtFunction", " ;#it{p}_{T} (GeV/#it{c});A#times#epsilon_{J/#psi}", 15, 0, 15);
    hPtFunction->GetYaxis()->SetRangeUser(0, 0.5);
    hPtFunction->Draw("SAME");

    TLatex * textPt = new TLatex (0.5,0.45,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    textPt->Draw("SAME");
    textPt->SetTextSizePixels(24);
    TLatex * text1Pt = new TLatex (0.5,0.41,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    text1Pt->SetTextSizePixels(22);
    text1Pt->Draw("SAME");

    //Graph vs Y
    TCanvas *cYFunction = new TCanvas("cYFunction", "", 900, 800);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TH1F* hYFunction=new TH1F("hYFunction", " ;#it{y};A#times#epsilon_{J/#psi}", 6, -4, -2.5);
    hYFunction->GetYaxis()->SetRangeUser(0, 0.3);
    hYFunction->Draw();

    TLatex * textY = new TLatex (-3.9,0.28,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    textY->Draw("SAME");
    textY->SetTextSizePixels(24);
    TLatex * text1Y = new TLatex (-3.9,0.26,"J/#psi #rightarrow #mu^{+}#mu^{-}, 0.3 < #it{p}_{T} < 15 GeV/#it{c}");
    text1Y->SetTextSizePixels(22);
    text1Y->Draw("SAME");


    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);

    TFile* results;
    
    for(int i=0; i<5; i++)
    {
        results = TFile::Open(Form("/Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/AcceptanceEfficiency/AccEff_Hadro_Finer_binning/Cent-%ito%i/AccEffi/iter-%i/AccEffiValues.root", centrality[i], centrality[i+1], iteration[i]));
        if(!results) 
        {
            Error("GetAccEff","Cannot open Acceptance Efficiency results file");
            return ;
        }

        cPtFunction->cd();
        TH1F* graphAccEffvsPt = (TH1F*)results->Get("histoAccEffiVsPt;1");
        graphAccEffvsPt->SetLineColor(lcolor[i]);
        graphAccEffvsPt->SetMarkerColor(lcolor[i]);
        graphAccEffvsPt->SetMarkerStyle(kFullCircle);
        graphAccEffvsPt->SetMarkerSize(0.9);
        graphAccEffvsPt->SetFillStyle(0);
        graphAccEffvsPt->Draw("epsame");


        cYFunction->cd();
        TH1F* graphAccEffvsY = (TH1F*)results->Get("histoAccEffiVsRap;1");
        graphAccEffvsY->SetLineColor(lcolor[i]);
        graphAccEffvsY->SetMarkerColor(lcolor[i]);
        graphAccEffvsY->SetMarkerStyle(kFullCircle);
        graphAccEffvsY->SetMarkerSize(0.9);
        graphAccEffvsY->SetFillStyle(0);
        graphAccEffvsY->Draw("epsame");
        legend->AddEntry(graphAccEffvsPt,Form("%i-%i %% : iteration %i",centrality[i], centrality[i+1], iteration[i]),"l");     
        
    }
    cPtFunction->cd();
    legend->Draw();

    cYFunction->cd();
    legend->Draw();
}