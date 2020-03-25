#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"

void DrawResults()
{
    Double_t x[5]={5, 20, 40, 60, 80};
    Double_t dx[5]={5, 10, 10, 10, 10};

    Double_t y1[5]={8117, 9257, 4252, 2705, 1739};
    Double_t y1_syst[5]={733, 575, 164, 86, 43};
    Double_t y1_stat[5]={668, 523, 215, 91, 52}; 

    Double_t y2[5]={66029, 64295, 20190, 5337, 1501};
    Double_t y2_syst[5]={3116, 2274, 519, 130, 30};
    Double_t y2_stat[5]={1830, 1523, 622, 190, 62};

    Double_t y3[5]={289291, 303150, 108777, 31839, 6322};
    Double_t y3_syst[5]={5681, 5610, 2041, 554, 111};
    Double_t y3_stat[5]={2655, 2238, 964, 348, 114};

    TCanvas *cCrossSection = new TCanvas("c1","Raw number of Jpsi");
    TH1F *fOption = new TH1F("Cross section distribution","", 9, 0, 90);
  	fOption->GetYaxis()->SetTitle("Nb of J/#psi");
  	fOption->GetXaxis()->SetTitle("cent (%)");
  	fOption->SetMaximum(350000);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
	fOption->Draw();

//0-0.3
    TGraphErrors* gr1_syst = new TGraphErrors(5,x,y1,dx,y1_syst);
    gr1_syst->SetLineColor(46);
  	gr1_syst->SetFillColor(46);
  	gr1_syst->SetFillStyle(0);
  	// gr1_syst->SetLineWidth(2);
  	gr1_syst->Draw("5same");

    TGraphErrors* gr1_stat = new TGraphErrors(5,x,y1,0,y1_stat);
    gr1_stat->SetMarkerColor(46);
  	gr1_stat->SetMarkerStyle(2);
  	// gr1_stat->SetMarkerSize(2);
  	gr1_stat->SetLineColor(46);
  	// gr1_stat->SetLineWidth(2);
  	gr1_stat->Draw("psame");
//0.3-1
    TGraphErrors* gr2_syst = new TGraphErrors(5,x,y2,dx,y2_syst);
    gr2_syst->SetLineColor(36);
  	gr2_syst->SetFillColor(36);
  	gr2_syst->SetFillStyle(0);
  	// gr2_syst->SetLineWidth(2);
  	gr2_syst->Draw("5same");

    TGraphErrors* gr2_stat = new TGraphErrors(5,x,y2,0,y2_stat);
    gr2_stat->SetMarkerColor(36);
  	gr2_stat->SetMarkerStyle(2);
  	// gr2_stat->SetMarkerSize(2);
  	gr2_stat->SetLineColor(36);
  	// gr2_stat->SetLineWidth(2);
  	gr2_stat->Draw("psame");
//1-8
TGraphErrors* gr3_syst = new TGraphErrors(5,x,y3,dx,y3_syst);
    gr3_syst->SetLineColor(9);
  	gr3_syst->SetFillColor(9);
  	gr3_syst->SetFillStyle(0);
  	// gr3_syst->SetLineWidth(2);
  	gr3_syst->Draw("5same");

    TGraphErrors* gr3_stat = new TGraphErrors(5,x,y3,0,y3_stat);
    gr3_stat->SetMarkerColor(9);
  	gr3_stat->SetMarkerStyle(2);
  	// gr3_stat->SetMarkerSize(2);
  	gr3_stat->SetLineColor(9);
  	// gr3_stat->SetLineWidth(2);
  	gr3_stat->Draw("psame");

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
    // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry((TObject*)0, "Pb-Pb collisions #sqrt{s_{NN} }= 5.02 TeV", "");
    legend->AddEntry(gr3_stat,"1 - 8 GeV/c","l");	
    legend->AddEntry(gr2_stat,"0.3 - 1 GeV/c","l");
    legend->AddEntry(gr1_stat,"0 - 0.3 GeV/c","l");
    legend->Draw();
}