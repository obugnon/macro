#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


void GetAccEffResults(Bool_t isDraw, Bool_t isExportToLatex)
{
    TFile* results;
    results = TFile::Open("/Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/AcceptanceEfficiency/AccEff_Coherent/AccEffi/AccEffiValues.root");
    if(!results) {
        Error("GetAccEffResults","No file found for acceptance efficiency results.")
        return;
    }    

    Double_t Cent, AccEff, dCent, dAccEff, dBinom;
    Double_t nRec, nGen;
    TString rangeName;
    Double_t minY=-4;
    Double_t maxY=-2.5;
    Double_t minPt, maxPt;
    Double_t minCent, maxCent; 

    TH1F* graphAccEff = (TH1F*)results->Get("histoAccEffiCent;1");
    TH1F* graphGen = (TH1F*)results->Get("histoGenCent;1");
    TH1F* graphRec = (TH1F*)results->Get("histoRecCent;1");

    Int_t nBins=graphAccEff->GetNbinsX();

    const Double_t tLowEdge[6]={0, 10, 30, 50, 70, 90};
    TH1D* graph = new TH1D("Coherent J/psi Acceptance Efficiency", "", nBins, tLowEdge);
    if(isDraw)
    {
        TCanvas *cAccEff = new TCanvas("","J/#psi Acceptance Efficiency");
        
        graph->GetXaxis()->SetTitle("centrality (%)");
        graph->GetYaxis()->SetTitle("A#times#epsilon");
        graph->SetMarkerColor(kBlue);
        graph->SetLineColor(kBlue);
        graph->SetFillColor(kBlue);
        graph->SetFillStyle(0);
    }


    for(int i=1; i<nBins+1; i++)
    {
        AccEff=graphAccEff->GetBinContent(i);
        Cent=graphAccEff->GetBinCenter(i);
        dCent=graphAccEff->GetBinWidth(i);

        nRec=graphRec->GetBinContent(i);
        nGen=graphGen->GetBinContent(i);

        if(nRec > nGen) 
        {
            dAccEff = 1/nGen;
            cout << " !!! nRec > nGen !!!" << endl;
        }    
        else 
        {
        dBinom = TMath::Sqrt(nRec/nGen * (1-nRec/nGen)/nGen);
        dAccEff = TMath::Max(dBinom,  1/nGen);
        }

        minPt=0;
        maxPt=0.3;
        minCent=Cent-dCent/2;
        maxCent=Cent+dCent/2;

        rangeName=SetRangeValue(kAccEff, minY, maxY, minPt, maxPt, minCent, maxCent);
        if(isExportToLatex) printf("%.f-%.f & $~ %.5f ~\\pm~ %.5f~(%.3f) ~$ \\\\\n", minCent, maxCent, AccEff, dAccEff, dAccEff/AccEff*100);
        else printf("{\"%s\", {%.9f, %.9f }},\n", rangeName.Data(), AccEff, dAccEff);
        if(isDraw)
        {
            graph->SetBinContent(i, AccEff);
            graph->SetBinError(i, dAccEff);
        }
        
    }
    if(isDraw)
    {
        graph->Draw();
        gStyle->SetOptStat(0);
    }
    else results->Close();
}
