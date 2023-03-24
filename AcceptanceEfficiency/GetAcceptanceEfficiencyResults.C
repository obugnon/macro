/*
 *  SignalExtraction.C
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
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultsAcceptanceEfficiency.C>



void GetAccEff(int minCent, int maxCent, int iter, Bool_t isOtherRange)
{
    TFile* results;
    if (isOtherRange) results = TFile::Open(Form("/Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/AcceptanceEfficiency/AccEff_Hadro_Finer_binning/Cent-%ito%i/AccEffi/iter-%i/AccEffiValues.root", minCent, maxCent, iter));
    else results = TFile::Open(Form("/Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/AcceptanceEfficiency/AccEff_Hadro/Cent-%ito%i/AccEffi/iter-%i/AccEffiValues.root", minCent, maxCent, iter));

    if(!results) 
    {
        Error("GetAccEff","Cannot open Acceptance Efficiency results file");
        return ;
    }

    TFile* results_syst;
    if (isOtherRange) results_syst = TFile::Open(Form("/Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/AcceptanceEfficiency/AccEff_Hadro_Finer_binning/Cent-%ito%i/AccEffSyst/AccEffiValues.root", minCent, maxCent));
    else results_syst = TFile::Open(Form("/Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/AcceptanceEfficiency/AccEff_Hadro/Cent-%ito%i/AccEffSyst/AccEffiValues.root", minCent, maxCent));

    if(!results_syst) 
    {
        Error("GetAccEff","Cannot open Acceptance Efficiency systematic uncertainty file");
        return ;
    }   

    Double_t meanPt, nAccEff, dPt, errStat, errSyst;
    TString rangeName;
    Double_t minY=-4;
    Double_t maxY=-2.5;
    Double_t minPt, maxPt;

    Double_t errStatLow, errStatUp, errStatBinom;
    Double_t nRec, nGen, dBinom;

    Double_t xSyst, ySyst;


    TH1F* graphAccEff = (TH1F*)results->Get("histoAccEffiVsPt;1");
    const Int_t nBins=graphAccEff->GetNbinsX();
    // TGraph* graphStat = (TGraph*)results_stat->Get("acceff_rms_pt");
    TGraph* graphSyst;

    Int_t nTests = 60;
    TH1F* graphGen = (TH1F*)results->Get("histoGenPt;1");
    TH1F* graphRec = (TH1F*)results->Get("histoRecPt;1");
    // TGraphAsymmErrors* graphAccEffForStat = new TGraphAsymmErrors(graphRec, graphGen);

    
    for(int i=1; i<nBins+1; i++)
    {
        errSyst=0;
        for (int j=1; j<nTests+1; j++)
        {
            graphSyst = (TGraph*) results_syst->Get(Form("acceff_variation_over_ref_pt;%i", j));
            graphSyst->GetPoint(i-1, xSyst, ySyst);
            if(errSyst != TMath::Max(errSyst, TMath::Abs(ySyst-1)))
            {
                errSyst = TMath::Max(errSyst, TMath::Abs(ySyst-1));
            }
        }

        nAccEff=graphAccEff->GetBinContent(i);
        meanPt=graphAccEff->GetBinCenter(i);
        dPt=graphAccEff->GetBinWidth(i);
        errStat = graphAccEff->GetBinError(i);

        // errStatLow = graphAccEffForStat->GetErrorYlow(i);
        // errStatUp = graphAccEffForStat->GetErrorYhigh(i);
        // errStatBinom = TMath::Max(errStatLow,  errStatUp);
        nRec=graphRec->GetBinContent(i);
        nGen=graphGen->GetBinContent(i);
        if(nRec > nGen) dBinom = 1/nGen;
        else 
        {
            dBinom = TMath::Sqrt(nRec/nGen * (1-nRec/nGen)/nGen);
            dBinom = TMath::Max(dBinom,  1/nGen);
        }
    
        minPt=meanPt-dPt/2;
        maxPt=meanPt+dPt/2;
        // cout << "For pt from " << minPt << " to " << maxPt << endl ;
        // cout << "Acc Eff = " << nAccEff << endl;
        // cout << "err stat = " << errStat << " represents " << errStat/nAccEff*100 << "%% at x = "<< meanPt << endl;
        // cout << "err stat with TAsymmErrors = " << errStatBinom << " represents " << errStatBinom/nAccEff*100 << "%% at x = "<< meanPt << endl;
        // cout << "err stat with binomial error = " << dBinom << " represents " << dBinom/nAccEff*100 << "%% at x = "<< meanPt << endl;
        // cout << "err syst = " << errSyst*nAccEff << " represents " << errSyst*100 << "%% at x = " << xSyst <<endl;

        rangeName=SetRangeValue(kAccEff, minY, maxY, minPt, maxPt, minCent, maxCent);
        rangeName=SetNameSystAccEff(minY, maxY, minPt, maxPt);
        printf("{\"%s\", {%.2f}},\n", rangeName.Data(), errSyst*100);

    }
    results->Close();
}

void GetAccEffFull(Bool_t isOtherRange)
{
    int tminCent[5]={0, 10, 30, 50, 70};
    int tmaxCent[5]={10, 30, 50, 70, 90};
    int tIter[5]={2,2,2,3,2};

    for(int j=0; j<5; j++)
    {
        printf("//%d-%d\n", tminCent[j], tmaxCent[j]);
        GetAccEff(tminCent[j], tmaxCent[j], tIter[j], isOtherRange);
        printf("\n");
    }
}

void DrawAccEffVScent()
{
    Double_t tminCent[6]={0, 10, 30, 50, 70, 90};
    int tmaxCent[5]={10, 30, 50, 70, 90};
    TH1D* graph = new TH1D("Hadronic J/psi Acceptance Efficiency", "", 5, tminCent);

    TCanvas *cAccEff = new TCanvas("","J/#psi Acceptance Efficiency");
    
    graph->GetXaxis()->SetTitle("centrality (%)");
    graph->GetYaxis()->SetTitle("A#times#epsilon");
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->SetFillColor(kBlue);
    graph->SetFillStyle(0);

    Double_t nAccEff, errStat;
    TString rangeName;

    for(int j=0; j<5; j++)
    {
        rangeName=SetRangeValue(kAccEff, -4, -2.5, 0, 0.3, tminCent[j], tminCent[j+1]);
        nAccEff = AccEffHadro[rangeName][kAccEffValue];
        errStat=AccEffHadro[rangeName][kAccEffStatError];
        graph->SetBinContent(j, nAccEff);
        graph->SetBinError(j, errStat);
    }
    graph->Draw();    
    
}
