/*
 *  DrawSignalExtraction.C
 *
 *  Created by Ophelie Bugnon on 15/10/19.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "SignalExtraction.C"

//const char* fileName = "AnalysisResults_15o_AOD229.root";
//const char* fileName = "AnalysisResults_18q_AOD225.root";
// const char* fileName = "AnalysisResults_18r_AOD225.root";
//const char* fileName = "AnalysisResults_18AOD225_Merged.root";
const char* fileName = "AnalysisResults_DataMerged.root";
// const char* fileName = "SubtractedHistos.root";
const char* fileLocation = "~/Documents/ALICE/AnalyseJPsi/AnalysisResults_5TeV";

Double_t min_y = -4;
Double_t max_y = -2.5;
Double_t min_pt = 0;
Double_t max_pt = 0.3;
Int_t min_cent = 70;
Int_t max_cent = 90;
Double_t min_mass = 0.;
Double_t max_mass = 10.;
Double_t min_fit = 2.2;
Double_t max_fit = 4.5;
Efunction f_BackGround = kVWGQuadratic;
Efunction f_Signal = kCBExtended;
Etails p_tails = kEMB;

void DrawSignalExtraction(const char* file=fileName, Double_t minPt=min_pt, Double_t maxPt=max_pt, Int_t minCent=min_cent, Int_t maxCent=max_cent, Efunction fBackground=f_BackGround, Efunction fSignal=f_Signal, Etails pTails=p_tails)
{
    std::vector<Double_t> resultsTest;
    resultsTest.clear();
    resultsTest = SignalExtraction(Form("%s/%s", fileLocation, file), min_y, max_y, minPt, maxPt, minCent, maxCent, min_mass, max_mass, fBackground, fSignal, pTails, min_fit, max_fit, kTRUE, kFALSE, kV0);
    if (resultsTest.size()==0) cout << "No signal extraction for this configuration" << endl;
}