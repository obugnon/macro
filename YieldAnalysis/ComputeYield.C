/*
 *  ComputeYield.C
 *
 *  Created by Ophelie Bugnon on 22/05/10.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsPPCrossSection.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsSignalExtraction.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsAcceptanceEfficiency.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


//nombre de pt bins pour chaque classe en centralité de 0-10% à 70-90%
// int nRanges[5]={13,13,13, 12, 11}; // 1 bin for incoherent
int nRanges[5]={14,14,14, 13, 12}; // 2 bin for incoherent

Double_t nBranchinRatio = 0.05961;
Double_t errBranchinRatio = 0.5; //globale (%)
Double_t drapidité = 1.5;

//-------------------------------------------------------------------------------------------------
//tableaux de la valeur moyenne du bin, de son ecart type et de sa valeur inférieure
std::vector<std::vector<Double_t>> xRange{
    {0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5},
    {0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5},
    {0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5},
    {0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11, 13.5},
    {0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11}
};

std::vector<std::vector<Double_t>> dxRange{
    {0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5},
    {0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5},
    {0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5},
    {0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1.5},
    {0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1}
};
//-------------------------------------------------------------------------------------------------
//number of MB events
//de 0-10% à 70-90%
    Double_t nMB[5]={581154046, 1162308093, 1162308093, 1162308093, 1162308093};
    Double_t errMB[5]={1523538, 3047077, 3047077, 3047077, 3047077}; //globale absolue

//-------------------------------------------------------------------------------------------------
void ComputeYield(Double_t minY, Double_t maxY, Double_t minPt, Double_t maxPt, Int_t minCent, Int_t maxCent)
{
    Int_t iCent=0;
    if (minCent==0) iCent=0;
    else if (minCent==10) iCent=1;
    else if (minCent==30) iCent=2;
    else if (minCent==50) iCent=3;
    else if (minCent==70) iCent=4;
    
    Double_t nYield;
    Double_t errYield_stat;
    Double_t errYield_syst;
    Double_t errYield_global;

    errYield_global =  TMath::Sqrt(TMath::Power((errMB[iCent]/nMB[iCent]), 2)  +  TMath::Power((errBranchinRatio/100), 2)  +  TMath::Power((EffLossTrack[iCent]/100),2) + TMath::Power((EffLossTrigg[iCent]/100),2) + TMath::Power((systCentralityLimits[iCent]/100),2))*100;

    printf("Global uncertainty is %.3f %% \n", errYield_global);

    std::vector<Double_t> vectNJpsi = SignalResultsHadro[SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent)];

    std::vector<Double_t> vectAccEff = AccEffHadro[SetRangeValue(kAccEff, minY, maxY, minPt, maxPt, minCent, maxCent)];

    Double_t EffResponseMTR = SystEffMTR[SetNameSystAccEff(minY, maxY, minPt, maxPt)][0];
    Double_t MCinput_correlations = SystMCinput_correlations[SetNameSystAccEff(minY, maxY, minPt, maxPt)][0];
    
    nYield=vectNJpsi[kJpsiMean]/(nBranchinRatio*nMB[iCent]*vectAccEff[kAccEffValue]*(drapidité*(maxPt-minPt)));
    errYield_stat=TMath::Sqrt(TMath::Power((vectNJpsi[kJpsiStatError]/vectNJpsi[kJpsiMean]), 2)  +  TMath::Power((vectAccEff[kAccEffStatError]/vectAccEff[kAccEffValue]), 2))*nYield;

    errYield_syst=TMath::Sqrt(TMath::Power((vectNJpsi[kJpsiSysError]/vectNJpsi[kJpsiMean]), 2)  + TMath::Power(MC_inputVSpt_stat/100, 2) + TMath::Power(MCinput_correlations/100, 2)  +  TMath::Power((EffMCH/100), 2)  + TMath::Power((EffResponseMTR/100), 2) + TMath::Power((EffMTR/100), 2)  + TMath::Power(0.01, 2)) * nYield;

    TString sName = SetRangeValue(kInvariantYield, minY, maxY, minPt, maxPt, minCent, maxCent);
    printf("{\"%s\", {%.8f, %.8f, %.8f }},\n", sName.Data(), nYield, errYield_stat, errYield_syst);
    
}

void ExportResultsYield(Double_t minY, Double_t maxY, Int_t minCent, Int_t maxCent)
{
    Int_t iCent=0;
    if (minCent==0) iCent=0;
    else if (minCent==10) iCent=1;
    else if (minCent==30) iCent=2;
    else if (minCent==50) iCent=3;
    else if (minCent==70) iCent=4;
    for(int i=0; i<nRanges[iCent]; i++)
    {
    
        ComputeYield(minY, maxY, xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i], minCent, maxCent);
 
    }
}