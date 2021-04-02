/*
 *  HadronicNumberYield.C
 *
 *  Created by Ophelie Bugnon on 14/09/20.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>
 
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultsYield.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultsAcceptanceEfficiency.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultsSignalExtraction.C>

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>

Double_t nBranchinRatio = 0.05961;
Double_t errBranchinRatio = 0.5; //globale (%)
Double_t drapidité = 1.5;

//ranges pT 2 bins for incoherent
Double_t x[3]={0.15, 0.475, 0.825};
Double_t dx[3]={0.15, 0.175, 0.175};

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
//number of MB events
//de 0-10% à 70-90%
Double_t nMB[5]={581154046, 1162308093, 1162308093, 1162308093, 1162308093};
Double_t errMB[5]={1523538, 3047077, 3047077, 3047077, 3047077}; //globale absolue

void ExtractNJpsi(Int_t minCent, Int_t maxCent, Double_t minPt, Double_t maxPt, Double_t minY=-4, Double_t maxY=-2.5)
{
    int centClass;
    if(minCent==0) centClass=0;
    else if(minCent==10) centClass=1;
    else if(minCent==30) centClass=2;
    else if(minCent==50) centClass=3;
    else if(minCent==70) centClass=4;
    
    Double_t nJPsiHadro;
    Double_t errJPsiHadro_stat;
    Double_t errJPsiHadro_syst;
    Double_t errJPsiHadro_global;
    Double_t errJPsiHadro_systTot;
    
    errJPsiHadro_global =  TMath::Sqrt(TMath::Power((errMB[centClass]/nMB[centClass]), 2)  +  TMath::Power((errBranchinRatio/100), 2)  +  TMath::Power((EffLossTrack[centClass]/100),2) + TMath::Power((EffLossTrigg[centClass]/100),2) + TMath::Power((systCentralityLimits[centClass]/100),2))*100;
    cout << "Error global = " << errJPsiHadro_global << endl;

    std::vector<Double_t> vectYield = IntegratedYieldResults[SetRangeValue(kInvariantYield, minY, maxY, minPt, maxPt, minCent, maxCent)];
    std::vector<Double_t> vectAccEff = AccEffHadro[SetRangeValue(kAccEff, minY, maxY, minPt, maxPt, minCent, maxCent)];

    Double_t EffResponseMTR = SystEffMTR[SetNameSystAccEff(minY, maxY, minPt, maxPt)][0];
    Double_t MCinput_correlations = SystMCinput_correlations[SetNameSystAccEff(minY, maxY, minPt, maxPt)][0];

    nJPsiHadro=vectYield[kYieldValue]*(nBranchinRatio*nMB[centClass]*vectAccEff[kAccEffValue]*(drapidité*(maxPt-minPt)));
    errJPsiHadro_stat=TMath::Sqrt(TMath::Power((vectYield[kYieldStat]/vectYield[kYieldValue]), 2)  +  TMath::Power((vectAccEff[kAccEffStatError]/vectAccEff[kAccEffValue]), 2))*nJPsiHadro;

    errJPsiHadro_syst=TMath::Sqrt(TMath::Power((vectYield[kYieldSyst]/vectYield[kYieldValue]), 2)  +  TMath::Power(MC_inputVSpt_stat/100, 2) + TMath::Power(MCinput_correlations/100, 2) +  TMath::Power((EffMCH/100), 2) + TMath::Power((EffResponseMTR/100), 2) + TMath::Power((EffMTR/100), 2)  + TMath::Power(0.01, 2)) * nJPsiHadro;
    // cout << "Error uncorrelated systematic  = " << errJPsiHadro_syst << endl;

    errJPsiHadro_systTot = TMath::Sqrt(TMath::Power(errJPsiHadro_syst/nJPsiHadro, 2) + TMath::Power(errJPsiHadro_global/nJPsiHadro, 2))*nJPsiHadro;
    // cout << "Error systematic total = " << errJPsiHadro_systTot << endl;

    nJPsiHadro = nJPsiHadro + 0.5 - (nJPsiHadro<0);
    errJPsiHadro_stat = errJPsiHadro_stat + 0.5 - (errJPsiHadro_stat<0);
    errJPsiHadro_syst = errJPsiHadro_syst + 0.5 - (errJPsiHadro_syst<0);
    errJPsiHadro_systTot = errJPsiHadro_systTot + 0.5 - (errJPsiHadro_systTot<0);

    TString sName = SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent);
    // printf("{\"%s\", {%.f, %.f, %.f }},\n", sName.Data(), nJPsiHadro, errJPsiHadro_stat, errJPsiHadro_syst); 
    printf("{\"%s\", {%.f, %.f, %.f }},\n", sName.Data(), nJPsiHadro, errJPsiHadro_stat, errJPsiHadro_syst); 

    
}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void ExportResultsYield(Double_t minY, Double_t maxY)
{
    Int_t centRange[6] = {0, 10, 30, 50, 70, 90};

    for(int j=0; j<5; j++)
    {
        printf("\n");
        printf("//%d-%d\n", centRange[j], centRange[j+1]);
        for(int i=0; i<3; i++)
        {
            ExtractNJpsi(centRange[j], centRange[j+1], x[i]-dx[i], x[i]+dx[i], minY, maxY);
        } 
    }
}