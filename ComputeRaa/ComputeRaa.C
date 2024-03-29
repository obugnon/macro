/*
 *  ComputeRaa.C
 *
 *  Created by Ophelie Bugnon on 22/06/20.
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

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsPPCrossSection.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsSignalExtraction.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsAcceptanceEfficiency.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsRaa.C>


//nombre de pt bins pour chaque classe en centralité de 0-10% à 70-90%
int nRanges[5]={15, 15, 15, 14, 13}; 

Double_t nBranchinRatio = 0.05961;
Double_t errBranchinRatio = 0.5; //globale (%)
Double_t drapidité = 1.5;

//-------------------------------------------------------------------------------------------------
//tableaux de la valeur moyenne du bin, de son ecart type et de sa valeur inférieure
std::vector<std::vector<Double_t>> xRange{
    {0.15, 0.65, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5},
    {0.15, 0.65, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5},
    {0.15, 0.65, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5},
    {0.15, 0.65, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11, 13.5},
    {0.15, 0.65, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11}
};

std::vector<std::vector<Double_t>> dxRange{
    {0.15, 0.35, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5},
    {0.15, 0.35, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5},
    {0.15, 0.35, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5},
    {0.15, 0.35, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1.5},
    {0.15, 0.35, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1}
};

//-------------------------------------------------------------------------------------------------
//number of MB events
//de 0-10% à 70-90%
    Double_t nMB[5]={581154046, 1162308093, 1162308093, 1162308093, 1162308093};
    Double_t errMB[5]={1523538, 3047077, 3047077, 3047077, 3047077}; //absolute value

//-------------------------------------------------------------------------------------------------
//Taa
//de 0-10% à 70-90%
    Double_t nTaa[5]={23.26, 11.5835, 3.9165, 0.9756, 0.161165}; 
    Double_t errTaa[5]={0.168, 0.1135, 0.065, 0.02335, 0.00365}; //absolute value

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
void ComputeRAAvsPt(Double_t minY, Double_t maxY, Double_t minPt, Double_t maxPt, Int_t minCent, Int_t maxCent)
{
    Int_t iCent=0;
    if (minCent==0) iCent=0;
    else if (minCent==10) iCent=1;
    else if (minCent==30) iCent=2;
    else if (minCent==50) iCent=3;
    else iCent=4;
    

    Double_t nRaa;
    Double_t errRaa_stat;
    Double_t errRaa_syst;
    Double_t errRaa_global;

    errRaa_global =  TMath::Power((errMB[iCent]/nMB[iCent]), 2) - TMath::Power((errBranchinRatio/100), 2) + TMath::Power((systCentralityLimits[iCent]/100),2) + TMath::Power((errTaa[iCent]/nTaa[iCent]), 2) +  TMath::Power((EffLossTrack[iCent]/100),2) + TMath::Power((EffLossTrigg[iCent]/100),2) + TMath::Power((systGlobCRpp/100), 2);
    errRaa_global = TMath::Sqrt(errRaa_global)*100;

    printf("Global error is %.3f %% \n", errRaa_global);
    // printf("Global error on CRpp is %.3f %% \n", systGlobCRpp);
    // printf("Global error on nMB is %.3f %% \n", errMB[iCent]/nMB[iCent]*100);
    // printf("Global error on Taa is %.3f %% \n", errTaa[iCent]/nTaa[iCent]*100);
    // printf("Global error on Eff loss track is %.3f %% \n", EffLossTrack[iCent]);
    // printf("Global error on Eff loss trig is %.3f %% \n", EffLossTrigg[iCent]/100);
    // printf("Global error on BR is %.3f %% \n", errBranchinRatio);
     
    std::vector<Double_t> vectNJpsi = SignalResultsHadro[SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent)];

    std::vector<Double_t> vectAccEff = AccEffHadro[SetRangeValue(kAccEff, minY, maxY, minPt, maxPt, minCent, maxCent)];
    Double_t effMTRresponse = SystEffMTR[SetNameSystAccEff(minY, maxY, minPt, maxPt)][0];
    Double_t MCinput_correlations = SystMCinput_correlations[SetNameSystAccEff(minY, maxY, minPt, maxPt)][0];

    std::vector<Double_t> vectCRpp = ppCrossSection[SetRangeValue(kCRpp, minY, maxY, minPt, maxPt)];

    nRaa=vectNJpsi[kJpsiMean]/(nBranchinRatio*nMB[iCent]*vectAccEff[kAccEffValue]*nTaa[iCent]*TMath::Power(10,-3)*(vectCRpp[kCRppValue]*drapidité*(maxPt-minPt)));
    
    // cout << "nJpsi = " << vectNJpsi[kJpsiMean] << " and nMB = " << nMB[iCent] << " and AccEff = " << vectAccEff[kAccEff] << " and Taa = " <<  nTaa[iCent] << " and ref pp = " << (vectCRpp[kCRpp]*drapidité*(maxPt-minPt)) << endl;

    errRaa_stat=TMath::Sqrt(TMath::Power((vectNJpsi[kJpsiStatError]/vectNJpsi[kJpsiMean]), 2)  +  TMath::Power((vectAccEff[kAccEffStatError]/vectAccEff[kAccEffValue]), 2)  +  TMath::Power((vectCRpp[kCRppStatError]/vectCRpp[kCRppValue]), 2))*nRaa;

    errRaa_syst=TMath::Sqrt(TMath::Power((vectNJpsi[kJpsiSysError]/vectNJpsi[kJpsiMean]), 2) + TMath::Power(MC_inputVSpt_stat/100, 2) + TMath::Power(MCinput_correlations/100, 2) +  TMath::Power((EffMCH/100), 2)  + TMath::Power((effMTRresponse/100), 2) + TMath::Power((EffMTR/100), 2)  + TMath::Power(0.01, 2) + TMath::Power((vectCRpp[kCRppSystError]/vectCRpp[kCRppValue]), 2)) * nRaa;

    // printf("Uncorrelated error is %.3f %% \n", errRaa_syst/nRaa*100);
    // printf("Uncorrelated error on Njpsi is %.3f %% \n", vectNJpsi[kJpsiSysError]/vectNJpsi[kJpsiMean]*100);
    // printf("Uncorrelated error on MC input is %.3f %% \n", MC_inputVSpt[iCent]);
    // printf("Uncorrelated error on Eff track is %.3f %% \n", EffMCH);
    // printf("Uncorrelated error on Eff trig is %.3f %% \n", TMath::Sqrt(effMTRresponse*effMTRresponse+EffMTR*EffMTR));
    // printf("Uncorrelated error on matching is %.3f %% \n", 0.01*100);
    // printf("Uncorrelated error on CRpp is %.3f %% \n", vectCRpp[kCRppSystError]/vectCRpp[kCRppValue]*100);

    
    TString sName = SetRangeValue(kRaa, minY, maxY, minPt, maxPt, minCent, maxCent);
    printf("{\"%s\", {%.8f, %.8f, %.8f }},\n", sName.Data(), nRaa, errRaa_stat, errRaa_syst);
}
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
void ExportResultsRAAvsPt(Double_t minY, Double_t maxY)
{
    Int_t cent[6]={0, 10, 30, 50, 70, 90};
    for(int j=0; j<5; j++)
    {
        printf("\n//%i-%i%%\n", cent[j], cent[j+1]);
        Int_t iCent=0;
        if (cent[j]==0) iCent=0;
        else if (cent[j]==10) iCent=1;
        else if (cent[j]==30) iCent=2;
        else if (cent[j]==50) iCent=3;
        else iCent=4;

        for(int i=0; i<nRanges[iCent]; i++)
        {
            ComputeRAAvsPt(minY, maxY, xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i], cent[j], cent[j+1]);
        }
    }     
}

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
void ComputeRAAvsCent(Double_t minY, Double_t maxY, Double_t minPt, Double_t maxPt, Int_t minCent, Int_t maxCent)
{
    Int_t iCent=0;
    if (minCent==0) iCent=0;
    else if (minCent==10) iCent=1;
    else if (minCent==30) iCent=2;
    else if (minCent==50) iCent=3;
    else iCent=4;


    Double_t nRaa;
    Double_t errRaa_stat;
    Double_t errRaa_syst;
    Double_t errRaa_global;

    std::vector<Double_t> vectNJpsi = SignalResultsLowPt[SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent)];
    std::vector<Double_t> vectAccEff = AccEffHadro[SetRangeValue(kAccEff, minY, maxY, minPt, maxPt, minCent, maxCent)];
    std::vector<Double_t> vectCRpp = ppCrossSection[SetRangeValue(kCRpp, minY, maxY, minPt, maxPt)];
    Double_t effMTRresponse = SystEffMTR[SetNameSystAccEff(minY, maxY, minPt, maxPt)][0];
    Double_t MCinput_correlations = SystMCinput_correlations[SetNameSystAccEff(minY, maxY, minPt, maxPt)][0];


    // printf("nJpsi = %.f and nMB = %.f and AccEff = %.2f and Taa = %.2f and ref pp = %.2f" ,vectNJpsi[kJpsiMean], nMB[iCent], vectAccEff[kAccEff], nTaa[iCent], (vectCRpp[kCRpp]*drapidité*(maxPt-minPt)));

    errRaa_global =  TMath::Power((errMB[iCent]/nMB[iCent]), 2) - TMath::Power((errBranchinRatio/100), 2) +  TMath::Power((EffMCH/100), 2)  + TMath::Power((effMTRresponse/100), 2) + TMath::Power((EffMTR/100), 2) + TMath::Power(0.01, 2) + TMath::Power((systGlobCRpp/100), 2) + TMath::Power((vectCRpp[kCRppSystError]/vectCRpp[kCRppValue]), 2);

    errRaa_global = TMath::Sqrt(errRaa_global)*100;

    printf("Global error is %.3f %% \n", errRaa_global);
    // printf("Global error on nMB is %.3f %% \n", errMB[iCent]/nMB[iCent]*100);
    // printf("Global error on BR is %.3f %% \n", errBranchinRatio);
    // printf("Global error on MC input is %.3f %% \n", MC_inputVScent[iPt]);
    // printf("Global error on Eff track is %.3f %% \n", EffMCH);
    // printf("Global error on Eff trig is %.3f %% \n", TMath::Sqrt(effMTRresponse*effMTRresponse+EffMTR*EffMTR));
    // printf("Global error on matching is %.3f %% \n", 0.01*100);
    // printf("Global error on CRpp is %.3f %% \n", TMath::Sqrt(TMath::Power((systGlobCRpp/100), 2) + TMath::Power((vectCRpp[kCRppSystError]/vectCRpp[kCRppValue]), 2))*100);        

    nRaa=vectNJpsi[kJpsiMean]/(nBranchinRatio*nMB[iCent]*vectAccEff[kAccEffValue]*nTaa[iCent]*TMath::Power(10,-3)*(vectCRpp[kCRppValue]*drapidité*(maxPt-minPt)));
    
    errRaa_stat=TMath::Sqrt(TMath::Power((vectNJpsi[kJpsiStatError]/vectNJpsi[kJpsiMean]), 2)  +  TMath::Power((vectAccEff[kAccEffStatError]/vectAccEff[kAccEffValue]), 2)  +  TMath::Power((vectCRpp[kCRppStatError]/vectCRpp[kCRppValue]), 2))*nRaa;

    errRaa_syst=TMath::Sqrt(TMath::Power((vectNJpsi[kJpsiSysError]/vectNJpsi[kJpsiMean]), 2) + TMath::Power(MC_inputVSpt_stat/100, 2) + TMath::Power(MCinput_correlations/100, 2) + TMath::Power((EffLossTrack[iCent]/100),2) + TMath::Power((EffLossTrigg[iCent]/100),2) + TMath::Power((errTaa[iCent]/nTaa[iCent]), 2) + TMath::Power((systCentralityLimits[iCent]/100),2)) * nRaa;

    // printf("Uncorrelated error is %.3f %% \n", errRaa_syst/nRaa*100);
    // printf("Uncorrelated error on Njpsi is %.3f %% \n", vectNJpsi[kJpsiSysError]/vectNJpsi[kJpsiMean]*100);
    // printf("Uncorrelated error on Eff loss track is %.3f %% \n", EffLossTrack[iCent]);
    // printf("Uncorrelated error on Eff loss trig is %.3f %% \n", EffLossTrigg[iCent]);
    // printf("Uncorrelated error on Taa is %.3f %% \n", errTaa[iCent]/nTaa[iCent]*100);


    TString sName = SetRangeValue(kRaa, minY, maxY, minPt, maxPt, minCent, maxCent);
    printf("{\"%s\", {%.8f, %.8f, %.8f }},\n", sName.Data(), nRaa, errRaa_stat, errRaa_syst);
}

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
void ExportResultsRAAvsCent(Double_t minY, Double_t maxY)
{
    Double_t ptRange[4]={0, 0.3, 1, 2};
    Int_t cent[6]={0, 10, 30, 50, 70, 90};
    Int_t iCent=0;

    for(int i=0; i<3; i++)
    {
        printf("\n//%.1f-%.1f\n", ptRange[i], ptRange[i+1]);
        for(int j=0; j<5; j++)
        {
            if (cent[j]==0) iCent=0;
            else if (cent[j]==10) iCent=1;
            else if (cent[j]==30) iCent=2;
            else if (cent[j]==50) iCent=3;
            else iCent=4;
            ComputeRAAvsCent(minY, maxY, ptRange[i], ptRange[i+1], cent[j], cent[j+1]);
        }
    }     
}

void extractSigmaDiff(Int_t minCent, Int_t maxCent, Double_t minPt1, Double_t maxPt1, Double_t minPt2, Double_t maxPt2, Double_t minY, Double_t maxY)
{
    Int_t iCent=0;
    if (minCent==0) iCent=0;
    else if (minCent==10) iCent=1;
    else if (minCent==30) iCent=2;
    else if (minCent==50) iCent=3;
    else iCent=4;


    Double_t nRaa1;
    Double_t errRaa1_stat;
    Double_t errRaa1_syst;
    Double_t errRaa1_global;
    
    Double_t nRaa2;
    Double_t errRaa2_stat;
    Double_t errRaa2_syst;
    Double_t errRaa2_global;

    std::vector<Double_t> vectNJpsi1 = SignalResultsLowPt[SetRangeValue(kRawYield, minY, maxY, minPt1, maxPt1, minCent, maxCent)];
    std::vector<Double_t> vectAccEff1 = AccEffHadro[SetRangeValue(kAccEff, minY, maxY, minPt1, maxPt1, minCent, maxCent)];
    std::vector<Double_t> vectCRpp1 = ppCrossSection[SetRangeValue(kCRpp, minY, maxY, minPt1, maxPt1)];
    Double_t MCinput1_correlations = SystMCinput_correlations[SetNameSystAccEff(minY, maxY, minPt1, maxPt1)][0];
    Double_t effMTRresponse1 = SystEffMTR[SetNameSystAccEff(minY, maxY, minPt1, maxPt1)][0];


    std::vector<Double_t> vectNJpsi2 = SignalResultsLowPt[SetRangeValue(kRawYield, minY, maxY, minPt2, maxPt2, minCent, maxCent)];
    std::vector<Double_t> vectAccEff2 = AccEffHadro[SetRangeValue(kAccEff, minY, maxY, minPt2, maxPt2, minCent, maxCent)];
    std::vector<Double_t> vectCRpp2 = ppCrossSection[SetRangeValue(kCRpp, minY, maxY, minPt2, maxPt2)];
    Double_t MCinput2_correlations = SystMCinput_correlations[SetNameSystAccEff(minY, maxY, minPt2, maxPt2)][0];
    Double_t effMTRresponse2 = SystEffMTR[SetNameSystAccEff(minY, maxY, minPt2, maxPt2)][0];

    
    nRaa1=vectNJpsi1[kJpsiMean]/(nBranchinRatio*nMB[iCent]*vectAccEff1[kAccEffValue]*nTaa[iCent]*TMath::Power(10,-3)*(vectCRpp1[kCRppValue]*drapidité*(maxPt1-minPt1)));
    errRaa1_stat=TMath::Sqrt(TMath::Power((vectNJpsi1[kJpsiStatError]/vectNJpsi1[kJpsiMean]), 2)  +  TMath::Power((vectAccEff1[kAccEffStatError]/vectAccEff1[kAccEffValue]), 2)  +  TMath::Power((vectCRpp1[kCRppStatError]/vectCRpp1[kCRppValue]), 2))*nRaa1;
    errRaa1_syst=TMath::Sqrt(TMath::Power((vectNJpsi1[kJpsiSysError]/vectNJpsi1[kJpsiMean]), 2) + TMath::Power(MC_inputVSpt_stat/100, 2) + TMath::Power(MCinput1_correlations/100, 2) +  TMath::Power((EffMCH/100), 2)  + TMath::Power((effMTRresponse1/100), 2) + TMath::Power((EffMTR/100), 2)  + TMath::Power(0.01, 2) + TMath::Power((vectCRpp1[kCRppSystError]/vectCRpp1[kCRppValue]), 2)) * nRaa1;

    nRaa2=vectNJpsi2[kJpsiMean]/(nBranchinRatio*nMB[iCent]*vectAccEff2[kAccEffValue]*nTaa[iCent]*TMath::Power(10,-3)*(vectCRpp2[kCRppValue]*drapidité*(maxPt2-minPt2)));
    errRaa2_stat=TMath::Sqrt(TMath::Power((vectNJpsi2[kJpsiStatError]/vectNJpsi2[kJpsiMean]), 2)  +  TMath::Power((vectAccEff2[kAccEffStatError]/vectAccEff2[kAccEffValue]), 2)  +  TMath::Power((vectCRpp2[kCRppStatError]/vectCRpp2[kCRppValue]), 2))*nRaa2;
    errRaa2_syst=TMath::Sqrt(TMath::Power((vectNJpsi2[kJpsiSysError]/vectNJpsi2[kJpsiMean]), 2) + TMath::Power(MC_inputVSpt_stat/100, 2) + TMath::Power(MCinput2_correlations/100, 2) +  TMath::Power((EffMCH/100), 2)  + TMath::Power((effMTRresponse2/100), 2) + TMath::Power((EffMTR/100), 2)  + TMath::Power(0.01, 2) + TMath::Power((vectCRpp2[kCRppSystError]/vectCRpp2[kCRppValue]), 2)) * nRaa2;


    TString sName1 = SetRangeValue(kRaa, minY, maxY, minPt1, maxPt1, minCent, maxCent);
    printf("{\"%s\", {%.8f, %.8f, %.8f }},\n", sName1.Data(), nRaa1, errRaa1_stat, errRaa1_syst);
    TString sName2 = SetRangeValue(kRaa, minY, maxY, minPt2, maxPt2, minCent, maxCent);
    printf("{\"%s\", {%.8f, %.8f, %.8f }},\n", sName2.Data(), nRaa2, errRaa2_stat, errRaa2_syst);

    Double_t sigma;
    sigma = TMath::Sqrt(TMath::Power(errRaa1_stat, 2)+TMath::Power(errRaa1_syst, 2)+TMath::Power(errRaa2_stat, 2)+TMath::Power(errRaa2_syst, 2));
    printf("Difference is %.2f, with sigma %.2f given a deviation of %.2f\n", TMath::Abs(nRaa1-nRaa2), sigma, TMath::Abs(nRaa1-nRaa2)/sigma);
}