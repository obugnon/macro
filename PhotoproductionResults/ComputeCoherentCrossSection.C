#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>


#include "TFile.h"
#include "TMath.h"

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsSignalExtraction.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsAcceptanceEfficiency.C>


enum kAnalysisType{kWithYield, kParam}; 
//Values from UPC paper
// Double_t fI=0.055;
// Double_t fD=0.055;
// Double_t errfD=0.01; //absolute uncertainty
// Double_t errfI=0.001; //absolute uncertainty
// Double_t errfD_stat=0.; 
// Double_t errfI_stat=0.;

//Values from with correct pT cut
Double_t fI=0.0888;
Double_t fD=0.0664;
Double_t errfD=0.0126; //systematic final value, absolute uncertainty
Double_t errfI=0.0338; //systematic final value, absolute uncertainty
// Double_t errfD=0.0114; //systematic from difference with the paper, absolute uncertainty
// Double_t errfI=0.0338; //systematic from difference with the paper, absolute uncertainty
// Double_t errfD_stat=0.0126; //absolute uncertainty
// Double_t errfI_stat=0.0027; //absolute uncertainty

//Values from paper on excess in peripheral collisions at 2.76TeV
// Double_t fI=0.14;
// Double_t fD=0.10;
// Double_t errfD=0.06; //absolute uncertainty 
// Double_t errfI=0.16; //absolute uncertainty 
// Double_t errfD_stat=0.; 
// Double_t errfI_stat=0.;

Double_t nBranchinRatio = 0.05961;
Double_t errBranchinRatio = 0.5; //globale (%)
Double_t drapidité = 1.5;

Double_t nMB = 5230386400;
Double_t errMB = 13711845;
Double_t inelCrossSection = 7.67; //(en b)
Double_t errinelCrossSection = 0.3835; //(en b) //changed with public note uncertainty
Double_t nLint = 756.3; //(inverse microbarn)
Double_t errLint = 18.9; //16.6 manuscrit et 18.9 papier absolute uncertainty (inverse microbarn)

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void ComputeNexcess(kAnalysisType aType=kParam, Int_t minCent=0, Int_t maxCent=10, Double_t minPt=0, Double_t maxPt=0.3, Double_t minY=-4, Double_t maxY=-2.5)
{
    Double_t nExcess;
    Double_t errExcess_stat;
    Double_t errExcess_syst;
    Double_t sigma;
    Double_t nSigma;

    std::vector<Double_t> vectHadro;
    std::vector<Double_t> vectTotal = SignalResultsLowPt[SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent)];

    switch (aType)
    {
        case kWithYield:
        vectHadro = HadroResultsWithYield[SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent)];
        break;
    
        case kParam: 
        vectHadro = HadroResultsWithParam[SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent)];
        break;
    }

    nExcess=vectTotal[kJpsiMean]-vectHadro[kJpsiMean];
    errExcess_stat=TMath::Sqrt(vectTotal[kJpsiStatError]*vectTotal[kJpsiStatError]+vectHadro[kJpsiStatError]*vectHadro[kJpsiStatError]);
    errExcess_syst=TMath::Sqrt(vectTotal[kJpsiSysError]*vectTotal[kJpsiSysError]+vectHadro[kJpsiSysError]*vectHadro[kJpsiSysError]);

    sigma = TMath::Sqrt(TMath::Power(errExcess_stat/nExcess, 2) + TMath::Power(errExcess_syst/nExcess, 2))*nExcess;
    nSigma = nExcess/sigma;
    
    TString sName = SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent);
    printf("{\"%s\", {%.f, %.f, %.f }}, // %.2fsigma\n", sName.Data(), nExcess, errExcess_stat, errExcess_syst, nSigma);
    printf("#sigma = %.2f giving a significance of %.2f\n", sigma, nSigma); 
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void ComputeCoherentCrossSection(Int_t minCent, Int_t maxCent, Double_t minPt=0, Double_t maxPt=0.3, Double_t minY=-4, Double_t maxY=-2.5)
{
    Int_t iCent=0;
    if (minCent==0) iCent=0;
    else if (minCent==10) iCent=1;
    else if (minCent==30) iCent=2;
    else if (minCent==50) iCent=3;
    else if (minCent==70) iCent=4;
    
    Double_t nExcess;
    Double_t errExcess_stat;
    Double_t errExcess_syst;
    Double_t nCoh;
    Double_t errCoh_stat;
    Double_t errCoh_syst;
    Double_t nCohCR;
    Double_t errCohCR_stat;
    Double_t errCohCR_syst;
    Double_t errCohCR_global;

    std::vector<Double_t> vectTotal = SignalResultsLowPt[SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent)];
    std::vector<Double_t> vectHadro = HadroResultsWithParam[SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent)];
    std::vector<Double_t> vectAccEff = AccEffCoherent[SetRangeValue(kAccEff, minY, maxY, minPt, maxPt, minCent, maxCent)];

    Double_t effMTRresponse = SystEffMTR_Starlight;

    //Number of J/psi in excess
    nExcess=vectTotal[kJpsiMean]-vectHadro[kJpsiMean];
    errExcess_stat=TMath::Sqrt(vectTotal[kJpsiStatError]*vectTotal[kJpsiStatError]+vectHadro[kJpsiStatError]*vectHadro[kJpsiStatError]);
    errExcess_syst=TMath::Sqrt(vectTotal[kJpsiSysError]*vectTotal[kJpsiSysError]+vectHadro[kJpsiSysError]*vectHadro[kJpsiSysError]);

    //Number of J/psi from coherent photoproduction
    nCoh = nExcess/(1+fI+fD);
    errCoh_stat = TMath::Sqrt(TMath::Power(errExcess_stat/nExcess, 2))*nCoh;
    errCoh_syst = TMath::Sqrt(TMath::Power(errExcess_syst/nExcess, 2) + TMath::Power(errfD/(1+fI+fD),2) + TMath::Power(errfI/(1+fI+fD),2))*nCoh;

    TString sName = SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent);
    // printf("J/psi from coherent photoproduction is %.1f \\pm %.1f \\ %.1f\n", nCoh, errCoh_stat, errCoh_syst);
    // printf("Err stat with fi and fd is %.3f and without is %.3f\n", errCoh_stat, errExcess_stat/nExcess*nCoh);
    // printf("Err stat on Ncoh is %.3f%% and err syst on Ncoh is %.3f%%\n", errCoh_stat/nCoh*100, errCoh_syst/nCoh*100);
    printf("{\"%s\", {%.1f, %.1f, %.1f }},\n\n", sName.Data(), nCoh, errCoh_stat, errCoh_syst);

    //Cross section for J/psi coherent photoproduction
    errCohCR_global = TMath::Power((errBranchinRatio/100), 2) + TMath::Power(errfD/(1+fI+fD),2) + TMath::Power(errfI/(1+fI+fD),2) + TMath::Power((EffMCH/100), 2)  + TMath::Power((effMTRresponse/100), 2) + TMath::Power((EffMTR/100), 2) + TMath::Power(0.01, 2) + TMath::Power(MC_inputStarlight/100, 2) + TMath::Power(errLint/nLint, 2) + TMath::Power((MC_pTcut_Starlight/100), 2);
    errCohCR_global = TMath::Sqrt(errCohCR_global)*100;

    printf("Global error is %.3f %% \n", errCohCR_global);

    nCohCR=nCoh/(nBranchinRatio*vectAccEff[kAccEffValue]*nLint*drapidité);

    errCohCR_stat=TMath::Sqrt(TMath::Power(errCoh_stat/nCoh, 2)  +  TMath::Power(vectAccEff[kAccEffStatError]/vectAccEff[kAccEffValue], 2))*nCohCR;

    errCohCR_syst=TMath::Sqrt(TMath::Power(errExcess_syst/nExcess, 2) + TMath::Power((EffLossTrack[iCent]/100),2) + TMath::Power((EffLossTrigg[iCent]/100),2) + TMath::Power(systCentralityLimits[iCent]/100, 2)) * nCohCR;

    printf("Coherent photoproduction cross section per unit of rapidity = %.3f \\pm %.3f (%.3f%%) \\pm %.3f (%.3f%%) \\pm %.3f (%.3f%%)\n", nCohCR, errCohCR_stat , errCohCR_stat/nCohCR*100, errCohCR_syst, errCohCR_syst/nCohCR*100, errCohCR_global/100*nCohCR, errCohCR_global);
    printf("{\"%s\", {%.3f, %.3f, %.3f , %.3f}},\n\n", sName.Data(), nCohCR, errCohCR_stat, errCohCR_syst, errCohCR_global/100*nCohCR);

}


//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void ComputeCoherentYield(Int_t minCent, Int_t maxCent, Double_t minPt=0, Double_t maxPt=0.3, Double_t minY=-4, Double_t maxY=-2.5)
{
    Int_t iCent=0;
    if (minCent==0) iCent=0;
    else if (minCent==10) iCent=1;
    else if (minCent==30) iCent=2;
    else if (minCent==50) iCent=3;
    else if (minCent==70) iCent=4;
    
    Double_t nExcess;
    Double_t errExcess_stat;
    Double_t errExcess_syst;
    Double_t nCoh;
    Double_t errCoh_stat;
    Double_t errCoh_syst;
    Double_t nCohYield;
    Double_t errCohYield_stat;
    Double_t errCohYield_syst;
    Double_t errCohYield_global;

    Double_t nEvMB, errEvMB;
    nEvMB=nMB*20/90;
    errEvMB=errMB*20/90;

    std::vector<Double_t> vectTotal = SignalResultsLowPt[SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent)];
    std::vector<Double_t> vectHadro = HadroResultsWithParam[SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent)];
    std::vector<Double_t> vectAccEff = AccEffCoherent[SetRangeValue(kAccEff, minY, maxY, minPt, maxPt, minCent, maxCent)];

    Double_t effMTRresponse = SystEffMTR_Starlight;

    //Number of J/psi in excess
    nExcess=vectTotal[kJpsiMean]-vectHadro[kJpsiMean];
    errExcess_stat=TMath::Sqrt(vectTotal[kJpsiStatError]*vectTotal[kJpsiStatError]+vectHadro[kJpsiStatError]*vectHadro[kJpsiStatError]);
    errExcess_syst=TMath::Sqrt(vectTotal[kJpsiSysError]*vectTotal[kJpsiSysError]+vectHadro[kJpsiSysError]*vectHadro[kJpsiSysError]);
    
    //Number of J/psi from coherent photoproduction
    nCoh = nExcess/(1+fI+fD);
    errCoh_stat = TMath::Sqrt(TMath::Power(errExcess_stat/nExcess, 2))*nCoh;
    errCoh_syst = TMath::Sqrt(TMath::Power(errExcess_syst/nExcess, 2) + TMath::Power(errfD/(1+fI+fD),2) + TMath::Power(errfI/(1+fI+fD),2))*nCoh;

    //Yield for J/psi coherent photoproduction
    errCohYield_global = TMath::Power((errBranchinRatio/100), 2) + TMath::Power(errfD/(1+fI+fD),2) + TMath::Power(errfI/(1+fI+fD),2) + TMath::Power((EffMCH/100), 2)  + TMath::Power((effMTRresponse/100), 2) + TMath::Power((EffMTR/100), 2) + TMath::Power(0.01, 2) + TMath::Power(MC_inputStarlight/100, 2) + TMath::Power(errEvMB/nEvMB, 2) + TMath::Power((MC_pTcut_Starlight/100), 2);
    errCohYield_global = TMath::Sqrt(errCohYield_global)*100;

    printf("Global error is %.3f %% \n", errCohYield_global);

    nCohYield=nCoh/(nBranchinRatio*vectAccEff[kAccEffValue]*nEvMB*drapidité);

    errCohYield_stat=TMath::Sqrt(TMath::Power(errCoh_stat/nCoh, 2)  +  TMath::Power(vectAccEff[kAccEffStatError]/vectAccEff[kAccEffValue], 2))*nCohYield;

    errCohYield_syst=TMath::Sqrt(TMath::Power(errExcess_syst/nExcess, 2) + TMath::Power((EffLossTrack[iCent]/100),2) + TMath::Power((EffLossTrigg[iCent]/100),2) + TMath::Power(systCentralityLimits[iCent]/100, 2)) * nCohYield;

    TString sName = SetRangeValue(kRawYield, minY, maxY, minPt, maxPt, minCent, maxCent);
    // printf("{\"%s\", {%.1f, %.1f, %.1f }},\n\n", sName.Data(), nCoh, errCoh_stat, errCoh_syst);
    // printf("Coherent photoproduction cross section per unit of rapidity = %.3f \\pm %.3f (%.3f%%) \\pm %.3f (%.3f%%) \\pm %.3f (%.3f%%)\n", nCohYield, errCohYield_stat , errCohYield_stat/nCohYield*100, errCohYield_syst, errCohYield_syst/nCohYield*100, errCohYield_global/100*nCohYield, errCohYield_global);
    printf("{\"%s\", {%.10f, %.10f, %.10f , %.10f}},\n\n", sName.Data(), nCohYield, errCohYield_stat, errCohYield_syst, errCohYield_global/100*nCohYield);

}