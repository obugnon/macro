/*
 *  RunAnalysis.C
 *
 *  Created by Ophelie Bugnon on 12/11/20.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TFile.h"
#include "TMatrixD.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TF1.h"

#include "Common.C"
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsSignalExtraction.C>

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
std::vector<Double_t> RunSingleAnalysis(eFunction fRaa=kPol1, eFunction fCRpp=kPowLawFree , eFunction fAccEff=kRatioLevy, Double_t minFitRaa=0.65, Double_t fixedParamRaa=0, Int_t minCent=70, Int_t maxCent=90)
{
    // std::cout << "------------------------ New Test ------------------------" << std::endl;

    TString sRange=SetRangeFunction(-4, -2.5, minCent, maxCent);
    TString sTestRaa=SetNameTest(kRaa, fRaa, kFALSE, kFALSE, minFitRaa, 15, fixedParamRaa);
    TString sTestCRpp=SetNameTest(kCRpp, fCRpp, kFALSE, kFALSE);
    TString sTestAccEff=SetNameTest(kAccEff, fAccEff, kTRUE, kFALSE, 0, 10);

    printf("\\item test \\Raa\\ with %s, test pp reference cross section with %s and test \\AccEff\\ with %s.\n", sTestRaa.Data(), sTestCRpp.Data(), sTestAccEff.Data());

    Int_t nPar=GetNPar(sRange, sTestRaa, sTestCRpp, sTestAccEff);
    // std::cout << "Number of parameters is " << nPar << std::endl;

    Double_t parameters[nPar];
    Double_t parameterErrors[nPar];

    GetFullFunctionParameters(sRange, sTestRaa, sTestCRpp, sTestAccEff, parameters, parameterErrors);
    
    TF1* function = GetFullFunction(sRange, sTestRaa, sTestCRpp, sTestAccEff);
    function->SetParameters(parameters);
    function->SetParErrors(parameterErrors);

    Double_t derPar[nPar];
    Double_t x[1] = {0.15};
    function->GradientPar(x, derPar);
    TMatrixD jac(1, nPar);

    for (int ipar = 0; ipar < nPar; ++ipar)
    {
        // std::cout << "Parameter " << ipar << std::endl;
        // std::cout << function->GetParameter(ipar) << std::endl;
        // std::cout << "Derivative: " << derPar[ipar] << std::endl;
        jac[0][ipar] = derPar[ipar];
    }
    // std::cout << "Jacobian Matrix" << std::endl;
    // jac.Print();

    TMatrixDSym covMatrix = GetFullCovarianceMatrix(sRange, sTestRaa, sTestCRpp, sTestAccEff);
    // std::cout << "Full covariance Matrix" << std::endl;
    // covMatrix.Print();

    TMatrixD tmp(jac, TMatrixD::kMult, covMatrix);
    TMatrixD Err2(tmp, TMatrixD::kMultTranspose, jac);
    // Err2.Print("V");

    //Results on ration R(pT)
    Double_t valueRatio = function->Eval(0.15);
    Double_t errRatio = TMath::Sqrt(Err2[0][0]);
    Double_t errRelatRation = errRatio/valueRatio*100;

    std::cout << "Value of the ratio is " << valueRatio << " with uncertainty " << errRatio  <<"(" << errRelatRation << "%)"<< std::endl;

    //Compute each contribution to sigma R
    TMatrixDSym covMatrixRaa = GetExtendedCovarianceMatrix(kRaa, sRange, sTestRaa, sTestCRpp, sTestAccEff);
    TMatrixD tmpRaa(jac, TMatrixD::kMult, covMatrixRaa);
    TMatrixD Err2Raa(tmpRaa, TMatrixD::kMultTranspose, jac);
    Double_t errRatioRaa = TMath::Sqrt(Err2Raa[0][0]);
    Double_t errRelatRationRaa = errRatioRaa/valueRatio*100;
    // covMatrixRaa.Print();
    std::cout << "Contribution from Raa is " << errRatioRaa << "    (" << errRelatRationRaa << "%)" << std::endl;


    TMatrixDSym covMatrixCRpp = GetExtendedCovarianceMatrix(kCRpp, sRange, sTestRaa, sTestCRpp, sTestAccEff);
    TMatrixD tmpCRpp(jac, TMatrixD::kMult, covMatrixCRpp);
    TMatrixD Err2CRpp(tmpCRpp, TMatrixD::kMultTranspose, jac);
    Double_t errRatioCRpp = TMath::Sqrt(Err2CRpp[0][0]);
    Double_t errRelatRationCRpp = errRatioCRpp/valueRatio*100;
    // covMatrixCRpp.Print();
    std::cout << "Contribution from CR pp is " << errRatioCRpp << "    (" << errRelatRationCRpp << "%)" << std::endl;


    TMatrixDSym covMatrixAccEff = GetExtendedCovarianceMatrix(kAccEff, sRange, sTestRaa, sTestCRpp, sTestAccEff);
    TMatrixD tmpAccEff(jac, TMatrixD::kMult, covMatrixAccEff);
    TMatrixD Err2AccEff(tmpAccEff, TMatrixD::kMultTranspose, jac);
    Double_t errRatioAccEff = TMath::Sqrt(Err2AccEff[0][0]);
    Double_t errRelatRationAccEff = errRatioAccEff/valueRatio*100;
    // covMatrixAccEff.Print();
    std::cout << "Contribution from Acc Eff is " << errRatioAccEff << "    (" << errRelatRationAccEff << "%)" << std::endl;


    //GetNumber of Jpsi in 1-8 GeV/c
    std::vector<Double_t> vectNJpsi = SignalResultsLowPt[SetRangeValue(kRawYield, -4, -2.5, 1, 8, minCent, maxCent)];
    Double_t NJpsi = vectNJpsi[kJpsiMean];
    Double_t errJpsi_stat = vectNJpsi[kJpsiStatError];
    Double_t errJpsi_syst = vectNJpsi[kJpsiSysError];

    Double_t nJpsiHadro = NJpsi*valueRatio;
    Double_t errSystTot = TMath::Sqrt(TMath::Power(errRatio/valueRatio, 2)+TMath::Power(errJpsi_syst/NJpsi, 2))*nJpsiHadro;
    Double_t errStatTot = errJpsi_stat/NJpsi*nJpsiHadro;

    //ExportResults 
    std::vector<Double_t> vResults;
    vResults.push_back(NJpsi);
    vResults.push_back(errJpsi_stat);
    vResults.push_back(errJpsi_syst);
    vResults.push_back(valueRatio);
    vResults.push_back(errRatio);
    
    vResults.push_back(nJpsiHadro);
    vResults.push_back(errStatTot);
    vResults.push_back(errSystTot);
    // std::cout << vResults[5]+ 0.5 - (vResults[5]<0) << std::endl;
    // printf("%.f\n", 1.5);

    return vResults;
} 


void RunAnalysis(Int_t minCent=10, Int_t maxCent=30)
{
    //General
    const Int_t nFunctionCRpp = 3;
    eFunction tfuncCRpp[nFunctionCRpp] = {kPowLawFree, kPowLawFixed, kLevy};

    std::vector<eFunction> tfuncAccEff;
    // tfuncAccEff = {kRatioLevy, kPol3};
    tfuncAccEff = {kRatioLevy};
    const Int_t nFunctionAccEff = tfuncAccEff.size();


    const Int_t nMinFitRaa = 2;
    Double_t minFitRaa[nMinFitRaa] = {0.65, 1};


    std::vector<eFunction> tfuncRaa;
    std::vector<Double_t> vectPt0;
    if(minCent>30) 
    {
        tfuncRaa = {kConst, kPol1};
        vectPt0 = {0, 0};
    }
    else 
    {
        tfuncRaa = {kWoodSaxonFree, kWoodSaxonFixed, kWoodSaxonFixed, kWoodSaxonFixed, kWoodSaxonFixed};
        if(minCent==0) vectPt0 = {0, 3.097, 1.549, 2.10, 4.20}; //0-10
        else if(minCent==10) vectPt0 = {0, 3.097, 1.549, 2.19, 4.38}; //10-30
        else vectPt0 = {0, 3.097, 1.549, 2.34, 4.68}; //30-50
    }
    const Int_t nFunctionRaa = tfuncRaa.size();

    Int_t nTestTot=nFunctionCRpp*nFunctionAccEff*nFunctionRaa*nMinFitRaa;

    //Output variables
    std::vector<Double_t> vectResults;
    Int_t testNumber=1;

    Double_t n_jpsi_allTest = 0;
    Double_t stat_jpsi_allTest = 0;
    Double_t syst_jpsi_allTest = 0;
    Double_t R_allTest = 0;
    Double_t syst_R_allTest = 0;
    Double_t standardDeviation_jpsi_allTest = 0;
    Double_t squareMean = 0;


    for(int i=0; i<nFunctionCRpp; i++)
    {
        for (int j=0; j<nFunctionAccEff; j++)
        {
            for(int k=0; k<nFunctionRaa; k++)
            {
                for(int l=0; l<nMinFitRaa; l++)
                {
                    vectResults.clear();
                    if(k==2 || (k==3 && l==1) || (k==4 && l==1) || (k==5 && l==1)) continue; //remove bad contribution from Raa parametrizations
                    vectResults = RunSingleAnalysis(tfuncRaa[k], tfuncCRpp[i] , tfuncAccEff[j], minFitRaa[l], vectPt0[k], minCent, maxCent);
                    // 0 = Njpsi from signal extaction
                    // 1 = ErrStat from signal extaction
                    // 2 = ErrSyst from signal extaction
                    // 3 = Ratio of Integral
                    // 4 = Error on Ratio
                    // 5 = nJpsi hadronic
                    // 6 = err stat on nJpsi hadronic
                    // 7 = err syst on nJpsi hadronic
                    printf("\\hline\n");
                    printf("%i & %.2f & %.2f & %.2f & $ %.f \\pm %.f \\pm %.f $ \\\\\n",testNumber, vectResults[4]/vectResults[3]*100, vectResults[2]/vectResults[0]*100, vectResults[7]/vectResults[5]*100, vectResults[5], vectResults[6], vectResults[7]);
                    
                    testNumber++;

                    n_jpsi_allTest += vectResults[5];
                    stat_jpsi_allTest += vectResults[6];
                    syst_jpsi_allTest += vectResults[7];
                    squareMean += TMath::Power(vectResults[5],2);
                    R_allTest += vectResults[3];
                    syst_R_allTest += vectResults[4];
                }
            }
        } 
    }
    nTestTot=testNumber-1;
    n_jpsi_allTest /= nTestTot;
    stat_jpsi_allTest /= nTestTot;
    syst_jpsi_allTest /= nTestTot;
    R_allTest /= nTestTot;
    syst_R_allTest /= nTestTot;

    squareMean = squareMean/nTestTot;
    standardDeviation_jpsi_allTest = TMath::Sqrt(nTestTot/(nTestTot - 1)*(squareMean-TMath::Power(n_jpsi_allTest,2)));

    Double_t fullSyst_jpsi = TMath::Sqrt(TMath::Power(syst_jpsi_allTest/n_jpsi_allTest, 2)+TMath::Power(standardDeviation_jpsi_allTest/n_jpsi_allTest, 2))*n_jpsi_allTest;
    Double_t fullSyst_jpsi_test = TMath::Sqrt(TMath::Power(vectResults[2]/vectResults[0], 2)+TMath::Power(syst_R_allTest/R_allTest, 2)+TMath::Power(standardDeviation_jpsi_allTest/n_jpsi_allTest, 2))*n_jpsi_allTest;
    // printf("\nContribution from Njpsi is %.3f, from R is %.3f and from standart deviation %.3f for total of %.f\n", vectResults[2]/vectResults[0]*100, syst_R_allTest/R_allTest*100, standardDeviation_jpsi_allTest/n_jpsi_allTest*100, fullSyst_jpsi_test);

    printf("{\"%s\", {%.1f, %.1f, %.1f }},\n", SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, minCent, maxCent).Data() , n_jpsi_allTest, stat_jpsi_allTest, fullSyst_jpsi);

}