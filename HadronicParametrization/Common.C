/*
 *  Common.C
 *
 *  Created by Ophelie Bugnon on 10/11/20.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>


#include "TFile.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/HadronicParametrization/FitFunctions.C>

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TDirectory* GetInputFunctionsDirectory(eVariable kVar, TString sRange, TString sTest)
{
    TString sFile;
    switch (kVar)
    {
    case kCRpp:
        sFile = "/Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/FitppCrossSection/FitFunctionsCRpp.root";
    break;
    
    case kAccEff:
        sFile = "/Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/FitAcceptanceEfficiency/FitFunctionsAcceptanceEfficiency.root";
    break;

    case kRaa:
        sFile = "/Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/FitRaa/FitFunctionsRaa.root";
    break;
    }
    
    
    TFile *inputFile = TFile::Open(Form("%s", sFile.Data()));
    if(!inputFile) printf("Could not find input file !\n");

    TDirectory *inputList;
    if(kVar==kCRpp) inputList = (TDirectory*) inputFile->Get(Form("%s", sTest.Data()));
    else  inputList = (TDirectory*) inputFile->Get(Form("%s_%s", sRange.Data(), sTest.Data()));

    if(!inputList) printf("Could not find input directory for %s range %s and test %s!\n", sFile.Data(), sRange.Data(), sTest.Data());
    return inputList;
    inputFile->Close();
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
std::vector<std::vector<Double_t>> GetFunctionParameters(eVariable kVar, TString sRange, TString sTest)
{
    TDirectory* inputDirectory = GetInputFunctionsDirectory(kVar, sRange, sTest);

    std::vector<std::vector<Double_t>> *vect;
    inputDirectory->GetObject("parameters", vect);
    std::vector<std::vector<Double_t>> parameters = *vect;
    return parameters;
}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
Int_t GetNPar(TString sRange, TString sTestRaa, TString sTestCRpp, TString sTestAccEff)
{
    std::vector<std::vector<Double_t>> vectRaa=GetFunctionParameters(kRaa, sRange, sTestRaa);
    std::vector<std::vector<Double_t>> vectCRpp=GetFunctionParameters(kCRpp, sRange, sTestCRpp);
    std::vector<std::vector<Double_t>> vectAccEff=GetFunctionParameters(kAccEff, sRange, sTestAccEff);

    Int_t nParFull = vectRaa[0].size()+vectCRpp[0].size()+vectAccEff[0].size()-3;
    return nParFull;
}    
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void GetFullFunctionParameters(TString sRange, TString sTestRaa, TString sTestCRpp, TString sTestAccEff, Double_t par[], Double_t parError[])
{
    std::vector<std::vector<Double_t>> vectRaa=GetFunctionParameters(kRaa, sRange, sTestRaa);
    std::vector<std::vector<Double_t>> vectCRpp=GetFunctionParameters(kCRpp, sRange, sTestCRpp);
    std::vector<std::vector<Double_t>> vectAccEff=GetFunctionParameters(kAccEff, sRange, sTestAccEff);

    for(int i=0; i<vectRaa[0].size()-1; i++)
    {
        par[i]=vectRaa[0][i];
        parError[i]=vectRaa[1][i];

    }
    for(int i=0; i<vectCRpp[0].size()-1; i++)
    {
        par[i+vectRaa[0].size()-1]=vectCRpp[0][i];
        parError[i+vectRaa[0].size()-1]=vectCRpp[1][i];    
    }
    for(int i=0; i<vectAccEff[0].size()-1; i++)
    {
        par[i+vectRaa[0].size()+vectCRpp[0].size()-2]=vectAccEff[0][i];
        parError[i+vectRaa[0].size()+vectCRpp[0].size()-2]=vectAccEff[1][i];
    }
    
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TMatrixDSym GetCovarianceMatrix(eVariable kVar, TString sRange, TString sTest)
{
    TDirectory* inputDirectory = GetInputFunctionsDirectory(kVar, sRange, sTest);
    TMatrixDSym *covarianceMatrix;
    inputDirectory->GetObject("covarianceMatrix", covarianceMatrix);
    TMatrixDSym covMatrix=*covarianceMatrix;
    return covMatrix;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TMatrixDSym GetFullCovarianceMatrix(TString sRange, TString sTestRaa, TString sTestCRpp, TString sTestAccEff)
{
    //input matrix
    TMatrixDSym covMatrixRaa = GetCovarianceMatrix(kRaa, sRange, sTestRaa);
    TMatrixDSym covMatrixCRpp = GetCovarianceMatrix(kCRpp, sRange, sTestCRpp);
    TMatrixDSym covMatrixAccEff = GetCovarianceMatrix(kAccEff, sRange, sTestAccEff);

    // covMatrixRaa.Print();
    // covMatrixCRpp.Print();
    // covMatrixAccEff.Print();

    const Int_t nRowRaa = covMatrixRaa.GetNrows();
    const Int_t nRowCRpp = covMatrixCRpp.GetNrows();
    const Int_t nRowAccEff = covMatrixAccEff.GetNrows();

    const Int_t nRowFull = nRowRaa+nRowCRpp+nRowAccEff;

    TMatrixDSym covMatrixFull(nRowFull);
    TMatrixDSub(covMatrixFull, 0, nRowRaa-1, 0, nRowRaa-1) = covMatrixRaa;
    TMatrixDSub(covMatrixFull, nRowRaa, nRowRaa+nRowCRpp-1, nRowRaa, nRowRaa+nRowCRpp-1) = covMatrixCRpp;
    TMatrixDSub(covMatrixFull, nRowRaa+nRowCRpp, nRowRaa+nRowCRpp+nRowAccEff-1, nRowRaa+nRowCRpp, nRowRaa+nRowCRpp+nRowAccEff-1) = covMatrixAccEff;

    return covMatrixFull;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TMatrixDSym GetExtendedCovarianceMatrix(eVariable kVar, TString sRange, TString sTestRaa, TString sTestCRpp, TString sTestAccEff)
{
    //input matrix
    TMatrixDSym covMatrixRaa = GetCovarianceMatrix(kRaa, sRange, sTestRaa);
    TMatrixDSym covMatrixCRpp = GetCovarianceMatrix(kCRpp, sRange, sTestCRpp);
    TMatrixDSym covMatrixAccEff = GetCovarianceMatrix(kAccEff, sRange, sTestAccEff);

    // covMatrixRaa.Print();
    // covMatrixCRpp.Print();
    // covMatrixAccEff.Print();

    const Int_t nRowRaa = covMatrixRaa.GetNrows();
    const Int_t nRowCRpp = covMatrixCRpp.GetNrows();
    const Int_t nRowAccEff = covMatrixAccEff.GetNrows();

    const Int_t nRowFull = nRowRaa+nRowCRpp+nRowAccEff;

    TMatrixDSym covMatrixFull(nRowFull);
    if(kVar==kRaa) TMatrixDSub(covMatrixFull, 0, nRowRaa-1, 0, nRowRaa-1) = covMatrixRaa;
    else if(kVar==kCRpp) TMatrixDSub(covMatrixFull, nRowRaa, nRowRaa+nRowCRpp-1, nRowRaa, nRowRaa+nRowCRpp-1) = covMatrixCRpp;
    else if(kVar==kAccEff) TMatrixDSub(covMatrixFull, nRowRaa+nRowCRpp, nRowRaa+nRowCRpp+nRowAccEff-1, nRowRaa+nRowCRpp, nRowRaa+nRowCRpp+nRowAccEff-1) = covMatrixAccEff;

    // covMatrixFull.Print();
    return covMatrixFull;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TF1* GetFullFunction(TString sRange, TString sTestRaa, TString sTestCRpp, TString sTestAccEff)
{
    TF1* fRatioIntegral;

    Int_t nParFull = GetNPar(sRange, sTestRaa, sTestCRpp, sTestAccEff);
    Double_t parameters[nParFull];
    Double_t parameterErrors[nParFull];
    GetFullFunctionParameters(sRange, sTestRaa, sTestCRpp, sTestAccEff, parameters, parameterErrors);

    

    //def functions 
    if(sTestRaa.Contains("WoodSaxonLike") && sTestCRpp.Contains("Levy") && sTestAccEff.Contains("RatioLevy"))
    {
        fRatioIntegral = new TF1("FunctionRatioIntegral", Ratio_WoodSaxon_Levy_RatioLevy, xMinParam, xMaxParam, nParFull);
    }

    if(sTestRaa.Contains("WoodSaxonLike") && sTestCRpp.Contains("Levy") && sTestAccEff.Contains("Pol3"))
    {
        fRatioIntegral = new TF1("FunctionRatioIntegral", Ratio_WoodSaxon_Levy_Polynomial, xMinParam, xMaxParam, nParFull);
    }

    if(sTestRaa.Contains("WoodSaxonLike") && sTestCRpp.Contains("PowerLaw") && sTestAccEff.Contains("RatioLevy"))
    {
        fRatioIntegral = new TF1("FunctionRatioIntegral", Ratio_WoodSaxon_PowerLaw_RatioLevy, xMinParam, xMaxParam, nParFull);
    }

    if(sTestRaa.Contains("WoodSaxonLike") && sTestCRpp.Contains("PowerLaw") && sTestAccEff.Contains("Pol3"))
    {
        fRatioIntegral = new TF1("FunctionRatioIntegral", Ratio_WoodSaxon_PowerLaw_Polynomial, xMinParam, xMaxParam, nParFull);
    }

    if((sTestRaa.Contains("Linear") || sTestRaa.Contains("Const")) && sTestCRpp.Contains("Levy") && sTestAccEff.Contains("RatioLevy"))
    {
        fRatioIntegral = new TF1("FunctionRatioIntegral", Ratio_Linear_Levy_RatioLevy, xMinParam, xMaxParam, nParFull);
    }

    if((sTestRaa.Contains("Linear") || sTestRaa.Contains("Const")) && sTestCRpp.Contains("Levy") && sTestAccEff.Contains("Pol3"))
    {
        fRatioIntegral = new TF1("FunctionRatioIntegral", Ratio_Linear_Levy_Polynomial, xMinParam, xMaxParam, nParFull);
    }

    if((sTestRaa.Contains("Linear") || sTestRaa.Contains("Const")) && sTestCRpp.Contains("PowerLaw") && sTestAccEff.Contains("RatioLevy"))
    {
        fRatioIntegral = new TF1("FunctionRatioIntegral", Ratio_Linear_PowerLaw_RatioLevy, xMinParam, xMaxParam, nParFull);
    }

    if((sTestRaa.Contains("Linear") || sTestRaa.Contains("Const")) && sTestCRpp.Contains("PowerLaw") && sTestAccEff.Contains("Pol3"))
    {
        fRatioIntegral = new TF1("FunctionRatioIntegral", Ratio_Linear_PowerLaw_Polynomial, xMinParam, xMaxParam, nParFull);
    }

    fRatioIntegral->SetParameters(parameters);
    fRatioIntegral->SetParErrors(parameterErrors);

    return fRatioIntegral;
}
