/*
 *  SystematicTest.C
 *
 *  Created by Ophelie Bugnon on 23/07/20.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"

#include "FitYield.C"
#include "SetAndGetFitResults.C"


const int numberFitFunction = 1;
eFunction test_Fit[numberFitFunction] = {kPowLawFree};
const int numberFitRanges = 4;
Double_t minFit[numberFitRanges] = {0, 0.3, 0.65, 1};
const int numberErrorTest = 3;

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void SystematicTest(Double_t minCent, Double_t maxCent, Double_t minPt, Double_t maxPt, Bool_t isSavedFunction=kFALSE, Bool_t isSavedValue=kFALSE, Double_t minY=-4, Double_t maxY=-2.5)
{
    TString nameTest;
    TString rangeName = SetRangeValue(kInvariantYield, minY, maxY, minPt, maxPt, minCent, maxCent);

    std::vector<Double_t> resultsTest;
    std::vector<Double_t> fitResult;
    std::vector<Double_t> fitError;
    std::vector<TString> badFit;

    Bool_t isStat, isSyst;
    for(int h=0; h<numberFitRanges; h++)
    {
        for(int i=0; i<numberFitFunction ; i++)
        {
            for(int j=0; j<numberErrorTest ; j++)
            {
                if(j==0)
                {
                    isStat=kTRUE;
                    isSyst=kFALSE;
                }
                else if(j==1)
                {
                    isStat=kFALSE;
                    isSyst=kTRUE;
                }
                else 
                {
                    isStat=kFALSE;
                    isSyst=kFALSE;
                }
                resultsTest.clear();
                fitResult.clear();
                fitError.clear();

                nameTest = SetNameTest(kInvariantYield, test_Fit[i], isStat, isSyst,  minFit[h], 15);
                cout << "------------------------------------------------------------------" << endl;
                cout << "------------------------------------------------------------------" << endl;
                cout << "Test for " << nameTest.Data() << endl;
                resultsTest = FitYield(test_Fit[i], minCent, maxCent, isStat, isSyst, minPt, maxPt, isSavedFunction, minFit[h], 15 );
                if (resultsTest.size()==0) continue;

                for (int r = 0; r < numberOfFitVariables*2; r = r+2)
                {
                    fitResult.push_back(resultsTest[r]);
                    fitError.push_back(resultsTest[r+1]);
                }
                if (resultsTest[6]>2.5 || resultsTest[8]!=0 || resultsTest[10]!=3)
                {
                    badFit.push_back(Form("%s",nameTest.Data()));
                    // continue;
                }
                if(isSavedValue) SetFitResults(rangeName,nameTest, fitResult, fitError);
            } 
        }
    }
    cout << "Number of bad fit : " << badFit.size() << endl;
    for(int m=0; m<badFit.size(); m++)
    {
        cout << "Test : " << badFit[m] << " failed " << endl;
    }    
}



//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void DrawSystematicTest(Double_t minCent, Double_t maxCent, Double_t minPt, Double_t maxPt, Double_t minY=-4, Double_t maxY=-2.5) 
{
    TString nameTest;
    TString rangeName = SetRangeValue(kInvariantYield, minY, maxY, minPt, maxPt, minCent, maxCent);

    int numberOfTests = numberFitFunction*numberErrorTest*numberFitRanges;

    TH1F* hNumSyst[numberOfFitVariables];
    for (int iVariable =0; iVariable<numberOfFitVariables; iVariable++) 
    {
        hNumSyst[iVariable] = new TH1F(Form("histo%s_%s", arrayFitVariableNames[iVariable].Data(), rangeName.Data()),"",numberOfTests,0,numberOfTests);
    }

    int counter=1;

    std::vector<Double_t> fitResult;
    std::vector<Double_t> fitError;
    Bool_t isStat, isSyst;

    for(int h=0; h<numberFitRanges; h++)
    {
        for(int i=0; i<numberFitFunction ; i++)
        {
            for(int j=0; j<numberErrorTest ; j++)
            {
                if(j==0)
                {
                        isStat=kTRUE;
                        isSyst=kFALSE;
                }
                else if(j==1)
                {
                    isStat=kFALSE;
                    isSyst=kTRUE;
                }
                else 
                {
                    isStat=kFALSE;
                    isSyst=kFALSE;
                }
                fitResult.clear();
                fitError.clear();

                nameTest = SetNameTest(kInvariantYield, test_Fit[i], isStat, isSyst,  minFit[h], 15);
                GetFitResults(rangeName,nameTest, fitResult, fitError);
                
                if (fitResult.size()==0 || fitResult[0]==0.0)
                        {
                            cout << "There is no results for the test " << nameTest.Data() << endl;
                            continue;
                        }
                cout << "Results for the test " << nameTest.Data() << endl;
                cout << fitResult[0] << " pm " << fitError[0] << endl;

                for(int iVariable=0;iVariable<numberOfFitVariables;iVariable++)
                {             
                    hNumSyst[iVariable]->SetBinContent(counter,fitResult[iVariable]);
                    hNumSyst[iVariable]->SetBinError(counter,fitError[iVariable]);
                    hNumSyst[iVariable]->GetXaxis()->SetBinLabel(counter,nameTest);
                }
                counter++;
                cout << "" << endl;
            } 
        }
    }    
    for (int iVariable =0; iVariable<numberOfFitVariables; iVariable++) 
    {
        TCanvas* c = new TCanvas(Form("systematic_%s_%s", arrayFitVariableNames[iVariable].Data(), rangeName.Data()), Form("%s in range %s",arrayFitVariableNames[iVariable].Data(),rangeName.Data()));
        hNumSyst[iVariable]->Draw();
    }

}


//compute the number of Jpsi and statistical, systematic uncertainties
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void GetIntegratedYield(Double_t minCent, Double_t maxCent, Double_t minPt, Double_t maxPt, Double_t minY=-4, Double_t maxY=-2.5)
{
    TString nameTest;
    TString rangeName = SetRangeValue(kInvariantYield, minY, maxY, minPt, maxPt, minCent, maxCent);
    std::vector<Double_t> fitResult;
    std::vector<Double_t> fitError;

    int numberOfTests = numberFitFunction*numberErrorTest*numberFitRanges;

    int counter =1; 
     
    std::vector<Double_t> result_yield;
    std::vector<Double_t> error_stat;
    std::vector<Double_t> error_syst;
    std::vector<Double_t> error_statAndSyst;

    Double_t n_yield = 0;
    Double_t errStat_yield = 0;
    Double_t errSyst_yield = 0;
    Double_t standardDeviation_yield = 0;

    Double_t sum_weight = numberFitFunction*numberFitRanges;
    Double_t a=0;


    Bool_t isStat, isSyst;
    for(int h=0; h<numberFitRanges; h++)
    {
        for(int i=0; i<numberFitFunction ; i++)
        {
            for(int j=0; j<numberErrorTest ; j++)    
            {
                fitResult.clear();
                fitError.clear();

                if(j==0)
                {
                    isStat=kTRUE;
                    isSyst=kFALSE;
                }
                else if(j==1)
                {
                isStat=kFALSE;
                isSyst=kTRUE;
                }
                else 
                {
                isStat=kFALSE;
                isSyst=kFALSE;
                }

                nameTest = SetNameTest(kInvariantYield, test_Fit[i], isStat, isSyst,  minFit[h], 15);
                GetFitResults(rangeName,nameTest, fitResult, fitError);
                if (fitResult.size()==0 || fitResult[0]==0.0)
                {
                    cout << "There is no results for the test " << nameTest.Data() << endl;
                    continue;
                }

                if(isStat && !isSyst)
                {
                    error_stat.push_back(fitError[0]);
                }

                else if(!isStat && isSyst)
                {
                    error_syst.push_back(fitError[0]);
                }

                else if(!isStat && !isSyst) 
                {
                    result_yield.push_back(fitResult[0]);
                    error_statAndSyst.push_back(fitError[0]);
                }
                counter++;
            } 
        }
    }



    Int_t nb_tests = result_yield.size();
        // cout << nb_tests << endl;

    for(int m=0; m<nb_tests; m++)
    {
        // cout << result_yield[m] << " pm " << error_statAndSyst[m] << " pm " << error_stat[m]<< " pm " << error_syst[m] << endl;
        n_yield += result_yield[m];
        errStat_yield += error_stat[m];
        errSyst_yield += error_syst[m];

        a += TMath::Power(result_yield[m],2);
    }
    n_yield = n_yield/sum_weight;
    errStat_yield = errStat_yield/sum_weight;
    errSyst_yield = errSyst_yield/sum_weight;
    standardDeviation_yield = TMath::Sqrt(sum_weight/(sum_weight - 1)*(a/sum_weight-TMath::Power(n_yield,2)));

    // cout << "Yield is " << n_yield << " pm " << errStat_yield << " pm " << errSyst_yield << " pm " << standardDeviation_yield << endl;
    Double_t fullerrSyst_yield = TMath::Sqrt(TMath::Power(errSyst_yield/n_yield, 2)+TMath::Power(standardDeviation_yield/n_yield, 2))*n_yield;

    printf("{\"%s\", {%.8f, %.8f, %.8f }},\n", rangeName.Data(), n_yield, errStat_yield, fullerrSyst_yield);

}   

