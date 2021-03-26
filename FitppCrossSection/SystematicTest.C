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

#include "FitCrossSection.C"
#include "SetAndGetFitResults.C"


const int numberFitFunction = 3;
// eFunction test_Fit[numberFitFunction] = {kPowLawFree, kPowLawFixed, kLevy, kUA1};
eFunction test_Fit[numberFitFunction] = {kPowLawFree, kPowLawFixed, kLevy};
const int numberErrorTest = 3;

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void SystematicTest(Double_t minPt, Double_t maxPt, Bool_t isSavedValue=kFALSE, Bool_t isSavedFunction=kFALSE, Double_t minY=-4, Double_t maxY=-2.5)
{
    TString nameTest;
    TString rangeName = SetRangeValue(kCRpp, minY, maxY, minPt, maxPt);

    std::vector<Double_t> resultsTest;
    std::vector<Double_t> fitResult;
    std::vector<Double_t> fitError;
    std::vector<TString> badFit;


    Bool_t isStat, isSyst;
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

            nameTest = SetNameTest(kCRpp, test_Fit[i], isStat, isSyst);
            cout << "------------------------------------------------------------------" << endl;
            cout << "------------------------------------------------------------------" << endl;
            cout << "Test for " << nameTest.Data() << endl;
            resultsTest = FitCrossSection(test_Fit[i], isStat, isSyst, minPt, maxPt, isSavedFunction);
            if (resultsTest.size()==0) continue; 

            if (resultsTest[2]>2.5 || resultsTest[4]!=0 || resultsTest[6]!=3)
            {
                badFit.push_back(Form("%s",nameTest.Data()));
                // continue;
            }

            for (int r = 0; r < numberOfFitVariables*2; r = r+2)
            {
                fitResult.push_back(resultsTest[r]);
                fitError.push_back(resultsTest[r+1]);
            }
            if(isSavedValue) SetFitResults(rangeName,nameTest, fitResult, fitError);
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
void DrawSystematicTest(Double_t minPt, Double_t maxPt, Double_t minY=-4, Double_t maxY=-2.5) 
{
    TString nameTest;
    TString rangeName = SetRangeValue(kCRpp, minY, maxY, minPt, maxPt);

    int numberOfTests = numberFitFunction*numberErrorTest;

    TH1F* hNumSyst[numberOfFitVariables];
    for (int iVariable =0; iVariable<numberOfFitVariables; iVariable++) 
    {
        hNumSyst[iVariable] = new TH1F(Form("histo%s_%s", arrayFitVariableNames[iVariable].Data(), rangeName.Data()),"",numberOfTests,0,numberOfTests);
    }

    int counter=1;

    std::vector<Double_t> fitResult;
    std::vector<Double_t> fitError;
    Bool_t isStat, isSyst;
    
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

            nameTest = SetNameTest(kCRpp, test_Fit[i], isStat, isSyst);
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
    for (int iVariable =0; iVariable<numberOfFitVariables; iVariable++) 
    {
        TCanvas* c = new TCanvas(Form("systematic_%s_%s", arrayFitVariableNames[iVariable].Data(), rangeName.Data()), Form("%s in range %s",arrayFitVariableNames[iVariable].Data(),rangeName.Data()));
        hNumSyst[iVariable]->Draw();
    }

}


//compute the number of Jpsi and statistical, systematic uncertainties
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void GetIntegratedCrossSection(Double_t minPt, Double_t maxPt, Bool_t isWeighted=kFALSE, Double_t minY=-4, Double_t maxY=-2.5)
{
    TString nameTest;
    TString rangeName = SetRangeValue(kCRpp, minY, maxY, minPt, maxPt);
    std::vector<Double_t> fitResult;
    std::vector<Double_t> fitError;

    int numberOfTests = numberFitFunction*numberErrorTest;;

    int counter =1; 
     
    std::vector<Double_t> result_cross_section;
    std::vector<Double_t> error_stat;
    std::vector<Double_t> error_syst;
    std::vector<Double_t> error_statAndSyst;

    Double_t n_cr = 0;
    Double_t stat_cr = 0;
    Double_t syst_cr = 0;
    Double_t standardDeviation_cr = 0;

    
    Double_t sum_weight = 0;
    Double_t a=0;
    Double_t weight = 0;


    Bool_t isStat, isSyst;
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

            nameTest = SetNameTest(kCRpp, test_Fit[i], isStat, isSyst );
            GetFitResults(rangeName,nameTest, fitResult, fitError);
            if (fitResult.size()==0 || fitResult[0]==0.0)
            {
                cout << "There is no results for the test " << nameTest.Data() << endl;
                continue;
            }

            if (isWeighted)
            {
                if(test_Fit[i]==kPowLawFree || test_Fit[i]==kPowLawFixed) weight = 1;
                else weight = 2;
            }
            else weight = 1;
            sum_weight += (weight/numberErrorTest);
            cout << sum_weight << endl;

            if(isStat && !isSyst)
            {
                error_stat.push_back(fitError[0]);
                error_stat.push_back(weight);
            }

            else if(!isStat && isSyst)
            {
                error_syst.push_back(fitError[0]);
                error_syst.push_back(weight);
            }

            else if(!isStat && !isSyst) 
            {
                result_cross_section.push_back(fitResult[0]);
                result_cross_section.push_back(weight);

                error_statAndSyst.push_back(fitError[0]);
                error_statAndSyst.push_back(weight);
            }
            counter++;
        } 
    }

    Int_t nb_tests = result_cross_section.size();
    cout << nb_tests << endl;
    for(int m=0; m<nb_tests; m+=2)
    {
        cout << result_cross_section[m] << " pm " << error_statAndSyst[m] << " pm " << error_stat[m]<< " pm " << error_syst[m] << " with weight " << result_cross_section[m+1] << endl;
        n_cr += result_cross_section[m+1]*result_cross_section[m];
        stat_cr += error_stat[m+1]*error_stat[m];
        syst_cr += error_syst[m+1]*error_syst[m];

        a += result_cross_section[m+1]*TMath::Power(result_cross_section[m],2);
    }
    n_cr = n_cr/sum_weight;
    stat_cr = stat_cr/sum_weight;
    syst_cr = syst_cr/sum_weight;
    // standardDeviation_cr = TMath::Sqrt(a/sum_weight-TMath::Power(n_cr,2));
    a = a/sum_weight;
    standardDeviation_cr = TMath::Sqrt(sum_weight/(sum_weight - 1)*(a-TMath::Power(n_cr,2)));
    

    cout << "Cross section is " << n_cr << " pm " << stat_cr << " pm " << syst_cr << " pm " << standardDeviation_cr << endl;
    Double_t fullSyst_cr = TMath::Sqrt(TMath::Power(syst_cr/n_cr, 2)+TMath::Power(standardDeviation_cr/n_cr, 2))*n_cr;

    printf("{\"%s\", {%.5f, %.5f, %.5f }},\n", rangeName.Data(), n_cr, stat_cr, fullSyst_cr);

}   