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
#include "SetAndGetFitResults.C"
#include "SignalExtraction.C"

//const char* fileName = "AnalysisResults_15o_AOD229.root";
//const char* fileName = "AnalysisResults_18q_AOD225.root";
//const char* fileName = "AnalysisResults_18r_AOD225.root";
const char* fileName = "AnalysisResults_DataMerged.root";
const char* fileLocation ="~/Documents/ALICE/AnalyseJPsi/AnalysisResults_5TeV";

Double_t min_y = -4.;
Double_t max_y = -2.5;
Double_t min_pt = 0.;
Double_t max_pt = 0.3;
Int_t min_cent = 0;
Int_t max_cent = 90;
const int numberFitRanges = 2;
Double_t min_fit[numberFitRanges] = {2.2, 2.4};
Double_t max_fit[numberFitRanges] = {4.5, 4.7} ;
const int numberBackgroundFunction = 2;
Efunction test_BackGround[numberBackgroundFunction] = {kVWGQuadratic, kPol2Pol3};
const int numberSignalFunction = 2;
Efunction test_Signal[numberSignalFunction] = {kCBExtended, kNA60};
const int numberTailsTest = 4;
Etails test_tails[numberTailsTest] = {kEMB, kPP, kSTARLIGHTcoh, kSTARLIGHTincoh};



//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void SystematicTest(const char* file=fileName, Double_t minY=min_y, Double_t maxY=max_y, Double_t minPt=min_pt, Double_t maxPt=max_pt, Int_t minCent=min_cent, Int_t maxCent=max_cent)
{
    TString nameTest;
    TString rangeName = SetRangeName(minY, maxY, minPt, maxPt, minCent, maxCent);

    std::vector<Double_t> resultsTest;
    std::vector<Double_t> fitResult;
    std::vector<Double_t> fitError;
    std::vector<TString> badFit;

    for(int i=0; i<numberFitRanges ; i++)
    {
       for(int j=0; j<numberBackgroundFunction ; j++)
        {
            for(int k=0; k<numberSignalFunction ; k++)
            {
                for(int l=0; l<numberTailsTest ; l++)
                {
                    resultsTest.clear();
                    fitResult.clear();
                    fitError.clear();
                    nameTest = SetNameTest(test_BackGround[j], test_Signal[k], test_tails[l], min_fit[i], max_fit[i]);
                    cout << "------------------------------------------------------------------" << endl;
                    cout << "------------------------------------------------------------------" << endl;
                    cout << "Test for " << nameTest.Data() << endl;
                    resultsTest = SignalExtraction(Form("%s/%s", fileLocation, fileName), minY, maxY, minPt, maxPt, minCent, maxCent,2., 5., test_BackGround[j], test_Signal[k], test_tails[l], min_fit[i], max_fit[i], kFALSE, kFALSE);
                    if (resultsTest.size()==0) continue; 

                    for (int r = 0; r < numberOfFitVariables*2; r = r+2)
                    {
                        fitResult.push_back(resultsTest[r]);
                        fitError.push_back(resultsTest[r+1]);
                    }
                    if (resultsTest[10] != 0 || resultsTest[11] != 3 || resultsTest[8]>2.5) 
                    {
                        badFit.push_back(Form("%s",nameTest.Data()));
                        continue;
                    }
                    SetFitResults(rangeName,nameTest, fitResult, fitError);
                } 
            } 
        } 
    }
    cout << "Number of bad fit : " << badFit.size() << endl;
    for(int m=0; m<badFit.size(); m++)
    {
        cout << "Test : " << badFit[m] << " failed " << endl;
    }

}
//Set manually a test result
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void SetTestResults( Double_t minRange, Double_t maxRange, Efunction fBackground, Efunction fSignal, Etails pTails, Bool_t isSaved, Bool_t isDraw, const char* file=fileName, Double_t minY=min_y, Double_t maxY=max_y, Double_t minPt=min_pt, Double_t maxPt=max_pt, Int_t minCent=min_cent, Int_t maxCent=max_cent)
{
    TString nameTest;
    TString rangeName = SetRangeName(minY, maxY, minPt, maxPt, minCent, maxCent);

    std::vector<Double_t> resultsTest;
    std::vector<Double_t> fitResult;
    std::vector<Double_t> fitError;
    std::vector<TString> badFit;

    resultsTest.clear();
    fitResult.clear();
    fitError.clear();
    nameTest = SetNameTest(fBackground, fSignal, pTails, minRange, maxRange);
    
    resultsTest = SignalExtraction(Form("%s/%s", fileLocation, fileName), minY, maxY, minPt, maxPt, minCent, maxCent,2., 5., fBackground, fSignal, pTails, minRange, maxRange, isDraw, kFALSE);
    if (resultsTest.size()==0) return; 

    for (int r = 0; r < numberOfFitVariables*2; r = r+2)
    {
        fitResult.push_back(resultsTest[r]);
        fitError.push_back(resultsTest[r+1]);
    }
    if (resultsTest[10] != 0 || resultsTest[11] != 3 || resultsTest[8]>2.5) 
    {
        cout << "Test failed with sit status " <<resultsTest[10]<< ", chi2 = "<< resultsTest[8] <<" and cov matrix status " << resultsTest[11]<< endl;
    }
    if(isSaved) SetFitResults(rangeName,nameTest, fitResult, fitError);

}

//Draw results for all variables as function of the test
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void DrawSystematicTest(Double_t minY=min_y, Double_t maxY=max_y, Double_t minPt=min_pt, Double_t maxPt=max_pt, Int_t minCent=min_cent, Int_t maxCent=max_cent)
{
    TString nameTest;
    TString rangeName = SetRangeName(minY, maxY, minPt, maxPt, minCent, maxCent);
    std::vector<Double_t> fitResult;
    std::vector<Double_t> fitError;
    int numberOfTests = numberFitRanges*numberBackgroundFunction*numberSignalFunction*numberTailsTest;

    TH1F* hNumSyst[numberOfFitVariables];
    for (int iVariable =0; iVariable<numberOfFitVariables; iVariable++) 
    {
        hNumSyst[iVariable] = new TH1F(Form("histo%s_%s", arrayFitVariableNames[iVariable].Data(), rangeName.Data()),"",numberOfTests,0,numberOfTests);
        hNumSyst[iVariable]->SetTitle(Form("%i-%i%",minCent, maxCent));
    }
   
    int counter = 1;
    for(int i=0; i<numberFitRanges ; i++)
    {
       for(int j=0; j<numberBackgroundFunction ; j++)
        {
            for(int k=0; k<numberSignalFunction ; k++)
            {
                for(int l=0; l<numberTailsTest ; l++)
                {
                    fitResult.clear();
                    fitError.clear();

                    nameTest = SetNameTest(test_BackGround[j], test_Signal[k], test_tails[l], min_fit[i], max_fit[i]);
                    GetFitResults(rangeName,nameTest, fitResult, fitError);
                    if (fitResult.size()==0 || fitResult[0]==0.0)
                    {
                        cout << "There is no results for the test " << nameTest.Data() << endl;
                        //counter++;
                        continue;
                    }
                    cout << "Results for the test " << nameTest.Data() << endl;

                    for(int iVariable=0;iVariable<numberOfFitVariables;iVariable++)
                    {             
                        hNumSyst[iVariable]->SetBinContent(counter,fitResult[iVariable]);
                        hNumSyst[iVariable]->SetBinError(counter,fitError[iVariable]);
                        hNumSyst[iVariable]->GetXaxis()->SetBinLabel(counter,nameTest);
                        cout << "result = " << fitResult[iVariable] << endl;
                        
                    }
                    counter++;
                    cout << "" << endl;
                    cout << "" << endl;
                } 
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
void GetNumberSignal(Double_t minY=min_y, Double_t maxY=max_y, Double_t minPt=min_pt, Double_t maxPt=max_pt, Int_t minCent=min_cent, Int_t maxCent= max_cent)
{
    TString nameTest;
    TString rangeName = SetRangeName(minY, maxY, minPt, maxPt, minCent, maxCent);
    std::vector<Double_t> fitResult;
    std::vector<Double_t> fitError;
    int numberOfTests = numberFitRanges*numberBackgroundFunction*numberSignalFunction*numberTailsTest;

    TH1F* hNumSyst = new TH1F(Form("histo_JPsi_%s", rangeName.Data()),"",numberOfTests,0,numberOfTests);
    hNumSyst->SetTitle(Form("p_{T} in [%.1f,%.1f] in %i-%i",minPt, maxPt,minCent, maxCent));
    int counter =1; 

    std::vector<Double_t> result_nb_signal;
    std::vector<Double_t> error_nb_signal;
    Double_t nb_signal = 0;
    Double_t stat_signal = 0;
    Double_t syst_signal = 0;
    Double_t a=0;


    for(int i=0; i<numberFitRanges ; i++)
    {
       for(int j=0; j<numberBackgroundFunction ; j++)
        {
            for(int k=0; k<numberSignalFunction ; k++)
            {
                for(int l=0; l<numberTailsTest ; l++)
                {
                    fitResult.clear();
                    fitError.clear();

                    nameTest = SetNameTest(test_BackGround[j], test_Signal[k], test_tails[l], min_fit[i], max_fit[i]);
                    GetFitResults(rangeName,nameTest, fitResult, fitError);
                    if (fitResult.size()==0 || fitResult[0]==0.0)
                    {
                        cout << "There is no results for the test " << nameTest.Data() << endl;
                        continue;
                    }
                    result_nb_signal.push_back(fitResult[0]);
                    error_nb_signal.push_back(fitError[0]);

                    hNumSyst->SetBinContent(counter,fitResult[0]);
                    hNumSyst->SetBinError(counter,fitError[0]);
                    hNumSyst->GetXaxis()->SetBinLabel(counter,nameTest);
                    hNumSyst->GetXaxis()->SetSizeLabel(2.5);
                    counter++;
                } 
            } 
        } 
    }
    Int_t nb_tests = result_nb_signal.size();
    for(int m=0; m<nb_tests; m++)
    {
        cout << result_nb_signal[m] << "pm" << error_nb_signal[m] << endl;
        nb_signal += result_nb_signal[m];
        stat_signal += error_nb_signal[m];
        a += TMath::Power(result_nb_signal[m],2);
    }
    nb_signal = nb_signal/nb_tests;
    stat_signal = stat_signal/nb_tests;
    syst_signal = TMath::Sqrt(a/nb_tests-TMath::Power(nb_signal,2));
    cout << "Number of Jpsi is " << nb_signal << " pm " << stat_signal << " pm " << syst_signal << endl;
    

        TCanvas* c = new TCanvas(Form("systematic_JPsi_%s", rangeName.Data()), Form("Number of Jpsi in range %s",rangeName.Data()));
        hNumSyst->Draw();
        TF1 *mean = new TF1("mean","[0]",-1,numberOfTests+1);
        mean->SetParameter(0, nb_signal);
        mean->SetLineColor(46);
        mean->Draw("SAME");
        TF1 *plussyst = new TF1("plussyst","[0]+[1]",-1,numberOfTests+1);
        plussyst->SetParameter(0, nb_signal);
        plussyst->SetParameter(1, syst_signal);
        plussyst->SetLineColor(9);
        plussyst->Draw("SAME");
        TF1 *minussyst = new TF1("minussyst","[0]-[1]",-1,numberOfTests+1);
        minussyst->SetParameter(0, nb_signal);
        minussyst->SetParameter(1, syst_signal);
        minussyst->SetLineColor(9);
        minussyst->Draw("SAME");

        TLatex* text = new TLatex();
        text->SetNDC();
        text->SetTextAlign(12);
        text->SetTextFont(43);
        text->SetTextSize(15);
        text->DrawLatex(0.4, 0.4, Form("N_{J/#psi} = %i #pm %i (stat) #pm %i (syst)", (Int_t)nb_signal, (Int_t)stat_signal, (Int_t)syst_signal));

}
