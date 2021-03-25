/*
 *  SystematicTest.C
 *
 *  Created by Ophelie Bugnon on 15/06/19.
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
#include "SetAndGetFitResults.C"
#include "SignalExtraction.C"

//const char* fileName = "AnalysisResults_15o_AOD229.root";
//const char* fileName = "AnalysisResults_18q_AOD225.root";
//const char* fileName = "AnalysisResults_18r_AOD225.root";
// const char* fileName = "AnalysisResults_DataMerged.root";
const char* fileName = "SubtractedHistos.root";
const char* fileLocation ="~/Documents/ALICE/AnalyseJPsi/AnalysisResults_5TeV";

Double_t min_y = -4;
Double_t max_y = -2.5;
Double_t min_pt = 0.;
Double_t max_pt = 0.3;
Int_t min_cent = 0;
Int_t max_cent = 10;
const int numberFitRanges = 2;
Double_t min_fit[numberFitRanges] = {2.2, 2.4};
Double_t max_fit[numberFitRanges] = {4.5, 4.7} ;
const int numberBackgroundFunction = 2;
Efunction test_BackGround[numberBackgroundFunction] = {kVWGQuadratic, kPol2Pol3};
const int numberSignalFunction = 2;
Efunction test_Signal[numberSignalFunction] = {kCBExtended, kNA60};
const int numberTailsTest = 4;
Etails test_tails[numberTailsTest] = {kEMB, kPP, kSTARLIGHTcoh, kSTARLIGHTincoh};
// Etails test_tails[numberTailsTest] = {kEMB,kPP};



//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void SystematicTest(const char* file=fileName, Double_t minY=min_y, Double_t maxY=max_y, Double_t minPt=min_pt, Double_t maxPt=max_pt, Int_t minCent=min_cent, Int_t maxCent=max_cent, Bool_t isSaved=kFALSE)
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
                    resultsTest = SignalExtraction(Form("%s/%s", fileLocation, fileName), minY, maxY, minPt, maxPt, minCent, maxCent,2., 5., test_BackGround[j], test_Signal[k], test_tails[l], min_fit[i], max_fit[i], kTRUE, isSaved, kV0);
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
                    if(isSaved) SetFitResults(rangeName,nameTest, fitResult, fitError);
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
    
    resultsTest = SignalExtraction(Form("%s/%s", fileLocation, fileName), minY, maxY, minPt, maxPt, minCent, maxCent,2., 5., fBackground, fSignal, pTails, minRange, maxRange, isDraw, isSaved, kV0);
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
void DrawSystematicTest(Double_t minY=min_y, Double_t maxY=max_y, Double_t minPt=min_pt, Double_t maxPt=max_pt, Int_t minCent=min_cent, Int_t maxCent=max_cent, Bool_t isSaved=kFALSE)
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
        hNumSyst[iVariable]->SetTitle(Form("%.1f-%.1f GeV/c in %i-%i%%",minPt, maxPt,minCent, maxCent));
    }
   
    int counter = 1;
    for(int j=0; j<numberBackgroundFunction ; j++)
    {
       for(int i=0; i<numberFitRanges ; i++)
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

    TString arrayYTitle[numberOfFitVariables] = {"N_{J/#psi}", "N_{#psi2S}", "#mu_{J/#psi} GeV/c", "#sigma_{J/#psi} GeV/c", "#chi^{2}/NDF"};
    Double_t arrayYRange[2*numberOfFitVariables] ={0, 0, 0, 0, 2.6, 3.4, 0.02, 0.09, 0, 2};
    for (int iVariable =0; iVariable<numberOfFitVariables; iVariable++) 
    {
        TCanvas* c = new TCanvas(Form("systematic_%s_%s", arrayFitVariableNames[iVariable].Data(), rangeName.Data()), Form("%s in range %s",arrayFitVariableNames[iVariable].Data(),rangeName.Data()), 800, 400);
        c->SetBottomMargin(0.16);
        hNumSyst[iVariable]->GetYaxis()->SetTitle(Form("%s",arrayYTitle[iVariable].Data()));
        hNumSyst[iVariable]->GetXaxis()->SetRange(1, counter+1);
        if(iVariable!=0 && iVariable!=1)
        hNumSyst[iVariable]->GetYaxis()->SetRangeUser(arrayYRange[2*iVariable], arrayYRange[2*iVariable+1]);
        hNumSyst[iVariable]->Draw();

        if(isSaved) c->SaveAs(".pdf");
    }

}

//compute the number of Jpsi and statistical, systematic uncertainties
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void GetNumberSignal(Double_t minY=min_y, Double_t maxY=max_y, Double_t minPt=min_pt, Double_t maxPt=max_pt, Int_t minCent=min_cent, Int_t maxCent= max_cent, Bool_t isWeighted=kTRUE, Bool_t isDraw=kTRUE)
{
    TString nameTest;
    TString rangeName = SetRangeName(minY, maxY, minPt, maxPt, minCent, maxCent);
    std::vector<Double_t> fitResult;
    std::vector<Double_t> fitError;
    int numberOfTests = numberFitRanges*numberBackgroundFunction*numberSignalFunction*numberTailsTest;

    TH1F* hNumSyst = new TH1F(Form("histo_JPsi_%s", rangeName.Data()),"",numberOfTests,0,numberOfTests);
    hNumSyst->SetTitle(Form("%.1f-%.1f GeV/c in %i-%i%%",minPt, maxPt,minCent, maxCent));
    hNumSyst->GetYaxis()->SetTitle("N_{J/#psi}");
    int counter =1; 

    std::vector<Double_t> result_nb_signal;
    std::vector<Double_t> error_nb_signal;
    std::vector<Double_t> weight_test;
    Double_t nb_signal = 0;
    Double_t stat_signal = 0;
    Double_t syst_signal = 0;
    Double_t sum_weight = 0;
    Double_t a=0;
    Double_t weight = 0;


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
                        // cout << "There is no results for the test " << nameTest.Data() << endl;
                        continue;
                    }
                    result_nb_signal.push_back(fitResult[0]);
                    error_nb_signal.push_back(fitError[0]);

                    if (isWeighted)
                    {
                    if(minPt == 0. && numberTailsTest == 4 && test_tails[l]==kPP) weight = 6;
                    else if(minPt == 0.3 && numberTailsTest == 4 && test_tails[l]==kPP) weight = 4;
                    else if(test_tails[l]==kPP) weight = 2;
                    else weight = 1;
                    }
                    else weight = 1;
                    weight_test.push_back(weight);

                    hNumSyst->SetBinContent(counter,fitResult[0]);
                    hNumSyst->SetBinError(counter,fitError[0]);
                    hNumSyst->GetXaxis()->SetBinLabel(counter,nameTest);
                    counter++;
                } 
            } 
        } 
    }

    Int_t nb_tests = result_nb_signal.size();
    for(int m=0; m<nb_tests; m++)
    {
        cout << result_nb_signal[m] << "pm" << error_nb_signal[m] << " with weight " << weight_test[m] << endl;
        nb_signal += weight_test[m]*result_nb_signal[m];
        stat_signal += weight_test[m]*error_nb_signal[m];
        sum_weight += weight_test[m];

        a += weight_test[m]*TMath::Power(result_nb_signal[m],2);
    }
    nb_signal = nb_signal/sum_weight;
    stat_signal = stat_signal/sum_weight;
    a = a/sum_weight;
    syst_signal = TMath::Sqrt((sum_weight/(sum_weight-1))*(a-TMath::Power(nb_signal,2)));

    nb_signal = nb_signal + 0.5 - (nb_signal<0);
    stat_signal = stat_signal + 0.5 - (stat_signal<0);
    syst_signal = syst_signal + 0.5 - (syst_signal<0);
    printf("{\"%s\", {%d., %d., %d. }},\n", rangeName.Data(), (Int_t)nb_signal, (Int_t)stat_signal, (Int_t)syst_signal);
    
    if(isDraw)
    {
        TCanvas* c = new TCanvas(Form("systematic_JPsi_%s", rangeName.Data()), Form("Number of Jpsi in range %s",rangeName.Data()), 800, 400);
        c->SetBottomMargin(0.16);
        hNumSyst->GetXaxis()->SetRange(1, counter+1);
        hNumSyst->GetYaxis()->SetRangeUser(nb_signal-20*syst_signal, nb_signal+10*syst_signal);
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
        c->SaveAs(".pdf");
    }
}


//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void ExportResults()
{
    Double_t minPtRange[16]={0, 0.3, 0.3, 0.65, 1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 10, 12 };
    Double_t maxPtRange[16]={0.3, 1, 0.65, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 12, 15};
    Double_t yRange[7]={-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5};
    // Double_t minPtRange[3]={0, 0.3, 1};
    // Double_t maxPtRange[3]={0.3, 1, 8};
    Int_t minCentRange[6]={0, 0, 10, 30, 50, 70};
    Int_t maxCentRange[6]={90, 10, 30, 50, 70, 90};

    for(int j=0; j<6; j++)
    {
        printf("//%d-%d\n",minCentRange[j],maxCentRange[j]);
        for(int i = 0; i<16; i++)
        {
            GetNumberSignal(-4, -2.5, minPtRange[i], maxPtRange[i], minCentRange[j],maxCentRange[j], kTRUE, kFALSE);
            // GetNumberSignal(yRange[i], yRange[i+1], 0.3, 15, minCentRange[j],maxCentRange[j], kTRUE, kFALSE);

        }
        printf("\n");
        printf("\n");

    }
}