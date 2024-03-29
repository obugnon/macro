/*
An interface to ease getting and setting the fit results/range/test
*/
//needed root classes
#include "Riostream.h"
#include <TFile.h>
#include "TROOT.h"
#include "TH3.h"
#include "TH1.h"
#include "TList.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"

const Int_t maxNumberOfRanges = 1000;
const Int_t maxNumberOfTests = 100;
const Int_t numberOfFitVariables = 2;
TString arrayFitVariableNames[numberOfFitVariables] = {"NCrossSection", "chi2"};

//-------------------------------------------------------------------------------------------------------------------------------------//
//Return the histogram which contains the values.
TH3F *GetFitCrossSectionHisto(Bool_t createIfNotFound = kFALSE){
  TH3F *histoFitCrossSection;
  
  TFile *outputFile = new TFile("$LOWPT/macro/ResultFiles/FitResultsCrossSection.root");
  histoFitCrossSection =  ((TH3F*) outputFile->Get("histoFitCrossSection"));

  if(!histoFitCrossSection && createIfNotFound){
    cout<<"fit results file not found: creating a new histogram"<<endl;
    histoFitCrossSection =  new TH3F("histoFitCrossSection","",maxNumberOfRanges,0,maxNumberOfRanges,maxNumberOfTests,0,maxNumberOfTests,numberOfFitVariables,0,numberOfFitVariables);
  }
  return histoFitCrossSection;
}
//-------------------------------------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------------------------------------//
void SetFitResults(TString rangeName, TString testName, std::vector<double> fitResults, std::vector<double> fitResults_err ){

  TH3F *histoFitResults = GetFitCrossSectionHisto(kTRUE);

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
  int xBinInHistoFitResults = histoFitResults->GetXaxis()->FindFixBin(rangeName.Data());
  if(xBinInHistoFitResults == -1)
  {
    //This used to work with FindLastBinAbove(), no clue on why it is not working! Replace with a loop over the x-axis labels.
    Int_t nonFilledBin;
    for(nonFilledBin=0;nonFilledBin<maxNumberOfRanges;nonFilledBin++){
      TString strTempoLabel = histoFitResults->GetXaxis()->GetBinLabel(nonFilledBin+1);
      if(strTempoLabel == ""){
        break;
      }
    }
    histoFitResults->GetXaxis()->SetBinLabel(nonFilledBin+1,rangeName);
    xBinInHistoFitResults = nonFilledBin+1;
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
  int yBinInHistoFitResults = histoFitResults->GetYaxis()->FindFixBin(testName.Data());
  if(yBinInHistoFitResults == -1)
  {
    //This used to work with FindLastBinAbove(), no clue on why it is not working! Replace with a loop over the x-axis labels.
    Int_t nonFilledBin;
    for(nonFilledBin=0;nonFilledBin<maxNumberOfTests;nonFilledBin++){
      TString strTempoLabel = histoFitResults->GetYaxis()->GetBinLabel(nonFilledBin+1);
      if(strTempoLabel == ""){
        break;
      }
    }
    histoFitResults->GetYaxis()->SetBinLabel(nonFilledBin+1,testName);
    yBinInHistoFitResults = nonFilledBin+1;
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//


  for(Int_t iVariable=0;iVariable<numberOfFitVariables;iVariable++){
    histoFitResults->SetBinContent(xBinInHistoFitResults, yBinInHistoFitResults, iVariable+1, fitResults[iVariable]);
    histoFitResults->SetBinError(xBinInHistoFitResults, yBinInHistoFitResults, iVariable+1, fitResults_err[iVariable]);
    cout<<fitResults[iVariable]<<endl;
    cout<<fitResults_err[iVariable]<<endl;
  
  }

  TFile *outputFile = new TFile("$LOWPT/macro/ResultFiles/FitResultsCrossSection.root", "RECREATE");
  histoFitResults->Write();
  outputFile->Close();
}
//-------------------------------------------------------------------------------------------------------------------------------------//


//-------------------------------------------------------------------------------------------------------------------------------------//
void GetFitResults(TString rangeName, TString testName, std::vector<double> &fitResults, std::vector<double> &fitResults_err){

  TH3F *histoFitResults = GetFitCrossSectionHisto(kFALSE);
  if(!histoFitResults){
    cout<<"Fit Results file not found"<<endl;
    return;
  }
  int xBinInHistoFitResults = histoFitResults->GetXaxis()->FindFixBin(rangeName.Data());
  if(xBinInHistoFitResults == -1)
  {
    cout<<"Range not found"<<endl;
    return;
  }
  int yBinInHistoFitResults = histoFitResults->GetYaxis()->FindFixBin(testName.Data());
  if(yBinInHistoFitResults == -1)
  {
    cout<<"test not found"<<endl;
    return;
  }

  for(Int_t iVariable=0;iVariable<numberOfFitVariables;iVariable++){
    fitResults.push_back(histoFitResults->GetBinContent(xBinInHistoFitResults, yBinInHistoFitResults, iVariable+1));
    fitResults_err.push_back(histoFitResults->GetBinError(xBinInHistoFitResults, yBinInHistoFitResults, iVariable+1));
  }
}
