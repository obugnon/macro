#include <vector>
#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TObjArray.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"


//__________________________________________________________________________________________________
std::vector<int> runListToVector(TString runList)
{
  std::vector<int> vectorRuns;
  TString strRunList = gSystem->GetFromPipe(Form("cat %s",runList.Data()));
  TObjArray *objRunList = strRunList.Tokenize("\n");
  TIter nextRun(objRunList);
  TObjString *strRun;
  while ((strRun=(TObjString*)nextRun())) 
  {
    TString currentRun = strRun->GetName();
    int runNumber = currentRun.Atoi();
    vectorRuns.push_back(runNumber);
  }
  return vectorRuns;
}

//__________________________________________________________________________________________________
const bool Contains( std::vector<int>& Vec, const int& Element ) 
{
    if (std::find(Vec.begin(), Vec.end(), Element) != Vec.end()) return kTRUE;
    else  return kFALSE;
}

//__________________________________________________________________________________________________
void Events_CMUL(const char* path)
{
    TFile* analysis = TFile::Open(Form("%s/AnalysisResults.root",path));
    if (!analysis) return;
    TList* eventHistos = (TList*)analysis->Get("EventHistos_CMUL7");
    TH1* hEvents = (TH1*)eventHistos->FindObject("fHistoPSEventsPerRun");
    TH1* hEventsbPS = (TH1*)eventHistos->FindObject("fHistoEventsBeforePSPerRun");

    std::vector<int> vectorRun = runListToVector(Form("%s/runList.txt",path));
    Int_t firstRun = vectorRun[0];
    Int_t lastRun = vectorRun[vectorRun.size()-1];
    cout << firstRun << " and " << lastRun << endl;

    Double_t nCMUL = 0.0;
    Double_t nTotCMUL = 0.0;
    Double_t nCMULbPS = 0.0;
    Double_t nTotbPS = 0.0;
    for (int runNumber = firstRun; runNumber < lastRun+1; runNumber++)
    {
        if (!Contains(vectorRun,runNumber)) continue;
        nCMUL = hEvents->GetBinContent(hEvents->FindBin(runNumber));
        nTotCMUL += nCMUL;

        nCMULbPS = hEventsbPS->GetBinContent(hEventsbPS->FindBin(runNumber));
        nTotbPS += nCMULbPS;
        
        cout << nCMUL << " CMUL triggers after PS in run " << runNumber << endl;
    }
    cout << nTotbPS << " CMUL triggers before PS the period " << endl;
    cout << nTotCMUL << " CMUL triggers after PS the period " << endl;

}
