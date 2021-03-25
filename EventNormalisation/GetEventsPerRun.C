/*
 *  GetEventsPerRun.C
 *
 *  Created by Ophelie Bugnon on 24/01/20.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TList.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TString.h"

//const char* PFile = "AnalysisResults_15o_AOD229.root";
//const char* PFile = "AnalysisResults_18q_AOD225.root";
//const char* PFile = "AnalysisResults_18r_AOD225.root";
//const char* PFile = "AnalysisResults_16e2_merged.root";
//const char* PFile = "AnalysisResults_19a2.root";
//const char* PFile = "AnalysisResults_18c2_incohJpsi.root";
//const char* PFile = "AnalysisResults_18c2_cohJpsi.root";
//const char* PFile = "AnalysisResults_20a6_incohJpsi.root";
//const char* PFile = "AnalysisResults_20a6_cohJpsi.root";


//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
std::vector<int> runListToVector(TString runList)
{
  std::vector<int> vectorRuns;
  TString strRunList = gSystem->GetFromPipe(Form("cat ~/Documents/ALICE/AnalyseJPsi/AnalysisResults_5TeV/%s",runList.Data()));
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
//__________________________________________________________________________________________________
const bool Contains( std::vector<int>& Vec, const int& Element ) 
{
    if (std::find(Vec.begin(), Vec.end(), Element) != Vec.end()) return kTRUE;
    else  return kFALSE;
}

//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
void GetEventPerRun(const char* file, TString runList )
{   
    //Opening Root file
    TFile* analysis = TFile::Open(Form("~/Documents/ALICE/AnalyseJPsi/AnalysisResults_5TeV/%s", file));
    if (!analysis) return;
    TList* eventHistos = (TList*)analysis->Get("EventHistos_CMUL7");
    TH1I* hPSeventsPerRun = (TH1I*)eventHistos->FindObject("fHistoPSEventsPerRun");

    //Creating output file
    FILE* fichier = NULL;
    fichier = fopen("CMULtot.txt", "w");

    //Get the run list
    std::vector<int> vectorRun = runListToVector(runList);
    Int_t firstRun = vectorRun[0];
    Int_t lastRun = vectorRun[vectorRun.size()-1];
    cout << firstRun << " and " << lastRun << endl;

    Int_t nCMUL = 0.0;
    Int_t nTotCMUL = 0.0;

    //write the output file
    for (int runNumber = firstRun; runNumber < lastRun+1; runNumber++)
    {
        if (!Contains(vectorRun,runNumber)) continue;
        nCMUL = hPSeventsPerRun->GetBinContent(hPSeventsPerRun->FindBin(runNumber));
        nTotCMUL += nCMUL;
        cout << nCMUL << " CMUL triggers after PS in run " << runNumber << endl;
        fprintf(fichier, "%d %d\n", runNumber, nCMUL);
    }
    cout << nTotCMUL << " CMUL triggers after PS the period " << endl;
    fclose(fichier);
}