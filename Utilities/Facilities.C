#include <vector>
#include <iostream>
#include <algorithm>

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