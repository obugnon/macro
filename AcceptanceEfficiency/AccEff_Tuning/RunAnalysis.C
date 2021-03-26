/*
This macro can run in two modes, either as a full analysis chain (with excutedFunction = kAllFunctions) or by excuting one single function (by provding the key of this function from the enum above). In the former case, this macro will call respectively the various functions from an iterationStart to iterationEnd steps. If earlyStop is true, the iterations will stop if the acceptance efficiency difference between steps is smaller than the provided threschold.
*/

/*
Run with (the b option is recommended if you are running with kAllFunctions):
root -l -b -q RunAnalysis.C+

*/

#include "RootIncluders.C"
#include "Common.C"
#include "Iteration/GetWeights.C"
#include "Iteration/FillHistoAndGetAccEffi.C"
#include "Iteration/FitYields.C"
#include "Iteration/FitGen.C"
#include "Iteration/DrawAndCompareAccEffi.C"
#include "Iteration/CompareIterationsShape.C"
#include "Iteration/CompareIterationsAccEffi.C"
#include "HtmlAndLatexRenderer.C"

using namespace std;

enum enumFunctions
{
  kFillHistoAndGetAccEffi,
  kFillHistoAndGetAccEffiSpecificRange,
  kDrawAndCompareAccEffi,
  kFitPtGen,
  kFitRapGen,
  kFitPtYield,
  kFitRapYield,
  kGetPtWeight,
  kGetRapWeight,
  kCompareIterationsEffi,
  kCompareIterationsShape,
  kAllFunctions
};

enum enumEarlyStopingMode
{
  kAccEffiDelta,
  kFitFuncDelta
};

void RunAnalysis(int excutedFunction = kAllFunctions, TString strCentralityBins = "0,10;10,30", int iterationStart = 0, int iterationEnd = 6, Bool_t earlyStop = kTRUE, int earlyStopMode = kAccEffiDelta, Double_t stopThreshold = 0.2)
{
  if (verbose == 0)
  {
    gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");
  }

  // Get the centrality bins, loop over them and excute the necessarly functions
  vector<vector<int>> centralityBins = StringToVector(strCentralityBins);
  int numberOfCentralityBins = (int)centralityBins.size();
  for (int iCent = 0; iCent < numberOfCentralityBins; iCent++)
  {
    PrintMessage(Form("Centrality %d-%d", centralityBins[iCent][0], centralityBins[iCent][1]), "^");
    if (iterationStart == 0)
    {
      gSystem->Exec(Form("mkdir -p Cent-%dto%d/AccEffi", centralityBins[iCent][0], centralityBins[iCent][1]));
      gSystem->Exec(Form("mkdir -p Cent-%dto%d/PtShapeIterations/iter-0", centralityBins[iCent][0], centralityBins[iCent][1]));
      gSystem->Exec(Form("mkdir -p Cent-%dto%d/RapShapeIterations/iter-0", centralityBins[iCent][0], centralityBins[iCent][1]));
      gSystem->Exec(Form("mkdir -p Cent-%dto%d/FinalPlots", centralityBins[iCent][0], centralityBins[iCent][1]));
    }
    std::vector<std::vector<TString>> accEffiVsRap;
    std::vector<std::vector<TString>> accEffiVsPt;
    //Loop over the iteration steps only if excutedFunction = kAllFunction, otherwise ignore the iterationEnd arguments
    if (excutedFunction == kAllFunctions)
    {
      int numberOfPerformedIterations = 0;
      Double_t maxDiffWRTPreviousStep_FitFunc = 100;
      Double_t maxDiffWRTPreviousStep_AccEffi = 100;
      for (int iteration = iterationStart; iteration <= iterationEnd; iteration++)
      {

        PrintMessage(Form("Iteration-step %d", iteration), "_");
        //Check if the necessary files from previous steps are available
        if (iteration > 0)
        {
          //Check the acceptance values
          if (gSystem->AccessPathName(Form("Cent-%dto%d/AccEffi/iter-%d/AccEffiValues.root", centralityBins[iCent][0], centralityBins[iCent][1], iteration - 1)))
          {
            cout << Form("AccEffi values from previous iteration step (%d) are not available. Execute root -l AccEffiIterator.C\\(%d\\)", iteration - 1, iteration - 1) << endl;
            return;
          }
        }

        if (iteration == 0)
        {
          FillHistoAndGetAccEffi(centralityBins[iCent][0], centralityBins[iCent][1], iteration);
          DrawAndCompareAccEffi(centralityBins[iCent][0], centralityBins[iCent][1], iteration);
          FitPtGen(centralityBins[iCent][0], centralityBins[iCent][1]);
          FitRapGen(centralityBins[iCent][0], centralityBins[iCent][1]);
        }
        else
        {
          FitPtYield(centralityBins[iCent][0], centralityBins[iCent][1], iteration);
          FitRapYield(centralityBins[iCent][0], centralityBins[iCent][1], iteration);
          if (earlyStop && earlyStopMode == kFitFuncDelta)
          {
            // Get the weights in rap and y
            TH1 *histoPtWeight = GetPtWeight(centralityBins[iCent][0], centralityBins[iCent][1], iteration, kFALSE);
            Double_t maxPtFitFuncDelta = TMath::Max(histoPtWeight->GetMaximum() - 1, 1 - histoPtWeight->GetMinimum());

            TH1 *histoRapWeight = GetRapWeight(centralityBins[iCent][0], centralityBins[iCent][1], iteration, kFALSE);
            Double_t maxRapFitFuncDelta = TMath::Max(histoRapWeight->GetMaximum() - 1, 1 - histoRapWeight->GetMinimum());

            maxDiffWRTPreviousStep_FitFunc = TMath::Max(maxPtFitFuncDelta, maxRapFitFuncDelta);
            PrintMessage(Form("Maximum fit-funtions diff between steps %d and %d is %2.2f %%", iteration, iteration - 1, 100 * maxDiffWRTPreviousStep_FitFunc), "-");
            if (100 * maxDiffWRTPreviousStep_FitFunc < stopThreshold)
            {
              PrintMessage(Form("Stopping at step %d", numberOfPerformedIterations), "=");
              break;
            }
          }
          FillHistoAndGetAccEffi(centralityBins[iCent][0], centralityBins[iCent][1], iteration);
          maxDiffWRTPreviousStep_AccEffi = DrawAndCompareAccEffi(centralityBins[iCent][0], centralityBins[iCent][1], iteration);
          numberOfPerformedIterations++;
          if (earlyStop && earlyStopMode == kAccEffiDelta)
          {
            PrintMessage(Form("Maximum Acc*Effi diff between steps %d and %d is %2.2f %%", iteration, iteration - 1, 100 * maxDiffWRTPreviousStep_AccEffi), "-");
            if (100 * maxDiffWRTPreviousStep_AccEffi < stopThreshold)
            {
              PrintMessage(Form("Stopping at step %d", numberOfPerformedIterations), "=");
              break;
            }
          }
        }
      } // End iteration for one centrality
      CompareIterationsShape(centralityBins[iCent][0], centralityBins[iCent][1], numberOfPerformedIterations);
      CompareIterationsAccEffi(centralityBins[iCent][0], centralityBins[iCent][1], numberOfPerformedIterations, accEffiVsRap, accEffiVsPt);
      MakeSummaryPerCentrality(centralityBins[iCent][0], centralityBins[iCent][1], numberOfPerformedIterations, 100 * maxDiffWRTPreviousStep_AccEffi, accEffiVsRap, accEffiVsPt);
    } // End condition if excutedFunction == kAllFunctions
    else
    {
      switch (excutedFunction)
      {
      case kFillHistoAndGetAccEffi:
        FillHistoAndGetAccEffi(centralityBins[iCent][0], centralityBins[iCent][1], iterationStart);
        break;
      case kFillHistoAndGetAccEffiSpecificRange:
        FillHistoAndGetAccEffi(centralityBins[iCent][0], centralityBins[iCent][1], iterationStart, kTRUE);
        break;
      case kDrawAndCompareAccEffi:
        DrawAndCompareAccEffi(centralityBins[iCent][0], centralityBins[iCent][1], iterationStart);
        break;
      case kFitPtGen:
        FitPtGen(centralityBins[iCent][0], centralityBins[iCent][1]);
        break;
      case kFitRapGen:
        FitRapGen(centralityBins[iCent][0], centralityBins[iCent][1]);
        break;
      case kFitPtYield:
        FitPtYield(centralityBins[iCent][0], centralityBins[iCent][1], iterationStart);
        break;
      case kFitRapYield:
        FitRapYield(centralityBins[iCent][0], centralityBins[iCent][1], iterationStart);
        break;
      case kGetPtWeight:
        GetPtWeight(centralityBins[iCent][0], centralityBins[iCent][1], iterationStart, kTRUE);
        break;
      case kGetRapWeight:
        GetRapWeight(centralityBins[iCent][0], centralityBins[iCent][1], iterationStart, kTRUE);
        break;
      case kCompareIterationsEffi:
        CompareIterationsAccEffi(centralityBins[iCent][0], centralityBins[iCent][1], -1, accEffiVsRap, accEffiVsPt);
        break;
      case kCompareIterationsShape:
        CompareIterationsShape(centralityBins[iCent][0], centralityBins[iCent][1], -1);
        break;
      }
    }
  } // End loop over centrality bins
  // collect html summaries and open html only if excutedfunction = all
  if (excutedFunction == kAllFunctions)
  {
    CollectCentSummaries(centralityBins);
  }
}