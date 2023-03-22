This is a collection of macros that can be used for J/psi MC tuning in different centrality bins in Pb-Pb at 5.023 TeV. There are two main sets of macros, the ones used for Acc*Effi wehigting via an iteration procedure (in ./Iteration) and the ones used for filling the uncorrected raw Jpsi yield (in ./RawYields). 


-------------------------------------------<./RawYields>--------------------------------------------
-This directory has only one macro that must be ran only once before starting the iteration pricedure.
-Its output is root file containing the two histograms for the raw yield distributions (vs pt and y) in different centrality bins. 



-------------------------------------------<./Iteration>--------------------------------------------
-This directory contains the main macros needed to perform an iteration procedure. The whole chain (detailed below), can be run using the RunAnalysis.C.
++++++++++++++++++++++++++++++++++++++RunAnalysis.C.++++++++++++++++++++++++++++++++++++++
+This macro is the one to run for performing the iteration and in turn it will excute the needed function from the macros below. It can run in two modes, either as a full analysis chain (with excutedFunction = kAllFunctions) or by excuting one single function (by provding the key of this function from the enum above). In the former case, this macro will call respectively the various functions from an iterationStart to iterationEnd steps. If earlyStop is true, 
+arguments:

-int excutedFunction: expected values: kAllFunctions: "to run all the functions in series", or one of the keys provided below to run the various functions

-TString strCentralityBins: A string containing the centrality bins for which the analysis will be performed (Make sure that their raw yields is filled in ./RawYields "TODO: Make that macro a part of the analysis chain"). Example: "0,10;10,30" will run on two centrality bins 0,10 and 10,30

-int iterationStart: usually 0 but can be set otherwise if the 0th step has been already performed.
******The following three parameters are only needed if excutedFunction == kAllFunctions******
-int iterationEnd: maximum number of of iterations perfromed for a given centrality.
-Bool_t earlyStop: If this is true (recommended), the iteration procedure will stop when the difference between two successive steps in terms of (acc*effi or fit-functions) is smaller than the provided threschold.
-int earlyStopMode (kAccEffiDelta or kFitFuncDelta)
-Double_t stopThreshold, in percentage (e.g 0.2 will stop the iteration if the difference is less than 0.2%)

+Running example:
root -l -b -q RunAnalysis.C+\(kAllFunctions,\"0,10;10,30\",0,3,kTRUE,kAccEffiDelta,0.2\)

-Macros in the directory ./Iteration:  FillHistoAndGetAccEffi.C, DrawAndCompareAccEffi.C, FitGen.C,FitYields.C, GetWeights.C, CompareIterationsShape.C, CompareIterationsAccEffi.C
+++++++++++++++++++++++++++++++++FillHistoAndGetAccEffi.C+++++++++++++++++++++++++++++++++
+standalone run option: kFillHistoAndGetAccEffi
+This macro takes as inputs a tree that contains dimuon information (generated and reconstructed) and fill them in histograms. The Acc*Effi is calculated as well.
This macros has been modified to weight the pt and rapidity shape via an iterative procedure. It also need the centrality bin as an argument.

+++++++++++++++++++++++++++++++++DrawAndCompareAccEffi.C.+++++++++++++++++++++++++++++++++
+standalone run option: kDrawAndCompareAccEffi
+This macro compare the pt and rapidity dependence of the Acc*Eff that correposnds to two different steps A and B.  

+++++++++++++++++++++++++++++++++++++++++FitGen.C+++++++++++++++++++++++++++++++++++++++++
+standalone run option: kFitPtGen | kFitRapGen
+For the zeroth iteration the fit is performed on the generated jpsi in a given centrality which must result in parameters similar to the one used for the MC generation (with no centrality dependence: THIS HAS TO BE CHECKED)
TODO: Incorporate this into a more general FitYields.C that can fit either corrected yield or generated distributions.

+++++++++++++++++++++++++++++++++++++++FitYields.C.+++++++++++++++++++++++++++++++++++++++
+standalone run option: kFitPtYield | kFitRapYield
+For a given iteration step, this macro calculate the yield of J/psi vs pt or y using the Acc*Eff values obtained from the previous step. It then fits the distribution and store the fit parameters to be used in the next step. The pt fit function is a power law plus an exponential for high pt. A zero-centered gaussian is used for the rapidity fit.

+++++++++++++++++++++++++++++++++++++++GetWeights.C+++++++++++++++++++++++++++++++++++++++
+standalone run option: kGetPtWeight | kGetRapWeight
+This macro return the weight as a function of y for a given iteration step. The weight is given by the ratio of the normalised y fit from the current step and the one before (for the step-0 take the generated input shape).
The macro is basically the one sent by Laure (weight.C).

+++++++++++++++++++++++++++++++++CompareIterationsShape.C+++++++++++++++++++++++++++++++++
+standalone run option: kCompareIterationsShape
+This macro compare the pt and rapidity shapes obtained for various iteration and centrality bins. It is not an intermediate step but rather used after the last iteration to produce a summary plot.

++++++++++++++++++++++++++++++++CompareIterationsAccEffi.C++++++++++++++++++++++++++++++++
+standalone run option: kCompareIterationsEffi
+This macro compare the pt and rapidity Acc$effi obtained for various iteration and centrality bins. It is not an intermediate step but rather used after the last iteration to produce a summary plot


##Needed inputs:
-input files for the MC. The macro FillHistoAndGetAccEffi.C suports TChain reading, just make sure to properly assign the variable "inputMCFilesDir".


##Output after each iteration step n > 0 and multiplicity bin x-y
-Cent-xtoy/AccEffi/iter-n/AccEffiValues.root: root file containing Acc*Eff histograms vs pt and y
-Cent-xtoy/AccEffi/iter-n/AccEffVsXX_ComparisonToIter-n-1.pdf: pdf files comparing the results to the previous step
-Cent-xtoy/RapShapeIterations/iter-n/RapYield.pdf: Rapidity yeild fitted
-Cent-xtoy/RapShapeIterations/iter-n/values.txt: rapidty fit results
-Cent-xtoy/RapShapeIterations/iter-n/RapWeight.pdf: rapidty weight (current shape divided by the previous)
-Same for pt


-------------------------------------------<other macros>----------------------------------------
-Common.C: Basic macro with some functions and flags to be used across diiferent macros.
-HTMLAndLatexRenderer.C: A collection of helper methods to make an HTML summary. (TODO: Implement similar functions for latex)

