/*
Basic macro with some functions and flags to be used across diiferent macros
*/

//-------------------------------------------------------------------------//
// Some settings and global variables
int verbose = 0; // > 0 to prit results of each step
TString ptAxisTitle = "#it{p}_{T} (GeV/#it{c})";
TString yAxisTitle = "#it{y}";
TString yieldAxisTitle = "";
TString accEffiAxisTitle = "Acc#times#varepsilon";
TString plusOrMinusSymbol = "&plusmn;";
Bool_t overWriteDirs = kTRUE;
TFile *inputFileRawYield = TFile::Open("./RawYields/NJpsi_Raw.root");
TString inputMCFilesDir = "./GridAnalysis/outputFiles/";

// Fit ranges (update, this was automatically taken as the raw yield extent but added manually for now to allow exclusion of photoproduction in the fit to the yield):
Double_t ptFitRangeMin = 0.3;
Double_t ptFitRangeMax = 15;
Double_t rapFitRangeMin = -4;
Double_t rapFitRangeMax = -2.5;
//Input files for CMUL weighting:
TH1F *histoNJpsiVsCentralityGen = nullptr;
TH1I *histoNeventPerRunGen = nullptr;
TFile *fileData = TFile::Open("./DataWeights/CMUL_Centrality_Weights.root");
TH1F *histoNJpsiVsCentrality = ((TH1F *)fileData->Get("histoNJpsiVsCentrality"));
TH1I *histoCMULEventPerRun = ((TH1I *)fileData->Get("histoCMULEventPerRun"));
Double_t totalCMULEvents = histoCMULEventPerRun->Integral();
Double_t totalNumberOfJpsi = histoNJpsiVsCentrality->Integral();
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
void SetHistoStyle(TH1 *histoSomething, Int_t color, Int_t markerStyle, Float_t markerSize, TString xAxisTitle, TString yAxisTitle, Int_t rebinFactor, Bool_t isRatio)
{
  histoSomething->Rebin(rebinFactor);
  histoSomething->SetLineColor(color);
  histoSomething->SetMarkerColor(color);
  histoSomething->SetMarkerStyle(markerStyle);
  histoSomething->SetMarkerSize(markerSize);
  histoSomething->GetXaxis()->SetTitle(xAxisTitle);
  histoSomething->GetYaxis()->SetTitle(yAxisTitle);

  if (!isRatio)
  {
    histoSomething->GetXaxis()->SetLabelSize(0.04);
    histoSomething->GetXaxis()->SetTitleSize(0.04);
    histoSomething->GetXaxis()->SetTitleOffset(1);

    histoSomething->GetYaxis()->SetLabelSize(0.04);
    histoSomething->GetYaxis()->SetTitleSize(0.04);
    histoSomething->GetYaxis()->SetTitleOffset(1.2);
  }
  else
  {
    histoSomething->GetXaxis()->SetLabelSize(0.15);
    histoSomething->GetXaxis()->SetTitleSize(0.12);
    histoSomething->GetXaxis()->SetTitleOffset(1.2);

    histoSomething->GetYaxis()->SetLabelSize(0.09);
    histoSomething->GetYaxis()->SetTitleSize(0.09);
    histoSomething->GetYaxis()->SetTitleOffset(0.45);
  }
}

//-------------------------------------------------------------------------//
vector<vector<int>> StringToVector(TString strBins = "0,90;20,30", TString outerSplitter = ";", TString innerSplitter = ",")
{
  vector<vector<int>> vectorBins;
  vector<int> vectorSingleBin;

  TObjArray *objBins = strBins.Tokenize(outerSplitter);
  TIter nextBin(objBins);
  TObjString *strBin;
  while ((strBin = (TObjString *)nextBin()))
  {
    TString currentBin = strBin->GetName();
    vectorSingleBin.clear();

    TObjArray *objValues = currentBin.Tokenize(innerSplitter);
    vectorSingleBin.push_back(((TObjString *)objValues->UncheckedAt(0))->String().Atoi());
    vectorSingleBin.push_back(((TObjString *)objValues->UncheckedAt(1))->String().Atoi());
    vectorBins.push_back(vectorSingleBin);
  }
  return vectorBins;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// Terminal printing
// Progress bar
void printProgress(double percentage)
{
  int PBWIDTH = 50;
  TString PBSymbol = "|";
  TString PBSTR = "";
  for (int i = 0; i < PBWIDTH; i++)
  {
    PBSTR += PBSymbol;
  }
  int val = (int)(percentage * 100);
  int lpad = (int)(percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR.Data(), rpad, "");
  fflush(stdout);
}
// Messages
void PrintMessage(TString message = "hello worlds", TString symbol = "-")
{
  // Get the number of chars in string and add "." if even:
  int maxMessageLength = 70;
  if (message.Length() < maxMessageLength && message.Length() % 2 != 0)
  {
    message += ".";
  }
  int numberOfAddedSymbols = maxMessageLength - message.Length();
  TString symbols = "";
  for (int i = 0; i < numberOfAddedSymbols / 2; i++)
  {
    symbols += symbol;
  }
  message = symbols + message + symbols;
  cout << message << endl
       << endl;
}
// Print acc*effi values
void PrintEffi(TString title, TH1 *histoAccEffi)
{
  TString strAccEffi = title;

  Int_t numberOfBins = histoAccEffi->GetNbinsX();
  for (int iBin = 1; iBin <= numberOfBins; iBin++)
  {
    Double_t accEffi = histoAccEffi->GetBinContent(iBin);
    Double_t accEffiErr = histoAccEffi->GetBinError(iBin);

    strAccEffi.Append(Form("\n %g,%g: %4.4f +/- %4.4f \n", histoAccEffi->GetXaxis()->GetBinLowEdge(iBin), histoAccEffi->GetXaxis()->GetBinUpEdge(iBin), accEffi, accEffiErr));
  }
  cout << strAccEffi << endl;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
//For a given centrality bin, this will retrieve the latest iteration.
TString GetLastIteration(int centMin = 0, int centMax = 10, TString identifier = "AccEffi")
{
  TSystemDirectory filesDir(Form("Cent-%dto%d/%s", centMin, centMax, identifier.Data()), Form("Cent-%dto%d/%s", centMin, centMax, identifier.Data()));
  TList *listOfSubDir = filesDir.GetListOfFiles();
  TString latestIter = "iter-0";
  if (listOfSubDir)
  {
    TSystemFile *subDirInList;
    TString subDirInListName;
    TIter next(listOfSubDir);
    while ((subDirInList = (TSystemFile *)next()))
    {
      subDirInListName = subDirInList->GetName();
      if (subDirInList->IsDirectory() && subDirInListName.Contains("iter-"))
      {
        if (subDirInListName > latestIter)
          latestIter = subDirInListName;
      }
    }
  }
  return latestIter;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
Double_t pT_shape(Double_t *x, Double_t *par)
{
  // power-law
  Double_t A = par[0];
  Double_t B = par[1];
  Double_t n1 = par[2];
  Double_t n2 = par[3];
  // pt
  Double_t pt = x[0];
  return (A * pt) / pow(1 + pow(pt / B, n1), n2);
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
const Double_t *GetBinsArray(int cerntMin, int centMax, TString ptOrRap, int &numberOfBins)
{

  TH1F *histo = ((TH1F *)inputFileRawYield->Get(Form("histo%sYield_cent%dto%d", ptOrRap.Data(), cerntMin, centMax)));
  numberOfBins = histo->GetNbinsX();
  const Double_t *arrayOfBins = histo->GetXaxis()->GetXbins()->GetArray();
  return arrayOfBins;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
std::vector<double> GetAVisuallyGoodYRangeUser(Double_t min, Double_t max, Double_t zoomFactor)
{

  // a list of options:
  Double_t rangeOptions[] = {0.0005, 0.001, 0.002, 0.005, 0.008, 0.01, 0.03, 0.05, 0.08, 0.1, 0.2, 0.3, 0.5, 0.8, 1};
  int numberOfRangeOptions = sizeof(rangeOptions) / sizeof(rangeOptions[0]);
  Double_t bestOption = 0.0005;
  for (int iOption = 0; iOption < numberOfRangeOptions; iOption++)
  {
    bestOption = rangeOptions[iOption];
    if (min > (1 - zoomFactor * rangeOptions[iOption]) && max < (1 + zoomFactor * rangeOptions[iOption]))
    {
      break;
    }
  }
  std::vector<double> bestRange;

  bestRange.push_back(1 - bestOption);
  bestRange.push_back(1 + bestOption * 1.15);
  return bestRange;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// custom colors
std::map<TString, Int_t> myColorsMap{
    {"red", TColor::GetColor("#d50000")},
    {"blue", TColor::GetColor("#4285F4")},
    {"cyan", TColor::GetColor("#00bbbb")},
    {"orange", TColor::GetColor("#FF5722")},
    {"yellow", TColor::GetColor("#F4B400")},
    {"green", TColor::GetColor("#0F9D58")},
    {"purple", TColor::GetColor("#823cff")},
    {"pink", TColor::GetColor("#FFABF5")},
    {"grey", TColor::GetColor("#596E80")},
};
// array of colors for iteration purpose
const int arrayOfColors[] = {myColorsMap["red"], myColorsMap["blue"], myColorsMap["yellow"], myColorsMap["grey"], myColorsMap["green"], myColorsMap["orange"], myColorsMap["cyan"], myColorsMap["purple"], myColorsMap["pink"]};
//hex to root-int-color
Int_t GetColorIndex(TString strHexIndex = "#000099")
{
  return TColor::GetColor(strHexIndex);
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
void SetCanvasStyle(TCanvas *can)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  int font = 42;
  gROOT->SetStyle("Plain");
  gStyle->SetFrameBorderMode(0);
  TGaxis::SetMaxDigits(4);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02, "y");
  gStyle->SetLabelSize(0.05, "xyz");
  gStyle->SetLabelFont(font, "xyz");
  gStyle->SetLabelOffset(0.01, "xyz");
  gStyle->SetTitleFont(font, "xyz");
  gStyle->SetTitleOffset(1.1, "xy");
  gStyle->SetTitleSize(0.05, "xyz");
  gStyle->SetMarkerSize(1.3);
  gStyle->SetPalette(1, 0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(10);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetEndErrorSize(8);
  can->SetFillColor(0);
  can->SetBorderMode(0);
  can->SetBorderSize(0);
  can->SetLeftMargin(0.18);
  can->SetRightMargin(0.1);
  can->SetBottomMargin(0.152);
  // can->SetTopMargin(0.);
  can->SetFrameBorderMode(0);
}
//-------------------------------------------------------------------------//

TLegend* SetLegend(Int_t centMin, Int_t centMax)
{
  TLegend* leg;
  leg = new TLegend(0.21, 0.70, 0.46, 0.85);

  leg->SetMargin(0.1);
  leg->SetFillStyle(0);
  leg->SetLineColorAlpha(0, 0);
  leg->SetTextColor(kBlack);
  leg->AddEntry("NULL", Form("centrality %i-%i %%", centMin, centMax), "");
  
  return leg;
}