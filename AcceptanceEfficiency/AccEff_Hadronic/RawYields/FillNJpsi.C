class RawDataVsCentrality
{
public:
  RawDataVsCentrality(double centMin, double centMax, std::vector<double> rapBins, std::vector<double> rapValues, std::vector<double> rapErrors, std::vector<double> ptBins, std::vector<double> ptValues, std::vector<double> ptErrors)
  {
    fCentMin = centMin;
    fCentMax = centMax;
    fRapBins = rapBins;
    fRapValues = rapValues;
    fRapErrors = rapErrors;
    fPtBins = ptBins;
    fPtValues = ptValues;
    fPtErrors = ptErrors;
  }
  RawDataVsCentrality()
  {
  }
  double fCentMin;
  double fCentMax;
  std::vector<double> fRapBins;
  std::vector<double> fPtBins;
  std::vector<double> fRapValues;
  std::vector<double> fRapErrors;
  std::vector<double> fPtValues;
  std::vector<double> fPtErrors;

  bool VerifyData()
  {
    if (fCentMax <= fCentMin)
    {
      cout << "Cent-min and cent-max seems to be switched" << endl;
      return false;
    }
    if (fRapBins.size() != fRapValues.size() + 1 || fRapBins.size() != fRapErrors.size() + 1 || fRapBins.size() == 0)
    {
      return false;
    }
    if (fPtBins.size() != fPtValues.size() + 1 || fPtBins.size() != fPtErrors.size() + 1 || fPtBins.size() == 0)
    {
      return false;
    }
    return true;
  }

  TH1F *MakeRapHisto()
  {
    TH1F *histoRapYield = new TH1F(Form("histoRapYield_cent%gto%g", fCentMin, fCentMax), "", fRapBins.size() - 1, &fRapBins[0]);
    for (int iRap = 0; iRap < fRapBins.size() - 1; iRap++)
    {
      Double_t rapidityBinWidth = (fRapBins[iRap + 1] - fRapBins[iRap]);
      histoRapYield->SetBinContent(iRap + 1, (fRapValues[iRap]) / rapidityBinWidth);
      Double_t error = fRapErrors[iRap]; // fRapValues[iRap];
      histoRapYield->SetBinError(iRap + 1, error / rapidityBinWidth);
    }
    return histoRapYield;
  }
  TH1F *MakePtHisto()
  {
    TH1F *histoPtYield = new TH1F(Form("histoPtYield_cent%gto%g", fCentMin, fCentMax), "", fPtBins.size() - 1, &fPtBins[0]);
    for (int iPt = 0; iPt < fPtBins.size() - 1; iPt++)
    {
      Double_t ptidityBinWidth = (fPtBins[iPt + 1] - fPtBins[iPt]);
      histoPtYield->SetBinContent(iPt + 1, (fPtValues[iPt]) / ptidityBinWidth);
      Double_t error = fPtErrors[iPt] ;// fPtValues[iPt];
      histoPtYield->SetBinError(iPt + 1, error / ptidityBinWidth);
    }
    return histoPtYield;
  }
};
void FillNJpsi()
{

  int const numberOfCentralityBins = 5;
  RawDataVsCentrality dataVsCentrality[numberOfCentralityBins];
  dataVsCentrality[0] = RawDataVsCentrality(0, 10, {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5}, {15186, 54056, 86443, 99776, 82142, 24417}, {499, 1138, 1605, 1756, 1675, 966}, {0, 0.3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15}, {8519, 68751, 115505, 82025, 44022, 22053, 11791, 5850, 3200, 1827, 1067, 1004, 542}, {771, 1854, 1971, 1430, 842, 573, 345, 204, 146, 97, 78, 69, 77});

  dataVsCentrality[1] = RawDataVsCentrality(10, 30, {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5}, {17145, 55092, 91821, 103394, 80290, 26871}, {472, 906, 1319, 1538, 1452, 895}, {0, 0.3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15}, {9840, 66379, 112970, 83586, 48732, 26806, 14570, 8152, 4325, 2443, 1324, 1504, 836}, {583, 1657, 1629, 1208, 725, 491, 308, 197, 136, 96, 93, 72, 78});

  dataVsCentrality[2] = RawDataVsCentrality(30, 50, {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5}, {6429, 20287, 31962, 35935, 28150, 8194}, {200, 387, 559, 633, 566, 326}, {0, 0.3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15}, {4372, 20591, 35643, 28588, 18378, 11340, 6608, 3855, 2197, 1223, 732, 738, 395}, {230, 644, 693, 482, 321, 223, 149, 91, 80, 60, 62, 70, 44});

  dataVsCentrality[3] = RawDataVsCentrality(50, 70, {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5}, {2074, 5936, 9226, 10323, 8054, 2475}, {93, 151, 198, 218, 194, 112}, {0, 0.3, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15}, {2820, 5580, 10008, 7863, 5628, 3592, 2193, 1383, 743, 698, 243, 159}, {99, 194, 244, 172, 129, 89, 69, 49, 35, 60, 33, 19});

  dataVsCentrality[4] = RawDataVsCentrality(70, 90, {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5}, {384., 1369., 1997., 2052., 1684., 471.}, {35., 56., 70., 70., 63., 56.}, {0, 0.3, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12}, {1785., 1518, 2108., 1508., 1083., 728., 414., 252., 155., 135., 50.}, {56., 65, 76., 55., 46., 36., 27., 20., 16., 15., 8.});

  
  TFile *outputFile = TFile::Open("NJpsi_Raw.root", "recreate");
  outputFile->cd();
  for (int iCent = 0; iCent < numberOfCentralityBins; iCent++)
  {
    if (dataVsCentrality[iCent].VerifyData())
    {
      cout << Form("Filling histograms for centrality %g-%g", dataVsCentrality[iCent].fCentMin, dataVsCentrality[iCent].fCentMax) << endl;

      TH1F *histoVsPt = dataVsCentrality[iCent].MakePtHisto();
      histoVsPt->Write();
      delete histoVsPt;

      TH1F *histoVsRap = dataVsCentrality[iCent].MakeRapHisto();
      histoVsRap->Write();
      delete histoVsRap;
    }
  }

  outputFile->Close();
  delete outputFile;
}
