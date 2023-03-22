#include "RootIncluders.C"
#include "Common.C"
#include "Iteration/FitYields.C"
#include "Iteration/GetWeights.C"
#include "TRandom.h"
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
int findLatestIteration(int centMin = 0, int centMax = 10)
{
    int foundIter = -1;
    TString baseDir = Form("Cent-%dto%d/AccEffi", centMin, centMax);
    for (int iter = 0; iter < 100; ++iter)
    {
        TString dirName = Form("%s/iter-%i", baseDir.Data(), iter);
        if (gSystem->AccessPathName(dirName.Data()))
        {
            break;
        }
        foundIter = iter;
    }
    if (foundIter < 0)
    {
        std::cout << "Error: cannot find last iteration for centrality class " << centMin << " - " << centMax << std::endl;
    }
    return foundIter;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
TH1 *getAccEffInput(int centMin = 0, int centMax = 10, bool isPt = true)
{
    int lastIter = findLatestIteration(centMin, centMax);
    if (lastIter < 0)
    {
        return nullptr;
    }
    TFile *inputAccEffFile = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-%d/AccEffiValues.root", centMin, centMax, lastIter));
    TString hname = isPt ? "histoAccEffiVsPt" : "histoAccEffiVsRap";
    TH1 *histoAccEff = static_cast<TH1 *>(inputAccEffFile->Get(hname.Data()));
    if (histoAccEff)
    {
        histoAccEff->SetDirectory(0);
    }
    else
    {
        std::cout << "Cannot find " << hname.Data() << std::endl;
    }

    delete inputAccEffFile;
    return histoAccEff;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
TH1 *getInput(int centMin = 0, int centMax = 10, bool isPt = true)
{
    int lastIter = findLatestIteration(centMin, centMax);
    if (lastIter < 0)
    {
        return nullptr;
    }

    //-----------
    // Gets generated distribution from the latest step
    // TFile *inputMCFile = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-%i/AccEffiValues.root", centMin, centMax, lastIter));
    // TString hname = Form("histoGen%s", (isPt ? "Pt" : "Rapidity"));
    // TH1 *histo = static_cast<TH1 *>(inputMCFile->Get(hname.Data()));
    // if (isPt)
    // {
    //     FitPtYield(histo, false);
    // }
    // else
    // {
    //     FitRapYield(histo, false);
    // }
    //-----------

    //-----------
    // Gets the generated distribution as the Raw yield divided by the latest AccxEff
    // The fit functions are designed to get the AccxEff of the asked iteration - 1
    // Since the original idea was to use it as the input of the current iteration.
    // So, if we want the very latest iteration, we must ask for the lastIter + 1
    if (isPt)
    {
        FitPtYield(centMin, centMax, lastIter + 1, true);
    }
    else
    {
        FitRapYield(centMin, centMax, lastIter + 1, true);
    }
    TString hname = Form("histo%sYield_cent%dto%d", (isPt ? "Pt" : "Rap"), centMin, centMax);
    TH1 *histo = static_cast<TH1 *>(gROOT->FindObject(hname.Data()));
    //-----------

    if (!histo)
    {
        std::cout << "Cannot find " << hname.Data() << std::endl;
        return nullptr;
    }
    double normFactor = 1. / histo->Integral("width");
    histo->Scale(normFactor);
    TF1 *func = static_cast<TF1 *>(histo->GetListOfFunctions()->At(0));
    func->SetParameter(0, func->GetParameter(0) * normFactor);
    return histo;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
void smearHisto(TH1 *histo)
{
    for (int ibin = 1; ibin <= histo->GetXaxis()->GetNbins(); ++ibin)
    {
        histo->SetBinContent(ibin, gRandom->Gaus(histo->GetBinContent(ibin), histo->GetBinError(ibin)));
    }
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
TF1 *getSmearedFunc(const TH1 *inputHisto, const char *smearedHistoName, bool isPt)
{
    TH1 *histo = static_cast<TH1 *>(inputHisto->Clone(smearedHistoName));
    smearHisto(histo);
    int customStepCounter;
    if (isPt)
    {
        customStepCounter = FitPtYield(histo, false);
    }
    else
    {
        customStepCounter = FitRapYield(histo, false);
    }
    if (customStepCounter >= 100)
    {
        std::cout << "Warning: fit to " << smearedHistoName << " did not converge" << std::endl;
        return nullptr;
    }
    TF1 *func = static_cast<TF1 *>(histo->GetListOfFunctions()->At(0));
    if (isPt && func->Eval(10.) < 1e-5)
    {
        std::cout << "Warning: fit to " << smearedHistoName << " has problems" << std::endl;
        return nullptr;
    }
    return func;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
std::vector<TF1 *> getSmearedFunctions(int centMin = 0, int centMax = 10, bool isPt = true, int nRndm = 50)
{
    std::vector<TF1 *> funcs;
    TH1 *histo = getInput(centMin, centMax, isPt);
    for (int irndm = 0; irndm < nRndm; ++irndm)
    {
        TString smearName = Form("%s_rndm_%i", histo->GetName(), irndm);
        TF1 *func = getSmearedFunc(histo, smearName.Data(), isPt);
        if (func)
        {
            funcs.emplace_back(func);
        }
        else
        {
            std::cout << "Sample again" << std::endl;
            --irndm;
        }
    }
    return funcs;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
TChain *getChain(const char *periods = "LHC16e2,LHC16e2plus,LHC19a2_18q,LHC19a2_18r", const char *baseDir = "GridAnalysis/outputFiles", const char *baseFilename = "AnalysisResultsNew.root")
{
    TChain *mergedTree = new TChain("eventsTree");
    TString speriods(periods);
    TObjArray *objArray = speriods.Tokenize(",");
    for (int iperiod = 0; iperiod < objArray->GetEntries(); ++iperiod)
    {
        TString filename = Form("%s/%s/%s", baseDir, objArray->At(iperiod)->GetName(), baseFilename);
        if (gSystem->AccessPathName(filename.Data()))
        {
            std::cout << "Cannot find " << filename.Data() << std::endl;
            continue;
        }
        std::cout << "Adding " << filename.Data() << std::endl;
        mergedTree->Add(filename.Data());
    }
    return mergedTree;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
std::vector<std::vector<TH1 *>> reweightProd(TChain *tree, std::vector<std::vector<TF1 *>> newPtRap, TF1 *refPt, TF1 *refRap, int centMin, int centMax)
{
    enum enumDimuon
    {
        kDimuonInvMass,
        kDimuonPt,
        kDimuonRapidity
    };

    size_t nfuncs = newPtRap.size();

    int numberOfRapidityBins;
    const Double_t *arrayRapidityBins = GetBinsArray(centMin, centMax, "Rap", numberOfRapidityBins);

    int numberOfPtBins;
    const Double_t *arrayPtBins = GetBinsArray(centMin, centMax, "Pt", numberOfPtBins);

    TH1 *histoPtWeight = nullptr;
    TH1 *histoRapWeight = nullptr;

    int iteration = findLatestIteration(centMin, centMax);
    histoPtWeight = GetPtWeight(centMin, centMax, iteration);
    histoRapWeight = GetRapWeight(centMin, centMax, iteration);
    for (int iIteration = iteration - 1; iIteration > 0; iIteration--)
    {
        histoPtWeight->Multiply(GetPtWeight(centMin, centMax, iIteration));
        histoRapWeight->Multiply(GetRapWeight(centMin, centMax, iIteration));
    }

    std::vector<TH1F> histoGenPt, histoGenRapidity, histoRecPt, histoRecRapidity;
    for (size_t ifunc = 0; ifunc < nfuncs; ++ifunc)
    {
        histoGenPt.emplace_back(Form("histoGenPt_%zu", ifunc), "", numberOfPtBins, arrayPtBins);
        histoGenPt.back().Sumw2();
        histoRecPt.emplace_back(Form("histoRecPt_%zu", ifunc), "", numberOfPtBins, arrayPtBins);
        histoRecPt.back().Sumw2();
        histoGenRapidity.emplace_back(Form("histoGenRapidity_%zu", ifunc), "", numberOfRapidityBins, arrayRapidityBins);
        histoGenRapidity.back().Sumw2();
        histoRecRapidity.emplace_back(Form("histoRecRapidity_%zu", ifunc), "", numberOfRapidityBins, arrayRapidityBins);
        histoRecRapidity.back().Sumw2();
    }

    std::vector<double> normPt, normRap;
    for (auto pair : newPtRap)
    {
        normPt.emplace_back(refPt->Integral(arrayPtBins[0], arrayPtBins[numberOfPtBins]) / pair[0]->Integral(arrayPtBins[0], arrayPtBins[numberOfPtBins]));
        normRap.emplace_back(refRap->Integral(arrayRapidityBins[0], arrayRapidityBins[numberOfRapidityBins]) / pair[1]->Integral(arrayRapidityBins[0], arrayRapidityBins[numberOfRapidityBins]));
    }

    Int_t nEntries = tree->GetEntries();

    std::vector<double> *tempoVectorGenDimuon = 0;
    std::vector<double> *tempoVectorRecDimuon = 0;
    Int_t runNumber;
    Float_t centrality;
    tree->SetBranchAddress("GenDimuon", &tempoVectorGenDimuon);
    tree->SetBranchAddress("RecDimuon", &tempoVectorRecDimuon);
    tree->SetBranchAddress("runNumber", &runNumber);
    tree->SetBranchAddress("centrality", &centrality);

    if (!histoNeventPerRunGen)
    {
        PrintMessage("Finding the number of simulated events per run and centrality", "-");
        histoNeventPerRunGen = static_cast<TH1I *>(histoCMULEventPerRun->Clone("histoNEventsPerRunGen"));
        histoNeventPerRunGen->Reset();
        histoNJpsiVsCentralityGen = static_cast<TH1F *>(histoNJpsiVsCentrality->Clone("histoNJpsiVsCentralityGen"));
        histoNJpsiVsCentralityGen->Reset();
        for (Int_t iEvent = 0; iEvent < nEntries; iEvent++)
        {
            if (iEvent % 100 == 0)
            {
                printProgress((double)iEvent / nEntries);
            }
            tree->GetEntry(iEvent);
            histoNeventPerRunGen->Fill(runNumber);
            histoNJpsiVsCentralityGen->Fill(centrality);
        }
        std::cout << std::endl;

        histoNeventPerRunGen->Scale(totalCMULEvents / histoNeventPerRunGen->Integral());
        histoCMULEventPerRun->Divide(histoNeventPerRunGen);
        histoNJpsiVsCentralityGen->Scale(totalNumberOfJpsi / histoNJpsiVsCentralityGen->Integral());
        histoNJpsiVsCentrality->Divide(histoNJpsiVsCentralityGen);
    }

    PrintMessage("Re-weight production", "-");

    for (Int_t iEvent = 0; iEvent < nEntries; iEvent++)
    {

        if (iEvent % 100 == 0)
        {
            printProgress((double)iEvent / nEntries);
        }

        tree->GetEntry(iEvent);
        if (centrality < centMin || centrality >= centMax)
        {
            continue;
        }

        //Weight CMUL per run and CMUL vs centrality

        double commonWeight = histoCMULEventPerRun->GetBinContent(histoCMULEventPerRun->GetXaxis()->FindBin(runNumber)) * histoNJpsiVsCentrality->GetBinContent((int)(centrality / 10.) + 1);
        double ptGen = tempoVectorGenDimuon->at(kDimuonPt);
        double rapGen = tempoVectorGenDimuon->at(kDimuonRapidity);

        bool hasGen = (rapGen < arrayRapidityBins[numberOfRapidityBins] && rapGen > arrayRapidityBins[0] && ptGen > arrayPtBins[0] && ptGen < arrayPtBins[numberOfPtBins]);

        bool hasRec = false;
        double ptRec = 0, rapRec = 0;
        if (tempoVectorRecDimuon->size() > 0)
        {
            ptRec = tempoVectorRecDimuon->at(kDimuonPt);
            rapRec = tempoVectorRecDimuon->at(kDimuonRapidity);
            hasRec = true;
        }

        if (hasGen)
        {
            // This is the original weight, that would produce the final AccxEff
            commonWeight *= histoPtWeight->GetBinContent(histoPtWeight->GetXaxis()->FindBin(tempoVectorGenDimuon->at(kDimuonPt))) * histoRapWeight->GetBinContent(histoRapWeight->GetXaxis()->FindBin(tempoVectorGenDimuon->at(kDimuonRapidity)));
            // This is the part of the new weight which is common to all functions
            commonWeight *= 1. / (refPt->Eval(ptGen) * refRap->Eval(rapGen));
        }

        for (size_t ifunc = 0; ifunc < nfuncs; ++ifunc)
        {
            double weightDimuon = commonWeight;
            if (hasGen)
            {
                // This is the modification w.r.t. the original weight which depends on the function
                weightDimuon *= newPtRap[ifunc][0]->Eval(ptGen) * newPtRap[ifunc][1]->Eval(rapGen) * normPt[ifunc] * normRap[ifunc];
                histoGenPt[ifunc].Fill(ptGen, weightDimuon);
                histoGenRapidity[ifunc].Fill(rapGen, weightDimuon);
            }
            if (hasRec && weightDimuon > 0)
            {
                histoRecPt[ifunc].Fill(ptRec, weightDimuon);
                histoRecRapidity[ifunc].Fill(rapRec, weightDimuon);
            }
        }
    }
    std::cout << std::endl;

    std::vector<std::vector<TH1 *>> accEff;
    for (size_t ifunc = 0; ifunc < nfuncs; ++ifunc)
    {
        TH1F *histoAccEffiVsPt = new TH1F(histoRecPt[ifunc]);
        histoAccEffiVsPt->SetName(Form("histoAccEffiVsPt_%zu", ifunc));
        histoAccEffiVsPt->Divide(&histoGenPt[ifunc]);

        TH1F *histoAccEffiVsRap = new TH1F(histoRecRapidity[ifunc]);
        histoAccEffiVsRap->SetName(Form("histoAccEffiVsRap_%zu", ifunc));
        histoAccEffiVsRap->Divide(&histoGenRapidity[ifunc]);
        accEff.push_back({histoAccEffiVsPt, histoAccEffiVsRap});
    }

    return accEff;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
std::vector<std::vector<TF1 *>> getPairs(const std::vector<TF1 *> ptFuncs, const std::vector<TF1 *> rapFuncs, bool allCombinations = false)
{
    std::vector<std::vector<TF1 *>> funcs;
    if (allCombinations)
    {
        for (auto ptFunc : ptFuncs)
        {
            for (auto rapFunc : rapFuncs)
            {
                funcs.push_back({ptFunc, rapFunc});
            }
        }
        return funcs;
    }

    if (ptFuncs.size() != rapFuncs.size())
    {
        std::cout << "Error: the number of pt functions and rapidity functions must be the same" << std::endl;
        return funcs;
    }

    for (size_t ifunc = 0; ifunc < ptFuncs.size(); ++ifunc)
    {
        funcs.push_back({ptFuncs[ifunc], rapFuncs[ifunc]});
    }
    return funcs;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
std::vector<TH1 *> getUncertainties(std::vector<std::vector<TH1 *>> accEff, TH1 *refAccEffPt, TH1 *refAccEffRap)
{
    std::vector<TH1 *> errHistos{static_cast<TH1 *>(refAccEffPt->Clone("errHisto_pt")), static_cast<TH1 *>(refAccEffRap->Clone("errHisto_y"))};
    for (int itype = 0; itype < 2; ++itype)
    {
        for (int ibin = 1; ibin <= errHistos[itype]->GetXaxis()->GetNbins(); ++ibin)
        {
            double sum = 0.;
            double squaresum = 0.;
            double nentries = 0.;
            for (auto histos : accEff)
            {
                double binCont = histos[itype]->GetBinContent(ibin);
                nentries += 1.;
                sum += binCont;
                squaresum += binCont * binCont;
            }
            double rms = TMath::Sqrt((squaresum - sum * sum / nentries) / (nentries - 1));
            errHistos[itype]->SetBinError(ibin, rms);
        }
    }
    return errHistos;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
void normalizeFunc(TF1 *func, double min, double max)
{
    func->SetParameter(0, func->GetParameter(0) / func->Integral(min, max));
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
void normalizeFunc(TF1 *func, TH1 *histo)
{
    normalizeFunc(func, histo->GetXaxis()->GetBinLowEdge(1), histo->GetXaxis()->GetBinUpEdge(histo->GetXaxis()->GetNbins()));
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
void drawUncertainties(std::vector<std::vector<TH1 *>> accEff, TH1 *refAccEffPt, TH1 *refAccEffRap, std::vector<std::vector<TF1 *>> newPtRap, TF1 *refPt, TF1 *refRap, Int_t centMin, Int_t centMax, TString typeUncert)
{
    std::vector<TF1 *> refPtRap{refPt, refRap};
    std::vector<TH1 *> refAccEffPtRap{refAccEffPt, refAccEffRap};
    auto rmsHistos = getUncertainties(accEff, refAccEffPt, refAccEffRap);
    gSystem->Exec(Form("mkdir -p Cent-%dto%d/AccEff%s/", centMin, centMax, typeUncert.Data()));
    TFile *outputFile = TFile::Open(Form("Cent-%dto%d/AccEff%s/AccEffiValues.root", centMin, centMax, typeUncert.Data()), "recreate");

    for (int itype = 0; itype < 2; ++itype)
    {
        TLegend* legend = SetLegend(centMin, centMax);
        //
        TString var = itype == 0 ? "pt" : "y";
        TString name = "input_variation_" + var;
        TCanvas *canFunc = new TCanvas(name, name);
        if (itype == 0)
        {
            canFunc->SetLogy();
        }
        SetCanvasStyle(canFunc);
        refPtRap[itype]->SetLineColor(1);
        normalizeFunc(refPtRap[itype], refAccEffPtRap[itype]);
        refPtRap[itype]->Draw();
        for (auto funcs : newPtRap)
        {
            funcs[itype]->Draw("same");
            normalizeFunc(funcs[itype], refAccEffPtRap[itype]);
        }
        legend->Draw();
        canFunc->SaveAs(Form("Cent-%dto%d/AccEff%s/InputVariationVs%s_%d-%d.png", centMin, centMax, typeUncert.Data(), var.Data(), centMin, centMax));

        //
        name = "acceff_variation_" + var;
        TCanvas *canAccEff = new TCanvas(name, name);
        SetCanvasStyle(canAccEff);
        TString drawOpt = "";
        for (auto histos : accEff)
        {
            SetHistoStyle(histos[itype], kBlue, kFullCircle, 0.8, (itype == 0 ? ptAxisTitle : yAxisTitle), accEffiAxisTitle, 1, kFALSE);
            histos[itype]->Draw(drawOpt.Data());
            drawOpt = "same";
        }
        legend->Draw();
        canAccEff->SaveAs(Form("Cent-%dto%d/AccEff%s/AccEffVariationVs%s_%d-%d.png", centMin, centMax, typeUncert.Data(), var.Data(), centMin, centMax));

        //
        name = "acceff_variation_over_ref_" + var;
        TCanvas *canRatio = new TCanvas(name, name);
        SetCanvasStyle(canRatio);
        drawOpt = "ap";
        for (auto histos : accEff)
        {
            TGraph *gr = new TGraph(histos[itype]);
            for (int ibin = 1; ibin <= histos[itype]->GetXaxis()->GetNbins(); ++ibin)
            {
                gr->SetPoint(ibin - 1, histos[itype]->GetBinCenter(ibin), histos[itype]->GetBinContent(ibin) / refAccEffPtRap[itype]->GetBinContent(ibin));
            }
            gr->GetYaxis()->SetRangeUser(0.9, 1.1);
            gr->GetYaxis()->SetTitle("AccEff ratio");
            gr->GetXaxis()->SetTitle(histos[itype]->GetXaxis()->GetTitle());
            gr->SetMarkerColor(kBlue);
            gr->Draw(drawOpt.Data());
            drawOpt = "p";
            if(typeUncert=="Syst")
            {
                gr->SetName(Form("%s", name.Data()));
                gr->Write();
            }
        }
        legend->Draw();
        canRatio->SaveAs(Form("Cent-%dto%d/AccEff%s/AccEffVariationOverRefVs%s_%d-%d.png", centMin, centMax, typeUncert.Data(), var.Data(), centMin, centMax));
        
        //
        name = "acceff_rms_" + var;
        TCanvas *canRMS = new TCanvas(name, name);
        SetCanvasStyle(canRMS);
        TGraph *gr = new TGraph(rmsHistos[itype]);
        // rmsHistos[itype]->Draw();
        // SetHistoStyle(rmsHistos[itype], kBlack, kFullCircle, 0.9, (itype == 0 ? ptAxisTitle : yAxisTitle), accEffiAxisTitle, 1, kFALSE);
        for (int ibin = 1; ibin <= rmsHistos[itype]->GetXaxis()->GetNbins(); ++ibin)
        {
            gr->SetPoint(ibin - 1, rmsHistos[itype]->GetBinCenter(ibin), rmsHistos[itype]->GetBinError(ibin) / rmsHistos[itype]->GetBinContent(ibin));
        }
        gr->SetMarkerStyle(kFullCircle);
        gr->SetMarkerColor(kBlue);
        gr->GetYaxis()->SetTitle("AccEff RMS/mean");
        gr->GetXaxis()->SetTitle(rmsHistos[itype]->GetXaxis()->GetTitle());
        gr->Draw("ap");
        legend->Draw();
        canRMS->SaveAs(Form("Cent-%dto%d/AccEff%s/AccEffRMSOverMeanVs%s_%d-%d.png", centMin, centMax, typeUncert.Data(), var.Data(), centMin, centMax));
        if(typeUncert=="Stat")
        {
            gr->SetName(Form("%s", name.Data()));
            gr->Write();
        }
    }
    delete outputFile;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
void accEffStat(int centMin = 0, int centMax = 10, int nRndms = 100)
{
    auto funcPt = getSmearedFunctions(centMin, centMax, true, nRndms);
    auto funcRap = getSmearedFunctions(centMin, centMax, false, nRndms);
    auto refPt = static_cast<TF1 *>(getInput(centMin, centMax, true)->GetListOfFunctions()->At(0));
    auto refRap = static_cast<TF1 *>(getInput(centMin, centMax, false)->GetListOfFunctions()->At(0));
    auto funcs = getPairs(funcPt, funcRap, false);

    TChain *tree = getChain();
    auto accEff = reweightProd(tree, funcs, refPt, refRap, centMin, centMax);

    auto refAccEffPt = getAccEffInput(centMin, centMax, true);
    auto refAccEffRap = getAccEffInput(centMin, centMax, false);

    drawUncertainties(accEff, refAccEffPt, refAccEffRap, funcs, refPt, refRap, centMin, centMax, "Stat");
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
std::vector<TF1 *> getSystematicFunctions(bool isPt)
{
    std::vector<TF1 *> funcs;
    if (isPt)
    {
        TString formula = "[0] * x / TMath::Power(1. + TMath::Power(x/[1],[2]), [3])";
        std::vector<std::vector<double>> parameters = {
            {1., 4.608, 1.6918, 4.2423},     // 2.5 < y < 2.8
            {1., 5.0297, 1.6365, 4.7689},    // 2.8 < y < 3
            {1., 4.4274, 1.7741, 4.026},     // 3 < y < 3.2
            {1., 4.7547, 1.6953, 4.4997},    // 3.2 < y < 3.5
            {1., 4.5144, 1.7925, 4.2223},    // 3.5 < y < 3.8
            {1., 5.448, 1.5637, 5.5774},     // 3.8 < y < 4
            {1., 4.75208, 1.69247, 4.49224}, // 2.5 < y < 4 => Reference must be last
        };

        for (auto &pars : parameters)
        {
            TF1 *func = new TF1("ptFunc", formula.Data(), 0., 100.);
            func->SetNpx(5000);
            func->SetParameters(pars.data());
            funcs.emplace_back(func);
        }
    }
    else
    {
        TString formula = "[0] * TMath::Exp(-0.5*x*x/[1]/[1])";
        std::vector<std::vector<double>> parameters = {
            {1.6018, 3.5004},   // 0 < pt < 1
            {3.0763, 3.4912},   // 1 < pt < 2
            {2.7625, 3.2579},   // 2 < pt < 3
            {1.832, 3.1353},    // 3 < pt < 4
            {1.1356, 2.9421},   // 4 < pt < 5
            {0.78027, 2.5091},  // 5 < pt < 6
            {0.38906, 2.7312},  // 6 < pt < 7
            {0.27124, 2.3781},  // 7 < pt < 8
            {0.10791, 2.5896},  // 8 < pt < 10
            {0.050845, 2.2532}, // 10 < pt < inf
            // {0.026694, 1.815},  // ??
            // {0.003918, 2.249}}; // ??
            {1., 2.98887} // 0 < pt < inf => Reference must be last
        };
        for (auto &pars : parameters)
        {
            TF1 *func = new TF1("yFunc", formula.Data(), -4., -2.5);
            func->SetNpx(5000);
            func->SetParameters(pars.data());
            funcs.emplace_back(func);
        }
    }
    return funcs;
}
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
void accEffSyst(int centMin = 0, int centMax = 10)
{
    auto funcPt = getSystematicFunctions(true);
    auto funcRap = getSystematicFunctions(false);
    auto refPt = funcPt.back();
    funcPt.pop_back();
    auto refRap = funcRap.back();
    funcRap.pop_back();

    auto funcs = getPairs(funcPt, funcRap, true);

    TChain *tree = getChain();
    auto accEff = reweightProd(tree, funcs, refPt, refRap, centMin, centMax);

    auto refAccEffPt = getAccEffInput(centMin, centMax, true);
    auto refAccEffRap = getAccEffInput(centMin, centMax, false);

    drawUncertainties(accEff, refAccEffPt, refAccEffRap, funcs, refPt, refRap, centMin, centMax, "Syst");
}
