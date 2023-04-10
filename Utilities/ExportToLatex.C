#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsRaa.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsPPCrossSection.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsSignalExtraction.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsAcceptanceEfficiency.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


//-------------------------------------------------------------------------------------------------
//tableaux de la valeur moyenne du bin, de son ecart type et de sa valeur inférieure
std::vector<std::vector<Double_t>> xRange{
    {0.15, 0.65, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5},
    {0.15, 0.65, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5},
    {0.15, 0.65, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5},
    {0.15, 0.65, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11, 13.5},
    {0.15, 0.65, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11}
};

std::vector<std::vector<Double_t>> dxRange{
    {0.15, 0.35, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5},
    {0.15, 0.35, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5},
    {0.15, 0.35, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5},
    {0.15, 0.35, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1.5},
    {0.15, 0.35, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1}
};

Double_t pTRangeVScent[4]={0, 0.3, 1, 2};

Int_t cent[6]={0, 10, 30, 50, 70, 90};

void ExportToLatex(eVariable varType)
{
    Int_t nCentRanges;
    if(varType==kCRpp)nCentRanges=1;
    else nCentRanges=5;

    TString sType;
    switch (varType)
    {
    case kCRpp:
        sType = "\\frac{\\mathrm{d}^{2}\\sigma}{\\mathrm{d}p_{\\mathrm{T}}\\mathrm{d}y}";
        break;

    case kRaa:
        sType = "R_{\\mathrm{AA}}^{\\jpsi}";
        break;

    case kAccEff:
        sType = "\\left(\\mathcal{A}\\times\\epsilon\\right)_{\\jpsi}";
        break;

    case kRawYield:
        sType = "N_{\\jpsi}^{\\mathrm{raw}}";
    }

    Double_t nValue, errStat, errSyst, ratioStat, ratioSyst;

    for(int i=0; i<nCentRanges; i++)
    {
        if(varType!=kCRpp) printf("\\multicolumn{2}{c}{%i-%i\\%%}\\\\ \n", cent[i], cent[i+1]);
        printf("\\hline \n");
        printf("\\pt\\ (\\gevc) & $ %s \\pm stat (\\%%) \\pm syst (\\%%)$  \\\\ \n", sType.Data());
        printf("\\hline \n");

        for(int j=0; j<xRange[i].size(); j++)
        {
            if(varType==kRaa)
            {
                // cout << SetRangeValue(kRaa, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1]) << endl;
                nValue=RaaResultsVSpT[SetRangeValue(kRaa, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1])][kRaaValue];
                // cout << nValue << endl;
                errStat=RaaResultsVSpT[SetRangeValue(kRaa, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1])][kRaaStat];
                errSyst=RaaResultsVSpT[SetRangeValue(kRaa, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1])][kRaaSyst];
            }

            // if(varType==kRawYield)
            // {
            //     // cout << SetRangeValue(kRaa, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1]) << endl;
            //     nValue=RaaResults[SetRangeValue(kRaa, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1])][kRaaValue];
            //     // cout << nValue << endl;
            //     errStat=RaaResults[SetRangeValue(kRaa, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1])][kRaaStat];
            //     errSyst=RaaResults[SetRangeValue(kRaa, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1])][kRaaSyst];
            // }
            
            else if(varType==kCRpp)
            {
                // cout << SetRangeValue(kRaa, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1]) << endl;
                nValue=ppCrossSection[SetRangeValue(kCRpp, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1])][kCRppValue];
                // cout << nValue << endl;
                errStat=ppCrossSection[SetRangeValue(kCRpp, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1])][kCRppStatError];
                errSyst=ppCrossSection[SetRangeValue(kCRpp, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1])][kCRppSystError];
            }

            else if(varType==kAccEff)
            {
                // cout << SetRangeValue(kRaa, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1]) << endl;
                nValue=AccEffHadro[SetRangeValue(kAccEff, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1])][kAccEffValue];
                // cout << nValue << endl;
                errStat=AccEffHadro[SetRangeValue(kAccEff, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1])][kAccEffStatError];
            }

            ratioStat=errStat/nValue*100;
            if(varType!=kAccEff) ratioSyst=errSyst/nValue*100;

            if(varType==kAccEff) printf("%.f-%.f & $~ %.4f ~\\pm~ %.4f~(%.1f) ~$ \\\\ \n", xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], nValue, errStat, ratioStat);

            else if(varType==kCRpp) printf("%.f-%.f & $~ %.2e ~\\pm~ %.2e~(%.1f) ~\\pm~ %.2e~(%.1f) ~$ \\\\ \n", xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], nValue, errStat, ratioStat, errSyst, ratioSyst);

            else printf("%.f-%.f & $~ %.2f ~\\pm~ %.2f~(%.1f) ~\\pm~ %.2f~(%.1f) ~$ \\\\ \n", xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], nValue, errStat, ratioStat, errSyst, ratioSyst);
   
        }
        printf("\\hline \n");
        printf("\n");
    }
}

void ExportToLatexVScent()
{
    TString sType = "R_{\\mathrm{AA}}^{\\jpsi}";
    Int_t nCentRanges=5;

    Double_t nValue, errStat, errSyst, ratioStat, ratioSyst;

    for(int i=0; i<3; i++)
    {
        printf("\\multicolumn{2}{c}{%.1f-%.1f%~(\\gevc)}\\\\ \n", pTRangeVScent[i], pTRangeVScent[i+1]);
        printf("\\hline \n");
        printf("cent (\\%%) & $ %s \\pm stat (\\%%) \\pm syst (\\%%)$  \\\\ \n", sType.Data());
        printf("\\hline \n");
    
        for(int j=0; j<5; j++)
        {

            // cout << SetRangeValue(kRaa, -4, -2.5, xRange[i][j]-dxRange[i][j], xRange[i][j]+dxRange[i][j], cent[i], cent[i+1]) << endl;
            nValue=RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, pTRangeVScent[i], pTRangeVScent[i+1], cent[j], cent[j+1])][kRaaValue];
            // cout << nValue << endl;
            errStat=RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, pTRangeVScent[i], pTRangeVScent[i+1], cent[j], cent[j+1])][kRaaStat];
            errSyst=RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, pTRangeVScent[i], pTRangeVScent[i+1], cent[j], cent[j+1])][kRaaSyst];

            ratioStat=errStat/nValue*100;
            ratioSyst=errSyst/nValue*100;
            printf("%i-%i & $~ %.2f ~\\pm~ %.2f~(%.1f) ~\\pm~ %.2f~(%.1f) ~$ \\\\ \n", cent[j], cent[j+1], nValue, errStat, ratioStat, errSyst, ratioSyst);
            
        }
        printf("\\hline \n");
        printf("\n");
    }
}