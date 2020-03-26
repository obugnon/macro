#include <iostream>
#include <stdio.h>
#include "TString.h"
//#include "FitFunctions.C"

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TString SetNameTest(Efunction fBackground, Efunction fSignal, Etails fTails, Double_t minFit, Double_t maxFit)
{
    TString sBackground, sSignal, sTails;
    switch (fBackground)
    {
        case kVWGQuadratic:
        sBackground = "qVWG";
        break;
        case kPol2Pol3:
        sBackground = "Pol2Pol3";
        break;
    }
    switch (fSignal)
    {
        case kCBExtended:
        sSignal = "CB2";
        break;
        case kNA60:
        sSignal = "NA60";
        break;
    }
    switch (fTails)
    {
        case kEMB:
        sTails = "Emb";
        break;
        case kPP:
        sTails = "ppData";
        break;
        case kSTARLIGHTcoh:
        sTails = "StarlightCoh";
        break;
        case kSTARLIGHTincoh:
        sTails = "StarlightIncoh";
        break;
    }
    TString s;
    s.Form("%s_%s_%s_m[%.1f-%.1f]",sBackground.Data(), sSignal.Data(), sTails.Data(), minFit, maxFit);
    return s;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TString SetRangeName(Double_t minY, Double_t maxY, Double_t minPt, Double_t maxPt, Int_t minCent, Int_t maxCent)
{
    TString s;
    s.Form("rapidity%.1f-%.1f_pT%.1f-%.1f_centrality%d-%d", minY, maxY, minPt, maxPt, minCent, maxCent);
    return s;

}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TString SetTailsName( Etails fTails, Efunction fSignal, Efunction fBackground, Double_t minY, Double_t maxY, Double_t minPt, Double_t maxPt, Double_t minFit, Double_t maxFit)
{
    TString sSignal, sTails, sBackground;
    switch (fBackground)
    {
        case kVWGQuadratic:
        sBackground = "qVWG";
        break;
        case kPol2Pol3:
        sBackground = "Pol2Pol3";
        break;
    }
    switch (fSignal)
    {
        case kCBExtended:
        sSignal = "CB2";
        break;
        case kNA60:
        sSignal = "NA60";
        break;
    }
    switch (fTails)
    {
        case kEMB:
        sTails = "Emb";
        break;
        case kPP:
        sTails = "ppData";
        break;
        case kSTARLIGHTcoh:
        sTails = "StarlightCoh";
        break;
        case kSTARLIGHTincoh:
        sTails = "StarlightIncoh";
        break;
    }    
    TString s;
    if(fTails == kPP)
    {
        if(minFit < 2.4 && maxFit < 4.7) s.Form("%s_%s_%s_firstRange_rapidity%.1f-%.1f_pT%.1f-%.1f", sSignal.Data(), sTails.Data(), sBackground.Data(), minY, maxY, minPt, maxPt);
        else s.Form("%s_%s_%s_secondRange_rapidity%.1f-%.1f_pT%.1f-%.1f", sSignal.Data(), sTails.Data(), sBackground.Data(), minY, maxY, minPt, maxPt);
    }
    else s.Form("%s_%s_rapidity%.1f-%.1f_pT%.1f-%.1f", sSignal.Data(), sTails.Data(), minY, maxY, minPt, maxPt);
    return s; 
}