#include <iostream>
#include <stdio.h>
#include "TString.h"

enum eVariable{kRawYield, kInvariantYield, kCRpp, kRaa, kAccEff};
enum eFunction{kPowLawFree, kPowLawFixed, kLevy, kUA1, kWoodSaxonFixed, kWoodSaxonFree, kPol1, kConst, kPol3, kRatioLevy, kRatioPowerLaw};

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TString SetNameTest(eVariable varType, eFunction fFit, Bool_t isStatistical, Bool_t isSystematic, Double_t minFit=0, Double_t maxFit=15, Double_t fixedParam=0)//Attention ordre a changé
{
    TString sFunction;
    switch (fFit)
    {
        case kPowLawFree:
        sFunction = "PowerLawFree";
        break;
        case kPowLawFixed:
        sFunction = "PowerLawFixed";
        break;
        case kLevy:
        sFunction = "Levy";
        break;
        case kUA1:
        sFunction = "UA1";
        break;
        case kWoodSaxonFixed:
        sFunction = Form("WoodSaxonLikeFixed%.3f", fixedParam);
        break;
        case kWoodSaxonFree:
        sFunction = "WoodSaxonLikeFree";
        break;
        case kPol1:
        sFunction = "Linear";
        break;
        case kPol3:
        sFunction = "Pol3";
        break;
        case kConst:
        sFunction = "Constant";
        break;
        case kRatioLevy:
        sFunction = "RatioLevy";
        break;
        case kRatioPowerLaw:
        sFunction = "RatioPowerLaw";
        break;
    }

    TString sError;
    if (isStatistical && !isSystematic) sError = "ErrStat";
    else if (isSystematic && !isStatistical) sError = "ErrSyst";
    else sError = "ErrStatAndSyst";

    TString s;
    switch (varType)
    {
        case kInvariantYield:
            if(minFit == 0.65) s.Form("%s_%s_%.2f-%.f",sFunction.Data(), sError.Data(), minFit, maxFit);
            else if(minFit == 0.3) s.Form("%s_%s_%.1f-%.f",sFunction.Data(), sError.Data(), minFit, maxFit);
            else s.Form("%s_%s_%.f-%.f",sFunction.Data(), sError.Data(), minFit, maxFit);
        break; 

        case kRaa:
            if(minFit == 0.65) s.Form("%s_%s_%.2f-%.f",sFunction.Data(), sError.Data(), minFit, maxFit);
            else s.Form("%s_%s_%.f-%.f",sFunction.Data(), sError.Data(), minFit, maxFit);
        break;

        case kAccEff:
            if(minFit == 0.65) s.Form("%s_%s_%.2f-%.f", sFunction.Data(), sError.Data(), minFit, maxFit);
            else s.Form("%s_%s_%.f-%.f", sFunction.Data(), sError.Data(), minFit, maxFit);
        break;

        case kCRpp:
            s.Form("%s_%s",sFunction.Data(), sError.Data());
        break;
    }
    return s;
}


//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TString SetRangeValue(eVariable varType, Double_t minY, Double_t maxY, Double_t minPt, Double_t maxPt, Int_t minCent=0, Int_t maxCent=0)
{
    TString s;
    switch (varType)
    {
        case kCRpp:
            s.Form("rapidity%.1f-%.1f_pT%.1f-%.1f", minY, maxY, minPt, maxPt); //NameSystAccEff + Cr pp name result + range value Crpp Name
        break;
   
        default:
            s.Form("rapidity%.1f-%.1f_pT%.1f-%.1f_centrality%d-%d", minY, maxY, minPt, maxPt, minCent, maxCent); //rawYieldResults + Acceptance efficiency results + range value name yield + RaaNameResults + range value name Raa
        break;
    }
     
    return s;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TString SetRangeFunction(Double_t minY, Double_t maxY, Int_t minCent, Int_t maxCent)
{
    TString s;
    s.Form("rapidity%.1f-%.1f_centrality%d-%d", minY, maxY, minCent, maxCent); //results Acc Eff (necessite changement sur test name) + range function name yield + range function name Raa
    return s;
}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TString SetFullnameTest(eFunction fCRpp, eFunction fAccEff, eFunction fRaa, Double_t minFitRaa=0.65, Double_t maxFitRaa=15, Double_t fixedParam=0)
{
    TString sFunctionCRpp;
    switch (fCRpp)
    {
        case kPowLawFree:
        sFunctionCRpp = "PowerLawFree";
        break;
        case kPowLawFixed:
        sFunctionCRpp = "PowerLawFixed";
        break;
        case kLevy:
        sFunctionCRpp = "Levy";
        break;
        case kUA1:
        sFunctionCRpp = "UA1";
        break;
    }
    TString sFunctionRaa;
    switch (fRaa)
    {    
        case kWoodSaxonFixed:
        sFunctionRaa = Form("WoodSaxonLikeFixed%.3f", fixedParam);
        break;
        case kWoodSaxonFree:
        sFunctionRaa = "WoodSaxonLikeFree";
        break;
        case kPol1:
        sFunctionRaa = "Linear";
        break;
        case kConst:
        sFunctionRaa = "Constant";
    }
    TString sFunctionAccEff;
    switch (fAccEff)
    {    
        case kPol3:
        sFunctionAccEff = "Pol3";
        break;
        case kRatioLevy:
        sFunctionAccEff = "RatioLevy";
        break;
        case kRatioPowerLaw:
        sFunctionAccEff = "RatioPowerLaw";
        break;
    }
    
    TString sRangeRaa;
    if(minFitRaa == 0.65) sRangeRaa.Form("%.2f-%.f", minFitRaa, maxFitRaa);
    else sRangeRaa.Form("%.f-%.f", minFitRaa, maxFitRaa);

    TString s;
    s.Form("%s_%s_%s_%s", sFunctionRaa.Data(), sRangeRaa.Data(), sFunctionCRpp.Data(), sFunctionAccEff.Data());
    return s;
}