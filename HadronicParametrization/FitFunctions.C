/*
 *  FitFunctions.C
 *
 *  Created by Ophelie Bugnon on 05/11/20.
 *
 */

Double_t xMinParam = 0;
Double_t xMaxParam = 8;
Double_t xMinIntegral = 0;
Double_t xMaxIntegral = 0.3;



//Function definitions
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
Double_t RatioLevyFunction(Double_t *x, Double_t *par)
{
	//par[0] = dNdy
	//par[1] = n1
	//par[2] = T1
	//par[3] = mass1

	Double_t bigCoef1=((par[1]-1)*(par[1]-2))/(par[1]*par[2]*(par[1]*par[2]+par[3]*(par[1]-2)));

	Double_t inPower1 = 1.0 + (TMath::Sqrt(x[0]*x[0]+par[3]*par[3])-par[3]) /(par[1]*par[2]);

	//par[4] = n2
	//par[5] = T2
	//par[6] = mass2

	Double_t bigCoef2 = ((par[4]-1)*(par[4]-2))/(par[4]*par[5]*(par[4]*par[5]+par[6]*(par[4]-2)));
	Double_t inPower2 = 1.0 + (TMath::Sqrt(x[0]*x[0]+par[6]*par[6])-par[6]) /(par[4]*par[5]);

	if((x[0] *bigCoef2 * TMath::Power(inPower2,(-1.0)*par[4]))==0) return 0;
	else return par[0] * (x[0] * bigCoef1 * TMath::Power(inPower1,(-1.0)*par[1]))/(x[0] *bigCoef2 * TMath::Power(inPower2,(-1.0)*par[4]));
}
//---------------------------------------------------------------------------------------------------
Double_t Pol3Function(Double_t *x, Double_t *par)
{
	//par[0] = Normalization
	//par[1] = slope
	return  par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}
//---------------------------------------------------------------------------------------------------
Double_t PowerLawFunction(Double_t *x, Double_t *par)
{
    //par[0] = Normalization
	//par[1] = pT0
	//par[2] = n1
	//par[3] = n2
    return par[0]*x[0]/TMath::Power((1+TMath::Power(x[0]/par[1],par[3])),par[2]);
}
//---------------------------------------------------------------------------------------------------
Double_t LevyFunction(Double_t *x, Double_t *par)
{
	//par[0] = Normalization
	//par[1] = n
	//par[2] = T
	//par[3] = mass
	return  x[0]*par[0]*(par[1]-1)*(par[1]-2)/( par[1]*par[2]*(par[1]*par[2]+par[3]*(par[1]-2) ) ) * TMath::Power(( 1 + (TMath::Sqrt(par[3]*par[3]+x[0]*x[0]) - par[3])/(par[1]*par[2])),-par[1]);
}
//---------------------------------------------------------------------------------------------------
Double_t WoodSaxonLikeFunction(Double_t *x, Double_t *par)
{
    //par[0] = Normalization
	//par[1] = DeltaRaa
	//par[2] = pT0
	//par[3] = sigmapT
    return par[0]*(1+ par[1]/(1+TMath::Exp((x[0]-par[2])/par[3])));
}
//---------------------------------------------------------------------------------------------------
Double_t Pol1Function(Double_t *x, Double_t *par)
{
	//par[0] = Normalization
	//par[1] = slope
	return  par[0] + par[1]*x[0];
}
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//Define Combined function
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
Double_t WoodSaxon_Levy_RatioLevy(Double_t *x, Double_t *par){
	//For Raa : Wood-Saxon = 4 parameters
	//For pp cross section : Levy = 4 parameters
	//For AccEff : ratio of Levy = 7 parameters

	return WoodSaxonLikeFunction(x, &par[0]) * LevyFunction(x, &par[4]) * RatioLevyFunction(x, &par[8]);
}
//---------------------------------------------------------------------------------------------------
Double_t WoodSaxon_Levy_Polynomial(Double_t *x, Double_t *par){
	//For Raa : Wood-Saxon = 4 parameters
	//For pp cross section : Levy = 4 parameters
	//For AccEff : 3rd degree polynomial = 4 parameters

	return WoodSaxonLikeFunction(x, &par[0]) * LevyFunction(x, &par[4]) * Pol3Function(x, &par[8]);
}
//---------------------------------------------------------------------------------------------------
Double_t WoodSaxon_PowerLaw_RatioLevy(Double_t *x, Double_t *par){
	//For Raa : Wood-Saxon = 4 parameters
	//For pp cross section : Power law = 4 parameters
	//For AccEff : ratio of Levy = 7 parameters

	return WoodSaxonLikeFunction(x, &par[0]) * PowerLawFunction(x, &par[4]) * RatioLevyFunction(x, &par[8]);
}
//---------------------------------------------------------------------------------------------------
Double_t WoodSaxon_PowerLaw_Polynomial(Double_t *x, Double_t *par){
	//For Raa : Wood-Saxon = 4 parameters
	//For pp cross section : Power law = 4 parameters
	//For AccEff : 3rd degree polynomial = 4 parameters

	return WoodSaxonLikeFunction(x, &par[0]) * PowerLawFunction(x, &par[4]) * Pol3Function(x, &par[8]);
}
//---------------------------------------------------------------------------------------------------
Double_t Linear_Levy_RatioLevy(Double_t *x, Double_t *par){
	//For Raa : Pol1 function = 2 parameters
	//For pp cross section : Levy = 4 parameters
	//For AccEff : ratio of Levy = 7 parameters

	return Pol1Function(x, &par[0]) * LevyFunction(x, &par[2]) * RatioLevyFunction(x, &par[6]);
}
//---------------------------------------------------------------------------------------------------
Double_t Linear_Levy_Polynomial(Double_t *x, Double_t *par){
	//For Raa : Pol1 function = 2 parameters
	//For pp cross section : Levy = 4 parameters
	//For AccEff : 3rd degree polynomial = 4 parameters

	return Pol1Function(x, &par[0]) * LevyFunction(x, &par[2]) * Pol3Function(x, &par[6]);
}
//---------------------------------------------------------------------------------------------------
Double_t Linear_PowerLaw_RatioLevy(Double_t *x, Double_t *par){
	//For Raa : Pol1 function = 2 parameters
	//For pp cross section : Power law = 4 parameters
	//For AccEff : ratio of Levy = 7 parameters

	return Pol1Function(x, &par[0]) * PowerLawFunction(x, &par[2]) * RatioLevyFunction(x, &par[6]);
}
//---------------------------------------------------------------------------------------------------
Double_t Linear_PowerLaw_Polynomial(Double_t *x, Double_t *par){
	//For Raa : Pol1 function = 2 parameters
	//For pp cross section : Power law = 4 parameters
	//For AccEff : 3rd degree polynomial = 4 parameters

	return Pol1Function(x, &par[0]) * PowerLawFunction(x, &par[2]) * Pol3Function(x, &par[6]);
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
TF1* func_WoodSaxon_Levy_RatioLevy = new TF1("WoodSaxon_Levy_RatioLevy", WoodSaxon_Levy_RatioLevy, xMinParam, xMaxParam, 15);
Double_t Ratio_WoodSaxon_Levy_RatioLevy(Double_t *x, Double_t *par){
	//For Raa : Wood-Saxon = 4 parameters
	//For pp cross section : Levy = 4 parameters
	//For AccEff : ratio of Levy = 7 parameters

    func_WoodSaxon_Levy_RatioLevy->SetParameters(par);

	return (func_WoodSaxon_Levy_RatioLevy->Integral(xMinIntegral, xMaxIntegral))/(func_WoodSaxon_Levy_RatioLevy->Integral(1., 8.));
}
//---------------------------------------------------------------------------------------------------
Double_t Ratio_WoodSaxon_Levy_Polynomial(Double_t *x, Double_t *par){
	//For Raa : Wood-Saxon = 4 parameters
	//For pp cross section : Levy = 4 parameters
	//For AccEff : 3rd degree polynomial = 4 parameters

	TF1* func = new TF1("WoodSaxon_Levy_Polynomial", WoodSaxon_Levy_Polynomial, xMinParam, xMaxParam, 12);

    func->SetParameters(par);
    
    return (func->Integral(xMinIntegral, xMaxIntegral))/(func->Integral(1., 8.));
}
//---------------------------------------------------------------------------------------------------
TF1* func_WoodSaxon_PowerLaw_RatioLevy = new TF1("WoodSaxon_PowerLaw_RatioLevy", WoodSaxon_PowerLaw_RatioLevy, xMinParam, xMaxParam, 15);
Double_t Ratio_WoodSaxon_PowerLaw_RatioLevy(Double_t *x, Double_t *par){
	//For Raa : Wood-Saxon = 4 parameters
	//For pp cross section : Power law = 4 parameters
	//For AccEff : ratio of Levy = 7 parameters

    func_WoodSaxon_PowerLaw_RatioLevy->SetParameters(par);

	return (func_WoodSaxon_PowerLaw_RatioLevy->Integral(xMinIntegral, xMaxIntegral))/(func_WoodSaxon_PowerLaw_RatioLevy->Integral(1., 8.));
}
//---------------------------------------------------------------------------------------------------
Double_t Ratio_WoodSaxon_PowerLaw_Polynomial(Double_t *x, Double_t *par){
	//For Raa : Wood-Saxon = 4 parameters
	//For pp cross section : Power law = 4 parameters
	//For AccEff : 3rd degree polynomial = 4 parameters

	TF1* func = new TF1("WoodSaxon_PowerLaw_Polynomial", WoodSaxon_PowerLaw_Polynomial, xMinParam, xMaxParam, 12);

    func->SetParameters(par);

	return (func->Integral(xMinIntegral, xMaxIntegral))/(func->Integral(1., 8.));
}
//---------------------------------------------------------------------------------------------------
TF1* func_Linear_Levy_RatioLevy = new TF1("Linear_Levy_RatioLevy", Linear_Levy_RatioLevy, xMinParam, xMaxParam, 13);
Double_t Ratio_Linear_Levy_RatioLevy(Double_t *x, Double_t *par){
	//For Raa : Pol1 function = 2 parameters
	//For pp cross section : Levy = 4 parameters
	//For AccEff : ratio of Levy = 7 parameters
    
	func_Linear_Levy_RatioLevy->SetParameters(par);

	return (func_Linear_Levy_RatioLevy->Integral(xMinIntegral, xMaxIntegral))/(func_Linear_Levy_RatioLevy->Integral(1., 8.));
}
//---------------------------------------------------------------------------------------------------
Double_t Ratio_Linear_Levy_Polynomial(Double_t *x, Double_t *par){
	//For Raa : Pol1 function = 2 parameters
	//For pp cross section : Levy = 4 parameters
	//For AccEff : 3rd degree polynomial = 4 parameters

	TF1* func = new TF1("Linear_Levy_Polynomial", Linear_Levy_Polynomial, xMinParam, xMaxParam, 10);

    func->SetParameters(par);

	return (func->Integral(xMinIntegral, xMaxIntegral))/(func->Integral(1., 8.));
}
//---------------------------------------------------------------------------------------------------
TF1* func_Linear_PowerLaw_RatioLevy = new TF1("Linear_PowerLaw_RatioLevy", Linear_PowerLaw_RatioLevy, xMinParam, xMaxParam, 13);
Double_t Ratio_Linear_PowerLaw_RatioLevy(Double_t *x, Double_t *par){
	//For Raa : Pol1 function = 2 parameters
	//For pp cross section : Power law = 4 parameters
	//For AccEff : ratio of Levy = 7 parameters

    func_Linear_PowerLaw_RatioLevy->SetParameters(par);
    
    return (func_Linear_PowerLaw_RatioLevy->Integral(xMinIntegral, xMaxIntegral))/(func_Linear_PowerLaw_RatioLevy->Integral(1., 8.));
}
//---------------------------------------------------------------------------------------------------
Double_t Ratio_Linear_PowerLaw_Polynomial(Double_t *x, Double_t *par){
	//For Raa : Pol1 function = 2 parameters
	//For pp cross section : Power law = 4 parameters
	//For AccEff : 3rd degree polynomial = 4 parameters

	TF1* func = new TF1("Linear_PowerLaw_Polynomial", Linear_PowerLaw_Polynomial, xMinParam, xMaxParam, 10);

    func->SetParameters(par);

	return (func->Integral(xMinIntegral, xMaxIntegral))/(func->Integral(1., 8.));
}
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------