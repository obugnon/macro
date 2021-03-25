/*
 *  FitFunctions.C
 *
 *  Created by Ophelie Bugnon on 16/05/19.
 *
 */
#include <vector>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include "TF1.h"


Bool_t reject;

enum Efunction
{
	kVWG = 100,	kVWGQuadratic = 200, kPol2Pol3 = 300, kDoubleExp = 400, kExp = 500,
	kCB = 10, kCBExtended = 20,	kNA60 = 30
};

enum Etails
{
	kEMB = 0,	kPP = 1, kSTARLIGHTcoh = 2, kSTARLIGHTincoh = 3 
};

enum Epart
{
	kJPsi = 0, kPsi2S = 1
};


//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
Double_t VWG(Double_t *x, Double_t *par) // Variable width Gaussian
{
	//par[0]=1/(sigma*sqrt(2*PI)  Normalization
	//par[1]=mu  Mean
	//par[2]=alpha
	//par[3]=beta
	
	Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
 	return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
}


//___________________________________________________________________________________________________________
Double_t VWGQuadratic(Double_t *x, Double_t *par) //Quadratic Variable width Gaussian
{
	//par[0]=1/(sigma*sqrt(2*PI)  Normalization
	//par[1]=mu  Mean
	//par[2]=alpha
	//par[3]=beta
  	//par[4]=gamma
	
	if (reject && x[0] > 2.8 && x[0] < 3.8) {
      TF1::RejectPoint();
      return 0;
   }
	
  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1])+par[4]*((x[0]-par[1])/par[1])*((x[0]-par[1])/par[1]);
 	return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
}

//___________________________________________________________________________________________________________
Double_t Pol2Pol3(Double_t *x, Double_t *par)
{
	//par[0]=Normalization
	//par[1]=a1
	//par[2]=a2
	//par[3]=b1
  	//par[4]=b2
	//par[5]=b3

	if (reject && x[0] > 2.8 && x[0] < 3.8) {
      TF1::RejectPoint();
      return 0;
   }

	return par[0]*(1+par[1]*x[0]+par[2]*x[0]*x[0])/(1+par[3]*x[0]+par[4]*x[0]*x[0]+par[5]*x[0]*x[0]*x[0]);

}

//___________________________________________________________________________________________________________
Double_t DoubleExp(Double_t *x, Double_t *par)
{
	//par[0]=Norm1
	//par[1]=alpha1
	//par[2]=Norm2
	//par[3]=alpha2

	if (reject && x[0] > 2.8 && x[0] < 3.8) {
      TF1::RejectPoint();
      return 0;
   }

   return par[0]*TMath::Exp(par[1]*x[0])+par[2]*TMath::Exp(par[3]*x[0]); 
}

//___________________________________________________________________________________________________________
Double_t Exp(Double_t *x, Double_t *par)
{
	//par[0]=Norm
	//par[1]=alpha


	if (reject && x[0] > 2.8 && x[0] < 3.8) {
      TF1::RejectPoint();
      return 0;
   }

   return par[0]*TMath::Exp(par[1]*x[0]); 
}

//___________________________________________________________________________________________________________
Double_t CrystalBall(Double_t *x, Double_t *par) 
{
	//par[0]=N Normalization
	//par[1]=mu Mean
	//par[2]=sigma Width
	//par[3]=alphaL alpha for left tail
	//par[4]=nL for left tail

	Double_t t = (x[0]-par[1])/par[2];
	Double_t absAlphaL = TMath::Abs(par[3]);
	
	if (par[2] < 0) t = -t;
	
	if (t > -absAlphaL)
	{
		return par[0]*(TMath::Exp(-0.5*t*t));	
	}
	
	if (t <= -absAlphaL)
	{
		Double_t A = TMath::Power(par[4]/absAlphaL,par[4])*TMath::Exp(-0.5*absAlphaL*absAlphaL);
		Double_t B = par[4]/absAlphaL - absAlphaL;
		
		return par[0]*A*TMath::Power(B-t,-par[4]);
	}
	
	return 0. ;
}

//___________________________________________________________________________________________________________
Double_t CrystalBallExtended(Double_t *x, Double_t *par)
{
	//par[0]=N  Normalization
  	//par[1]=mu  Mean
  	//par[2]=sigma  Width
  	//par[3]=alphaL  Alpha for the left tail
	//par[4]=nL for the left tail
	//par[5]=alphaR  Alpha for the right tail
	//par[6]=nR for the right tail


	Double_t t = (x[0]-par[1])/par[2];
	if (par[2] < 0) t = -t;
	
	Double_t absAlphaL = TMath::Abs(par[3]);
  	Double_t absAlphaR = TMath::Abs(par[5]);

  	if (t > -absAlphaL && t < absAlphaR) // gaussian core
  	{
   		return par[0]*(TMath::Exp(-0.5*t*t));
 	}
		
	
	if (t <= -absAlphaL) //left tail
  	{
    		Double_t A =  TMath::Power(par[4]/absAlphaL,par[4])*TMath::Exp(-0.5*absAlphaL*absAlphaL);
    		Double_t B = par[4]/absAlphaL - absAlphaL;
	
   		return par[0]*(A/TMath::Power(B - t, par[4]));
 	}
  
  	
  	
  	if (t >= absAlphaR) //right tail
 	{
 		Double_t C =  TMath::Power(par[6]/absAlphaR,par[6])*TMath::Exp(-0.5*absAlphaR*absAlphaR);
    		Double_t D = par[6]/absAlphaR - absAlphaR;
    		
 		return par[0]*(C/TMath::Power(D + t, par[6]));
  	}
  
  	return 0. ; 
} 

//___________________________________________________________________________________________________________
Double_t DoubleCrystalBallExtended(Double_t *x, Double_t *par)
{
	//par[0]=N  Normalization
  	//par[1]=mu  Mean
  	//par[2]=sigma  Width
  	//par[3]=alphaL  Alpha for the left tail
	//par[4]=nL for the left tail
	//par[5]=alphaR  Alpha for the right tail
	//par[6]=nR for the right tail
	//par[7]=Npsi'
	
	Double_t absAlphaL = TMath::Abs(par[3]);
  	Double_t absAlphaR = TMath::Abs(par[5]);
	Int_t signAlphaL = TMath::Sign(1,par[3]);
  	Int_t signAlphaR = TMath::Sign(1,par[5]); 

	Double_t t1 = signAlphaL*signAlphaR*(x[0]-par[1])/par[2];
	//if (par[2] < 0) t1 = -t1;

	Double_t jpsi;
	Double_t psi2S;
  	if (t1 > -absAlphaL && t1 < absAlphaR) // gaussian core
  	{
   		jpsi = par[0]*(TMath::Exp(-0.5*t1*t1));
 	}
		
	if (t1 <= -absAlphaL) //left tail
  	{
    		Double_t A =  TMath::Power(par[4]/absAlphaL,par[4])*TMath::Exp(-0.5*absAlphaL*absAlphaL);
    		Double_t B = par[4]/absAlphaL - absAlphaL;
	
   		jpsi = par[0]*(A/TMath::Power(B - t1, par[4]));
 	}
  	
  	if (t1 >= absAlphaR) //right tail
 	{
 		Double_t C =  TMath::Power(par[6]/absAlphaR,par[6])*TMath::Exp(-0.5*absAlphaR*absAlphaR);
    		Double_t D = par[6]/absAlphaR - absAlphaR;
    		
 		jpsi = par[0]*(C/TMath::Power(D + t1, par[6]));
  	}
  
	//psi2S
	Double_t t2 = signAlphaL*signAlphaR*(x[0]-(par[1]+3.686109 - 3.096916))/(par[2]*1.05);
	if (t2 > -absAlphaL && t2 < absAlphaR) // gaussian core
  	{
   		psi2S = par[7]*(TMath::Exp(-0.5*t2*t2));
 	}
		
	if (t2 <= -absAlphaL) //left tail
  	{
    		Double_t A =  TMath::Power(par[4]/absAlphaL,par[4])*TMath::Exp(-0.5*absAlphaL*absAlphaL);
    		Double_t B = par[4]/absAlphaL - absAlphaL;
	
   		psi2S = par[7]*(A/TMath::Power(B - t2, par[4]));
 	}
  	
  	if (t2 >= absAlphaR) //right tail
 	{
 		Double_t C =  TMath::Power(par[6]/absAlphaR,par[6])*TMath::Exp(-0.5*absAlphaR*absAlphaR);
    		Double_t D = par[6]/absAlphaR - absAlphaR;
    		
 		psi2S = par[7]*(C/TMath::Power(D + t2, par[6]));
  	}
  	return jpsi+psi2S ; 
} 

//___________________________________________________________________________________________________________
Double_t NA60(Double_t *x, Double_t*par)
{
	//par[0]=N  Normalization
  	//par[1]=mu  Mean
  	//par[2]=sigma  Width
  	//par[3]=alphaL  Alpha for the left tail
	//par[4]=p1 for the left tail
	//par[5]=p2 for the left tail
	//par[6]=p3 for the left tail
	//par[7]=alphaR  Alpha for the right tail
	//par[8]=p1 for the right tail
	//par[9]=p2 for the right tail
	//par[10]=p3 for the right tail

	Double_t t = (x[0]-par[1])/par[2];
	if (par[2] < 0) t = -t;
	Double_t t0;

	if (t < par[3])
	{
		t0 = 1 + TMath::Power(par[4]*(par[3]-t), par[5]-par[6]*TMath::Sqrt(par[3]-t));
	}

	if (t >= par[3] && t < par[7])
	{
		t0 = 1;
	}

	if (t >= par[7])
	{
		t0 = 1 + TMath::Power(par[8]*(t-par[7]), par[9]-par[10]*TMath::Sqrt(t-par[7]));
	}

	return par[0]*TMath::Exp(-0.5*(t/t0)*(t/t0));
}

//___________________________________________________________________________________________________________
Double_t DoubleNA60(Double_t *x, Double_t*par)
{
	//par[0]=N  Normalization
  	//par[1]=mu  Mean
  	//par[2]=sigma  Width
  	//par[3]=alphaL  Alpha for the left tail
	//par[4]=p1 for the left tail
	//par[5]=p2 for the left tail
	//par[6]=p3 for the left tail
	//par[7]=alphaR  Alpha for the right tail
	//par[8]=p1 for the right tail
	//par[9]=p2 for the right tail
	//par[10]=p3 for the right tail
	//par[11]=N psi2s

	Double_t t1 = (x[0]-par[1])/par[2];
	if (par[2] < 0) t1 = -t1;
	Double_t tjpsi;
	Double_t tpsi2S;
	Double_t jpsi;
	Double_t psi2S;

	if (t1 < par[3])
	{
		tjpsi = 1 + TMath::Power(par[4]*(par[3]-t1), par[5]-par[6]*TMath::Sqrt(par[3]-t1));
	}

	if (t1 >= par[3] && t1 < par[7])
	{
		tjpsi = 1;
	}

	if (t1 >= par[7])
	{
		tjpsi = 1 + TMath::Power(par[8]*(t1-par[7]), par[9]-par[10]*TMath::Sqrt(t1-par[7]));
	}

	//psi2S

	Double_t t2 = (x[0]-(par[1]+3.686109 - 3.096916))/(par[2]*1.05);
	if (par[2] < 0) t2 = -t2;

		if (t2 < par[3])
	{
		tpsi2S = 1 + TMath::Power(par[4]*(par[3]-t2), par[5]-par[6]*TMath::Sqrt(par[3]-t2));
	}

	if (t2 >= par[3] && t2 < par[7])
	{
		tpsi2S = 1;
	}

	if (t2 >= par[7])
	{
		tpsi2S = 1 + TMath::Power(par[8]*(t2-par[7]), par[9]-par[10]*TMath::Sqrt(t2-par[7]));
	}

	return (par[0]*TMath::Exp(-0.5*(t1/tjpsi)*(t1/tjpsi)))+(par[11]*TMath::Exp(-0.5*(t2/tpsi2S)*(t2/tpsi2S)));
}

//___________________________________________________________________________________________________________
Double_t VWG_CBext(Double_t *x, Double_t *par)//4 parmeters for BG and 7 parameters for signal
{
	return VWG(x,par)+CrystalBallExtended(x,&(par[4]));
}

//___________________________________________________________________________________________________________
Double_t VWGquad_CBext(Double_t *x, Double_t *par)//5 parmeters for BG and 7 parameters for signal
{
	return VWGQuadratic(x,par)+CrystalBallExtended(x,&(par[5]));
}

//___________________________________________________________________________________________________________
Double_t VWGquad_DoubleCBext(Double_t *x, Double_t *par)//5 parmeters for BG and 7+1 parameters for signal
{
  	//BackGround
  	//par[0]=1/(sigma*sqrt(2*PI)  Normalization
	//par[1]=mu  Mean
	//par[2]=alpha
	//par[3]=beta
  	//par[4]=gamma

	//J/psi signal
  	//par[5]=N  Normalization
  	//par[6]=mu  Mean
  	//par[7]=sigma  Width
  	//par[8]=alphaL  Alpha for the left tail
	//par[9]=nL  n for the left tail
	//par[10]=alphaR  Alpha for the right tail
	//par[11]=nR  n for the right tail

  	//psi' signal
  	//par[12]=N  Normalization

	return VWGQuadratic(x,par)+DoubleCrystalBallExtended(x,&(par[5]));
}

//___________________________________________________________________________________________________________
Double_t VWGquad_DoubleNA60(Double_t *x, Double_t *par)//5 parmeters for BG and 11+1 parameters for signal
{
  	//BackGround
  	//par[0]=1/(sigma*sqrt(2*PI)  Normalization
	//par[1]=mu  Mean
	//par[2]=alpha
	//par[3]=beta
  	//par[4]=gamma

	//J/psi signal
  	//par[5]=N  Normalization
  	//par[6]=mu  Mean
  	//par[7]=sigma  Width
  	//par[8]=alphaL  Alpha for the left tail
	//par[9]=p1 for the left tail
	//par[10]=p2 for the left tail
	//par[11]=p3 for the left tail
	//par[12]=alphaR  Alpha for the right tail
	//par[13]=p1 for the right tail
	//par[14]=p2 for the right tail
	//par[15]=p3 for the right tail

  	//psi' signal
  	//par[16]=N  Normalization
  	
	return VWGQuadratic(x,par)+DoubleNA60(x,&(par[5]));
}

//___________________________________________________________________________________________________________
Double_t Pol2Pol3_DoubleCBext(Double_t *x, Double_t *par)//6 parmeters for BG and 7+1 parameters for signal
{
  	//BackGround
  	//par[0]=Normalization
	//par[1]=a1
	//par[2]=a2
	//par[3]=b1
  	//par[4]=b2
	//par[5]=b3

	//J/psi signal
  	//par[6]=N  Normalization
  	//par[7]=mu  Mean
  	//par[8]=sigma  Width
  	//par[9]=alphaL  Alpha for the left tail
	//par[10]=nL  n for the left tail
	//par[11]=alphaR  Alpha for the right tail
	//par[12]=nR  n for the right tail

  	//psi' signal
  	//par[13]=N  Normalization

	return Pol2Pol3(x,par)+DoubleCrystalBallExtended(x,&(par[6]));
}


//___________________________________________________________________________________________________________
Double_t Pol2Pol3_DoubleNA60(Double_t *x, Double_t *par)//6 parmeters for BG and 11+1 parameters for signal
{
  	//BackGround
  	//par[0]=Normalization
	//par[1]=a1
	//par[2]=a2
	//par[3]=b1
  	//par[4]=b2
	//par[5]=b3

	//J/psi signal
  	//par[6]=N  Normalization
  	//par[7]=mu  Mean
  	//par[8]=sigma  Width
  	//par[9]=alphaL  Alpha for the left tail
	//par[10]=p1 for the left tail
	//par[11]=p2 for the left tail
	//par[12]=p3 for the left tail
	//par[13]=alphaR  Alpha for the right tail
	//par[14]=p1 for the right tail
	//par[15]=p2 for the right tail
	//par[16]=p3 for the right tail

  	//psi' signal
  	//par[17]=N  Normalization

	return Pol2Pol3(x,par)+DoubleNA60(x,&(par[6]));
}

//___________________________________________________________________________________________________________
Double_t DoubleExp_DoubleCBext(Double_t *x, Double_t *par)//4 parameters for BG ans 7+1 for signals
{
	//BackGround
	//par[0]=Norm1
	//par[1]=alpha1
	//par[2]=Norm2
	//par[3]=alpha2

	//J/psi signal
  	//par[4]=N  Normalization
  	//par[5]=mu  Mean
  	//par[6]=sigma  Width
  	//par[7]=alphaL  Alpha for the left tail
	//par[8]=nL  n for the left tail
	//par[9]=alphaR  Alpha for the right tail
	//par[10]=nR  n for the right tail

  	//psi' signal
  	//par[11]=N  Normalization

	return DoubleExp(x,par)+DoubleCrystalBallExtended(x,&(par[4]));

}

//___________________________________________________________________________________________________________
Double_t DoubleExp_DoubleNA60(Double_t *x, Double_t *par)//4 parameters for BG ans 11+1 for signals
{
	//BackGround
	//par[0]=Norm1
	//par[1]=alpha1
	//par[2]=Norm2
	//par[3]=alpha2

	//J/psi signal
  	//par[4]=N  Normalization
  	//par[5]=mu  Mean
  	//par[6]=sigma  Width
  	//par[7]=alphaL  Alpha for the left tail
	//par[8]=p1 for the left tail
	//par[9]=p2 for the left tail
	//par[10]=p3 for the left tail
	//par[11]=alphaR  Alpha for the right tail
	//par[12]=p1 for the right tail
	//par[13]=p2 for the right tail
	//par[14]=p3 for the right tail

  	//psi' signal
  	//par[15]=N  Normalization

	return DoubleExp(x,par)+DoubleNA60(x,&(par[4]));

}
//___________________________________________________________________________________________________________
Double_t Exp_DoubleCBext(Double_t *x, Double_t *par)//4 parameters for BG ans 7+1 for signals
{
	//BackGround
	//par[0]=Norm
	//par[1]=alpha

	//J/psi signal
  	//par[2]=N  Normalization
  	//par[3]=mu  Mean
  	//par[4]=sigma  Width
  	//par[5]=alphaL  Alpha for the left tail
	//par[6]=nL  n for the left tail
	//par[7]=alphaR  Alpha for the right tail
	//par[8]=nR  n for the right tail

  	//psi' signal
  	//par[9]=N  Normalization

	return DoubleExp(x,par)+DoubleCrystalBallExtended(x,&(par[2]));

}

//___________________________________________________________________________________________________________
Double_t Exp_DoubleNA60(Double_t *x, Double_t *par)//4 parameters for BG ans 11+1 for signals
{
	//BackGround
	//par[0]=Norm
	//par[1]=alpha

	//J/psi signal
  	//par[2]=N  Normalization
  	//par[3]=mu  Mean
  	//par[4]=sigma  Width
  	//par[5]=alphaL  Alpha for the left tail
	//par[6]=p1 for the left tail
	//par[7]=p2 for the left tail
	//par[8]=p3 for the left tail
	//par[9]=alphaR  Alpha for the right tail
	//par[10]=p1 for the right tail
	//par[11]=p2 for the right tail
	//par[12]=p3 for the right tail

  	//psi' signal
  	//par[13]=N  Normalization

	return DoubleExp(x,par)+DoubleNA60(x,&(par[2]));

}
//___________________________________________________________________________________________________________
Int_t GetNPar(Efunction fName)
{
	switch (fName)
	{
		case kVWG:
		case kDoubleExp:
			return 4;
		break;

		case kExp:
			return 2;
		break;
	
		case kVWGQuadratic:
		case kCB:
			return 5;
		break;

		case kPol2Pol3:
			return 6;
		break;
	
		case kCBExtended:
			return 7;
		break;

		case kNA60:
			return 11;
		break;	
	}
}

