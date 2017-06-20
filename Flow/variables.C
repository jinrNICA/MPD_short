#include "variables.h"

Double_t ResEventPlane(Double_t chi, Int_t harm) //harm = 1 or 2;
{
  Double_t con = TMath::Sqrt(TMath::Pi()/2)/2 ;
  Double_t arg = chi * chi / 4.;
  if (harm == 1) { Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg)); return res; }
  else if (harm == 2) { Double_t res = con * chi * exp(-arg) * (ROOT::Math::cyl_bessel_i(0.5,arg) + ROOT::Math::cyl_bessel_i(1.5,arg)); return res; }
  else cout << "Wrong harmonic in resolution extrapolation! " <<endl;
}

Double_t ResEventPlane(Double_t chi, Int_t k, Int_t m) //harm = 1 or 2;
{
  Double_t con = TMath::Sqrt(TMath::Pi()/2)/2 ;
  Double_t arg = chi * chi / 4.;
  Double_t res = con * chi * exp(-arg) * (ROOT::Math::cyl_bessel_i((k-1)/2,arg) + ROOT::Math::cyl_bessel_i((k+1)/2,arg)); 
  return res;
}

Double_t Chi(Double_t res, Int_t harm, Int_t accuracy) //harm = 1 or 2; 
{

  Double_t chi   = 2.0;
  Double_t delta = 1.0;
  for(int i = 0; i < accuracy; i++)
  {
   if(ResEventPlane(chi, harm) < res) { chi = chi + delta ; }
   else                               { chi = chi - delta ; }
   delta = delta / 2.;
  }

  return chi ;
}

Double_t Chi(Double_t res, Int_t k, Int_t m, Int_t accuracy) //harm = 1 or 2; 
{

  Double_t chi   = 2.0;
  Double_t delta = 1.0;
  for(int i = 0; i < accuracy; i++)
  {
   if(ResEventPlane(chi, k, m) < res) { chi = chi + delta ; }
   else                               { chi = chi - delta ; }
   delta = delta / 2.;
  }

  return chi ;
}

//double cyl_bessel_i(double nu, double x) { return gsl_sf_bessel_Inu(nu, x); }