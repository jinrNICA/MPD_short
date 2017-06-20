#define compile

#define _N_ARM 2
#define _N_HARM 2
#define _N_SORTS 3
#define _N_CENTRALITY_BINS 2
#define _N_METHOD 4

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#ifdef compile
//This is for compilation of this file
#include <iostream>
#include <TMath.h>
#include <TProfile.h>
#include <TVector3.h>
//#include <TROOT.h>
//#include <./Math.h>
//#include "./gsl_sf_bessel.h"
#endif


using std::cout;
using std::endl;
using TMath::Sin;
using TMath::Cos;
using TMath::ATan2;
using TMath::Abs;
using TMath::LocMin;
using TMath::Pi;
using TMath::Abs;
using TMath::Sqrt;

const Float_t Cut_Pt_Min = 0.2;
const Float_t Cut_Eta_Min = 0.7;
const Float_t Cut_Eta_Max = 1.5;
const Int_t Cut_No_Of_hits_min = 15;

const float b_bins_range[] = {2.,7.,12.};
const int NimpactBins = 2;

//const float centralityBinsFlow[] = {20.,60.,100.};
const float centralityBinsFlow[] = {2.,7.,12.};
const int NcentralityBinsFlow = 2;

const int multiplicityBinsFlow[] = {18.,121.,10000.};
const int NmultiplicityBinsFlow = _N_CENTRALITY_BINS;

//const float centralityBinsRes[] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100};
const float centralityBinsRes[] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.};
const int NcentralityBinsRes = 16;

const int multiplicityBinsRes[] = {0.,8.,18.,38.,71.,121.,194.,293.,428.,614.,10000.};
const int NmultiplicityBinsRes = 10;

const float rapidityBins[] = {-2., -1.8, -1.6, -1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.};
const int NrapidityBins = 20;

const float etaBins[] = {-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,3.0,3.2};
const int NetaBins = 32;

const float ptbin[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.7,2.0,2.3,2.6,3.0};
const int   Nptbins = 18;

const TString arm_names[_N_ARM] = {TString("L"),TString("R")};
const TString side_names[_N_ARM] = {TString("Left"),TString("Right")};
const TString proj_name[2]      = {TString("x"),TString("y")};
const TString sorts_of_particles[_N_SORTS+1] = {TString("protons (2212)") , TString("kaons (321)"), TString("pions (211)"), TString("all sorts")};
const TString methods_names[_N_METHOD] = {TString("TPC |#eta|>0."), TString("TPC |#eta|>0.2"), TString("TPC |#eta|>0.5"), TString("FHCal")};
const TString eta_method[_N_METHOD-1] = {TString("|#eta|>0."), TString("|#eta|>0.2"), TString("|#eta|>0.5")};
const TString etaGap_name[] = {TString("#eta-gap=0"),TString("#eta-gap=0.4"),TString("#eta-gap=1.0")};

const Double_t Res1Fit[] = {0.510239,0.658662,0.744103,0.810397,0.848877,0.865456,0.864284,0.845521,0.788119,0.664607};
const Double_t Res1Flow[2] = {0.768456,0.807093};
const Double_t Res2Flow[2] = {0.37533,0.523963};

const Double_t EtaGap[] = {0.,0.2,0.5};
const Int_t    NetaGap  = _N_METHOD-1;
