#define _MAX_TRACKS 5000
#define _MAX_TRACKS_MC 200000
#define _N2_ZDC 225
#define _N_ZDC 15
#define _N_MODULES_TOTAL 90
#define _N_ARM 2
#define _N_HARM 2
#define _N_METHOD 2
#define _N_QCOMP 2

#include <iostream>

#include <TMath.h>
#include <TSystem.h>
//#include "TRoot.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
//#include "FairMCEventHeader.h"
//#include "MCHeader.h"
//#include "MpdEvent.h"
//#include "MPDEvent.h"
//#include "MpdZdcDigi.h"
//#include "ZDCHit.h"
//#include "TClonesArray.h"
//#include "MpdTrack.h"
//#include "FairMCTrack.h"
 

using std::cout;
using std::endl;

using TMath::ATan2;
using TMath::Sqrt;

void create_reduced_tree(TString inFiles,TString outFiles)
{
	gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
	mpdloadlibs(kTRUE,kTRUE); 
	
	TFile *inFile = new TFile(inFiles.Data(),"READ");
	TTree *inTree = (TTree*) inFile->Get("cbmsim");
	
	TFile  *outFile = new TFile(outFiles.Data(),"RECREATE");
	TTree  *outTree = new TTree("cbmsim_reduced","cbmsim_reduced");
	
	long int mc_side[_MAX_TRACKS_MC][10];
	long int mpd_side[_MAX_TRACKS];

	Float_t b_mc;
    	outTree->Branch("b_mc",&b_mc,"b_mc/F");
    	Float_t phiEP_mc;
    	outTree->Branch("phiEP_mc",&phiEP_mc,"phiEP_mc/F");
    	Float_t x_vertex_mc;
    	outTree->Branch("x_vertex_mc",&x_vertex_mc,"x_vertex_mc/F");
    	Float_t y_vertex_mc;
    	outTree->Branch("y_vertex_mc",&y_vertex_mc,"y_vertex_mc/F");
    	Float_t z_vertex_mc;
    	outTree->Branch("z_vertex_mc",&z_vertex_mc,"z_vertex_mc/F");

	Long_t n_tracks_mc;
    	outTree->Branch("n_tracks_mc",&n_tracks_mc,"n_tracks_mc/L"); 
    	Float_t eta_mc[_MAX_TRACKS];
    	outTree->Branch("eta_mc",eta_mc,"eta_mc[n_tracks_mc]/F");
    	Float_t pt_mc[_MAX_TRACKS];
    	outTree->Branch("pt_mc",pt_mc,"pt_mc[n_tracks_mc]/F");  
    	Int_t mother_ID_mc[_MAX_TRACKS];
    	outTree->Branch("mother_ID_mc",mother_ID_mc,"mother_ID_mc[n_tracks_mc]/I");
    	Int_t PDG_code_mc[_MAX_TRACKS];
    	outTree->Branch("PDG_code_mc",PDG_code_mc,"PDG_code_mc[n_tracks_mc]/I");
    	Float_t px_mc[_MAX_TRACKS];
    	outTree->Branch("px_mc",px_mc,"px_mc[n_tracks_mc]");
    	Float_t py_mc[_MAX_TRACKS];
    	outTree->Branch("py_mc",py_mc,"py_mc[n_tracks_mc]");
    	Float_t pz_mc[_MAX_TRACKS];
    	outTree->Branch("pz_mc",pz_mc,"pz_mc[n_tracks_mc]");
    	Float_t start_x_mc[_MAX_TRACKS];
    	outTree->Branch("start_x_mc",start_x_mc,"start_x_mc[n_tracks_mc]");
   	Float_t start_y_mc[_MAX_TRACKS];
    	outTree->Branch("start_y_mc",start_y_mc,"start_y_mc[n_tracks_mc]");
    	Float_t start_z_mc[_MAX_TRACKS];
    	outTree->Branch("start_z_mc",start_z_mc,"start_z_mc[n_tracks_mc]");
    	Float_t mass_mc[_MAX_TRACKS];
    	outTree->Branch("mass_mc",mass_mc,"mass_mc[n_tracks_mc]");
    	Float_t energy_mc[_MAX_TRACKS];
    	outTree->Branch("energy_mc",energy_mc,"energy_mc[n_tracks_mc]");
    
    	Long_t n_tracks_mpd;
    	outTree->Branch("n_tracks_mpd",&n_tracks_mpd,"n_tracks_mpd/L");
    	Float_t eta_mpd[_MAX_TRACKS];
    	outTree->Branch("eta_mpd",eta_mpd,"eta_mpd[n_tracks_mpd]/F");
    	Float_t phi_mpd[_MAX_TRACKS];
    	outTree->Branch("phi_mpd",phi_mpd,"phi_mpd[n_tracks_mpd]/F");
    	Float_t theta_mpd[_MAX_TRACKS];
    	outTree->Branch("theta_mpd",theta_mpd,"theta_mpd[n_tracks_mpd]/F");
    	Int_t TOF_flag_mpd[_MAX_TRACKS];
    	outTree->Branch("TOF_flag_mpd",TOF_flag_mpd,"TOF_flag_mpd[n_tracks_mpd]/I");
    	Float_t ZDC_energy_mpd[_N_MODULES_TOTAL];
    	outTree->Branch("ZDC_energy_mpd",ZDC_energy_mpd,"ZDC_energy_mpd[90]/F"); /////////////////////////////////
    	Float_t pid_tpc_prob_electron_mpd[_MAX_TRACKS];
    	outTree->Branch("pid_tpc_prob_electron_mpd",pid_tpc_prob_electron_mpd,"pid_tpc_prob_electron_mpd[n_tracks_mpd]/F");
    	Float_t pid_tpc_prob_pion_mpd[_MAX_TRACKS];
    	outTree->Branch("pid_tpc_prob_pion_mpd",pid_tpc_prob_pion_mpd,"pid_tpc_prob_pion_mpd[n_tracks_mpd]/F");
    	Float_t pid_tpc_prob_kaon_mpd[_MAX_TRACKS];
    	outTree->Branch("pid_tpc_prob_kaon_mpd",pid_tpc_prob_kaon_mpd,"pid_tpc_prob_kaon_mpd[n_tracks_mpd]/F");
    	Float_t pid_tpc_prob_proton_mpd[_MAX_TRACKS];
    	outTree->Branch("pid_tpc_prob_proton_mpd",pid_tpc_prob_proton_mpd,"pid_tpc_prob_proton_mpd[n_tracks_mpd]/F");
    	Float_t pid_tof_prob_electron_mpd[_MAX_TRACKS];
    	outTree->Branch("pid_tof_prob_electron_mpd",pid_tof_prob_electron_mpd,"pid_tof_prob_electron_mpd[n_tracks_mpd]/F");
    	Float_t pid_tof_prob_pion_mpd[_MAX_TRACKS];
    	outTree->Branch("pid_tof_prob_pion_mpd",pid_tof_prob_pion_mpd,"pid_tof_prob_pion_mpd[n_tracks_mpd]/F");
    	Float_t pid_tof_prob_kaon_mpd[_MAX_TRACKS];
    	outTree->Branch("pid_tof_prob_kaon_mpd",pid_tof_prob_kaon_mpd,"pid_tof_prob_kaon_mpd[n_tracks_mpd]/F");
    	Float_t pid_tof_prob_proton_mpd[_MAX_TRACKS];
    	outTree->Branch("pid_tof_prob_proton_mpd",pid_tof_prob_proton_mpd,"pid_tof_prob_proton_mpd[n_tracks_mpd]/F");
    	Float_t tof_beta_mpd[_MAX_TRACKS];
    	outTree->Branch("tof_beta_mpd",tof_beta_mpd,"tof_beta_mpd[n_tracks_mpd]/F");
    	Float_t tof_mass2_mpd[_MAX_TRACKS];
    	outTree->Branch("tof_mass2_mpd",tof_mass2_mpd,"tof_mass2_mpd[n_tracks_mpd]/F");
    	Float_t dEdx_tpc_mpd[_MAX_TRACKS];
    	outTree->Branch("dEdx_tpc_mpd",dEdx_tpc_mpd,"dEdx_tpc_mpd[n_tracks_mpd]/F");
    	Float_t chi2_mpd[_MAX_TRACKS];
    outTree->Branch("chi2_mpd",chi2_mpd,"chi2_mpd[n_tracks_mpd]/F");
    Float_t pt_error_mpd[_MAX_TRACKS];
    outTree->Branch("pt_error_mpd",pt_error_mpd,"pt_error_mpd[n_tracks_mpd]/F");
    Float_t theta_error_mpd[_MAX_TRACKS];
    outTree->Branch("theta_error_mpd",theta_error_mpd,"theta_error_mpd[n_tracks_mpd]/F");
    Float_t phi_error_mpd[_MAX_TRACKS];
    outTree->Branch("phi_error_mpd",phi_error_mpd,"phi_error_mpd[n_tracks_mpd]/F");
    Float_t DCA_x_mpd[_MAX_TRACKS];
    outTree->Branch("DCA_x_mpd",DCA_x_mpd,"DCA_x_mpd[n_tracks_mpd]/F");
    Float_t DCA_y_mpd[_MAX_TRACKS];
    outTree->Branch("DCA_y_mpd",DCA_y_mpd,"DCA_y_mpd[n_tracks_mpd]/F");
    Float_t DCA_z_mpd[_MAX_TRACKS];
    outTree->Branch("DCA_z_mpd",DCA_z_mpd,"DCA_z_mpd[n_tracks_mpd]/F");
    Float_t DCA_global_x_mpd[_MAX_TRACKS];
    outTree->Branch("DCA_global_x_mpd",DCA_global_x_mpd,"DCA_global_x_mpd[n_tracks_mpd]/F");
    Float_t DCA_global_y_mpd[_MAX_TRACKS];
    outTree->Branch("DCA_global_y_mpd",DCA_global_y_mpd,"DCA_global_y_mpd[n_tracks_mpd]/F");
    Float_t DCA_global_z_mpd[_MAX_TRACKS];
    outTree->Branch("DCA_global_z_mpd",DCA_global_z_mpd,"DCA_global_z_mpd[n_tracks_mpd]/F");
    Int_t n_hits_mpd[_MAX_TRACKS];
    outTree->Branch("n_hits_mpd",n_hits_mpd,"n_hits_mpd[n_tracks_mpd]/I");
    Int_t n_hits_poss_mpd[_MAX_TRACKS];
    outTree->Branch("n_hits_poss_mpd",n_hits_poss_mpd,"n_hits_poss_mpd[n_tracks_mpd]/I");
    Float_t signed_pt_mpd[_MAX_TRACKS];
    outTree->Branch("signed_pt_mpd",signed_pt_mpd,"signed_pt_mpd[n_tracks_mpd]/F");
    Float_t p_mpd[_MAX_TRACKS];
    outTree->Branch("p_mpd",p_mpd,"p_mpd[n_tracks_mpd]/F");

    outTree->Branch("id_from_mc_mpd",mpd_side,"id_from_mc_mpd[n_tracks_mpd]/L");
	
    FairMCEventHeader *MCHeader=0;
    inTree->SetBranchAddress("MCEventHeader.", &MCHeader);
    TClonesArray *MCTracks=0;
    inTree->SetBranchAddress("MCTrack", &MCTracks);
    MpdEvent *MPDEvent=0;
    inTree->SetBranchAddress("MPDEvent.", &MPDEvent);
    TClonesArray *ZDCHits = 0;
    inTree->SetBranchAddress("ZdcDigi",&ZDCHits);
    TClonesArray *MpdGlobalTracks=0;
    MpdZdcDigi* ZDCHit = 0;
    
    Int_t n_entries = inTree->GetEntries();
    //Int_t n_entries = 100;
    for (int i = 0; i < n_entries; ++i)
    {

	for (long int j = 0; j < _MAX_TRACKS_MC; ++j)
        {
		for (int k = 0; k < 10; ++k)
		{
			mc_side[j][k] = -1;
		}
	}
	for (long int j = 0; j < _MAX_TRACKS; ++j)
	{
		mpd_side[j] = -1;
	}
		
	for (long int j = 0; j < 90; ++j)
        {
		ZDC_energy_mpd[j] = 0;
	}
	    
	cout << "EVENT N "<< i <<endl;
	inTree->GetEntry(i);
	MpdGlobalTracks = MPDEvent->GetGlobalTracks();
	    
	Int_t number_of_zdchits = ZDCHits->GetEntries();
	for (Int_t zdchit_number = 0; zdchit_number < number_of_zdchits; ++zdchit_number)
	{
		ZDCHit = (MpdZdcDigi*) ZDCHits->At(zdchit_number);
		Int_t detector_ID = ZDCHit->GetDetectorID();//1,2
	        Int_t module_ID = ZDCHit->GetModuleID();//1-45
		Double_t energy_deposit_per_hit = ZDCHit->GetELoss();
			
		ZDC_energy_mpd[ (detector_ID - 1)*45 + module_ID - 1] += energy_deposit_per_hit;
	}
        
        b_mc = MCHeader->GetB();
        phiEP_mc = MCHeader->GetEP();
	//phi_event_mc = MCHeader->GetPhi();
	x_vertex_mc = MCHeader->GetX();
	y_vertex_mc = MCHeader->GetY();
	z_vertex_mc = MCHeader->GetZ();	
	    
	n_tracks_mpd = MpdGlobalTracks->GetEntriesFast();
	long int m = 0;
	for (long int j = 0; j < n_tracks_mpd; ++j)
	{ 
			MpdTrack *mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(j); 
			FairMCTrack *mctrack = (FairMCTrack*)MCTracks->At(mpdtrack->GetID());
			if (mctrack->GetMotherId() > -1)  continue; ///////////////////////////////////////Mother ID cut is here
			for (int k = 0; k < 10; ++k) 
			{
				if ((mc_side[mpdtrack->GetID()][k]) == -1) { mc_side[mpdtrack->GetID()][k] = m; break;} 
		    	}
			
			eta_mpd[m] = mpdtrack->GetEta();
			n_hits_mpd[m] = mpdtrack->GetNofHits();
			n_hits_poss_mpd[m] = mpdtrack->GetNofHitsPossTpc();
			pid_tpc_prob_electron_mpd[m] = mpdtrack->GetTPCPidProbElectron();
			pid_tpc_prob_pion_mpd[m] = mpdtrack->GetTPCPidProbPion();
			pid_tpc_prob_kaon_mpd[m] = mpdtrack->GetTPCPidProbKaon();
			pid_tpc_prob_proton_mpd[m] = mpdtrack->GetTPCPidProbProton();
			pid_tof_prob_electron_mpd[m] = mpdtrack->GetTOFPidProbElectron();
			pid_tof_prob_pion_mpd[m] = mpdtrack->GetTOFPidProbPion();
			pid_tof_prob_kaon_mpd[m] = mpdtrack->GetTOFPidProbKaon();
			pid_tof_prob_proton_mpd[m] = mpdtrack->GetTOFPidProbProton();
			tof_beta_mpd[m] = mpdtrack->GetTofBeta();
			tof_mass2_mpd[m] = mpdtrack->GetTofMass2();
			dEdx_tpc_mpd[m] = mpdtrack->GetdEdXTPC();
			chi2_mpd[m] = mpdtrack->GetChi2();
			pt_error_mpd[m] = mpdtrack->GetPtError();
			theta_error_mpd[m] = mpdtrack->GetThetaError();
			phi_error_mpd[m] = mpdtrack->GetPhiError();
			DCA_x_mpd[m] = mpdtrack->GetDCAX();
			DCA_y_mpd[m] = mpdtrack->GetDCAY();
			DCA_z_mpd[m] = mpdtrack->GetDCAZ();
			DCA_global_x_mpd[m] = mpdtrack->GetDCAGlobalX();
			DCA_global_y_mpd[m] = mpdtrack->GetDCAGlobalY();
			DCA_global_z_mpd[m] = mpdtrack->GetDCAGlobalZ(); 
			signed_pt_mpd[m] = mpdtrack->GetPt();
			p_mpd[m] = TMath::Sqrt( mpdtrack->GetPx()*mpdtrack->GetPx() + mpdtrack->GetPy()*mpdtrack->GetPy() + mpdtrack->GetPz()*mpdtrack->GetPz());
			phi_mpd[m] = mpdtrack->GetPhi();
            		theta_mpd[m] = mpdtrack->GetTheta();
            		TOF_flag_mpd[m] = mpdtrack->GetTofFlag();	
			
			++m;
	}
	n_tracks_mpd = m;
		
	n_tracks_mc = MCTracks->GetEntries();
	m = 0;
	for(long int j = 0; j < n_tracks_mc; ++j)
	{
		FairMCTrack *mctrack = (FairMCTrack*)MCTracks->At(j);
		if (mctrack->GetMotherId() > -1) continue; /////////////////////////Mother ID cut is here
		for (int k = 0; k < 10; ++k) 
		{
			if ( mc_side[j][k] != -1) mpd_side[mc_side[j][k]] = m;
		}
			
		Float_t theta_mc = TMath::ATan2(mctrack->GetPt(),mctrack->GetPz());
		eta_mc[m] = -TMath::Log(TMath::Tan(theta_mc/2.));
		pt_mc[m] = mctrack->GetPt();
		mother_ID_mc[m] = mctrack->GetMotherId();
            	PDG_code_mc[m] = mctrack->GetPdgCode();
            	px_mc[m] = mctrack->GetPx();
            	py_mc[m] = mctrack->GetPy();
            	pz_mc[m] = mctrack->GetPz();
            	start_x_mc[m] = mctrack->GetStartX();
            	start_y_mc[m] = mctrack->GetStartY();
            	start_z_mc[m] = mctrack->GetStartZ();
            	mass_mc[m] = mctrack->GetMass();
            	energy_mc[m] = mctrack->GetEnergy();   
			
		++m;
	}
	n_tracks_mc = m;
		
    	outTree->Fill();
    }
    
    outFile->cd();
    outTree->Write();
    outFile->Close();
}
