///////////////////////////////////////////////////////////////////////////////////////////////////
#define real_flow_cxx
#include "real_flow.h"

//#define debug
//#define make_tree

//#define check

void real_flow::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L real_flow.C
//      Root > real_flow t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   InitHisto();
   	#ifndef debug
   		nentries = fChain->GetEntriesFast();
   		//nentries = 250000;
   		//nentries = 10000;
   	#elif defined debug
   		nentries = 10000;
   		cout << nentries << " Events:" << endl;
   	#endif

   Long64_t nbytes = 0, nb = 0;
   Double_t Q[2][_N_ARM][_N_HARM][_N_METHOD][2]; //Q[weight=pT,1][R,L][harm][TPC,FHCal][x,y];
   Double_t QxTree[_N_HARM][_N_ARM][_N_METHOD];
   Double_t QyTree[_N_HARM][_N_ARM][_N_METHOD];
   Int_t    QwTree[_N_HARM][_N_ARM][_N_METHOD];
   for (Int_t harm=0;harm<_N_HARM;harm++){
   	for(Int_t side=0;side<_N_ARM;side++){
   		for (Int_t method=0;method<_N_METHOD;method++){
   			QxTree[harm][side][method] = 0.;
			QyTree[harm][side][method] = 0.;
			QwTree[harm][side][method] = 0;
   		}
   	}
   }
   #ifdef make_tree
   out_file = new TFile(Form("tree_%s",outFile.Data()),"recreate");
   out_tree = new TTree(Form("tree_EP"),"tree_EP");
   Double_t pt_tree[_MAX_TRACKS];

   	out_tree->Branch("b_mc", &b_mc, "b_mc/F");
	out_tree->Branch("phiEP_mc", &phiEP_mc, "phiEP_mc/F");
	out_tree->Branch("n_tracks_mpd", &n_tracks_mpd,"n_tracks_mpd/I");
	out_tree->Branch("eta_mpd", eta_mpd, "eta_mpd[n_tracks_mpd]/F");
	out_tree->Branch("phi_mpd",phi_mpd,"phi_mpd[n_tracks_mpd]/F");
	out_tree->Branch("pt_tree",pt_tree,"pt_tree[n_tracks_mpd]/D");
	out_tree->Branch("Qx_1_Left_FHCal",&QxTree[0][0][0],"Qx_1_Left_FHCal/D");
	out_tree->Branch("Qx_1_Left_TPC_eta_0",&QxTree[0][0][1],"Qx_1_Left_TPC_eta_0/D");
	out_tree->Branch("Qx_1_Left_TPC_eta_0_2",&QxTree[0][0][2],"Qx_1_Left_TPC_eta_0_2/D");
	out_tree->Branch("Qx_1_Left_TPC_eta_0_5",&QxTree[0][0][3],"Qx_1_Left_TPC_eta_0_5/D");
	out_tree->Branch("Qx_1_Right_FHCal",&QxTree[0][1][0],"Qx_1_Right_FHCal/D");
	out_tree->Branch("Qx_1_Right_TPC_eta_0",&QxTree[0][1][1],"Qx_1_Right_TPC_eta_0/D");
	out_tree->Branch("Qx_1_Right_TPC_eta_0_2",&QxTree[0][1][2],"Qx_1_Right_TPC_eta_0_2/D");
	out_tree->Branch("Qx_1_Right_TPC_eta_0_5",&QxTree[0][1][3],"Qx_1_Right_TPC_eta_0_5/D");
	out_tree->Branch("Qx_2_Left_FHCal",&QxTree[1][0][0],"Qx_2_Left_FHCal/D");
	out_tree->Branch("Qx_2_Left_TPC_eta_0",&QxTree[1][0][1],"Qx_2_Left_TPC_eta_0/D");
	out_tree->Branch("Qx_2_Left_TPC_eta_0_2",&QxTree[1][0][2],"Qx_2_Left_TPC_eta_0_2/D");
	out_tree->Branch("Qx_2_Left_TPC_eta_0_5",&QxTree[1][0][3],"Qx_2_Left_TPC_eta_0_5/D");
	out_tree->Branch("Qx_2_Right_FHCal",&QxTree[1][1][0],"Qx_2_Right_FHCal/D");
	out_tree->Branch("Qx_2_Right_TPC_eta_0",&QxTree[1][1][1],"Qx_2_Right_TPC_eta_0/D");
	out_tree->Branch("Qx_2_Right_TPC_eta_0_2",&QxTree[1][1][2],"Qx_2_Right_TPC_eta_0_2/D");
	out_tree->Branch("Qx_2_Right_TPC_eta_0_5",&QxTree[1][1][3],"Qx_2_Right_TPC_eta_0_5/D");
	out_tree->Branch("Qy_1_Left_FHCal",&QyTree[0][0][0],"Qy_1_Left_FHCal/D");
	out_tree->Branch("Qy_1_Left_TPC_eta_0",&QyTree[0][0][1],"Qy_1_Left_TPC_eta_0/D");
	out_tree->Branch("Qy_1_Left_TPC_eta_0_2",&QyTree[0][0][2],"Qy_1_Left_TPC_eta_0_2/D");
	out_tree->Branch("Qy_1_Left_TPC_eta_0_5",&QyTree[0][0][3],"Qy_1_Left_TPC_eta_0_5/D");
	out_tree->Branch("Qy_1_Right_FHCal",&QyTree[0][1][0],"Qy_1_Right_FHCal/D");
	out_tree->Branch("Qy_1_Right_TPC_eta_0",&QyTree[0][1][1],"Qy_1_Right_TPC_eta_0/D");
	out_tree->Branch("Qy_1_Right_TPC_eta_0_2",&QyTree[0][1][2],"Qy_1_Right_TPC_eta_0_2/D");
	out_tree->Branch("Qy_1_Right_TPC_eta_0_5",&QyTree[0][1][3],"Qy_1_Right_TPC_eta_0_5/D");
	out_tree->Branch("Qy_2_Left_FHCal",&QyTree[1][0][0],"Qy_2_Left_FHCal/D");
	out_tree->Branch("Qy_2_Left_TPC_eta_0",&QyTree[1][0][1],"Qy_2_Left_TPC_eta_0/D");
	out_tree->Branch("Qy_2_Left_TPC_eta_0_2",&QyTree[1][0][2],"Qy_2_Left_TPC_eta_0_2/D");
	out_tree->Branch("Qy_2_Left_TPC_eta_0_5",&QyTree[1][0][3],"Qy_2_Left_TPC_eta_0_5/D");
	out_tree->Branch("Qy_2_Right_FHCal",&QyTree[1][1][0],"Qy_2_Right_FHCal/D");
	out_tree->Branch("Qy_2_Right_TPC_eta_0",&QyTree[1][1][1],"Qy_2_Right_TPC_eta_0/D");
	out_tree->Branch("Qy_2_Right_TPC_eta_0_2",&QyTree[1][1][2],"Qy_2_Right_TPC_eta_0_2/D");
	out_tree->Branch("Qy_2_Right_TPC_eta_0_5",&QyTree[1][1][3],"Qy_2_Right_TPC_eta_0_5/D");
	out_tree->Branch("Qw_1_Left_FHCal",&QwTree[0][0][0],"Qw_1_Left_FHCal/I");
	out_tree->Branch("Qw_1_Left_TPC_eta_0",&QwTree[0][0][1],"Qw_1_Left_TPC_eta_0/I");
	out_tree->Branch("Qw_1_Left_TPC_eta_0_2",&QwTree[0][0][2],"Qw_1_Left_TPC_eta_0_2/I");
	out_tree->Branch("Qw_1_Left_TPC_eta_0_5",&QwTree[0][0][3],"Qw_1_Left_TPC_eta_0_5/I");
	out_tree->Branch("Qw_1_Right_FHCal",&QwTree[0][1][0],"Qw_1_Right_FHCal/I");
	out_tree->Branch("Qw_1_Right_TPC_eta_0",&QwTree[0][1][1],"Qw_1_Right_TPC_eta_0/I");
	out_tree->Branch("Qw_1_Right_TPC_eta_0_2",&QwTree[0][1][2],"Qw_1_Right_TPC_eta_0_2/I");
	out_tree->Branch("Qw_1_Right_TPC_eta_0_5",&QwTree[0][1][3],"Qw_1_Right_TPC_eta_0_5/I");
	out_tree->Branch("Qw_2_Left_FHCal",&QwTree[1][0][0],"Qw_2_Left_FHCal/I");
	out_tree->Branch("Qw_2_Left_TPC_eta_0",&QwTree[1][0][1],"Qw_2_Left_TPC_eta_0/I");
	out_tree->Branch("Qw_2_Left_TPC_eta_0_2",&QwTree[1][0][2],"Qw_2_Left_TPC_eta_0_2/I");
	out_tree->Branch("Qw_2_Left_TPC_eta_0_5",&QwTree[1][0][3],"Qw_2_Left_TPC_eta_0_5/I");
	out_tree->Branch("Qw_2_Right_FHCal",&QwTree[1][1][0],"Qw_2_Right_FHCal/I");
	out_tree->Branch("Qw_2_Right_TPC_eta_0",&QwTree[1][1][1],"Qw_2_Right_TPC_eta_0/I");
	out_tree->Branch("Qw_2_Right_TPC_eta_0_2",&QwTree[1][1][2],"Qw_2_Right_TPC_eta_0_2/I");
	out_tree->Branch("Qw_2_Right_TPC_eta_0_5",&QwTree[1][1][3],"Qw_2_Right_TPC_eta_0_5/I");
   #endif

   //MakeTTree("tree_EP", QxTree, QyTree, QwTree);
   //const Int_t N_det[_N_METHOD]={3,0,1,2};
   //const Int_t N_side[_N_ARM]={1,0};
   Double_t total_mom[_N_METHOD-1][_N_ARM];
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      #ifndef debug
      cout << "Event # " << jentry << endl;
      #endif

      TVector3* particle = new TVector3[n_tracks_mpd];
      Int_t seltracks=0;
      for (Long64_t track = 0; track < n_tracks_mpd; ++track){
      	if (signed_pt_mpd[track]<0.2 && signed_pt_mpd[track]>3.) continue;
      	particle[seltracks].SetPtEtaPhi(Abs(signed_pt_mpd[track]),eta_mpd[track],ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])));
      	//etaTree[track] = eta_mpd[track];
      	#ifdef make_tree
      	pt_tree[seltracks] = Abs(signed_pt_mpd[track]);
      	#endif
      	seltracks++;
      }

      //Qx,y_1_Left_FHCal
	  	GetQsZdc(ZDC_energy_mpd,1,1,QxTree[0][0][0],QyTree[0][0][0]);
	  	//Qx,y_1_Right_FHCal
	  	GetQsZdc(ZDC_energy_mpd,0,1,QxTree[0][1][0],QyTree[0][1][0]);
	  	//Qx,y_2_Left_FHCal
	  	GetQsZdc(ZDC_energy_mpd,1,2,QxTree[1][0][0],QyTree[1][0][0]);
	  	//Qx,y_2_Right_FHCal
	  	GetQsZdc(ZDC_energy_mpd,0,2,QxTree[1][1][0],QyTree[1][1][0]);
	  	//Qx,y_1_Left_TPC-eta-0
	  	GetQsTpc(0,1,particle,seltracks, -1, 1, QxTree[0][0][1],QyTree[0][0][1]);
	  	//Qx,y_1_Right_TPC-eta-0
	  	GetQsTpc(0,1,particle,seltracks,  1, 1, QxTree[0][1][1],QyTree[0][1][1]);
	  	//Qx,y_1_Left_TPC-eta-0_2
	  	GetQsTpc(1,1,particle,seltracks, -1, 1, QxTree[0][0][2],QyTree[0][0][2]);
	  	//Qx,y_1_Right_TPC-eta-0_2
	  	GetQsTpc(1,1,particle,seltracks,  1, 1, QxTree[0][1][2],QyTree[0][1][2]);
	  	//Qx,y_1_Left_TPC-eta-0_5
	  	GetQsTpc(2,1,particle,seltracks, -1, 1, QxTree[0][0][3],QyTree[0][0][3]);
	  	//Qx,y_1_Right_TPC-eta-0_5
	  	GetQsTpc(2,1,particle,seltracks,  1, 1, QxTree[0][1][3],QyTree[0][1][3]);
	  	//Qx,y_2_Left_TPC-eta-0
	  	GetQsTpc(0,1,particle,seltracks, -1, 1, QxTree[1][0][1],QyTree[1][0][1]);
	  	//Qx,y_2_Right_TPC-eta-0
	  	GetQsTpc(0,1,particle,seltracks,  1, 2, QxTree[1][1][1],QyTree[1][1][1]);
	  	//Qx,y_2_Left_TPC-eta-0_2
	  	GetQsTpc(1,1,particle,seltracks, -1, 2, QxTree[1][0][2],QyTree[1][0][2]);
	  	//Qx,y_2_Right_TPC-eta-0_2
	  	GetQsTpc(1,1,particle,seltracks,  1, 2, QxTree[1][1][2],QyTree[1][1][2]);
	  	//Qx,y_2_Left_TPC-eta-0_5
	  	GetQsTpc(2,1,particle,seltracks, -1, 2, QxTree[1][0][3],QyTree[1][0][3]);
	  	//Qx,y_2_Right_TPC-eta-0_5
	  	GetQsTpc(2,1,particle,seltracks,  1, 2, QxTree[1][1][3],QyTree[1][1][3]);
	  	//Qw_1_Left_FHCal
	  	QwTree[0][0][0] = GetTotalEnergy(ZDC_energy_mpd, 1);
	  	//Qw_1_Right_FHCal
	  	QwTree[0][1][0] = GetTotalEnergy(ZDC_energy_mpd, 0);
	  	//Qw_2_Left_FHCal
	  	QwTree[1][0][0] = QwTree[0][0][0];
	  	//Qw_2_Right_FHCal
	  	QwTree[1][1][0] = QwTree[0][1][0];
	  	//Qw_1_Left_TPC-eta-0
	  	QwTree[0][0][1] = GetMultiplicityTPC(0, particle, seltracks, -1);
	  	//Qw_1_Right_TPC-eta-0
	  	QwTree[0][1][1] = GetMultiplicityTPC(0, particle, seltracks,  1);
	  	//Qw_1_Left_TPC-eta-0_2
	  	QwTree[0][0][2] = GetMultiplicityTPC(1, particle, seltracks, -1);
	  	//Qw_1_Right_TPC-eta-0_2
	  	QwTree[0][1][2] = GetMultiplicityTPC(1, particle, seltracks,  1);
	  	//Qw_1_Left_TPC-eta-0_5
	  	QwTree[0][0][3] = GetMultiplicityTPC(2, particle, seltracks, -1);
	  	//Qw_1_Right_TPC-eta-0_5
	  	QwTree[0][1][3] = GetMultiplicityTPC(2, particle, seltracks,  1);
	  	//Qw_2_Left_TPC-eta-0
	  	QwTree[1][0][1] = QwTree[0][0][1];
	  	//Qw_2_Right_TPC-eta-0
	  	QwTree[1][1][1] = QwTree[0][1][1];
	  	//Qw_2_Left_TPC-eta-0_2
	  	QwTree[1][0][2] = QwTree[0][0][2];
	  	//Qw_2_Right_TPC-eta-0_2
	  	QwTree[1][1][2] = QwTree[0][1][2];
	  	//Qw_2_Left_TPC-eta-0_5
	  	QwTree[1][0][3] = QwTree[0][0][3];
	  	//Qw_2_Right_TPC-eta-0_5
	  	QwTree[1][1][3] = QwTree[0][1][3];
	  #ifdef make_tree
	  	out_tree->Fill();
	  #endif

      //Centrality selection;

      Int_t multiplicity = GetMultiplicityTPC();
      if (multiplicity==0) continue;
      /*Int_t cenrality_bin= GetCentralityBinFlow(multiplicity); //getting b bin of a given event for flow measurements
      if (cenrality_bin == -1) continue;
      Int_t cenrality_bin_res = GetCentralityBinRes(multiplicity);//getting the resolution bin in whic the event is 
	  if (cenrality_bin_res == -1) continue;*/

      Int_t cenrality_bin = -1;
      for(Int_t iC=0;iC<NimpactBins;iC++){
        if(b_mc>=b_bins_range[iC] && b_mc<b_bins_range[iC+1]){
          cenrality_bin = iC;
        }
      }
      //if (cenrality_bin==-1) continue;
      Int_t cenrality_bin_res = cenrality_bin;

      #ifdef debug
      Int_t w_sign = 0;
      Double_t Q1R = 0., Qy1R = 0., Q1L = 0., Qy1L = 0.;

      if (jentry%1000==0){
      	cout << "---------------------------------------------------------------------------------------------------------" << endl;
      	cout << "| Event " << jentry << "; Ntracks = " << n_tracks_mpd << "; Total momenta = " << GetTotalMomenta(0,particle,seltracks,1) << " + " << GetTotalMomenta(0,particle,seltracks,-1) << "|" << endl;
      	cout << "| Event " << jentry << "; Ntracks = " << n_tracks_mpd << "; Total momenta = " << GetTotalMomenta(1,particle,seltracks,1) << " + " << GetTotalMomenta(1,particle,seltracks,-1) << "|" << endl;
      	cout << "| Event " << jentry << "; Ntracks = " << n_tracks_mpd << "; Total momenta = " << GetTotalMomenta(2,particle,seltracks,1) << " + " << GetTotalMomenta(2,particle,seltracks,-1) << "|" << endl;
      	cout << "---------------------------------------------------------------------------------------------------------" << endl;
      	for (Int_t i_track=0;i_track<10;i_track++){
      		cout << "part[" << i_track << "].Pt() = " << particle[i_track].Pt() << "; part[" << i_track << "].Eta() = " << particle[i_track].Eta() << "; part[" << i_track << "].Phi() = " << particle[i_track].Phi() << "; ";
      		if (particle[i_track].Eta()>0) w_sign = 1;
      		else w_sign = -1;
      		cout << "Q1(track = " << i_track << ") TPC = " << w_sign*particle[i_track].Pt() / GetTotalMomenta(0,particle,seltracks,w_sign) * Cos(1.*particle[i_track].Phi()) << "; " << endl;
      	}
      	//for (Long64_t it = 0; it<n_tracks_mpd;it++){
//////////////////////////////////FIXIT/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      		////////////////////////////////////////
      		////////////////////////////////////////
      		////////////////////////////////////////
      		//if (particle[it].Eta()>0) Q1R += particle[it].Pt() / GetTotalMomenta(0,particle,n_tracks_mpd,1) * Cos(1.*particle[it].Phi());
      		//if (particle[it].Eta()>0) Qy1R += particle[it].Pt() / GetTotalMomenta(0,particle,n_tracks_mpd,1) * Sin(1.*particle[it].Phi());
      		//if (particle[it].Eta()<0) Q1L += -1*particle[it].Pt() / GetTotalMomenta(0,particle,n_tracks_mpd,1) * Cos(1.*particle[it].Phi());
      		//if (particle[it].Eta()<0) Qy1L += -1*particle[it].Pt() / GetTotalMomenta(0,particle,n_tracks_mpd,1) * Sin(1.*particle[it].Phi());
      	//}
      	/*cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
      	cout << "|Qx1R = " << Q1R << "; Qy1R = " << Qy1R << "; Psi_1_R = " << ATan2(Qy1R,Q1R)*180./Pi() << "; |" << endl;
      	cout << "|Qx1L = " << Q1L << "; Qy1L = " << Qy1L << "; Psi_1_L = " << ATan2(Qy1L,Q1L)*180./Pi() << "; |" << endl;
      	cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
      	*/
      }
      #endif
      for (Int_t igap=0;igap<_N_METHOD-1;igap++){
      	if (GetTotalMomenta(igap,particle,seltracks, 1)==0 || GetTotalMomenta(igap,particle,seltracks,-1)==0) continue;

      	Double_t Qq1x,Qq1y;
	  	Int_t bin = b_mc - 2;
	  	for (Int_t iC=0;iC<NcentralityBinsRes;iC++){
	  		if ((Int_t) b_mc >= centralityBinsRes[iC] && (Int_t) b_mc < centralityBinsRes[iC+1]){
	  			for (Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
					Double_t psi_N_R = GetPsiHalfTpc(igap,1,particle,seltracks, 1, i_harm + 1, Qq1x, Qq1y);
					Double_t psi_N_L = GetPsiHalfTpc(igap,1,particle,seltracks,-1, i_harm + 1, Qq1x, Qq1y);
					h_delta_fit[iC][i_harm][igap]->Fill(Abs(Unfold(psi_N_R, psi_N_L,i_harm+1)*180./Pi()));
				}
	  		}
	  	}
      }

      //skip the event if the ZDC signal is zero
	  if ((GetTotalEnergy(ZDC_energy_mpd, 0) == 0) || (GetTotalEnergy(ZDC_energy_mpd, 1) == 0)) continue;


	  Double_t q1x,q1y;
	  Int_t bin = b_mc - 2;
	  for (Int_t iC=0;iC<NcentralityBinsRes;iC++){
	  	if ((Int_t) b_mc >= centralityBinsRes[iC] && (Int_t) b_mc < centralityBinsRes[iC+1]){
	  		for (Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
				Double_t psi_N_R = GetPsiHalfZdc(ZDC_energy_mpd, 0, i_harm + 1, q1x, q1y);
				Double_t psi_N_L = GetPsiHalfZdc(ZDC_energy_mpd, 1, i_harm + 1, q1x, q1y);
				h_delta_fit[iC][i_harm][3]->Fill(Abs(Unfold(psi_N_R, psi_N_L,i_harm+1)*180./Pi()));
			}
	  	}
	  }
	  /*for(Int_t i_weigth=0;i_weigth<2;i_weigth++){
	  	for (Int_t n=0;n<_N_HARM;n++){
	  		GetQsZdc(ZDC_energy_mpd,0,n+1, Q[i_weigth][0][n][3][0], Q[i_weigth][0][n][3][1]);
	  		GetQsZdc(ZDC_energy_mpd,1,n+1, Q[i_weigth][1][n][3][0], Q[i_weigth][1][n][3][1]);
	  		for (Int_t igap=0;igap<_N_METHOD-1;igap++){
	  			GetQsTpc(igap,i_weigth,particle, n_tracks_mpd, 1, n+1, Q[i_weigth][0][n][igap][0], Q[i_weigth][0][n][igap][1]);
	  			GetQsTpc(igap,i_weigth,particle, n_tracks_mpd,-1, n+1, Q[i_weigth][1][n][igap][0], Q[i_weigth][1][n][igap][1]);
	  		}
	  	}
	  }*/

	  for (Int_t igap=0;igap<_N_METHOD-1;igap++){
      	if (GetTotalMomenta(igap,particle,seltracks, 1)==0 || GetTotalMomenta(igap,particle,seltracks,-1)==0) continue;
      	FillTPC(igap,1,cenrality_bin_res, particle, seltracks);
      }
      FillZDC(cenrality_bin_res,ZDC_energy_mpd);

	  for (Int_t igap=0;igap<_N_METHOD-1;igap++){
	  	//Fill3Sub(igap, particle, seltracks,ZDC_energy_mpd);
	  }

	  if (cenrality_bin==-1) continue;

	  //////////////////////////////////////////////////////////////////////////
	  //Int_t sign_side[2]={-1,1};
	  //const Int_t N_det[_N_METHOD]={3,0,1,2};
   	  //const Int_t N_side[_N_ARM]={1,0};
	  /*for (Int_t harm=0;harm<_N_HARM;harm++){
	  	for (Int_t side=0;side<_N_ARM;side++){
	  		GetQsZdc(ZDC_energy_mpd,N_side[side],harm+1, QxTree[harm][side][0], QyTree[harm][side][0]);
	  		QwTree[harm][side][0] = GetTotalEnergy(ZDC_energy_mpd, N_side[side]);
	  		for (Int_t det=1;det<_N_METHOD-1;det++){
	  			GetQsTpc(det-1,1,particle, n_tracks_mpd, sign_side[side], harm+1, QxTree[harm][side][det], QyTree[harm][side][det]);
	  			QwTree[harm][side][det] = GetMultiplicityTPC(det-1, particle, sign_side[side]);
	  		}
	  	}
	  }*/
	  	
	  //FillTTree();
	  //////////////////////////////////////////////////////////////////////////

	  //MakeTTree("tree_EP", QxTree, QyTree, QwTree, phiEP_mc, b_mc);

	  #ifdef debug
	  	if (jentry%1000==0){
	  		/*cout << "*********************************************************************************************************" << endl;
	  		//cout << "Event " << jentry << ": Q1x(R,TPC) = " << Q[0][0][0][0] << "; Q1x(L,TPC) = " << Q[1][0][0][0] << "|| Q1x(R,FHCal) = " << Q[0][0][1][0] << "; Q1x(L,FHCal) = " << Q[1][0][1][0] << ";" << endl;
	  		printf("Event %4i: Q1x(R,TPC) = %+1.4f; Q1x(L,TPC) = %+1.4f || Q1x(R,FHCal) = %+1.4f; Q1x(L,FHCal) = %+1.4f;\n",(Int_t) jentry,Q[0][0][0][0][0],Q[0][1][0][0][0],Q[0][0][0][3][0],Q[0][1][0][3][0]);
	  		//cout << "Event " << jentry << ": Q1y(R,TPC) = " << Q[0][0][0][1] << "; Q1y(L,TPC) = " << Q[1][0][0][1] << "|| Q1y(R,FHCal) = " << Q[0][0][1][1] << "; Q1y(L,FHCal) = " << Q[1][0][1][1] << ";" << endl;
	  		printf("Event %4i: Q1y(R,TPC) = %+1.4f; Q1y(L,TPC) = %+1.4f || Q1y(R,FHCal) = %+1.4f; Q1y(L,FHCal) = %+1.4f;\n",(Int_t) jentry,Q[0][0][0][0][1],Q[0][1][0][0][1],Q[0][0][0][3][1],Q[0][1][0][3][1]);
	  		//cout << "Event " << jentry << ": Q2x(R,TPC) = " << Q[0][1][0][0] << "; Q2x(L,TPC) = " << Q[1][1][0][0] << "|| Q2x(R,FHCal) = " << Q[0][1][1][0] << "; Q2x(L,FHCal) = " << Q[1][1][1][0] << ";" << endl;
	  		printf("Event %4i: Q2x(R,TPC) = %+1.4f; Q2x(L,TPC) = %+1.4f || Q2x(R,FHCal) = %+1.4f; Q2x(L,FHCal) = %+1.4f;\n",(Int_t) jentry,Q[0][0][1][0][0],Q[0][1][1][0][0],Q[0][0][1][3][0],Q[0][1][1][3][0]);
	  		//cout << "Event " << jentry << ": Q2y(R,TPC) = " << Q[0][1][0][1] << "; Q2y(L,TPC) = " << Q[1][1][0][1] << "|| Q2y(R,FHCal) = " << Q[0][1][1][1] << "; Q2y(L,FHCal) = " << Q[1][1][1][1] << ";" << endl;
	  		printf("Event %4i: Q1y(R,TPC) = %+1.4f; Q1y(L,TPC) = %+1.4f || Q1y(R,FHCal) = %+1.4f; Q1y(L,FHCal) = %+1.4f;\n",(Int_t) jentry,Q[0][0][1][0][1],Q[0][1][1][0][1],Q[0][0][1][3][1],Q[0][1][1][3][1]);
	  		cout << "*********************************************************************************************************" << endl;
	  		*/
	  		Double_t Q1xR=0., Q1yR=0., Q1xL=0., Q1yL=0.;
	  		Double_t Q1xRFHCal=0., Q1yRFHCal=0., Q1xLFHCal=0., Q1yLFHCal=0.;
			
	  		//GetQsTpc(0,1,particle, seltracks, 1, 1, Q1xR, Q1yR);
	  		//GetQsTpc(0,1,particle, seltracks, 1,-1, Q1xL, Q1yL);
	  		//cout << "Q1x(R,TPC) = " << Q1xR << "; Q1y(R,TPC) = " << Q1yR << "; Psi_1_R = " << GetPsiHalfTpc(0,1,particle,seltracks, 1,1,Q1xR,Q1yR)*180./Pi() << "; Psi_1_R(FHCal) = " << GetPsiHalfZdc(ZDC_energy_mpd,0,1,Q1xRFHCal,Q1yRFHCal)*180./Pi() << ";" << endl;
			//cout << "Q1x(L,TPC) = " << Q1xL << "; Q1y(L,TPC) = " << Q1yL << "; Psi_1_L = " << GetPsiHalfTpc(0,1,particle,seltracks,-1,1,Q1xL,Q1yL)*180./Pi() << "; Psi_1_L(FHCal) = " << GetPsiHalfZdc(ZDC_energy_mpd,1,1,Q1xLFHCal,Q1yLFHCal)*180./Pi() << ";" << endl;
	  		cout << "*********************************************************************************************************" << endl;
	  		cout << "Psi_RP = " << ATan2(Sin(phiEP_mc),Cos(phiEP_mc))*180./Pi() << "; Psi_1_Full(TPC) = " << GetPsiFullTpc(0,1,particle,seltracks,1)*180./Pi() << "; Psi_1_Full(FHCal) = " << GetPsiFullZdc(ZDC_energy_mpd,1)*180./Pi() << ";" << endl;
	  		cout << endl;
	  		Int_t count_L=0., count_R=0.;
	  		for (Long64_t it=0; it < seltracks; it++){
	  			if (particle[it].Eta()<0 && count_L==0){
	  			 cout << "//////////////////////////////////////////////////////////////////////////////" << endl;
	  			 cout << "phi_mc Left  =" << ATan2(py_mc[id_from_mc_mpd[it]],px_mc[id_from_mc_mpd[it]]) << "; phi_mpd L = " << particle[it].Phi() << ";" << endl;
	  			 count_L++;
	  			}
	  			if (particle[it].Eta()>0 && count_R==0){
	  			 cout << "phi_mc Right =" << ATan2(py_mc[id_from_mc_mpd[it]],px_mc[id_from_mc_mpd[it]]) << "; phi_mpd R = " << particle[it].Phi() << ";" << endl;
	  			 count_R++;
	  			 cout << "//////////////////////////////////////////////////////////////////////////////" << endl;
	  			}
	  			if (count_R>0 && count_L>0) break;
	  		}
	  	}
	  	if (jentry==nentries-1) cout << "Done." << endl;
	  #endif
	  for (Int_t i_arm=0;i_arm<_N_ARM;i_arm++){
	  	for (Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
	  		for (Int_t method=0;method<_N_METHOD;method++){
	  			for (Int_t i_proj=0;i_proj<2;i_proj++){
	  				h_Qvsb[i_arm][i_harm][method][i_proj]   ->Fill(b_mc,Q[0][i_arm][i_harm][method][i_proj]);
	  				h_QvsMult[i_arm][i_harm][method][i_proj]->Fill(multiplicity,Q[0][i_arm][i_harm][method][i_proj]);
	  			}
	  		}
	  	}
	  }

	  Double_t phi_EP = ATan2(Sin(phiEP_mc),Cos(phiEP_mc)); //unfold the generated event plane
	  Double_t Psi_1_FULL_ZDC = GetPsiFullZdc(ZDC_energy_mpd,1); //getting event plane using ZDC for 1st harmonic
	  Double_t Psi_2_FULL_ZDC = GetPsiFullZdc(ZDC_energy_mpd,2);
	  /*Double_t Psi_1_FULL_TPC[3];
	  Double_t Psi_2_FULL_TPC[3];
	  for (Int_t method=0;method<_N_METHOD-1;method++){
	  	Psi_1_FULL_TPC[method] = GetPsiFullTpc(method,0,particle,n_tracks_mpd,1);
	  	Psi_2_FULL_TPC[method] = GetPsiFullTpc(method,0,particle,n_tracks_mpd,2);	
	  	total_mom[method][0] = GetTotalMomenta(method,particle,n_tracks_mpd, 1);
	  	total_mom[method][1] = GetTotalMomenta(method,particle,n_tracks_mpd,-1);
	  }*/

	  //Loop over tracks (MPD);
	  //---------------------------------------------------------------------------------------------------------
	  for (Long64_t track = 0; track < n_tracks_mpd; ++track){
	  	if (id_from_mc_mpd[track] == -1) continue; //equivalent to mother ID cut
		if (n_hits_mpd[track] < Cut_No_Of_hits_min) continue; //n hits in TPC cut
		Int_t p_sort=-1;
		if (PDG_code_mc[id_from_mc_mpd[track]] == 2212)     p_sort = 0; //proton selection;
		else if (PDG_code_mc[id_from_mc_mpd[track]] == 321) p_sort = 1; //kaon selection;
		else if (PDG_code_mc[id_from_mc_mpd[track]] == 211) p_sort = 2; //pion selection;
		if (p_sort ==-1) continue;
		Int_t PDG = PDG_code_mc[id_from_mc_mpd[track]];
		Float_t Pt = Abs(signed_pt_mpd[track]);
		Float_t Eta = eta_mpd[track];
		Float_t Phi = ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track]));
		Float_t Rapidity = 0.5*TMath::Log((energy_mc[id_from_mc_mpd[track]] + pz_mc[id_from_mc_mpd[track]])/(energy_mc[id_from_mc_mpd[track]] - pz_mc[id_from_mc_mpd[track]]));

		Int_t sign;
		if (Eta < 0) sign = -1;
		else sign = 1;

		for (Int_t igap=0;igap<_N_METHOD-1;igap++){
			//To avoid autocorrelations:
			Double_t QRc=0.,QRs=0.,QR2c=0.,QR2s=0.,QLc=0.,QLs=0.,QL2c=0.,QL2s=0.;
			if (Eta>EtaGap[igap]){
				QRc = QxTree[0][1][igap+1]+ /*Q[0][0][0][igap][0] - Pt/total_mom[igap][0] * */Cos(Phi);
				QRs = QyTree[0][1][igap+1]+ /*Q[0][0][0][igap][1] - Pt/total_mom[igap][0] * */Sin(Phi);
				QR2c= QxTree[1][1][igap+1]+ /*Q[0][0][1][igap][0] - Pt/total_mom[igap][0] * */Cos(2*Phi);
				QR2s= QyTree[1][1][igap+1]+ /*Q[0][0][1][igap][1] - Pt/total_mom[igap][0] * */Sin(2*Phi);
			}
			else{
				if (Eta<-EtaGap[igap]){
					QLc = QxTree[0][0][igap+1]-  /*Q[0][1][0][igap][0] - Pt/total_mom[igap][1] * */Cos(Phi);
					QLs = QyTree[0][0][igap+1]-  /*Q[0][1][0][igap][1] - Pt/total_mom[igap][1] * */Sin(Phi);
					QL2c= QxTree[1][0][igap+1]-  /*Q[0][1][1][igap][0] - Pt/total_mom[igap][1] * */Cos(2*Phi);
					QL2s= QyTree[1][0][igap+1]-  /*Q[0][1][1][igap][1] - Pt/total_mom[igap][1] * */Sin(2*Phi);
				}
				else{
					QRc = QxTree[0][1][igap+1] /*Q[0][0][0][igap][0]*/; QLc = QxTree[0][0][igap+1] /*Q[0][1][0][igap][0]*/;
					QRs = QyTree[0][1][igap+1] /*Q[0][0][0][igap][1]*/; QLs = QyTree[0][0][igap+1] /*Q[0][1][0][igap][1]*/;
					QR2c= QxTree[1][1][igap+1] /*Q[0][0][1][igap][0]*/; QL2c= QxTree[1][0][igap+1] /*Q[0][1][1][igap][0]*/;
					QR2s= QyTree[1][1][igap+1] /*Q[0][0][1][igap][1]*/; QL2s= QyTree[1][0][igap+1] /*Q[0][1][1][igap][1]*/;
				}
			}
			Double_t Psi_1_F_TPC = ATan2(QRs+QLs,QRc+QLc);
			Double_t Psi_2_F_TPC = ATan2(QR2s+QL2s,QR2c+QL2c)/2;
			if ((Abs(Eta) > 0.2) && (Abs(Eta) < Cut_Eta_Max))
				p_flow_wrt_full_vs_pt[p_sort][cenrality_bin][0][0][igap]->Fill(Pt,sign*Cos(Psi_1_F_TPC - Phi));
			if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min))
				p_flow_wrt_full_vs_rapidity[p_sort][cenrality_bin][0][0][igap]->Fill(Rapidity,Cos(Psi_1_F_TPC - Phi));
			if (Pt > Cut_Pt_Min)
				p_flow_wrt_full_vs_eta[p_sort][cenrality_bin][0][0][igap]->Fill(Eta,Cos(Psi_1_F_TPC /*phi_EP*/ - Phi));
			if (Abs(Eta) < Cut_Eta_Max){
				p_flow_wrt_full_vs_pt[p_sort][cenrality_bin][1][0][igap]->Fill(Pt,Cos(2.*(Psi_1_F_TPC - Phi)));
				p_flow_wrt_full_vs_pt[p_sort][cenrality_bin][1][1][igap]->Fill(Pt,Cos(2.*(Psi_2_F_TPC - Phi)));
			}
			if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min)){
				p_flow_wrt_full_vs_rapidity[p_sort][cenrality_bin][1][0][igap]->Fill(Rapidity,Cos(2*(Psi_1_F_TPC - Phi)));
				p_flow_wrt_full_vs_rapidity[p_sort][cenrality_bin][1][1][igap]->Fill(Rapidity,Cos(2*(Psi_2_F_TPC - Phi)));
			}
			if (Pt > Cut_Pt_Min){	
				p_flow_wrt_full_vs_eta[p_sort][cenrality_bin][1][0][igap]->Fill(Eta,Cos(2.*(Psi_1_F_TPC - Phi)));
				p_flow_wrt_full_vs_eta[p_sort][cenrality_bin][1][1][igap]->Fill(Eta,Cos(2.*(Psi_2_F_TPC - Phi)));
			}
		}

		if ((Abs(Eta) > 0.2) && (Abs(Eta) < Cut_Eta_Max))
			p_flow_wrt_full_vs_pt[p_sort][cenrality_bin][0][0][3]->Fill(Pt,sign*Cos(Psi_1_FULL_ZDC - Phi)/Res1Flow[cenrality_bin]);
		if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min))
			p_flow_wrt_full_vs_rapidity[p_sort][cenrality_bin][0][0][3]->Fill(Rapidity,Cos(Psi_1_FULL_ZDC - Phi)/Res1Flow[cenrality_bin]);
		if (Pt > Cut_Pt_Min)
				p_flow_wrt_full_vs_eta[p_sort][cenrality_bin][0][0][3]->Fill(Eta,Cos(Psi_1_FULL_ZDC - Phi)/Res2Flow[cenrality_bin]);
		if (Abs(Eta) < Cut_Eta_Max){
			p_flow_wrt_full_vs_pt[p_sort][cenrality_bin][1][0][3]->Fill(Pt,Cos(2.*(Psi_1_FULL_ZDC - Phi))/Res2Flow[cenrality_bin]);
			p_flow_wrt_full_vs_pt[p_sort][cenrality_bin][1][1][3]->Fill(Pt,Cos(2.*(Psi_2_FULL_ZDC - Phi))/Res2Flow[cenrality_bin]);
		}
		if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min))
			p_flow_wrt_full_vs_rapidity[p_sort][cenrality_bin][1][0][3]->Fill(Rapidity,Cos(2*(Psi_1_FULL_ZDC - Phi))/Res1Flow[cenrality_bin]);
		if (Pt > Cut_Pt_Min)	
				p_flow_wrt_full_vs_eta[p_sort][cenrality_bin][1][0][3]->Fill(Eta,Cos(2.*(Psi_1_FULL_ZDC - Phi))/Res1Flow[cenrality_bin]);
		/////////////////////////////////////////////////////////////////////////////////////
		/*for (Int_t method=0;method<_N_METHOD-1;method++){
			if (Eta <-EtaGap[method]){
			 	h_Mult_vs_eta[0][method]->Fill(-Eta);
			 	p_Mult_vs_eta[0][method]->Fill(-Eta,QwTree[0][0][method+1],1);
			}
			if (Eta > EtaGap[method]){
			 	h_Mult_vs_eta[1][method]->Fill(Eta);
			 	p_Mult_vs_eta[1][method]->Fill(Eta,QwTree[0][1][method+1],1);
			}
			h_Ratio_eta[method]->Divide(h_Mult_vs_eta[0][method],h_Mult_vs_eta[1][method]);
		}*/

	  }//end of loop over tracks (MPD);
	  //---------------------------------------------------------------------------------------------------------

	  //Loop over tracks (MC);
	  //---------------------------------------------------------------------------------------------------------
	  for (Long64_t track = 0; track < n_tracks_mc; ++track){
	  	Int_t p_sort=-1;
		if (PDG_code_mc[track] == 2212)     p_sort = 0; //proton selection;
		else if (PDG_code_mc[track] == 321) p_sort = 1; //kaon selection;
		else if (PDG_code_mc[track] == 211) p_sort = 2; //pion selection;
		if (p_sort ==-1) continue;

		Int_t PDG = PDG_code_mc[id_from_mc_mpd[track]];
		Float_t Pt = Abs(pt_mc[track]);
		Float_t Eta = eta_mc[track];
		Float_t Phi = ATan2(py_mc[track],px_mc[track]);
		Float_t Rapidity = .5*TMath::Log((energy_mc[track] + pz_mc[track])/(energy_mc[track] - pz_mc[track]));

		Int_t sign;
		if (Eta < 0) sign = -1;
		else sign = 1;

		if ((Abs(Eta) > 0.2) && (Abs(Eta) < Cut_Eta_Max))
			p_flow_wrt_RP_vs_pt[p_sort][cenrality_bin][0]->Fill(Pt , sign*Cos(phi_EP - Phi));
		if (Abs(Eta) < Cut_Eta_Max)	
			p_flow_wrt_RP_vs_pt[p_sort][cenrality_bin][1]->Fill(Pt , Cos(2.*(phi_EP-Phi)));
		if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min))
			p_flow_wrt_RP_vs_rapidity[p_sort][cenrality_bin][0]->Fill(Rapidity , Cos(phi_EP	- Phi));
		if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min))	
			p_flow_wrt_RP_vs_rapidity[p_sort][cenrality_bin][1]->Fill(Rapidity , Cos(2.*(phi_EP - Phi)));
		if (Pt > Cut_Pt_Min)	
			p_flow_wrt_RP_vs_eta[p_sort][cenrality_bin][0]->Fill(Eta , Cos(phi_EP - Phi));
		if (Pt > Cut_Pt_Min)	
			p_flow_wrt_RP_vs_eta[p_sort][cenrality_bin][1]->Fill(Eta , Cos(2.*(phi_EP - Phi)));
	  }//end of loop over tracks (MC);
	  //---------------------------------------------------------------------------------------------------------

	  delete [] particle;

   }//end of loop over events;

	Double_t Res[_N_HARM][_N_HARM][_N_METHOD][NcentralityBinsRes];
	Double_t eRes[_N_HARM][_N_HARM][_N_METHOD][NcentralityBinsRes];
	Double_t Res3[_N_HARM][_N_HARM][_N_METHOD][NcentralityBinsRes];
	Double_t eRes3[_N_HARM][_N_HARM][_N_METHOD][NcentralityBinsRes];
	Double_t a,b,c,ea,eb,ec;
	for (Int_t harm=0; harm<_N_HARM; harm++){
		for (Int_t _harm=0;_harm<_N_HARM;_harm++){
			for (Int_t method=0;method<_N_METHOD;method++){
				for (Int_t iC=0;iC<NcentralityBinsRes;iC++){
					Res[harm][_harm][method][iC] = p_Res2Psi_vs_b[harm][_harm][method]->GetBinContent(iC+1);
					eRes[harm][_harm][method][iC]= p_Res2Psi_vs_b[harm][_harm][method]->GetBinError(iC+1);
					if(Res[harm][_harm][method][iC]>0){
						Res[harm][_harm][method][iC] = Sqrt(Res[harm][_harm][method][iC]);
						eRes[harm][_harm][method][iC]= Abs(eRes[harm][_harm][method][iC]/(2*Res[harm][_harm][method][iC]));
						h_Res_true[harm][_harm][method]->SetBinContent(iC+1,Res[harm][_harm][method][iC]);
						h_Res_true[harm][_harm][method]->SetBinError(iC+1,eRes[harm][_harm][method][iC]);
					}
					a = p_Cos_A_vs_b[harm][_harm][method]->GetBinContent(iC+1);
					b = p_Cos_B_vs_b[harm][_harm][method]->GetBinContent(iC+1);
					c = p_Cos_C_vs_b[harm][_harm][method]->GetBinContent(iC+1);
					ea= p_Cos_A_vs_b[harm][_harm][method]->GetBinError(iC+1);
					eb= p_Cos_B_vs_b[harm][_harm][method]->GetBinError(iC+1);
					ec= p_Cos_C_vs_b[harm][_harm][method]->GetBinError(iC+1);
					if(a*b/c >0){
						Res3[harm][_harm][method][iC] = Sqrt(a*b/c);
						eRes3[harm][_harm][method][iC]= 0.5*Sqrt(b*ea*ea/(a*c)+a*eb*eb/(b*c)+a*b*ec*ec/(c*c*c));
						h_Res3_true[harm][_harm][method]->SetBinContent(iC+1,Res3[harm][_harm][method][iC]);
						h_Res3_true[harm][_harm][method]->SetBinError(iC+1,eRes3[harm][_harm][method][iC]);
					}
				}
			}
		}
	}

   Write();
   //SaveTTree("Tree.root");
   #ifdef make_tree
   out_file->cd();
   out_tree->Write();
   out_file->Close();
   #endif
}









Int_t real_flow::GetMultiplicityTPC() //should be called in loop over events
{
	Int_t multiplicity = 0;
	for (Long64_t track = 0; track < n_tracks_mpd; ++track)//loop over mpdtracks, cut on mother ID will be everywhere
	{
		if (id_from_mc_mpd[track] == -1) continue; //equivalent to mother id cut	
		multiplicity++; //multiplicity in TPC before eta, nhits and pt cuts, but after mother ID cut
	}
	return multiplicity;
}

Int_t real_flow::GetMultiplicityTPC(Int_t gap, TVector3* part, Long64_t ntracks, Int_t sign) //should be called in loop over events
{
	Int_t multiplicity = 0;
	for (Long64_t track = 0; track < ntracks; ++track)//loop over mpdtracks, cut on mother ID will be everywhere
	{
		if (id_from_mc_mpd[track] == -1) continue; //equivalent to mother id cut
		if ( !(part[track].Eta()*sign > EtaGap[gap]) ) continue;//multiplicity in each side and eta-gap
		multiplicity++; //multiplicity in TPC before eta, nhits and pt cuts, but after mother ID cut
	}
	return multiplicity;
}

Int_t real_flow::GetCentralityBinFlow(Int_t multiplicity)
{
	Int_t centrality_bin = -1;
	Int_t  multiplicityBinsFlow[_N_CENTRALITY_BINS+1] = { 0,1000,2000};
	for (Int_t multiplicityBin = 0; multiplicityBin < _N_CENTRALITY_BINS; ++multiplicityBin)
	if (multiplicity > multiplicityBinsFlow[multiplicityBin] && multiplicity <= multiplicityBinsFlow[multiplicityBin+1]) 
		centrality_bin = multiplicityBin;
	return centrality_bin;
}

Int_t real_flow::GetCentralityBinRes(Int_t multiplicity)
{
	int centrality_bin = -1;
	for (int multiplicityBin = 0; multiplicityBin < NmultiplicityBinsRes; ++multiplicityBin)
	if (multiplicity > multiplicityBinsRes[multiplicityBin] && multiplicity <= multiplicityBinsRes[multiplicityBin+1]) 
		centrality_bin = multiplicityBin;
	return centrality_bin;
}

Double_t real_flow::GetTotalEnergy(Float_t* zdc_energy, Int_t zdc_ID)
{
	Double_t total_energy = 0.;
	for (int i = _MAX_MODULE/2 * zdc_ID; i < _MAX_MODULE/2 * (zdc_ID + 1); ++i)
	{
		if ((i==22) || (i==67)) continue;
	//	if ((i ==15)||(i ==21)||(i==23)||(i ==29)||(i==60)||
	//			(i==66)||(i==68)||(i==74)) continue;
		total_energy += zdc_energy[i];
	}
	return total_energy;
}

Double_t real_flow::GetPsiFullTpc(Int_t gap, Int_t weight, TVector3* part, Long64_t ntracks, Int_t harm)
{
	Double_t QcosR=0., QsinR=0.;
	GetQsTpc(gap, weight, part,ntracks,1,harm,QcosR,QsinR);

	Double_t QcosL=0., QsinL=0.;
	GetQsTpc(gap, weight, part,ntracks,-1,harm,QcosL,QsinL);

	Double_t psiEP = ATan2(QsinR + QsinL,QcosR + QcosL)/(Double_t) harm; // (-pi/n,pi/n]
	return psiEP;
}

Double_t real_flow::GetPsiFullZdc(Float_t* zdc_energy, Int_t n)
{
	Double_t QcosR=0., QsinR=0.;
	GetQsZdc(zdc_energy,0,n,QcosR, QsinR);

	Double_t QcosL=0., QsinL=0.;
	GetQsZdc(zdc_energy,1,n,QcosL,QsinL);

	Double_t psiEP = ATan2(QsinR + QsinL,QcosR + QcosL)/(Double_t) n; // (-pi/n,pi/n]
    return psiEP;
}

Double_t real_flow::GetTotalMomenta(Int_t gap, TVector3* part, Long64_t ntracks, Int_t sign)
{
	Double_t total_momenta = 0.;
	for (int ep_particle = 0; ep_particle < ntracks; ++ep_particle)
	{
		if (!(part[ep_particle].Eta()*sign > EtaGap[gap])) continue;
		total_momenta += part[ep_particle].Pt();
	}
	return total_momenta;
}

void real_flow::GetQsTpc(Int_t gap, Int_t weight, TVector3* part, Long64_t ntracks, Int_t sign, Int_t harm, Double_t &Qx, Double_t &Qy)
{
	Double_t Qcos=0., Qsin=0. , total_weight = 0.;


	if (weight==0){
	 	//total_weight = GetTotalMomenta(gap,part,sign);
	}
	if (weight==1){
		//total_weight = GetMultiplicityTPC(gap,part,sign);
	}

	Int_t side=0.;
	if (sign== 1) side=1;
	Float_t w_part=0.;
	for (Long64_t ep_particle = 0; ep_particle < ntracks; ++ep_particle)
	{
		/////////////////////////////////////////////////////////////
		if (!(part[ep_particle].Eta()*sign > EtaGap[gap])) continue;
		//w_part = sign;
		w_part = 1;
		//w_part = -sign;
		//if (harm==2) w_part = 1;
		if (weight==0) w_part = sign*part[ep_particle].Pt();
		Qcos += w_part /*/ total_weight*/ * Cos(harm*part[ep_particle].Phi());
		Qsin += w_part /*/ total_weight*/ * Sin(harm*part[ep_particle].Phi());
		h_Mult_vs_eta[harm-1][side][gap]->Fill(w_part*part[ep_particle].Eta());
	}

	Qx = Qcos;
	Qy = Qsin;
}

void real_flow::GetQsZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t harm, Double_t &Qx, Double_t &Qy)
{
	Double_t *phi_angle_of_modules = GetAngles();
	Double_t Qcos = 0 , Qsin = 0;
	Double_t w_sign = 0;
	Double_t total_energy = GetTotalEnergy(zdc_energy,zdc_ID);

	if (harm == 2) w_sign = 1;
	else if (harm == 1)
	{
		if (zdc_ID == 0) w_sign = 1;
		else if (zdc_ID == 1) w_sign = -1;
	}


	for (int module = _MAX_MODULE/2 * zdc_ID; module < _MAX_MODULE/2 * (zdc_ID + 1); ++module)
	{
		if ((module==22) || (module==67)) continue;
	    //if ((i ==15)||(i ==21)||(i==23)||(i ==29)||(i==60)||
	    //		(i==66)||(i==68)||(i==74)) continue;
		Qcos += w_sign*zdc_energy[module] / total_energy * Cos(harm*phi_angle_of_modules[module]);
        Qsin += w_sign*zdc_energy[module] / total_energy * Sin(harm*phi_angle_of_modules[module]);
	}

	Qx = Qcos; Qy = Qsin;
}

Double_t* real_flow::GetAngles()
{
	Double_t *phi_angle_of_module = new Double_t[_MAX_MODULE];

	for (int i = 0; i < _N_ARM; ++i)
	{
		Int_t x_axis_switch;
		if (i == 0) x_axis_switch = 1;
		else if (i == 1) x_axis_switch = -1;

		for (Int_t j = 0; j < _MAX_MODULE/2; ++j)
		{
			Double_t y = 0, x = 0;

			if ((j>=0) && (j<=4))
			{
				y = 45., x = (j-2)*15.;
				phi_angle_of_module[j + i*_MAX_MODULE/2] = ATan2(y,x_axis_switch*x);
			}
			else if ((j>=5) && (j<=39))
			{
				y = (3-(j+2)/7)*15, x = (3-(j+2)%7)*15;
				phi_angle_of_module[j + i*_MAX_MODULE/2] = ATan2(y,x_axis_switch*x);
			}
			else if ((j>=40) && (j<=44))
			{
				y = -45. , x = (j-42)*15.;
				phi_angle_of_module[j + i*_MAX_MODULE/2] = ATan2(y,x_axis_switch*x);
			}
		}
	}

	return phi_angle_of_module;
}
Float_t real_flow::Unfold(Float_t phiEP_mc, Float_t psi_N_FULL, Int_t harm)
{
	Float_t values[10] , absvalues[10];

	values[0] = phiEP_mc - psi_N_FULL;
	absvalues[0] = Abs(values[0]);
	values[1] = phiEP_mc + 2. * Pi() / harm - psi_N_FULL;
	absvalues[1] = Abs(values[1]);
	values[2] = phiEP_mc - 2. * Pi() / harm - psi_N_FULL;
	absvalues[2] = Abs(values[2]);
    return values[LocMin(3, absvalues)];
}

void real_flow::FillTPC(Int_t gap, Int_t weight, Int_t cenrality_bin, TVector3 *part, Long64_t ntracks)
{
	Bool_t centGood;
	if (cenrality_bin==-1) 
		centGood=kFALSE;
	else 
		centGood=kTRUE;
	if (gap<0 || gap>2){ 
		cout <<"gap index should be from 0 to 2 !" << endl;
	}
	else{
		Double_t qx, qy;
		Double_t psi_N_HALF;
		phiEP_mc = ATan2(Sin(phiEP_mc), Cos(phiEP_mc));
		if (centGood) h_psi_RP[cenrality_bin]->Fill(phiEP_mc*180./Pi());

		for (Int_t harm = 0; harm < _N_HARM; ++harm)
		{
			Double_t psi_N_R = GetPsiHalfTpc(gap,weight,part,ntracks, 1, harm+1, qx, qy);
			Double_t psi_N_L = GetPsiHalfTpc(gap,weight,part,ntracks, -1, harm+1, qx, qy);

			if (centGood) h_delta_psi[cenrality_bin][harm][gap]->Fill(Abs(Unfold(psi_N_R, psi_N_L,harm+1)*180./Pi()));
			//if (harm == 0) h_psi_RP[cenrality_bin]->Fill(psi_N_L*180./Pi());

			Double_t psi_N_FULL = GetPsiFullTpc(gap,weight,part,ntracks, harm+1);
			//if (harm==0) psi_N_FULL-=Pi();
			if (centGood) h_delta_EP_full_RP[cenrality_bin][harm][gap]->Fill(Abs(Unfold(phiEP_mc,psi_N_FULL,harm+1)*180./Pi()));

			for (Int_t _harm = 0; _harm < _N_HARM; ++_harm)
			{
				Double_t psi_N_R = GetPsiHalfTpc(gap,weight,part,ntracks, 1, _harm+1, qx, qy);
				//psi_N_R -=Pi();
				Double_t psi_N_L = GetPsiHalfTpc(gap,weight,part,ntracks, -1, _harm+1, qx, qy);
				//psi_N_L -=Pi();
				//cout << "Filling with centrality = " << centralityBinsRes[cenrality_bin]+0.1 << "Res2 = " << Cos((harm+1)*(psi_N_R - psi_N_L)) << endl;
				p_Res2Psi_vs_b[harm][_harm][gap]->Fill(b_mc, Cos((harm+1)*(psi_N_R - psi_N_L )));
				p_ResAPsi_vs_b[harm][_harm][gap]->Fill(b_mc, Cos((harm+1)*(psi_N_R - phiEP_mc)));
				p_ResBPsi_vs_b[harm][_harm][gap]->Fill(b_mc, Cos((harm+1)*(psi_N_L - phiEP_mc)));

				Double_t psi_N_FULL = GetPsiFullTpc(gap,weight,part,ntracks, harm+1);
				//psi_N_FULL-=Pi();///////////////////////////////check cosine distribution;
				p_true_Res_vs_b[harm][_harm][gap]->Fill(b_mc,Cos((harm+1)*(psi_N_FULL - phiEP_mc)));

				for (Int_t arm = 0; arm < _N_ARM; ++arm)
				{
					Double_t psi_N_HALF = GetPsiHalfTpc(gap,weight,part,ntracks, -2*arm + 1, _harm+1, qx, qy);
					p_true_Res_half_vs_b[arm][harm][_harm][gap]->Fill(b_mc,Cos((harm+1)*Abs(Unfold(psi_N_HALF,phiEP_mc,_harm+1))));
				}
			}

			for (Int_t arm = 0; arm < _N_ARM; ++arm)
			{
				psi_N_HALF = GetPsiHalfTpc(gap,weight,part,ntracks, -2*arm + 1, harm+1, qx, qy);
				if (centGood) h_qx[cenrality_bin][arm][harm][gap]->Fill(qx);
				if (centGood) h_qy[cenrality_bin][arm][harm][gap]->Fill(qy);
				p_qx_vs_b[arm][harm][gap]->Fill(b_mc, qx);
				p_qy_vs_b[arm][harm][gap]->Fill(b_mc, qy);
				if (centGood) h_delta_EP_half_RP[cenrality_bin][arm][harm][gap]->Fill(psi_N_HALF - phiEP_mc);
				if (centGood) h_cor_EP_RP[cenrality_bin][arm][harm][gap]->Fill(psi_N_HALF,phiEP_mc);

				if (centGood) h_psi[cenrality_bin][arm][harm][gap]->Fill(psi_N_HALF);
			}
		}
	}
}

void real_flow::FillZDC(Int_t cenrality_bin, Float_t *ZDC_energy_mpd)
{
	Bool_t centGood;
	if (cenrality_bin==-1) 
		centGood=kFALSE;
	else 
		centGood=kTRUE;
	Double_t qx,qy;
	for (Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
		Double_t psi_N_R = GetPsiHalfZdc(ZDC_energy_mpd, 0, i_harm + 1, qx, qy);
		Double_t psi_N_L = GetPsiHalfZdc(ZDC_energy_mpd, 1, i_harm + 1, qx, qy);
		//if(i_harm == 0) h_psi_RP[cenrality_bin]->Fill(psi_N_R*180./Pi());
		if (centGood) h_delta_psi[cenrality_bin][i_harm][3]->Fill(Abs(Unfold(psi_N_R, psi_N_L,i_harm+1)*180./Pi()));

		Double_t psi_N_FULL = GetPsiFullZdc(ZDC_energy_mpd, i_harm+1);
		if (centGood) h_delta_EP_full_RP[cenrality_bin][i_harm][3]->Fill(Abs(Unfold(phiEP_mc,psi_N_FULL,i_harm+1)*180./Pi()));

		for (Int_t _harm=0;_harm<_N_HARM;_harm++){
			Double_t psi_N_R = GetPsiHalfZdc(ZDC_energy_mpd, 0, _harm + 1, qx, qy);
			Double_t psi_N_L = GetPsiHalfZdc(ZDC_energy_mpd, 1, _harm + 1, qx, qy);
			p_Res2Psi_vs_b[i_harm][_harm][3]->Fill(b_mc, Cos((i_harm+1)*(psi_N_R - psi_N_L)));
			p_ResAPsi_vs_b[i_harm][_harm][3]->Fill(b_mc, Cos((i_harm+1)*(psi_N_R - phiEP_mc)));
			p_ResBPsi_vs_b[i_harm][_harm][3]->Fill(b_mc, Cos((i_harm+1)*(psi_N_L - phiEP_mc)));

			Double_t psi_N_FULL = GetPsiFullZdc(ZDC_energy_mpd, _harm+1);
			p_true_Res_vs_b[i_harm][_harm][3]->Fill(b_mc,Cos((i_harm+1)*(psi_N_FULL - phiEP_mc)));
			for (Int_t arm = 0; arm < _N_ARM; ++arm)
			{
				Double_t psi_N_HALF = GetPsiHalfZdc(ZDC_energy_mpd, arm, _harm + 1, qx, qy);
				p_true_Res_half_vs_b[arm][i_harm][_harm][3]->Fill(b_mc,Cos((i_harm+1)*(psi_N_HALF - phiEP_mc)));
			}
		}
		for (Int_t arm = 0; arm < _N_ARM; ++arm)
		{
			Double_t psi_N_HALF = GetPsiHalfZdc(ZDC_energy_mpd, arm, i_harm + 1, qx, qy);
			if (centGood) h_qx[cenrality_bin][arm][i_harm][3]->Fill(qx);
			if (centGood) h_qy[cenrality_bin][arm][i_harm][3]->Fill(qy);
			if (centGood) p_qx_vs_b[arm][i_harm][3]->Fill(centralityBinsRes[cenrality_bin] + 0.1, qx);
			if (centGood) p_qy_vs_b[arm][i_harm][3]->Fill(centralityBinsRes[cenrality_bin], qy);
			if (centGood) h_delta_EP_half_RP[cenrality_bin][arm][i_harm][3]->Fill(psi_N_HALF - phiEP_mc);
			if (centGood) h_cor_EP_RP[cenrality_bin][arm][i_harm][3]->Fill(psi_N_HALF,phiEP_mc);

			if (centGood) h_psi[cenrality_bin][arm][i_harm][3]->Fill(psi_N_HALF);
		}
	}
}

void real_flow::Fill3Sub(Int_t gap, TVector3 *part, Long64_t ntracks,Float_t *ZDC_energy_mpd)
{
	if (gap<0 || gap>2){ 
		cout <<"gap index should be from 0 to 2 !" << endl;
	}
	else{
		Double_t qx, qy;
		Double_t psi_N_HALF;
		phiEP_mc = ATan2(Sin(phiEP_mc), Cos(phiEP_mc));

		for (Int_t harm = 0; harm < _N_HARM; ++harm)
		{
			for (Int_t _harm = 0; _harm < _N_HARM; ++_harm)
			{
				Double_t psi_N_R = GetPsiHalfTpc(gap,0,part,ntracks, 1, _harm+1, qx, qy);
				Double_t psi_N_L = GetPsiHalfTpc(gap,0,part,ntracks, -1, _harm+1, qx, qy);
				Double_t psi_N_F = GetPsiFullZdc(ZDC_energy_mpd, _harm+1);
				//cout << "Filling with centrality = " << centralityBinsRes[cenrality_bin]+0.1 << "Res2 = " << Cos((harm+1)*(psi_N_R - psi_N_L)) << endl;
				p_Cos_A_vs_b[harm][_harm][gap]->Fill(b_mc, Cos((harm+1)*Abs(Unfold(psi_N_F, psi_N_R,_harm+1))));
				p_Cos_B_vs_b[harm][_harm][gap]->Fill(b_mc, Cos((harm+1)*Abs(Unfold(psi_N_F, psi_N_L,_harm+1))));
				p_Cos_C_vs_b[harm][_harm][gap]->Fill(b_mc, Cos((harm+1)*Abs(Unfold(psi_N_R, psi_N_L,_harm+1))));

			}
		}
	}
}

Double_t real_flow::GetPsiHalfTpc(Int_t gap, Int_t weight, TVector3* part, Long64_t ntracks, Int_t sign, Int_t harm, Double_t &Qx, Double_t &Qy)
{
	Double_t Qcos=0., Qsin=0.;
	GetQsTpc(gap,weight,part,ntracks,sign,harm,Qcos,Qsin);
	Qx = Qcos;
	Qy = Qsin;
	Double_t PsiEP = (1/(Double_t)harm) * ATan2(Qsin,Qcos);
	if (harm==1 && sign==1) PsiEP-=Pi();
	return PsiEP;
}

Double_t real_flow::GetPsiHalfZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t n, Double_t &qx, Double_t &qy)
{
	Double_t Qcos=0., Qsin=0.;
	GetQsZdc(zdc_energy,zdc_ID,n,Qcos,Qsin);

	qx = Qcos;
	qy = Qsin;
	Double_t PsiEP = (1/(Double_t)n) * ATan2(Qsin,Qcos);
	return PsiEP;
}

void real_flow::InitHisto(){
	for (Int_t harm=0;harm<_N_HARM;harm++){
	for (Int_t side=0; side<_N_ARM;side++){
		for (Int_t method=0;method<_N_METHOD-1;method++){
			sprintf(name,"h_Mult_vs_eta%i%i%i",harm,side,method);
			sprintf(title,"(%i) #eta (%s) distribution for %s;|#eta|;",harm+1,arm_names[side].Data(),etaGap_name[method].Data());
			h_Mult_vs_eta[harm][side][method] = new TH1F(name,title,44,-2.2,2.2); h_Mult_vs_eta[harm][side][method]->Sumw2();
		}
	}
	}
	/*for (Int_t side=0; side<_N_ARM;side++){
		for (Int_t method=0;method<_N_METHOD-1;method++){
			sprintf(name,"p_Mult_vs_eta%i%i",side,method);
			sprintf(title,"<#eta> vs Multiplicity(%s) for %s;N_{mult};",arm_names[side].Data(),etaGap_name[method].Data());
			p_Mult_vs_eta[side][method] = new TProfile(name,title,44,-2.2,2.2);
		}
	}*/
	for (Int_t method=0;method<_N_METHOD-1;method++){
		sprintf(name,"h_Ratio_eta%i",method);
		sprintf(title,"Ratio #eta distribution Left/Right for %s;|#eta|;",etaGap_name[method].Data());
		h_Ratio_eta[method] = new TH1F(name,title,44,-2.2,2.2); h_Ratio_eta[method]->Sumw2();
	}

	for (Int_t sort = 0; sort < _N_SORTS; ++sort){
		for (int harm = 0; harm < _N_HARM; ++harm){
			for (int centralityBin = 0; centralityBin < _N_CENTRALITY_BINS; ++centralityBin){
				sprintf(name,"p_flow_wrt_RP_vs_pt%i%i%i",sort,centralityBin,harm);
				sprintf(title,"v_{%i} wrt RP for %s at %.2f < b < %.2f;p_{t}, GeV/c;",harm+1,sorts_of_particles[sort].Data(),centralityBinsFlow[centralityBin],centralityBinsFlow[centralityBin+1]);
				p_flow_wrt_RP_vs_pt[sort][centralityBin][harm] = new TProfile(name,title,Nptbins,ptbin);
				sprintf(name,"p_flow_wrt_RP_vs_rapidity%i%i%i",sort,centralityBin,harm);
				sprintf(title,"v_{%i} wrt RP for %s at %.2f < b < %.2f;y;",harm+1,sorts_of_particles[sort].Data(),centralityBinsFlow[centralityBin],centralityBinsFlow[centralityBin+1]);
				p_flow_wrt_RP_vs_rapidity[sort][centralityBin][harm] = new TProfile(name,title,NrapidityBins,rapidityBins);
				sprintf(name,"p_flow_wrt_RP_vs_eta%i%i%i",sort,centralityBin,harm);
				sprintf(title,"v_{%i} wrt RP for %s at %.2f < b < %.2f;#eta;",harm+1,sorts_of_particles[sort].Data(),centralityBinsFlow[centralityBin],centralityBinsFlow[centralityBin+1]);
				p_flow_wrt_RP_vs_eta[sort][centralityBin][harm] = new TProfile(name,title,NetaBins,etaBins);
			}
		}
	}

   	for (Int_t sort = 0; sort < _N_SORTS; ++sort){
      for (int harm = 0; harm < _N_HARM; ++harm){
         for (int _harm = 0; _harm < _N_HARM; ++ _harm){
            for (int method = 0; method < _N_METHOD; ++method){
               for (int centralityBin = 0; centralityBin < _N_CENTRALITY_BINS; ++centralityBin){
                  	sprintf(name,"p_flow_wrt_full_vs_pt%i%i%i%i%i",sort,centralityBin,harm,_harm,method);
               		sprintf(title,"v_{%i} wrt #Psi_{%i,%s}^{FULL} for %s at %.2f < b < %.2f;p_{t}, GeV/c;",harm+1,_harm+1,methods_names[method].Data(),sorts_of_particles[sort].Data(),centralityBinsFlow[centralityBin], centralityBinsFlow[centralityBin+1]);
               		p_flow_wrt_full_vs_pt[sort][centralityBin][harm][_harm][method] = new TProfile(name,title,Nptbins,ptbin);
               		sprintf(name,"p_flow_wrt_full_vs_rapidity%i%i%i%i%i",sort,centralityBin,harm,_harm,method);
               		sprintf(title,"v_{%i} wrt #Psi_{%i,%s}^{FULL} for %s at %.2f < b < %.2f;y;",harm+1,_harm+1,methods_names[method].Data(),sorts_of_particles[sort].Data(),centralityBinsFlow[centralityBin], centralityBinsFlow[centralityBin+1]);
               		p_flow_wrt_full_vs_rapidity[sort][centralityBin][harm][_harm][method] = new TProfile(name,title,NrapidityBins,rapidityBins);
               		sprintf(name,"p_flow_wrt_full_vs_eta%i%i%i%i%i",sort,centralityBin,harm,_harm,method);
               		sprintf(title,"v_{%i} wrt #Psi_{%i,%s}^{FULL} for %s at %.2f < b < %.2f;#eta;",harm+1,_harm+1,methods_names[method].Data(),sorts_of_particles[sort].Data(),centralityBinsFlow[centralityBin], centralityBinsFlow[centralityBin+1]);
               		p_flow_wrt_full_vs_eta[sort][centralityBin][harm][_harm][method] = new TProfile(name,title,NetaBins,etaBins);

               }

            }
         }
      }
    }
    for (Int_t harm=0;harm<_N_HARM;harm++){
    	for (Int_t _harm=0;_harm<_N_HARM;_harm++){
    		for (Int_t method=0;method<_N_METHOD;method++){
    			sprintf(name,"p_true_Res_vs_b%i%i%i",harm,_harm,method);
				sprintf(title,"%s%i%s%i%s%s%s","<cos(",harm+1,"(#Psi_{",_harm + 1,"}^{FULL} - #Psi_{RP}))> using ",
						methods_names[method].Data(),";b,fm;");
				p_true_Res_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);

				sprintf(name,"%s%i%i%i","p_Res2Psi_vs_b",harm,_harm,method);
				sprintf(title,"%s%i%s%i%s%s%s%i%s%s%s","cos(",harm + 1,"(#Psi_{",_harm + 1,",",methods_names[method].Data(),"}^{R} - #Psi_{",
						_harm + 1,",",methods_names[method].Data(),"}^{L}));b,fm;");
				p_Res2Psi_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
				sprintf(name,"%s%i%i%i","p_ResAPsi_vs_b",harm,_harm,method);
				sprintf(title,"%s%i%s%i%s%s%s%i%s%s%s","cos(",harm + 1,"(#Psi_{",_harm + 1,",",methods_names[method].Data(),"}^{R} - #Psi_{",
						_harm + 1,",",methods_names[method].Data(),"}^{L}));b,fm;");
				p_ResAPsi_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
				sprintf(name,"%s%i%i%i","p_ResBPsi_vs_b",harm,_harm,method);
				sprintf(title,"%s%i%s%i%s%s%s%i%s%s%s","cos(",harm + 1,"(#Psi_{",_harm + 1,",",methods_names[method].Data(),"}^{R} - #Psi_{",
						_harm + 1,",",methods_names[method].Data(),"}^{L}));b,fm;");
				p_ResBPsi_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);

				//sprintf(name,"p_ResPsi_vs_b%i%i%i",harm,_harm,method);
				//sprintf(title,"Res(#Psi_{%i,%s}) for v_{%i};b,fm;",_harm+1,methods_names[method].Data(),harm+1);
				//p_ResPsi_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
				for (Int_t arm=0;arm<_N_ARM;arm++){
					sprintf(name,"p_true_Res_half_vs_b%i%i%i%i",arm,harm,_harm,method);
					sprintf(title,"<cos(%i(#Psi_{%i,%s}^{%s} - #Psi_{RP} ))> ;b,fm;",harm+1,_harm+1,methods_names[method].Data(),arm_names[arm].Data());
					p_true_Res_half_vs_b[arm][harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
				}
    		}
    		for (Int_t method=0;method<_N_METHOD-1;method++){
    			sprintf(name,"p_Cos_A_vs_b%i%i%i",harm,_harm,method);
				sprintf(title,"%s%i%s%i%s%i%s%s%s","<cos(",harm+1,"(#Psi_{",_harm + 1,"}^{A} - #Psi_{",_harm+1,"}^{B}))> (",eta_method[method].Data(),");b,fm;");
				p_Cos_A_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
				sprintf(name,"p_Cos_B_vs_b%i%i%i",harm,_harm,method);
				sprintf(title,"%s%i%s%i%s%i%s%s%s","<cos(",harm+1,"(#Psi_{",_harm + 1,"}^{A} - #Psi_{",_harm+1,"}^{C}))> (",eta_method[method].Data(),");b,fm;");
				p_Cos_B_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
				sprintf(name,"p_Cos_C_vs_b%i%i%i",harm,_harm,method);
				sprintf(title,"%s%i%s%i%s%i%s%s%s","<cos(",harm+1,"(#Psi_{",_harm + 1,"}^{B} - #Psi_{",_harm+1,"}^{C}))> (",eta_method[method].Data(),");b,fm;");
				p_Cos_C_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
			}
    	}
    }
    for(Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
    	for (Int_t i_method=0;i_method<_N_METHOD;i_method++){
    		for (int centralityBin = 0; centralityBin < _N_CENTRALITY_BINS; ++centralityBin){
    			sprintf(name,"%s%i%i%i","h_delta_psi",centralityBin,i_harm,i_method);
				sprintf(title,"%s%i%s%i%s%.2f%s%.2f%s%s","#Psi_{",i_harm+1,"}^{R} - #Psi_{",i_harm+1,"}^{L} for ",centralityBinsFlow[centralityBin],"< b <",centralityBinsFlow[centralityBin+1]," using ",methods_names[i_method].Data());
				h_delta_psi[centralityBin][i_harm][i_method] = new TH1F(name,title,36,0.,180./(i_harm+1)); h_delta_psi[centralityBin][i_harm][i_method]->Sumw2();

				sprintf(name,"%s%i%i%i","h_delta_EP_full_RP",centralityBin,i_harm,i_method);
				sprintf(title,"%s%i%s%.2f%s%.2f%s%s","#Psi_{",i_harm+1,"}^{Full} - #Psi_{RP} for ",centralityBinsFlow[centralityBin],"< b <",centralityBinsFlow[centralityBin+1]," using ",methods_names[i_method].Data());
				h_delta_EP_full_RP[centralityBin][i_harm][i_method] = new TH1F(name,title,36,0.,180./(i_harm+1)); h_delta_EP_full_RP[centralityBin][i_harm][i_method]->Sumw2();
    		}
    	}
    }
    for (int centralityBin = 0; centralityBin < _N_CENTRALITY_BINS; ++centralityBin){
    	sprintf(name,"%s%i","h_psi_RP",centralityBin);
		sprintf(title,"%s%.2f%s%.2f%s","#Psi_{RP} for ",centralityBinsFlow[centralityBin],"< b <",centralityBinsFlow[centralityBin+1]," using TPC");
		h_psi_RP[centralityBin] = new TH1F(name,title,36,-180.,180.); h_psi_RP[centralityBin]->Sumw2();
    }
    for(Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
    	for (Int_t i_method=0;i_method<_N_METHOD;i_method++){
    		for (int centralityBin = 0; centralityBin < NcentralityBinsRes; ++centralityBin){
    			sprintf(name,"%s%i%i%i","h_delta_fit",centralityBin,i_harm,i_method);
				sprintf(title,"%s%i%s%i%s%.2f%s%.2f%s%s","#Psi_{",i_harm+1,"}^{R} - #Psi_{",i_harm+1,"}^{L} for ",centralityBinsRes[centralityBin],"< b <",centralityBinsRes[centralityBin+1]," using ",methods_names[i_method].Data());
				h_delta_fit[centralityBin][i_harm][i_method] = new TH1F(name,title,36,0.,180./(i_harm+1)); h_delta_fit[centralityBin][i_harm][i_method]->Sumw2();
			}
		}
	}

    for (Int_t i_arm=0;i_arm<_N_ARM;i_arm++){
	  	for (Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
	  		for (Int_t method=0;method<_N_METHOD;method++){
	  			for (Int_t i_proj=0;i_proj<2;i_proj++){
	  				sprintf(name,"h_Qvsb%i%i%i%i",i_arm,i_harm,method,i_proj);
               		sprintf(title,"<Q_{%i,%s}^{%s}> (%s) vs b, fm",i_harm+1,proj_name[i_proj].Data(),arm_names[i_arm].Data(),methods_names[method].Data());
	  				h_Qvsb[i_arm][i_harm][method][i_proj]    = new TProfile(name,title,5,centralityBinsRes);
	  				sprintf(name,"h_QvsMult%i%i%i%i",i_arm,i_harm,method,i_proj);
               		sprintf(title,"<Q_{%i,%s}^{%s}> (%s) vs multiplicity from TPC;b, fm;",i_harm+1,proj_name[i_proj].Data(),arm_names[i_arm].Data(),methods_names[method].Data());
	  				h_QvsMult[i_arm][i_harm][method][i_proj] = new TProfile(name,title,30,18.,900.);
	  			}
	  		}
	  	}
	}
	for (int centralityBin = 0; centralityBin < NcentralityBinsRes; ++centralityBin){
		for (int arm = 0; arm < _N_ARM; ++arm){
			for (int harm = 0; harm < _N_HARM; ++harm){
				for (int method = 0; method < _N_METHOD; ++method){
					sprintf (name,"%s%i%i%i%i","h_qx",centralityBin,harm,arm,method);
					sprintf (title,"%s%i%s%s%s%.2f%s%.2f%s%s","Q_{x,",harm+1,"}^{",arm_names[arm].Data(),"} for ",centralityBinsRes[centralityBin],"< b <",centralityBinsRes[centralityBin+1]," using ",methods_names[method].Data());
					h_qx[centralityBin][arm][harm][method] = new TH1F(name,title,100,-1.5,1.5);

					sprintf (name,"%s%i%i%i%i","h_qy",centralityBin,harm,arm,method);
					sprintf (title,"%s%i%s%s%s%.2f%s%.2f%s%s","Q_{y,",harm+1,"}^{",arm_names[arm].Data(),"} for ",centralityBinsRes[centralityBin],"< b <",centralityBinsRes[centralityBin+1]," using ",methods_names[method].Data());
					h_qy[centralityBin][arm][harm][method] = new TH1F(name,title,100,-1.5,1.5);

					sprintf(name,"%s%i%i%i%i","h_delta_EP_half_RP",centralityBin,arm,harm,method);
					sprintf(title,"%s%i%s%s%s%.2f%s%.2f%s%s","#Psi_{",harm+1,"}^{",arm_names[arm].Data(),"} - #Psi_{RP} for ",centralityBinsRes[centralityBin],
							" < b < ",centralityBinsRes[centralityBin+1], " using ",methods_names[method].Data());
					h_delta_EP_half_RP[centralityBin][arm][harm][method] = new TH1F(name,title,100,-4.,4.);

					sprintf(name,"%s%i%i%i%i","h_cor_EP_RP",centralityBin,arm,harm,method);
					sprintf(title,"%s%i%s%s%s%.2f%s%.2f%s%s","#Psi_{",harm+1,"}^{",arm_names[arm].Data(),"}, #Psi_{RP} corr. for ",centralityBinsRes[centralityBin],
							" < b < ",centralityBinsRes[centralityBin+1], " using ",methods_names[method].Data());
					h_cor_EP_RP[centralityBin][arm][harm][method] = new TH2F(name,title,100,-4.,4.,100,-4.,4.);

					sprintf (name,"%s%i%i%i%i","h_psi",centralityBin,harm,arm,method);
					sprintf (title,"%s%i%s%s%s%.2f%s%.2f%s%s","#Psi_{",harm+1,"}^{",arm_names[arm].Data(),"} for ",centralityBinsRes[centralityBin],"< b <",centralityBinsRes[centralityBin+1]," using ",methods_names[method].Data());
					h_psi[centralityBin][arm][harm][method] = new TH1F(name,title,100,-4.,4.);
				}
			}
		}
	}
	for (Int_t arm=0;arm<_N_ARM;arm++){
		for(Int_t harm=0;harm<_N_HARM;harm++){
			for(Int_t method=0;method<_N_METHOD;method++){
				sprintf(name,"%s%i%i%i","p_qx_vs_b",arm,harm,method);
				sprintf(title,"%s%i%s%s%s%s","Q_{x}^{",harm+1,",",arm_names[arm].Data(),"} using ",methods_names[method].Data());
				p_qx_vs_b[arm][harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);

				sprintf(name,"%s%i%i%i","p_qy_vs_b",arm,harm,method);
				sprintf(title,"%s%i%s%s%s%s","Q_{y}^{",harm+1,",",arm_names[arm].Data(),"} using ",methods_names[method].Data());
				p_qy_vs_b[arm][harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
			}
		}
	}
	for (Int_t harm=0; harm<_N_HARM; harm++){
		for (Int_t _harm=0;_harm<_N_HARM;_harm++){
			for (Int_t method=0;method<_N_METHOD;method++){
				sprintf(name,"%s%i%i%i","h_Res_true",harm,_harm,method);
				sprintf(title,"%s_{%i}%s%i%s%s%s%i%s%s%s","Res",harm + 1,"{#Psi_{",_harm + 1,",",methods_names[method].Data(),"}^{R},#Psi_{",
						_harm + 1,",",methods_names[method].Data(),"}^{L}};b,fm;");
				h_Res_true[harm][_harm][method] = new TH1F(name,title,NcentralityBinsRes,centralityBinsRes);
				sprintf(name,"%s%i%i%i","h_Res3_true",harm,_harm,method);
				sprintf(title,"%s_{%i}%s%i%s%s%s","Res ",harm + 1," from 3 sub-events of #Psi_{",_harm + 1,"}",methods_names[method].Data(),";b,fm;");
				h_Res3_true[harm][_harm][method] = new TH1F(name,title,NcentralityBinsRes,centralityBinsRes);
			}
		}
	}
	for (Int_t iC=0;iC<NcentralityBinsRes;iC++){
		h_Res1_fit[iC] = new TH1F(Form("h_Res1_fit%i",iC),Form("Resolution Res_{%i}^{FHCal} given from fit",iC),NcentralityBinsRes,centralityBinsRes);
	}
}

void real_flow::Write(){
	f_out = new TFile(outFile.Data(),"recreate");
   	if(!f_out->cd()) cout << "Can not open output file" << endl;
	else{
		f_out->cd();
		for (Int_t i_sort=0;i_sort<_N_SORTS;i_sort++){	
			for (int harm = 0; harm < _N_HARM; ++harm){
   				for (int _harm = 0; _harm < _N_HARM; ++ _harm){
   					for (int method = 0; method < _N_METHOD; ++method){
   						for (int cenralityBin = 0; cenralityBin < _N_CENTRALITY_BINS; ++cenralityBin){
   							p_flow_wrt_full_vs_pt[i_sort][cenralityBin][harm][_harm][method]->Write();
   							p_flow_wrt_full_vs_rapidity[i_sort][cenralityBin][harm][_harm][method]->Write();
   							p_flow_wrt_full_vs_eta[i_sort][cenralityBin][harm][_harm][method]->Write();
   						}
   					}
   				}
   			}
		}
		for (Int_t i_arm=0;i_arm<_N_ARM;i_arm++){
	  		for (Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
	  			for (Int_t method=0;method<_N_METHOD;method++){
	  				for (Int_t i_proj=0;i_proj<2;i_proj++){
	  					h_Qvsb[i_arm][i_harm][method][i_proj]   ->Write();
	  					h_QvsMult[i_arm][i_harm][method][i_proj]->Write();
	  				}
	  			}
	  		}
	  	}
	  	for (Int_t sort = 0; sort < _N_SORTS; ++sort){
			for (int harm = 0; harm < _N_HARM; ++harm){
				for (int cenralityBin = 0; cenralityBin < _N_CENTRALITY_BINS; ++cenralityBin){
					p_flow_wrt_RP_vs_pt[sort][cenralityBin][harm]->Write();
					p_flow_wrt_RP_vs_rapidity[sort][cenralityBin][harm]->Write();
					p_flow_wrt_RP_vs_eta[sort][cenralityBin][harm]->Write();
				}
			}
		}
		for (int harm = 0; harm < _N_HARM; ++harm){
   			for (int _harm = 0; _harm < _N_HARM; ++ _harm){
   				for (int method = 0; method < _N_METHOD; ++method){
   					p_Res2Psi_vs_b[harm][_harm][method]->Write();
   					p_ResAPsi_vs_b[harm][_harm][method]->Write();
   					p_ResBPsi_vs_b[harm][_harm][method]->Write();
   					p_true_Res_vs_b[harm][_harm][method]->Write();
   				}
   				for (Int_t method=0;method<_N_METHOD-1;method++){
    			p_Cos_A_vs_b[harm][_harm][method]->Write();
				p_Cos_B_vs_b[harm][_harm][method]->Write();
				p_Cos_C_vs_b[harm][_harm][method]->Write();
			}
   			}
   		}
   		for (int centralityBin = 0; centralityBin < NcentralityBinsRes; ++centralityBin){
			for (int arm = 0; arm < _N_ARM; ++arm){
				for (int harm = 0; harm < _N_HARM; ++harm){
					for (int method = 0; method < _N_METHOD; ++method){
						h_qx[centralityBin][arm][harm][method]->Write();
						h_qy[centralityBin][arm][harm][method]->Write();
						h_delta_EP_half_RP[centralityBin][arm][harm][method]->Write();
						h_cor_EP_RP[centralityBin][arm][harm][method]->Write();
						h_psi[centralityBin][arm][harm][method]->Write();
					}
				}
			}
		}
		for (Int_t arm=0;arm<_N_ARM;arm++){
			for(Int_t harm=0;harm<_N_HARM;harm++){
				for(Int_t method=0;method<_N_METHOD;method++){
					p_qx_vs_b[arm][harm][method]->Write();
					p_qy_vs_b[arm][harm][method]->Write();
				}
			}
		}
		for(Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
    		for (Int_t i_method=0;i_method<_N_METHOD;i_method++){
    			for (int centralityBin = 0; centralityBin < _N_CENTRALITY_BINS; ++centralityBin){
    				h_delta_psi[centralityBin][i_harm][i_method]->Write();
    				h_delta_EP_full_RP[centralityBin][i_harm][i_method]->Write();
    			}
    		}
    	}
    	for (int centralityBin = 0; centralityBin < _N_CENTRALITY_BINS; ++centralityBin){
    		h_psi_RP[centralityBin]->Write();
    	}
    	for (Int_t harm=0; harm<_N_HARM; harm++){
			for (Int_t _harm=0;_harm<_N_HARM;_harm++){
				for (Int_t method=0;method<_N_METHOD;method++){
					h_Res_true[harm][_harm][method]->Write();
					h_Res3_true[harm][_harm][method]->Write();
				}
			}
		}
		for(Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
    		for (Int_t i_method=0;i_method<_N_METHOD;i_method++){
    			for (int centralityBin = 0; centralityBin < NcentralityBinsRes; ++centralityBin){
    				h_delta_fit[centralityBin][i_harm][i_method]->Write();
				}
			}
		}
		for (Int_t iC=0;iC<NcentralityBinsRes;iC++){
			h_Res1_fit[iC]->Write();
		}
		for (Int_t harm=0;harm<_N_HARM;harm++){
		for (Int_t side=0; side<_N_ARM;side++){
			for (Int_t method=0;method<_N_METHOD-1;method++){
				h_Mult_vs_eta[harm][side][method]->Write();
				//p_Mult_vs_eta[side][method]->Write();
			}
		}
		}
		for (Int_t method=0;method<_N_METHOD-1;method++){
			h_Ratio_eta[method]->Write();
		}
	}
	f_out->Close();
}
void real_flow::DrawResFit(){
	for (Int_t iC=0;iC<NcentralityBinsRes;iC++){
		//h_Res1_fit[iC]->SetBinContent(iC+1,Res1Fit[iC]);
	}
}
void real_flow::MakeTTree(TString outputTree, Double_t Qx[_N_HARM][_N_ARM][_N_METHOD], Double_t Qy[_N_HARM][_N_ARM][_N_METHOD], Int_t Qw[_N_HARM][_N_ARM][_N_METHOD])
//Do before event loop
{	
	/*out_tree = new TTree(outputTree.Data(),outputTree.Data());
	
	//Int_t Ndet = 4;
	//Int_t Nside = 2;
	//Int_t Nharm = 2;

	//const TString side_name[Nside]={TString("Left"),TString("Right")};
	//const TString det_name[Ndet]  ={TString("FHCal"),TString("TPC-eta-0"),TString("TPC-eta-0_2"),TString("TPC-eta-0_5")};

	out_tree->Branch("b_mc", &b_mc, "b_mc/F");
	out_tree->Branch("phiEP_mc", &phiEP_mc, "phiEP_mc/F");
	out_tree->Branch("Qx_1_Left_FHCal",&Qx[0][0][0],"Qx_1_Left_FHCal/D");
	out_tree->Branch("Qx_1_Left_TPC-eta-0",&Qx[0][0][1],"Qx_1_Left_TPC-eta-0/D");
	out_tree->Branch("Qx_1_Left_TPC-eta-0_2",&Qx[0][0][2],"Qx_1_Left_TPC-eta-0_2/D");
	out_tree->Branch("Qx_1_Left_TPC-eta-0_5",&Qx[0][0][3],"Qx_1_Left_TPC-eta-0_5/D");
	out_tree->Branch("Qx_1_Right_FHCal",&Qx[0][1][0],"Qx_1_Right_FHCal/D");
	out_tree->Branch("Qx_1_Right_TPC-eta-0",&Qx[0][1][1],"Qx_1_Right_TPC-eta-0/D");
	out_tree->Branch("Qx_1_Right_TPC-eta-0_2",&Qx[0][1][2],"Qx_1_Right_TPC-eta-0_2/D");
	out_tree->Branch("Qx_1_Right_TPC-eta-0_5",&Qx[0][1][3],"Qx_1_Right_TPC-eta-0_5/D");
	out_tree->Branch("Qx_2_Left_FHCal",&Qx[1][0][0],"Qx_2_Left_FHCal/D");
	out_tree->Branch("Qx_2_Left_TPC-eta-0",&Qx[1][0][1],"Qx_2_Left_TPC-eta-0/D");
	out_tree->Branch("Qx_2_Left_TPC-eta-0_2",&Qx[1][0][2],"Qx_2_Left_TPC-eta-0_2/D");
	out_tree->Branch("Qx_2_Left_TPC-eta-0_5",&Qx[1][0][3],"Qx_2_Left_TPC-eta-0_5/D");
	out_tree->Branch("Qx_2_Right_FHCal",&Qx[1][1][0],"Qx_2_Right_FHCal/D");
	out_tree->Branch("Qx_2_Right_TPC-eta-0",&Qx[1][1][1],"Qx_2_Right_TPC-eta-0/D");
	out_tree->Branch("Qx_2_Right_TPC-eta-0_2",&Qx[1][1][2],"Qx_2_Right_TPC-eta-0_2/D");
	out_tree->Branch("Qx_2_Right_TPC-eta-0_5",&Qx[1][1][3],"Qx_2_Right_TPC-eta-0_5/D");
	out_tree->Branch("Qy_1_Left_FHCal",&Qy[0][0][0],"Qy_1_Left_FHCal/D");
	out_tree->Branch("Qy_1_Left_TPC-eta-0",&Qy[0][0][1],"Qy_1_Left_TPC-eta-0/D");
	out_tree->Branch("Qy_1_Left_TPC-eta-0_2",&Qy[0][0][2],"Qy_1_Left_TPC-eta-0_2/D");
	out_tree->Branch("Qy_1_Left_TPC-eta-0_5",&Qy[0][0][3],"Qy_1_Left_TPC-eta-0_5/D");
	out_tree->Branch("Qy_1_Right_FHCal",&Qy[0][1][0],"Qy_1_Right_FHCal/D");
	out_tree->Branch("Qy_1_Right_TPC-eta-0",&Qy[0][1][1],"Qy_1_Right_TPC-eta-0/D");
	out_tree->Branch("Qy_1_Right_TPC-eta-0_2",&Qy[0][1][2],"Qy_1_Right_TPC-eta-0_2/D");
	out_tree->Branch("Qy_1_Right_TPC-eta-0_5",&Qy[0][1][3],"Qy_1_Right_TPC-eta-0_5/D");
	out_tree->Branch("Qy_2_Left_FHCal",&Qy[1][0][0],"Qy_2_Left_FHCal/D");
	out_tree->Branch("Qy_2_Left_TPC-eta-0",&Qy[1][0][1],"Qy_2_Left_TPC-eta-0/D");
	out_tree->Branch("Qy_2_Left_TPC-eta-0_2",&Qy[1][0][2],"Qy_2_Left_TPC-eta-0_2/D");
	out_tree->Branch("Qy_2_Left_TPC-eta-0_5",&Qy[1][0][3],"Qy_2_Left_TPC-eta-0_5/D");
	out_tree->Branch("Qy_2_Right_FHCal",&Qy[1][1][0],"Qy_2_Right_FHCal/D");
	out_tree->Branch("Qy_2_Right_TPC-eta-0",&Qy[1][1][1],"Qy_2_Right_TPC-eta-0/D");
	out_tree->Branch("Qy_2_Right_TPC-eta-0_2",&Qy[1][1][2],"Qy_2_Right_TPC-eta-0_2/D");
	out_tree->Branch("Qy_2_Right_TPC-eta-0_5",&Qy[1][1][3],"Qy_2_Right_TPC-eta-0_5/D");
	out_tree->Branch("Qw_1_Left_FHCal",&Qw[0][0][0],"Qw_1_Left_FHCal/I");
	out_tree->Branch("Qw_1_Left_TPC-eta-0",&Qw[0][0][1],"Qw_1_Left_TPC-eta-0/I");
	out_tree->Branch("Qw_1_Left_TPC-eta-0_2",&Qw[0][0][2],"Qw_1_Left_TPC-eta-0_2/I");
	out_tree->Branch("Qw_1_Left_TPC-eta-0_5",&Qw[0][0][3],"Qw_1_Left_TPC-eta-0_5/I");
	out_tree->Branch("Qw_1_Right_FHCal",&Qw[0][1][0],"Qw_1_Right_FHCal/I");
	out_tree->Branch("Qw_1_Right_TPC-eta-0",&Qw[0][1][1],"Qw_1_Right_TPC-eta-0/I");
	out_tree->Branch("Qw_1_Right_TPC-eta-0_2",&Qw[0][1][2],"Qw_1_Right_TPC-eta-0_2/I");
	out_tree->Branch("Qw_1_Right_TPC-eta-0_5",&Qw[0][1][3],"Qw_1_Right_TPC-eta-0_5/I");
	out_tree->Branch("Qw_2_Left_FHCal",&Qw[1][0][0],"Qw_2_Left_FHCal/I");
	out_tree->Branch("Qw_2_Left_TPC-eta-0",&Qw[1][0][1],"Qw_2_Left_TPC-eta-0/I");
	out_tree->Branch("Qw_2_Left_TPC-eta-0_2",&Qw[1][0][2],"Qw_2_Left_TPC-eta-0_2/I");
	out_tree->Branch("Qw_2_Left_TPC-eta-0_5",&Qw[1][0][3],"Qw_2_Left_TPC-eta-0_5/I");
	out_tree->Branch("Qw_2_Right_FHCal",&Qw[1][1][0],"Qw_2_Right_FHCal/I");
	out_tree->Branch("Qw_2_Right_TPC-eta-0",&Qw[1][1][1],"Qw_2_Right_TPC-eta-0/I");
	out_tree->Branch("Qw_2_Right_TPC-eta-0_2",&Qw[1][1][2],"Qw_2_Right_TPC-eta-0_2/I");
	out_tree->Branch("Qw_2_Right_TPC-eta-0_5",&Qw[1][1][3],"Qw_2_Right_TPC-eta-0_5/I");
	


	for (Int_t harm=0;harm<Nharm;harm++){
		for (Int_t side=0;side<Nside;side++){
			for (Int_t det=0;det<Ndet;det++){
				out_tree->Branch(Form("Qx_%i_%s_%s",harm+1,side_name[side].Data(),det_name[det].Data()), &Qx[harm][side][det],Form("Qx_%i_%s_%s/D",harm+1,side_name[side].Data(),det_name[det].Data()));
				out_tree->Branch(Form("Qy_%i_%s_%s",harm+1,side_name[side].Data(),det_name[det].Data()), &Qy[harm][side][det],Form("Qx_%i_%s_%s/D",harm+1,side_name[side].Data(),det_name[det].Data()));
				out_tree->Branch(Form("Qw_%i_%s_%s",harm+1,side_name[side].Data(),det_name[det].Data()), &Qw[harm][side][det],Form("Qx_%i_%s_%s/D",harm+1,side_name[side].Data(),det_name[det].Data()));
			}
		}
	}
		*/
}
void real_flow::FillTTree()//Do in event loop, after MakeTTree(...)
{
	out_tree->Fill();
}
void real_flow::SaveTTree(TString outputFile)//Do after event loop
{
	out_file = new TFile(outputFile.Data(),"recreate");
	out_file->cd();
	out_tree->Write();
	out_file->Close();
}
