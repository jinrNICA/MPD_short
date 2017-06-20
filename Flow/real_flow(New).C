
#define real_flow_cxx
#include "real_flow.h"

#define debug

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
   Long64_t nentries = fChain->GetEntriesFast();
   #elif defined debug
   Long64_t nentries = 10000;
   #endif

   Long64_t nbytes = 0, nb = 0;
   Double_t Q[_N_ARM][_N_HARM][2]; //Q[R,L][harm][x,y];
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      cout << "Event # " << jentry << endl;

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
      if (cenrality_bin==-1) continue;
      Int_t cenrality_bin_res = cenrality_bin;

      //skip the event if the ZDC signal is zero
	  if ((GetTotalEnergy(ZDC_energy_mpd, 0) == 0) || (GetTotalEnergy(ZDC_energy_mpd, 1) == 0)) continue;

	  FillZDC(cenrality_bin_res,ZDC_energy_mpd);

	  for (Int_t n=0;n<_N_HARM;n++){
	  	GetQsZdc(ZDC_energy_mpd,0,n+1,Q[0][n][0], Q[0][n][1]);
	  	GetQsZdc(ZDC_energy_mpd,1,n+1,Q[1][n][0], Q[1][n][1]);
	  }
	  for (Int_t i_arm=0;i_arm<_N_ARM;i_arm++){
	  	for (Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
	  		for (Int_t i_proj=0;i_proj<2;i_proj++){
	  			h_Qvsb[i_arm][i_harm][i_proj]   ->Fill(b_mc,Q[i_arm][i_harm][i_proj]);
	  			h_QvsMult[i_arm][i_harm][i_proj]->Fill(multiplicity,Q[i_arm][i_harm][i_proj]);
	  		}
	  	}
	  }

	  Double_t phi_EP = ATan2(Sin(phiEP_mc),Cos(phiEP_mc)); //unfold the generated event plane
	  Double_t Psi_1_FULL_ZDC = GetPsiFullZdc(ZDC_energy_mpd,1); //getting event plane using ZDC for 1st harmonic

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
		if ((Abs(Eta) > 0.2) && (Abs(Eta) < Cut_Eta_Max))
			p_flow_wrt_full_vs_pt[p_sort][cenrality_bin][0][0][1]->Fill(Pt,sign*Cos(Psi_1_FULL_ZDC - Phi));
		if (Abs(Eta) < Cut_Eta_Max)
			p_flow_wrt_full_vs_pt[p_sort][cenrality_bin][1][0][1]->Fill(Pt,Cos(2.*(Psi_1_FULL_ZDC - Phi)));

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
	  }//end of loop over tracks (MC);
	  //---------------------------------------------------------------------------------------------------------



   }//end of loop over events;

	Double_t Res[_N_HARM][_N_HARM][_N_METHOD][NcentralityBinsRes];
	Double_t eRes[_N_HARM][_N_HARM][_N_METHOD][NcentralityBinsRes];
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
				}
			}
		}
	}  

   Write("output.root");

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

Double_t real_flow::GetPsiFullZdc(Float_t* zdc_energy, Int_t n)
{
	Double_t QcosR=0., QsinR=0.;
	GetQsZdc(zdc_energy,0,n,QcosR, QsinR);

	Double_t QcosL=0., QsinL=0.;
	GetQsZdc(zdc_energy,1,n,QcosL,QsinL);

	Double_t psiEP = ATan2(QsinR + QsinL,QcosR + QcosL)/(Double_t) n; // (-pi/n,pi/n]
    return psiEP;
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

void real_flow::FillZDC(Int_t cenrality_bin, Float_t *ZDC_energy_mpd)
{
	Double_t qx,qy;
	for (Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
		Double_t psi_N_R = GetPsiHalfZdc(ZDC_energy_mpd, 0, i_harm + 1, qx, qy);
		Double_t psi_N_L = GetPsiHalfZdc(ZDC_energy_mpd, 1, i_harm + 1, qx, qy);
		h_delta_psi[cenrality_bin][i_harm][1]->Fill(Abs(Unfold(psi_N_R, psi_N_L,i_harm+1)));

		Double_t psi_N_FULL = GetPsiFullZdc(ZDC_energy_mpd, i_harm+1);

		for (Int_t _harm=0;_harm<_N_HARM;_harm++){
			Double_t psi_N_R = GetPsiHalfZdc(ZDC_energy_mpd, 0, _harm + 1, qx, qy);
			Double_t psi_N_L = GetPsiHalfZdc(ZDC_energy_mpd, 1, _harm + 1, qx, qy);
			p_Res2Psi_vs_b[i_harm][_harm][1]->Fill(b_mc, Cos((i_harm+1)*(psi_N_R - psi_N_L)));

			Double_t psi_N_FULL = GetPsiFullZdc(ZDC_energy_mpd, _harm+1);
			p_true_Res_vs_b[i_harm][_harm][1]->Fill(b_mc,Cos((i_harm+1)*(psi_N_FULL - phiEP_mc)));
			for (Int_t arm = 0; arm < _N_ARM; ++arm)
			{
				Double_t psi_N_HALF = GetPsiHalfZdc(ZDC_energy_mpd, arm, _harm + 1, qx, qy);
				p_true_Res_half_vs_b[arm][i_harm][_harm][1]->Fill(b_mc,Cos((i_harm+1)*(psi_N_HALF - phiEP_mc)));
			}
		}
		for (Int_t arm = 0; arm < _N_ARM; ++arm)
		{
			Double_t psi_N_HALF = GetPsiHalfZdc(ZDC_energy_mpd, arm, i_harm + 1, qx, qy);
			h_qx[cenrality_bin][arm][i_harm][1]->Fill(qx);
			h_qy[cenrality_bin][arm][i_harm][1]->Fill(qy);
			p_qx_vs_b[arm][i_harm][1]->Fill(centralityBinsRes[cenrality_bin] + 0.1, qx);
			p_qy_vs_b[arm][i_harm][1]->Fill(centralityBinsRes[cenrality_bin], qy);
			h_delta_EP_half_RP[cenrality_bin][arm][i_harm][1]->Fill(psi_N_HALF - phiEP_mc);

			h_psi[cenrality_bin][arm][i_harm][1]->Fill(psi_N_HALF);
		}
	}
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
	for (Int_t sort = 0; sort < _N_SORTS; ++sort){
		for (int harm = 0; harm < _N_HARM; ++harm){
			for (int centralityBin = 0; centralityBin < _N_CENTRALITY_BINS; ++centralityBin){
				sprintf(name,"p_flow_wrt_RP_vs_pt%i%i%i",sort,centralityBin,harm);
				sprintf(title,"v_{%i} wrt RP for %s at %.2f < b < %.2f;p_{t}, GeV/c;",harm+1,sorts_of_particles[sort].Data(),centralityBinsFlow[centralityBin],centralityBinsFlow[centralityBin+1]);
				p_flow_wrt_RP_vs_pt[sort][centralityBin][harm] = new TProfile(name,title,20,0.2,3.);
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
               		p_flow_wrt_full_vs_pt[sort][centralityBin][harm][_harm][method] = new TProfile(name,title,20,0.2,3.);
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
				p_true_Res_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,2.,12.);

				sprintf(name,"%s%i%i%i","p_Res2Psi_vs_b",harm,_harm,method);
				sprintf(title,"%s%i%s%i%s%s%s%i%s%s%s","cos(",harm + 1,"(#Psi_{",_harm + 1,",",methods_names[method].Data(),"}^{R} - #Psi_{",
						_harm + 1,",",methods_names[method].Data(),"}^{L});b,fm;");
				p_Res2Psi_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,2.,12.);

				//sprintf(name,"p_ResPsi_vs_b%i%i%i",harm,_harm,method);
				//sprintf(title,"Res(#Psi_{%i,%s}) for v_{%i};b,fm;",_harm+1,methods_names[method].Data(),harm+1);
				//p_ResPsi_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
				for (Int_t arm=0;arm<_N_ARM;arm++){
					sprintf(name,"p_true_Res_half_vs_b%i%i%i%i",arm,harm,_harm,method);
					sprintf(title,"<cos(%i(#Psi_{%i,%s}^{%s} - #Psi_{RP} ))> ;b,fm;",harm+1,_harm+1,methods_names[method].Data(),arm_names[arm].Data());
					p_true_Res_half_vs_b[arm][harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,2.,12.);
				}
    		}
    	}
    }
    for(Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
    	for (Int_t i_method=0;i_method<_N_METHOD;i_method++){
    		for (int centralityBin = 0; centralityBin < _N_CENTRALITY_BINS; ++centralityBin){
    			sprintf(name,"%s%i%i%i","h_delta_psi",centralityBin,i_harm,i_method);
				sprintf(title,"%s%i%s%i%s%.2f%s%.2f%s%s","#Psi_{",i_harm+1,"}^{R} - #Psi_{",i_harm+1,"}^{L} for ",centralityBinsFlow[centralityBin],"< b <",centralityBinsFlow[centralityBin+1]," using ",methods_names[i_method].Data());
				h_delta_psi[centralityBin][i_harm][i_method] = new TH1F(name,title,36,0.,Pi()/(i_harm+1)); h_delta_psi[centralityBin][i_harm][i_method]->Sumw2();
    		}
    	}
    }
    
    for (Int_t i_arm=0;i_arm<_N_ARM;i_arm++){
	  	for (Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
	  		for (Int_t i_proj=0;i_proj<2;i_proj++){
	  			sprintf(name,"h_Qvsb%i%i%i",i_arm,i_harm,i_proj);
               	sprintf(title,"<Q_{%i,%s}^{%s}> (FHCal) vs b, fm",i_harm+1,proj_name[i_proj].Data(),arm_names[i_arm].Data());
	  			h_Qvsb[i_arm][i_harm][i_proj]    = new TProfile(name,title,5,2.,12.);
	  			sprintf(name,"h_QvsMult%i%i%i",i_arm,i_harm,i_proj);
               	sprintf(title,"<Q_{%i,%s}^{%s}> (FHCal) vs multiplicity from TPC;b, fm;",i_harm+1,proj_name[i_proj].Data(),arm_names[i_arm].Data());
	  			h_QvsMult[i_arm][i_harm][i_proj] = new TProfile(name,title,30,18.,900.);
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
				sprintf(name,"%s%i%s%i%s%i%s","p_qx_vs_b[",arm,"][",harm,"][",method,"]");
				sprintf(title,"%s%i%s%s%s%s","Q_{x}^{",harm+1,",",arm_names[arm].Data(),"} using ",methods_names[method].Data());
				p_qx_vs_b[arm][harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);

				sprintf(name,"%s%i%s%i%s%i%s","p_qy_vs_b[",arm,"][",harm,"][",method,"]");
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
				h_Res_true[harm][_harm][method] = new TH1F(name,title,NcentralityBinsRes,2.,12.);
			}
		}
	}
}

void real_flow::Write(TString outputName){
	f_out = new TFile(outputName.Data(),"recreate");
   	if(!f_out->cd()) cout << "Can not open output file" << endl;
	else{
		f_out->cd();
		for (Int_t i_sort=0;i_sort<_N_SORTS;i_sort++){	
			for (int harm = 0; harm < _N_HARM; ++harm){
   				for (int _harm = 0; _harm < _N_HARM; ++ _harm){
   					for (int method = 0; method < _N_METHOD; ++method){
   						for (int cenralityBin = 0; cenralityBin < _N_CENTRALITY_BINS; ++cenralityBin){
   							p_flow_wrt_full_vs_pt[i_sort][cenralityBin][harm][_harm][method]->Write();
   						}
   					}
   				}
   			}
		}
		for (Int_t i_arm=0;i_arm<_N_ARM;i_arm++){
	  		for (Int_t i_harm=0;i_harm<_N_HARM;i_harm++){
	  			for (Int_t i_proj=0;i_proj<2;i_proj++){
	  				h_Qvsb[i_arm][i_harm][i_proj]   ->Write();
	  				h_QvsMult[i_arm][i_harm][i_proj]->Write();
	  			}
	  		}
	  	}
	  	for (Int_t sort = 0; sort < _N_SORTS; ++sort){
			for (int harm = 0; harm < _N_HARM; ++harm){
				for (int cenralityBin = 0; cenralityBin < _N_CENTRALITY_BINS; ++cenralityBin){
					p_flow_wrt_RP_vs_pt[sort][cenralityBin][harm]->Write();
				}
			}
		}
		for (int harm = 0; harm < _N_HARM; ++harm){
   			for (int _harm = 0; _harm < _N_HARM; ++ _harm){
   				for (int method = 0; method < _N_METHOD; ++method){
   					p_Res2Psi_vs_b[harm][_harm][method]->Write();
   					p_true_Res_vs_b[harm][_harm][method]->Write();
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
    			}
    		}
    	}
    	for (Int_t harm=0; harm<_N_HARM; harm++){
			for (Int_t _harm=0;_harm<_N_HARM;_harm++){
				for (Int_t method=0;method<_N_METHOD;method++){
					h_Res_true[harm][_harm][method]->Write();
				}
			}
		}
	}
}