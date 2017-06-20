/*#define _N_ARM 2
#define _N_HARM 2
#define _N_SORTS 3
#define _N_CENTRALITY_BINS 2
#define _N_METHOD 2*/

#include "variables.C"

void GetCalcRes(TString inFileName, TString inFitFileName, TString outFileName, Int_t accuracy){
	TFile* f_in = new TFile(inFileName.Data(),"read");
	TFile* f_fit= new TFile(inFitFileName.Data(),"read");
	TFile* f_out= new TFile(outFileName.Data(),"recreate");
	TProfile* p_Res2Psi_vs_b[_N_HARM][_N_HARM][_N_METHOD];
	TProfile* p_ResAPsi_vs_b[_N_HARM][_N_HARM][_N_METHOD];
	TProfile* p_ResBPsi_vs_b[_N_HARM][_N_HARM][_N_METHOD];
	TH1F*     h_Res_calc[_N_HARM][_N_HARM][_N_METHOD];
	TH1F*     h_Res3_calc[_N_HARM][_N_HARM][_N_METHOD-1];
	TH1F*     h_Res_out[_N_HARM][_N_HARM][_N_METHOD];
	TH1F*     h_Res3_out[_N_HARM][_N_HARM][_N_METHOD-1];
   	TH1F*     h_Res_true[_N_HARM][_N_HARM][_N_METHOD];
   	TH1F* 	  fitRes[_N_HARM][_N_METHOD];
   	TH1F*     ResRatio[_N_HARM][_N_METHOD];
   	TH1D*     h_proj_A[_N_HARM][_N_HARM][_N_METHOD];
   	TH1D*     h_proj_B[_N_HARM][_N_HARM][_N_METHOD];
   	TH1D*     h_proj_2[_N_HARM][_N_HARM][_N_METHOD];
   	TH1D*     h_rat[_N_HARM][_N_HARM][_N_METHOD];

	char name[200];
	char title[200];
	for (Int_t harm=0; harm<_N_HARM; harm++){
		fitRes[harm][3] = (TH1F*) f_fit->Get(Form("fitRes%i",harm));
		for (Int_t method=0;method<_N_METHOD;method++){
			fitRes[harm][method] = (TH1F*) f_fit->Get(Form("fitRes%i%i",harm,method));
			ResRatio[harm][method] = new TH1F(Form("ResRatio%i%i",harm,method),Form("Ratio <Cos(%i(#Psi_{%i}-#Psi_{RP}))> for %s",harm+1,harm+1,methods_names[method].Data()),NcentralityBinsRes,centralityBinsRes);
		}
		for (Int_t _harm=0;_harm<_N_HARM;_harm++){
			for (Int_t method=0;method<_N_METHOD;method++){
				h_Res3_calc[harm][_harm][method] = (TH1F*) f_in->Get(Form("h_Res3_true%i%i%i",harm,_harm,method));
				p_Res2Psi_vs_b[harm][_harm][method] = (TProfile*) f_in->Get(Form("p_Res2Psi_vs_b%i%i%i",harm,_harm,method));
				p_ResAPsi_vs_b[harm][_harm][method] = (TProfile*) f_in->Get(Form("p_ResAPsi_vs_b%i%i%i",harm,_harm,method));
				p_ResBPsi_vs_b[harm][_harm][method] = (TProfile*) f_in->Get(Form("p_ResBPsi_vs_b%i%i%i",harm,_harm,method));

				h_proj_A[harm][_harm][method] = new TH1D(Form("h_proj_A%i%i%i",harm,_harm,method),Form("h_proj_A%i%i%i",harm,_harm,method),10,2.,12.);
				h_proj_B[harm][_harm][method] = new TH1D(Form("h_proj_B%i%i%i",harm,_harm,method),Form("h_proj_B%i%i%i",harm,_harm,method),10,2.,12.);
				h_proj_2[harm][_harm][method] = new TH1D(Form("h_proj_2%i%i%i",harm,_harm,method),Form("<cos(%i(#Psi_{%i}^{R}-#Psi_{RP}))><cos(%i(#Psi_{%i}^{L}-#Psi_{RP}))> using %s",harm+1,_harm+1,harm+1,_harm+1,methods_names[method].Data()),NcentralityBinsRes,centralityBinsRes);
				h_rat[harm][_harm][method]    = new TH1D(Form("h_rat%i%i%i",harm,_harm,method),Form("<cos(%i(#Psi_{%i}^{R}-#Psi_{%i}^{L}))>/<cos(%i(#Psi_{%i}^{R}-#Psi_{RP}))><cos(%i(#Psi_{%i}^{L}-#Psi_{RP}))> using %s",harm+1,_harm+1,_harm+1,harm+1,_harm+1,harm+1,_harm+1,methods_names[method].Data()),NcentralityBinsRes,centralityBinsRes);

				sprintf(name,"%s%i%i%i","h_Res3_out",harm,_harm,method);
				sprintf(title,"%s_{%i}%s%i%s%s%s%i%s%s%s","(Calc 3 Sub) Res",harm + 1,"{#Psi_{",_harm + 1,",",methods_names[method].Data(),"}^{R},#Psi_{",_harm + 1,",",methods_names[method].Data(),"}^{L}};b,fm;");
				h_Res3_out[harm][_harm][method] = new TH1F(name,title,NcentralityBinsRes,centralityBinsRes);
      		}
			for (Int_t method=0;method<_N_METHOD;method++){
				h_Res_calc[harm][_harm][method] = (TH1F*) f_in->Get(Form("h_Res_true%i%i%i",harm,_harm,method));
                h_Res_true[harm][_harm][method] = (TH1F*) f_in->Get(Form("p_true_Res_vs_b%i%i%i",harm,_harm,method));

				sprintf(name,"%s%i%i%i","h_Res_out",harm,_harm,method);
				sprintf(title,"%s_{%i}%s%i%s%s%s%i%s%s%s","(Calc) Res",harm + 1,"{#Psi_{",_harm + 1,",",methods_names[method].Data(),"}^{R},#Psi_{",_harm + 1,",",methods_names[method].Data(),"}^{L}};b,fm;");
				h_Res_out[harm][_harm][method] = new TH1F(name,title,NcentralityBinsRes,centralityBinsRes);
			}
		}
	}

	for (Int_t harm=0;harm<_N_HARM;harm++){
		for (Int_t _harm=0;_harm<_N_HARM;_harm++){
			for (Int_t method=0;method<_N_METHOD;method++){
				h_proj_A[harm][_harm][method] = p_ResAPsi_vs_b[harm][_harm][method]->ProjectionX();
				h_proj_B[harm][_harm][method] = p_ResBPsi_vs_b[harm][_harm][method]->ProjectionX();
				h_proj_2[harm][_harm][method] ->Multiply(h_proj_A[harm][_harm][method],h_proj_B[harm][_harm][method]);
				h_rat[harm][_harm][method]    = p_Res2Psi_vs_b[harm][_harm][method]->ProjectionX();
				h_rat[harm][_harm][method]    ->Divide(h_proj_2[harm][_harm][method]);
			}
		}
	}

	// Fill calculated full resolution;
	for (Int_t method=0;method<_N_METHOD;method++){
		for (Int_t harm=0;harm<_N_HARM;harm++){
			//for (Int_t _harm=0;_harm<_N_HARM;_harm++){
				for (Int_t iC=0;iC<NcentralityBinsRes;iC++){
					Double_t res = h_Res_calc[harm][0][method]->GetBinContent(iC+1);
					Double_t chi = Chi(res,harm+1,accuracy);
					chi      *= TMath::Sqrt(2);
					Double_t true_res = ResEventPlane(chi,harm+1);
					cout << "harm = " << harm << "; _harm = " << 0 << "; method = " << method << "; b = " << iC << ": res2 = " << res << "; approx = " << true_res << "; fitted = " << fitRes[harm][method]->GetBinContent(iC+1) << "; true = " << h_Res_true[harm][0][method]->GetBinContent(iC+1) << "." << endl;
					h_Res_out[harm][0][method]->SetBinContent(iC+1,true_res);
				}
			//}
		}
	}
	for (Int_t method=0;method<_N_METHOD;method++){
		for (Int_t iC=0;iC<10;iC++){
			Double_t chi = Chi(h_Res_calc[1][1][method]->GetBinContent(iC+1),1,accuracy);
			chi      *= TMath::Sqrt(2);
			h_Res_out[1][1][method]->SetBinContent(iC+1,ResEventPlane(chi,1));
		}
	}
	/*for (Int_t method=0;method<_N_METHOD-1;method++){
		for (Int_t harm=0;harm<_N_HARM;harm++){
			for (Int_t _harm=0;_harm<_N_HARM;_harm++){
				for (Int_t iC=0;iC<10;iC++){
					Double_t chi = Chi(h_Res3_calc[harm][0][method]->GetBinContent(iC+1),harm+1,accuracy);
					chi      *= TMath::Sqrt(2);
					h_Res3_out[harm][0][method]->SetBinContent(iC+1,ResEventPlane(chi,harm+1));
				}
			}
		}
	}*/
	
	for (Int_t harm=0; harm<_N_HARM; harm++){
		for (Int_t method=0;method<_N_METHOD;method++){
			ResRatio[harm][method]->Divide(fitRes[harm][method],h_Res_out[harm][harm][method]);
		}
	}

	f_out->cd();
	for (Int_t harm=0; harm<_N_HARM; harm++){
		for (Int_t _harm=0;_harm<_N_HARM;_harm++){
			for (Int_t method=0;method<_N_METHOD;method++){
				if(_harm==0){
					fitRes[harm][method]->Write();
					ResRatio[harm][method]->Write();
				}
				h_Res_out[harm][_harm][method]->Write();
				h_Res3_calc[harm][_harm][method]->Write();
                h_Res_true[harm][_harm][method]->Write();
                h_proj_2[harm][_harm][method]->Write();
                p_Res2Psi_vs_b[harm][_harm][method]->Write();
                h_rat[harm][_harm][method]->Write();
			}
		}
	}
}