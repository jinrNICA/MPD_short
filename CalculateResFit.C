#include "variables.h"

void CalculateResFit(Int_t ResHarm, Int_t type, TString output_name)
{
	gROOT->LoadMacro("rxnpl.C");

	const TString methods_names[_N_METHOD] = {TString("TPC |#eta|>0."), TString("TPC |#eta|>0.2"), TString("TPC |#eta|>0.5"), TString("FHCal")};
	Double_t cent_bin;
	if (type == 0) cent_bin=_N_CENTRALITY_BINS;
	if (type == 1) cent_bin=NcentralityBinsRes;
	Double_t ResFit[_N_HARM][16][_N_METHOD];
	for (Int_t iharm=0; iharm<_N_HARM;iharm++){
		for (Int_t iC=0;iC<cent_bin;iC++){
			for (Int_t iM=0;iM<_N_METHOD;iM++){
				if (iM!=3) ResFit[iharm][iC][iM] = (Double_t) rxnpl(1,type,"output.root",iharm,iC,iM);
				if (iM==3) ResFit[iharm][iC][iM] = (Double_t) rxnpl(iharm+1,type,"output.root",0,iC,iM);
			}
		}
	}
	TFile* f;
	TH1F* fitRes[_N_HARM][_N_METHOD];
	if (type == 1){
		f = new TFile(output_name.Data(),"recreate");
		f->cd();
		for (Int_t iharm=0; iharm<_N_HARM; iharm++){
			for (Int_t iM=0;iM<_N_METHOD;iM++){
				fitRes[iharm][iM] = new TH1F(Form("fitRes%i%i",iharm,iM),Form("Fitted Resolution of %i harmonic for %s;b fm;",iharm+1,methods_names[iM].Data()),NcentralityBinsRes,centralityBinsRes);
			}
		}
		for (Int_t iharm=0; iharm<_N_HARM;iharm++){
			for (Int_t iM=0;iM<_N_METHOD;iM++){
				for (Int_t iC=0;iC<cent_bin;iC++){
					fitRes[iharm][iM]->SetBinContent(iC+1,ResFit[iharm][iC][iM]);
				}
			}
		}
		for (Int_t iharm=0; iharm<_N_HARM;iharm++){
			for (Int_t iM=0;iM<_N_METHOD;iM++){
				fitRes[iharm][iM]->Write();
			}
		}
		f->Close();
	}

}