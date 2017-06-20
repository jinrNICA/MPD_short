//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 11 15:15:19 2016 by ROOT version 5.34/19
// from TTree cbmsim_reduced/cbmsim_reduced
// found on file: ../mpd_data/7.7gev/reco/merged.root
//////////////////////////////////////////////////////////

#define _MAX_TRACKS 3000
#define _MAX_MODULE 90

#define GsiLustre

#ifndef real_flow_h
#define real_flow_h

#include "variables.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class real_flow {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TChain         *inChain;

   Long64_t        nentries;

   // Declaration of leaf types
   Float_t         b_mc;
   Float_t         phiEP_mc;
   Float_t         x_vertex_mc;
   Float_t         y_vertex_mc;
   Float_t         z_vertex_mc;
   Long64_t        n_tracks_mc;
   Float_t         eta_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Float_t         pt_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Int_t           mother_ID_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Int_t           PDG_code_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Float_t         px_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Float_t         py_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Float_t         pz_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Float_t         start_x_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Float_t         start_y_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Float_t         start_z_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Float_t         mass_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Float_t         energy_mc[_MAX_TRACKS];   //[n_tracks_mc]
   Long64_t        n_tracks_mpd;
   Float_t         eta_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         phi_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         theta_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Int_t           TOF_flag_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         ZDC_energy_mpd[_MAX_MODULE];
   Float_t         pid_tpc_prob_electron_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         pid_tpc_prob_pion_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         pid_tpc_prob_kaon_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         pid_tpc_prob_proton_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         pid_tof_prob_electron_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         pid_tof_prob_pion_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         pid_tof_prob_kaon_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         pid_tof_prob_proton_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         tof_beta_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         tof_mass2_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         dEdx_tpc_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         chi2_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         pt_error_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         theta_error_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         phi_error_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         DCA_x_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         DCA_y_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         DCA_z_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         DCA_global_x_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         DCA_global_y_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         DCA_global_z_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Int_t           n_hits_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Int_t           n_hits_poss_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Float_t         signed_pt_mpd[_MAX_TRACKS];   //[n_tracks_mpd]
   Long64_t        id_from_mc_mpd[_MAX_TRACKS];   //[n_tracks_mpd]

   TString         outFile;
   TFile*          f_out;
   TFile*          out_file;
   TTree*          out_tree;

   #ifdef GsiLustre
   TString in_tree=TString("/lustre/nyx/hades/user/parfenov/merged7GeV.root");
   #else 
   TString in_tree=TString("/mnt/pool/rhic/4/parfenovpeter/mpd_data/7.7gev/reco/merged.root");
   #endif

   TProfile* p_flow_wrt_full_vs_pt[_N_SORTS][_N_CENTRALITY_BINS][_N_HARM][_N_HARM][_N_METHOD];
   TProfile* p_flow_wrt_full_vs_rapidity[_N_SORTS][_N_CENTRALITY_BINS][_N_HARM][_N_HARM][_N_METHOD];
   TProfile* p_flow_wrt_full_vs_eta[_N_SORTS][_N_CENTRALITY_BINS][_N_HARM][_N_HARM][_N_METHOD];
   TProfile* p_flow_wrt_RP_vs_pt[_N_SORTS][_N_CENTRALITY_BINS][_N_HARM];
   TProfile* p_flow_wrt_RP_vs_rapidity[_N_SORTS][_N_CENTRALITY_BINS][_N_HARM];
   TProfile* p_flow_wrt_RP_vs_eta[_N_SORTS][_N_CENTRALITY_BINS][_N_HARM];
   TProfile* h_Qvsb[_N_ARM][_N_HARM][_N_METHOD][2];
   TProfile* h_QvsMult[_N_ARM][_N_HARM][_N_METHOD][2];
   TProfile* p_true_Res_vs_b[_N_HARM][_N_HARM][_N_METHOD];
   TProfile* p_Res2Psi_vs_b[_N_HARM][_N_HARM][_N_METHOD];
   TProfile* p_ResAPsi_vs_b[_N_HARM][_N_HARM][_N_METHOD];
   TProfile* p_ResBPsi_vs_b[_N_HARM][_N_HARM][_N_METHOD];
   TProfile* p_ResPsi_vs_b[_N_HARM][_N_HARM][_N_METHOD];
   TProfile* p_true_Res_half_vs_b[_N_ARM][_N_HARM][_N_HARM][_N_METHOD];
   TProfile* p_qx_vs_b[_N_ARM][_N_HARM][_N_METHOD];
   TProfile* p_qy_vs_b[_N_ARM][_N_HARM][_N_METHOD];
   TProfile* p_Cos_A_vs_b[_N_HARM][_N_HARM][_N_METHOD-1];
   TProfile* p_Cos_B_vs_b[_N_HARM][_N_HARM][_N_METHOD-1];
   TProfile* p_Cos_C_vs_b[_N_HARM][_N_HARM][_N_METHOD-1];
   TH1F*     h_delta_psi[_N_CENTRALITY_BINS][_N_HARM][_N_METHOD];
   TH1F*     h_delta_EP_full_RP[_N_CENTRALITY_BINS][_N_HARM][_N_METHOD];
   TH1F*     h_delta_fit[NcentralityBinsRes][_N_HARM][_N_METHOD];
   TH1F*     h_d_psi_fit[_N_CENTRALITY_BINS][_N_HARM][_N_METHOD];
   TH1F*     h_qx[NcentralityBinsRes][_N_ARM][_N_HARM][_N_METHOD];
   TH1F*     h_qy[NcentralityBinsRes][_N_ARM][_N_HARM][_N_METHOD];
   TH1F*     h_delta_EP_half_RP[NcentralityBinsRes][_N_ARM][_N_HARM][_N_METHOD];
   TH2F*     h_cor_EP_RP[NcentralityBinsRes][_N_ARM][_N_HARM][_N_METHOD];
   TH1F*     h_psi[NcentralityBinsRes][_N_ARM][_N_HARM][_N_METHOD];
   TH1F*     h_Res_true[_N_HARM][_N_HARM][_N_METHOD];
   TH1F*     h_Res3_true[_N_HARM][_N_HARM][_N_METHOD];
   TH1F*     h_Res1_fit[NcentralityBinsRes];
   TH1F*     h_psi_RP[NcentralityBinsRes];
   TH1F*     h_Mult_vs_eta[_N_HARM][_N_ARM][_N_METHOD-1];
   //TProfile* p_Mult_vs_eta[_N_ARM][_N_METHOD-1];
   TH1F*     h_Ratio_eta[_N_METHOD-1];

   char name[200];
   char title[200];

   // List of branches
   TBranch        *b_b_mc;   //!
   TBranch        *b_phiEP_mc;   //!
   TBranch        *b_x_vertex_mc;   //!
   TBranch        *b_y_vertex_mc;   //!
   TBranch        *b_z_vertex_mc;   //!
   TBranch        *b_n_tracks_mc;   //!
   TBranch        *b_eta_mc;   //!
   TBranch        *b_pt_mc;   //!
   TBranch        *b_mother_ID_mc;   //!
   TBranch        *b_PDG_code_mc;   //!
   TBranch        *b_px_mc;   //!
   TBranch        *b_py_mc;   //!
   TBranch        *b_pz_mc;   //!
   TBranch        *b_start_x_mc;   //!
   TBranch        *b_start_y_mc;   //!
   TBranch        *b_start_z_mc;   //!
   TBranch        *b_mass_mc;   //!
   TBranch        *b_energy_mc;   //!
   TBranch        *b_n_tracks_mpd;   //!
   TBranch        *b_eta_mpd;   //!
   TBranch        *b_phi_mpd;   //!
   TBranch        *b_theta_mpd;   //!
   TBranch        *b_TOF_flag_mpd;   //!
   TBranch        *b_ZDC_energy_mpd;   //!
   TBranch        *b_pid_tpc_prob_electron_mpd;   //!
   TBranch        *b_pid_tpc_prob_pion_mpd;   //!
   TBranch        *b_pid_tpc_prob_kaon_mpd;   //!
   TBranch        *b_pid_tpc_prob_proton_mpd;   //!
   TBranch        *b_pid_tof_prob_electron_mpd;   //!
   TBranch        *b_pid_tof_prob_pion_mpd;   //!
   TBranch        *b_pid_tof_prob_kaon_mpd;   //!
   TBranch        *b_pid_tof_prob_proton_mpd;   //!
   TBranch        *b_tof_beta_mpd;   //!
   TBranch        *b_tof_mass2_mpd;   //!
   TBranch        *b_dEdx_tpc_mpd;   //!
   TBranch        *b_chi2_mpd;   //!
   TBranch        *b_pt_error_mpd;   //!
   TBranch        *b_theta_error_mpd;   //!
   TBranch        *b_phi_error_mpd;   //!
   TBranch        *b_DCA_x_mpd;   //!
   TBranch        *b_DCA_y_mpd;   //!
   TBranch        *b_DCA_z_mpd;   //!
   TBranch        *b_DCA_global_x_mpd;   //!
   TBranch        *b_DCA_global_y_mpd;   //!
   TBranch        *b_DCA_global_z_mpd;   //!
   TBranch        *b_n_hits_mpd;   //!
   TBranch        *b_n_hits_poss_mpd;   //!
   TBranch        *b_signed_pt_mpd;   //!
   TBranch        *b_id_from_mc_mpd;   //!

   real_flow(TTree *tree=0);
   real_flow(TString inFile_name, TString outFile_name);
   virtual ~real_flow();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Init(TChain *chain);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Int_t    GetMultiplicityTPC();
   virtual Int_t    GetMultiplicityTPC(Int_t gap, TVector3* part, Long64_t ntracks, Int_t sign);
   virtual Int_t    GetCentralityBinFlow(Int_t multiplicity);
   virtual Int_t    GetCentralityBinRes(Int_t multiplicity);
   virtual Double_t GetTotalEnergy(Float_t* zdc_energy, Int_t zdc_ID);
   virtual Double_t GetPsiFullZdc(Float_t* zdc_energy, Int_t n);
   virtual void     GetQsZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t harm, Double_t &Qx, Double_t &Qy);
   virtual Double_t* GetAngles();
   virtual void     Write();
   virtual void     InitHisto();
   virtual Float_t  Unfold(Float_t phiEP_mc, Float_t psi_N_FULL, Int_t harm);
   virtual void     FillZDC(Int_t cenrality_bin, Float_t *ZDC_energy_mpd);
   virtual Double_t GetPsiHalfZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t n, Double_t &qx, Double_t &qy);
   virtual void     DrawResFit();
   virtual void     GetQsTpc(Int_t gap, Int_t weight, TVector3* part, Long64_t ntracks, Int_t sign, Int_t harm, Double_t &Qx, Double_t &Qy);
   virtual Double_t GetTotalMomenta(Int_t gap, TVector3* part, Long64_t ntracks, Int_t sign);
   virtual Double_t GetPsiHalfTpc(Int_t gap, Int_t weight, TVector3* part, Long64_t ntracks, Int_t sign, Int_t harm, Double_t &Qx, Double_t &Qy);
   virtual Double_t GetPsiFullTpc(Int_t gap, Int_t weight, TVector3* part, Long64_t ntracks, Int_t harm);
   virtual void     FillTPC(Int_t gap, Int_t weight, Int_t cenrality_bin, TVector3 *part, Long64_t ntracks);
   virtual void     Fill3Sub(Int_t gap, TVector3 *part, Long64_t ntracks,Float_t *ZDC_energy_mpd);
   virtual void     MakeTTree(TString outputTree, Double_t Qx[_N_HARM][_N_ARM][_N_METHOD], Double_t Qy[_N_HARM][_N_ARM][_N_METHOD], Int_t Qw[_N_HARM][_N_ARM][_N_METHOD]);
   virtual void     FillTTree();
   virtual void     SaveTTree(TString outputFile);
   //virtual void     GetFitHistZDC(Int_t cenrality_bin, Float_t *ZDC_energy_mpd);
};

#endif

#ifdef real_flow_cxx
real_flow::real_flow(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(in_tree.Data());
      if (!f || !f->IsOpen()) {
         f = new TFile(in_tree.Data());
      }
      f->GetObject("cbmsim_reduced",tree);

   }
   Init(tree);
   outFile = "output.root";

}
real_flow::real_flow(TString inFile_name, TString outFile_name)
{
   //inChain = new TChain("cbmsim_reduced");
   //inChain->Add(inFile_name.Data(),"read");
   TTree* tree;

   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile_name.Data());
   if (!f || !f->IsOpen()) {
      cout << "Open file: " << inFile_name.Data() << endl;
      f = new TFile(inFile_name.Data(),"read");
      f->GetObject("cbmsim_reduced",tree);

      Init(tree);
   }
   outFile = outFile_name;
   //f_out = new TFile(outFile_name.Data(),"recreate");
}

real_flow::~real_flow()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t real_flow::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t real_flow::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void real_flow::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("b_mc", &b_mc, &b_b_mc);
   fChain->SetBranchAddress("phiEP_mc", &phiEP_mc, &b_phiEP_mc);
   fChain->SetBranchAddress("x_vertex_mc", &x_vertex_mc, &b_x_vertex_mc);
   fChain->SetBranchAddress("y_vertex_mc", &y_vertex_mc, &b_y_vertex_mc);
   fChain->SetBranchAddress("z_vertex_mc", &z_vertex_mc, &b_z_vertex_mc);
   fChain->SetBranchAddress("n_tracks_mc", &n_tracks_mc, &b_n_tracks_mc);
   fChain->SetBranchAddress("eta_mc", eta_mc, &b_eta_mc);
   fChain->SetBranchAddress("pt_mc", pt_mc, &b_pt_mc);
   fChain->SetBranchAddress("mother_ID_mc", mother_ID_mc, &b_mother_ID_mc);
   fChain->SetBranchAddress("PDG_code_mc", PDG_code_mc, &b_PDG_code_mc);
   fChain->SetBranchAddress("px_mc", px_mc, &b_px_mc);
   fChain->SetBranchAddress("py_mc", py_mc, &b_py_mc);
   fChain->SetBranchAddress("pz_mc", pz_mc, &b_pz_mc);
   fChain->SetBranchAddress("start_x_mc", start_x_mc, &b_start_x_mc);
   fChain->SetBranchAddress("start_y_mc", start_y_mc, &b_start_y_mc);
   fChain->SetBranchAddress("start_z_mc", start_z_mc, &b_start_z_mc);
   fChain->SetBranchAddress("mass_mc", mass_mc, &b_mass_mc);
   fChain->SetBranchAddress("energy_mc", energy_mc, &b_energy_mc);
   fChain->SetBranchAddress("n_tracks_mpd", &n_tracks_mpd, &b_n_tracks_mpd);
   fChain->SetBranchAddress("eta_mpd", eta_mpd, &b_eta_mpd);
   fChain->SetBranchAddress("phi_mpd", phi_mpd, &b_phi_mpd);
   fChain->SetBranchAddress("theta_mpd", theta_mpd, &b_theta_mpd);
   fChain->SetBranchAddress("TOF_flag_mpd", TOF_flag_mpd, &b_TOF_flag_mpd);
   fChain->SetBranchAddress("ZDC_energy_mpd", ZDC_energy_mpd, &b_ZDC_energy_mpd);
   fChain->SetBranchAddress("pid_tpc_prob_electron_mpd", pid_tpc_prob_electron_mpd, &b_pid_tpc_prob_electron_mpd);
   fChain->SetBranchAddress("pid_tpc_prob_pion_mpd", pid_tpc_prob_pion_mpd, &b_pid_tpc_prob_pion_mpd);
   fChain->SetBranchAddress("pid_tpc_prob_kaon_mpd", pid_tpc_prob_kaon_mpd, &b_pid_tpc_prob_kaon_mpd);
   fChain->SetBranchAddress("pid_tpc_prob_proton_mpd", pid_tpc_prob_proton_mpd, &b_pid_tpc_prob_proton_mpd);
   fChain->SetBranchAddress("pid_tof_prob_electron_mpd", pid_tof_prob_electron_mpd, &b_pid_tof_prob_electron_mpd);
   fChain->SetBranchAddress("pid_tof_prob_pion_mpd", pid_tof_prob_pion_mpd, &b_pid_tof_prob_pion_mpd);
   fChain->SetBranchAddress("pid_tof_prob_kaon_mpd", pid_tof_prob_kaon_mpd, &b_pid_tof_prob_kaon_mpd);
   fChain->SetBranchAddress("pid_tof_prob_proton_mpd", pid_tof_prob_proton_mpd, &b_pid_tof_prob_proton_mpd);
   fChain->SetBranchAddress("tof_beta_mpd", tof_beta_mpd, &b_tof_beta_mpd);
   fChain->SetBranchAddress("tof_mass2_mpd", tof_mass2_mpd, &b_tof_mass2_mpd);
   fChain->SetBranchAddress("dEdx_tpc_mpd", dEdx_tpc_mpd, &b_dEdx_tpc_mpd);
   fChain->SetBranchAddress("chi2_mpd", chi2_mpd, &b_chi2_mpd);
   fChain->SetBranchAddress("pt_error_mpd", pt_error_mpd, &b_pt_error_mpd);
   fChain->SetBranchAddress("theta_error_mpd", theta_error_mpd, &b_theta_error_mpd);
   fChain->SetBranchAddress("phi_error_mpd", phi_error_mpd, &b_phi_error_mpd);
   fChain->SetBranchAddress("DCA_x_mpd", DCA_x_mpd, &b_DCA_x_mpd);
   fChain->SetBranchAddress("DCA_y_mpd", DCA_y_mpd, &b_DCA_y_mpd);
   fChain->SetBranchAddress("DCA_z_mpd", DCA_z_mpd, &b_DCA_z_mpd);
   fChain->SetBranchAddress("DCA_global_x_mpd", DCA_global_x_mpd, &b_DCA_global_x_mpd);
   fChain->SetBranchAddress("DCA_global_y_mpd", DCA_global_y_mpd, &b_DCA_global_y_mpd);
   fChain->SetBranchAddress("DCA_global_z_mpd", DCA_global_z_mpd, &b_DCA_global_z_mpd);
   fChain->SetBranchAddress("n_hits_mpd", n_hits_mpd, &b_n_hits_mpd);
   fChain->SetBranchAddress("n_hits_poss_mpd", n_hits_poss_mpd, &b_n_hits_poss_mpd);
   fChain->SetBranchAddress("signed_pt_mpd", signed_pt_mpd, &b_signed_pt_mpd);
   fChain->SetBranchAddress("id_from_mc_mpd", id_from_mc_mpd, &b_id_from_mc_mpd);
   Notify();
}

void real_flow::Init(TChain *chain)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers

   chain->SetBranchAddress("b_mc", &b_mc, &b_b_mc);
   chain->SetBranchAddress("phiEP_mc", &phiEP_mc, &b_phiEP_mc);
   chain->SetBranchAddress("x_vertex_mc", &x_vertex_mc, &b_x_vertex_mc);
   chain->SetBranchAddress("y_vertex_mc", &y_vertex_mc, &b_y_vertex_mc);
   chain->SetBranchAddress("z_vertex_mc", &z_vertex_mc, &b_z_vertex_mc);
   chain->SetBranchAddress("n_tracks_mc", &n_tracks_mc, &b_n_tracks_mc);
   chain->SetBranchAddress("eta_mc", eta_mc, &b_eta_mc);
   chain->SetBranchAddress("pt_mc", pt_mc, &b_pt_mc);
   chain->SetBranchAddress("mother_ID_mc", mother_ID_mc, &b_mother_ID_mc);
   chain->SetBranchAddress("PDG_code_mc", PDG_code_mc, &b_PDG_code_mc);
   chain->SetBranchAddress("px_mc", px_mc, &b_px_mc);
   chain->SetBranchAddress("py_mc", py_mc, &b_py_mc);
   chain->SetBranchAddress("pz_mc", pz_mc, &b_pz_mc);
   chain->SetBranchAddress("start_x_mc", start_x_mc, &b_start_x_mc);
   chain->SetBranchAddress("start_y_mc", start_y_mc, &b_start_y_mc);
   chain->SetBranchAddress("start_z_mc", start_z_mc, &b_start_z_mc);
   chain->SetBranchAddress("mass_mc", mass_mc, &b_mass_mc);
   chain->SetBranchAddress("energy_mc", energy_mc, &b_energy_mc);
   chain->SetBranchAddress("n_tracks_mpd", &n_tracks_mpd, &b_n_tracks_mpd);
   chain->SetBranchAddress("eta_mpd", eta_mpd, &b_eta_mpd);
   chain->SetBranchAddress("phi_mpd", phi_mpd, &b_phi_mpd);
   chain->SetBranchAddress("theta_mpd", theta_mpd, &b_theta_mpd);
   chain->SetBranchAddress("TOF_flag_mpd", TOF_flag_mpd, &b_TOF_flag_mpd);
   chain->SetBranchAddress("ZDC_energy_mpd", ZDC_energy_mpd, &b_ZDC_energy_mpd);
   chain->SetBranchAddress("pid_tpc_prob_electron_mpd", pid_tpc_prob_electron_mpd, &b_pid_tpc_prob_electron_mpd);
   chain->SetBranchAddress("pid_tpc_prob_pion_mpd", pid_tpc_prob_pion_mpd, &b_pid_tpc_prob_pion_mpd);
   chain->SetBranchAddress("pid_tpc_prob_kaon_mpd", pid_tpc_prob_kaon_mpd, &b_pid_tpc_prob_kaon_mpd);
   chain->SetBranchAddress("pid_tpc_prob_proton_mpd", pid_tpc_prob_proton_mpd, &b_pid_tpc_prob_proton_mpd);
   chain->SetBranchAddress("pid_tof_prob_electron_mpd", pid_tof_prob_electron_mpd, &b_pid_tof_prob_electron_mpd);
   chain->SetBranchAddress("pid_tof_prob_pion_mpd", pid_tof_prob_pion_mpd, &b_pid_tof_prob_pion_mpd);
   chain->SetBranchAddress("pid_tof_prob_kaon_mpd", pid_tof_prob_kaon_mpd, &b_pid_tof_prob_kaon_mpd);
   chain->SetBranchAddress("pid_tof_prob_proton_mpd", pid_tof_prob_proton_mpd, &b_pid_tof_prob_proton_mpd);
   chain->SetBranchAddress("tof_beta_mpd", tof_beta_mpd, &b_tof_beta_mpd);
   chain->SetBranchAddress("tof_mass2_mpd", tof_mass2_mpd, &b_tof_mass2_mpd);
   chain->SetBranchAddress("dEdx_tpc_mpd", dEdx_tpc_mpd, &b_dEdx_tpc_mpd);
   chain->SetBranchAddress("chi2_mpd", chi2_mpd, &b_chi2_mpd);
   chain->SetBranchAddress("pt_error_mpd", pt_error_mpd, &b_pt_error_mpd);
   chain->SetBranchAddress("theta_error_mpd", theta_error_mpd, &b_theta_error_mpd);
   chain->SetBranchAddress("phi_error_mpd", phi_error_mpd, &b_phi_error_mpd);
   chain->SetBranchAddress("DCA_x_mpd", DCA_x_mpd, &b_DCA_x_mpd);
   chain->SetBranchAddress("DCA_y_mpd", DCA_y_mpd, &b_DCA_y_mpd);
   chain->SetBranchAddress("DCA_z_mpd", DCA_z_mpd, &b_DCA_z_mpd);
   chain->SetBranchAddress("DCA_global_x_mpd", DCA_global_x_mpd, &b_DCA_global_x_mpd);
   chain->SetBranchAddress("DCA_global_y_mpd", DCA_global_y_mpd, &b_DCA_global_y_mpd);
   chain->SetBranchAddress("DCA_global_z_mpd", DCA_global_z_mpd, &b_DCA_global_z_mpd);
   chain->SetBranchAddress("n_hits_mpd", n_hits_mpd, &b_n_hits_mpd);
   chain->SetBranchAddress("n_hits_poss_mpd", n_hits_poss_mpd, &b_n_hits_poss_mpd);
   chain->SetBranchAddress("signed_pt_mpd", signed_pt_mpd, &b_signed_pt_mpd);
   chain->SetBranchAddress("id_from_mc_mpd", id_from_mc_mpd, &b_id_from_mc_mpd);
   Notify();
}

Bool_t real_flow::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void real_flow::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t real_flow::Cut(Long64_t entry)
{
   // This function may be called from Loop.
   // returns  1 if entry is accepted.
   // returns -1 otherwise.
   return 1;
}
#endif // #ifdef real_flow_cxx
