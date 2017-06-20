void MakeStyle(TString output){
    
    const Int_t NRGBs = 5;
    const Int_t NCont = 256;
    
    TStyle* style = new TStyle("style","my_style");
    
    style->SetTitleSize(0.08,"X");
    style->SetTitleSize(0.08,"Y");
    style->SetTitleOffset(1.1,"Y");
    style->SetTitleOffset(0.9,"X");
    style->SetFrameLineWidth(2);
    style->SetLabelSize(0.06,"X");
    style->SetLabelSize(0.06,"Y");
    style->SetPadTopMargin(0.1);
    style->SetPadBottomMargin(0.18);
    style->SetPadLeftMargin(0.18);
    style->SetPadTopMargin(0.1);
    style->SetPadRightMargin(0.1);
    //style->SetMarkerSize(2.);
    style->SetErrorX(0);
    style->SetOptStat(0);
    style->SetCanvasColor(0);
    style->SetTitleFillColor(0);
    style->SetTitleBorderSize(0);
    //style->SetMarkerStyle(20);
    style->SetHistLineWidth(3);
    style->SetLineWidth(2);
    style->SetCanvasBorderMode(0);
    style->SetLegendBorderSize(0);
    style->SetPadBorderMode(0);
    
    Double_t stops[NRGBs] = { 0.00, 0.30, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs] = { 0.00, 0.00, 0.57, 0.90, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.65, 0.95, 0.20, 0.00 };
    Double_t blue[NRGBs] = { 0.51, 0.55, 0.15, 0.00, 0.10 };
    TColor* color = new TColor();
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    style->SetNumberContours(NCont); 
    
    TFile* file = new TFile(output.Data(),"recreate");
    file->cd();
    style->Write();
    file->Close();
}
