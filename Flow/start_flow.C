void start_flow(TString input, TString output){
  gSystem->Load("real_flow_C.so");
  real_flow *t = new real_flow(input.Data(),output.Data());
  t->Loop();
}
