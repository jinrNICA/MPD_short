/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: rootlogon.C,v 1.4 2007/05/04 11:18:51 ivana Exp $ */

/// By Laurent Aphecetche

{
  //cout << "Loading MPD libraries ..." << endl;
  //gROOT->LoadMacro("${ALICE_ROOT}/MUON/loadlibs.C");
  //gInterpreter->ProcessLine("loadlibs()");
    
  cout << "Setting include path ..." << endl;
  TString includePath = "-I${ROOT_INCLUDE_DIR} ";
  includePath        += "-I${VMCWORKDIR}/base ";
  includePath        += "-I${VMCWORKDIR}/geobase ";
  includePath        += "-I${VMCWORKDIR}/tpc ";
  includePath        += "-I${VMCWORKDIR}/kalman ";
  includePath        += "-I${VMCWORKDIR}/lhetrack ";
  includePath        += "-I${VMCWORKDIR}/mcstack ";
  includePath        += "-I${VMCWORKDIR}/strawendcap ";
  includePath        += "-I${VMCWORKDIR}/etof ";
  includePath        += "-I${VMCWORKDIR}/tof ";
  includePath        += "-I${VMCWORKDIR}/sft ";
  includePath        += "-I${VMCWORKDIR}/parbase ";
  includePath        += "-I${VMCWORKDIR}/mpddata ";
  includePath        += "-I${VMCWORKDIR}/mpdbase ";
  includePath        += "-I${VMCWORKDIR}/fairtools ";
  includePath        += "-I${VMCWORKDIR}/mpdpid ";
  /*
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/FASTSIM ";
  includePath        += "-I${ALICE_ROOT}/EVGEN ";
  includePath        += "-I${ALICE_ROOT}/SHUTTLE/TestShuttle ";
  includePath        += "-I${ALICE_ROOT}/ITS ";
  includePath        += "-I${ALICE_ROOT}/MUON ";
  includePath        += "-I${ALICE_ROOT}/MUON/mapping";
  */
  gSystem->SetIncludePath(includePath.Data());
  
  gSystem->Load("libRIO");
  gSystem->Load("libGeom");
  gSystem->Load("libGeomPainter");
  gSystem->Load("libVMC");
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libPythia6");
  gSystem->Load("libPluto");
  gSystem->Load("libPhysics");
  gSystem->Load("libNet");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libMathMore");

  gSystem->Load("libProof");
  gSystem->Load("libProofPlayer");
  gSystem->Load("libGX11TTF");
  gSystem->Load("libGX11");

  gSystem->Load("libGdml");
  
  // Load other libraries
  gSystem->Load("libFairTools");
  gSystem->Load("libGeoBase");
  gSystem->Load("libParBase");
  gSystem->Load("libBase");
  gSystem->Load("libMCStack");
  gSystem->Load("libMpdField");
  gSystem->Load("libPassive");
  gSystem->Load("libGen");
  gSystem->Load("libTrkBase");
  gSystem->Load("libMpdBase");
  gSystem->Load("libMpdData");
  gSystem->Load("libMpdPid");
  gSystem->Load("libMpdgenerators");

  // HADGEN
  gSystem->Load("libHADGEN.so");
  gSystem->Load("libTHadgen.so");
  gSystem->Load("libMpdGeneralGenerator.so");

  // FEMTOSCOPY
  //gSystem->Load("libMpdFemto.so");
  gSystem->Load("libMpdPid.so");
  
  gSystem->Load("libKalman");
  gSystem->Load("libtpc");
  gSystem->Load("libLHETrack");
  gSystem->Load("libGeane");
  
  gSystem->Load("libtpc");
  gSystem->Load("libTof");
  //gSystem->Load("libEtof");
  gSystem->Load("libEmc");
  gSystem->Load("libStrawendcap");
  //gSystem->Load("libStt");
  //gSystem->Load("libSts");
  //gSystem->Load("libBbc");
  gSystem->Load("libZdc");
  //gSystem->Load("libFsa");
  gSystem->Load("libFfd");
  //gSystem->Load("libCpc");
  //gSystem->Load("libNDet");
  //gSystem->Load("libSft");
  gSystem->Load("libStrawECT");
  
  //gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
  //mpdloadlibs(kTRUE,kTRUE); // all libs
}
