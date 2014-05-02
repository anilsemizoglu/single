// Usage:
// > root -b doAll.C
// C++
#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <iomanip>

// ROOT
#include "TBenchmark.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TLegend.h"
#include "TPad.h"
#include "TH2F.h"
#include "TH1.h"
#include "THStack.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "Math/Vector4D.h"
#include "TRandom3.h"

// CMS2
//#include "/home/users/sanil/CORE/CMS2.h"
#include "/home/users/sanil/myheaders/SL.h"
#include "/home/users/sanil/CORE/MT2/MT2.h"
#include "/home/users/sanil/CORE/MT2/MT2Utility.h"

#include "/home/users/sanil/CORE/lostlepton/ETHLeptonSelections.h"
#include "trackSelections.h"
#include "muonSelections.h"
#include "susySelections.h"
#include "electronSelections.h"
#include "eventSelections.h"


#include "../myheaders/myheader.h"

using namespace std;
using namespace sl;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");
////////////////////////////////////////////////////////////////////////////////////////////
int mt2_l_bin = 0;
int mt2_h_bin = 200;
int mt2_n_bin = 20;

int mt_l_bin = 0;
int mt_h_bin = 300;
int mt_n_bin = 30;

int met_l_bin = 0;
int met_h_bin = 225;
int met_n_bin = 45;

	TH1F *mt2_h_data = new TH1F("mt2_data_","MT2_data",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_data = new  TH1F("mt_data_", "MT_data",mt_n_bin,0,mt_h_bin);
	TH1F   *met_data = new TH1F("met_data_","met_data",met_n_bin,0,met_h_bin);

	TH1F *mt2_h_wjets = new TH1F("mt2_wjets_","MT2_wjets",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_wjets = new  TH1F("mt_wjets_", "MT_wjets",mt_n_bin,0,mt_h_bin);
	TH1F   *met_wjets = new TH1F("met_wjets_","met_wjets",met_n_bin,0,met_h_bin);

	TH1F *mt2_h_wwjets = new TH1F("mt2_wwjets_","MT2_wwjets",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_wwjets = new  TH1F("mt_wwjets_", "MT_wwjets",mt_n_bin,0,mt_h_bin);
	TH1F   *met_wwjets = new TH1F("met_wwjets_","met_wwjets",met_n_bin,0,met_h_bin);

	TH1F *mt2_h_ttwjets = new TH1F("mt2_ttwjets_","MT2_ttwjets",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_ttwjets = new  TH1F("mt_ttwjets_", "MT_ttwjets",mt_n_bin,0,mt_h_bin);
	TH1F   *met_ttwjets = new TH1F("met_ttwjets_","met_ttwjets",met_n_bin,0,met_h_bin);

	TH1F *mt2_h_ttzjets = new TH1F("mt2_ttzjets_","MT2_ttzjets",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_ttzjets = new  TH1F("mt_ttzjets_", "MT_ttzjets",mt_n_bin,0,mt_h_bin);
	TH1F   *met_ttzjets = new TH1F("met_ttzjets_","met_ttzjets",met_n_bin,0,met_h_bin);

	TH1F *mt2_h_dy = new TH1F("mt2_dy_","MT2_dy",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_dy = new  TH1F("mt_dy_", "MT_dy",mt_n_bin,0,mt_h_bin);
	TH1F   *met_dy = new TH1F("met_dy_","met_dy",met_n_bin,0,met_h_bin);

	TH1F *mt2_h_ttbar = new TH1F("mt2_ttbar_","MT2_ttbar",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_ttbar = new  TH1F("mt_ttbar_", "MT_ttbar",mt_n_bin,0,mt_h_bin);
	TH1F   *met_ttbar = new TH1F("met_ttbar_","met_ttbar",met_n_bin,0,met_h_bin);





TRandom  *d = new TRandom();
TRandom3 rx;
TRandom3 ry;

//mt2_h_data->Sumw2();
//mt_h_data ->Sumw2();
//met_data  ->Sumw2();

/////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
//GLOBAL COUNTERS//
int btag1counter= 0;
int jet2counter=0;
int met100counter=0;
int vtx10counter=0; 	    
int file_count = 1;

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("tree");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    single.Init(tree);
    
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      single.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      SL::progress( nEventsTotal, nEventsChain );

      // Analysis Code
//event based variables
	float mu_max = 0.0, el_max = 0.0, xc = 0 , yc = 0, zc = 0;
	//COUNTERS
	int el_count = 0, mu_count = 0, nJets = 0, nBtags = 0, good_index = -1, mu_index = -1, el_index = -1;

     //jet looper
      for (unsigned int k = 0; k < jets_p4().size(); k++){
        if (jets_p4().at(k).pt()*jets_p4Correction().at(k) < 30) continue;
        if (fabs(jets_p4().at(k).eta()) > 2.5) continue;
        //if (ROOT::Math::VectorUtil::DeltaR(jets_p4().at(k), lr_p4()) < 0.4) continue;
        nJets++;
        if (btagDiscriminant().at(k) < 0.244) continue;
        nBtags++;
      }
	if(lr_p4().pt() < 30) continue;

//////////// SELECTION ///////
if(file_count == 1){
met_data->Fill(met());
}
else if(file_count == 2){
met_ttzjets->Fill(met(),scale_1fb()*5.2);
}
else if(file_count == 3){
met_wjets->Fill(met(),scale_1fb()*5.2);
}
else if(file_count == 4){
met_dy->Fill(met(),scale_1fb()*5.2);
}
else if(file_count == 5){
met_ttbar->Fill(met(),scale_1fb()*5.2);
}
else if(file_count == 6){
met_ttwjets->Fill(met(),scale_1fb()*5.2);
}
else if(file_count == 7){
met_wwjets->Fill(met(),scale_1fb()*5.2);
}

if(met() < 100 ) met100counter++;
if(met() < 100 ) continue;
/*
if (nJets < 1) jet2counter++;
if (nJets < 1) continue;
if (nBtags < 1) btag1counter++;
if (nBtags < 1) continue;
if(nvtxs() < 10) vtx10counter++;
if(nvtxs() < 10) continue;
*/
//////////// SELECTION ///////

/////////////    ADD FAKE    ////////////

//magnitude of the vector
	float r = 40, x = 0, y = 0, phi =0 , theta = 0;
	float a[2];
//Create Lepton fake
	LorentzVector fake_lep;
	do{
//random numbers x,y
	rx.RndmArray(2,a);
	x = a[0];
	y = a[1];
//phi 2*pi*x(0to1)
	phi = (x-0.5)*2*3.1415;
//arccos(2y(0to1)-1)
	theta = acos(2*y - 1);

	 xc = r*sin(theta)*cos(phi);
	 yc = r*sin(theta)*sin(phi);
	 zc = r*cos(theta);
//set px py pz
	fake_lep.SetPxPyPzE(xc,yc,zc,40); 
		} while  (fake_lep.Eta() > 2.4 );
//add neutrino by changing metPhi and magnitude of met
//calculate the phi of the fake neutrino from the fake lepton phi
//phi is -pi to pi
	float fake_metPhi = -fake_lep.Phi();
//now vector addition using fake_metPhi, metPhi, met, and 40(momentum of the new neutrino)
	float metx = cos(metPhi())*met();
	float mety = sin(metPhi())*met();

	float fake_metx = cos(fake_metPhi)*40;
	float fake_mety = sin(fake_metPhi)*40;
//fake met vector
	LorentzVector fake_met;
	fake_met.SetPxPyPzE(fake_metx,fake_mety,0,0); 
//real met vector
	LorentzVector real_met;
	real_met.SetPxPyPzE(metx,mety,0,0); 
//new met vector
	LorentzVector new_met;
	new_met = real_met + fake_met;

///////////////  END FAKE   ////////////////////////////////////

////////  MT //////
float dphi = lr_p4().Phi() - metPhi();
if(dphi > 3.1415) dphi = 6.28 - dphi;
float MT = sqrt(2 * lr_p4().pt() * met() *(1 - cos(dphi))); 
////////  MT //////

////////  MT2 //////
double mt2_event = MT2(new_met.pt(),new_met.Phi(), lr_p4(), fake_lep);
////////  MT2 //////

////////  FILL //////
if(file_count == 1){
mt2_h_data->Fill(mt2_event);
mt_h_data->Fill(MT);
}
else if(file_count == 2){
mt2_h_ttzjets->Fill(mt2_event,scale_1fb()*5.2);
mt_h_ttzjets->Fill(MT,scale_1fb()*5.2);
}
else if(file_count == 3){
mt2_h_wjets->Fill(mt2_event,scale_1fb()*5.2);
mt_h_wjets->Fill(MT,scale_1fb()*5.2);
}
else if(file_count == 4){
mt2_h_dy->Fill(mt2_event,scale_1fb()*5.2);
mt_h_dy->Fill(MT,scale_1fb()*5.2);
}
else if(file_count == 5){
mt2_h_ttbar->Fill(mt2_event,scale_1fb()*5.2);
mt_h_ttbar->Fill(MT,scale_1fb()*5.2);
}
else if(file_count == 6){
mt2_h_ttwjets->Fill(mt2_event,scale_1fb()*5.2);
mt_h_ttwjets->Fill(MT,scale_1fb()*5.2);
}
else if(file_count == 7){
mt2_h_wwjets->Fill(mt2_event,scale_1fb()*5.2);
mt_h_wwjets->Fill(MT,scale_1fb()*5.2);
}

}//event 
//if(met() < 40 ) ;
//if (nJets < 2) mt_h->Fill(MT);
////////  FILL //////

	/*
	TCanvas *c73 = new TCanvas("c73","c73");
	mt_h_data->Draw();
	TCanvas *c74 = new TCanvas("c74","c74");
	mt2_h_data->Draw();
	TCanvas *c75 = new TCanvas("c75","c75");
	met_data->Draw();
*/
if(file_count == 1){
     TFile* fout = new TFile("/home/users/sanil/single/data_hists_met100.root","RECREATE");
	mt2_h_data->Write();
	mt_h_data->Write();
	met_data->Write();
  fout->Close();
}
else if(file_count == 2){
  TFile* fout = new TFile("/home/users/sanil/single/ttzjets_hists_met100.root","RECREATE");
	mt2_h_ttzjets->Write();
	 mt_h_ttzjets->Write();
	  met_ttzjets->Write();
  fout->Close();
}
else if(file_count == 3){
    TFile* fout = new TFile("/home/users/sanil/single/wjets_hists_met100.root","RECREATE");
	mt2_h_wjets->Write();
	 mt_h_wjets->Write();
	  met_wjets->Write();
  fout->Close();
}
else if(file_count == 4){
       TFile* fout = new TFile("/home/users/sanil/single/dy_hists_met100.root","RECREATE");
	mt2_h_dy->Write();
	 mt_h_dy->Write();
	  met_dy->Write();
  fout->Close();
}
else if(file_count == 5){
    TFile* fout = new TFile("/home/users/sanil/single/ttbar_hists_met100.root","RECREATE");
	mt2_h_ttbar->Write();
	 mt_h_ttbar->Write();
	  met_ttbar->Write();
  fout->Close();
}
else if(file_count == 6){
  TFile* fout = new TFile("/home/users/sanil/single/ttwjets_hists_met100.root","RECREATE");
	mt2_h_ttwjets->Write();
	 mt_h_ttwjets->Write();
	  met_ttwjets->Write();
  fout->Close();
}
else if(file_count == 7){
   TFile* fout = new TFile("/home/users/sanil/single/wwjets_hists_met100.root","RECREATE");
	mt2_h_wwjets->Write();
	 mt_h_wwjets->Write();
	  met_wwjets->Write();
  fout->Close();
}


    // Clean Up
    delete tree;
    file->Close();
    delete file;

	file_count++;
  }	//file_loop
met_data->Draw("P+E");
met_wjets->Draw("same");
//Normalize(mt_h_data);
//Normalize(mt2_h_data);
//getbin(mt2_h_data,"mt2_0_data");
//getbin(mt_h_data,"mt_0_data");

  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
	

  
  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return 0;
	}
