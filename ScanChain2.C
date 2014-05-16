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

int jet_n_bin = 10;
int jet_h_bin = 10;

	TH1F *mt2_h_data = new TH1F("mt2_data_","MT2_data",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_data = new  TH1F("mt_data_", "MT_data",mt_n_bin,0,mt_h_bin);
	TH1F   *met_data = new TH1F("met_data_","met_data",met_n_bin,0,met_h_bin);
	TH1F   *jets_data= new TH1F("jets_data_","Jets_data",jet_n_bin,0,jet_h_bin);

	TH1F *mt2_h_ttbar = new TH1F("mt2_ttbar_","MT2_ttbar",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_ttbar = new  TH1F("mt_ttbar_", "MT_ttbar",mt_n_bin,0,mt_h_bin);
	TH1F   *met_ttbar = new TH1F("met_ttbar_","met_ttbar",met_n_bin,0,met_h_bin);
	TH1F   *jets_ttbar  = new TH1F("jets_ttbar_","Jets_ttbar",jet_n_bin,0,jet_h_bin);

	TH1F *mt2_h_wjets = new TH1F("mt2_wjets_","MT2_wjets",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_wjets = new  TH1F("mt_wjets_", "MT_wjets",mt_n_bin,0,mt_h_bin);
	TH1F   *met_wjets = new TH1F("met_wjets_","met_wjets",met_n_bin,0,met_h_bin);
	TH1F   *jets_wjets  = new TH1F("jets_wjets_","Jets_wjets",jet_n_bin,0,jet_h_bin);

////MUONS////


	TH1F *mt2_h_data_mu = new TH1F("mt2_data_mu_","MT2_data_mu",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_data_mu = new  TH1F("mt_data_mu_", "MT_data_mu",mt_n_bin,0,mt_h_bin);
	TH1F   *met_data_mu = new TH1F("met_data_mu_","met_data_mu",met_n_bin,0,met_h_bin);
	TH1F   *jets_data_mu= new TH1F("jets_data_mu_","Jets_data_mu",jet_n_bin,0,jet_h_bin);

	TH1F *mt2_h_ttbar_mu = new TH1F("mt2_ttbar_mu_","MT2_ttbar_mu",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_ttbar_mu = new  TH1F("mt_ttbar_mu_", "MT_ttbar_mu",mt_n_bin,0,mt_h_bin);
	TH1F   *met_ttbar_mu = new TH1F("met_ttbar_mu_","met_ttbar_mu",met_n_bin,0,met_h_bin);
	TH1F   *jets_ttbar_mu  = new TH1F("jets_ttbar_mu_","Jets_ttbar_mu",jet_n_bin,0,jet_h_bin);

	TH1F *mt2_h_wjets_mu = new TH1F("mt2_wjets_mu_","MT2_wjets_mu",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_wjets_mu = new  TH1F("mt_wjets_mu_", "MT_wjets_mu",mt_n_bin,0,mt_h_bin);
	TH1F   *met_wjets_mu = new TH1F("met_wjets_mu_","met_wjets_mu",met_n_bin,0,met_h_bin);
	TH1F   *jets_wjets_mu  = new TH1F("jets_wjets_mu_","Jets_wjets_mu",jet_n_bin,0,jet_h_bin);

/*
	TH1F *mt2_h_wwjets = new TH1F("mt2_wwjets_","MT2_wwjets",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_wwjets = new  TH1F("mt_wwjets_", "MT_wwjets",mt_n_bin,0,mt_h_bin);
	TH1F   *met_wwjets = new TH1F("met_wwjets_","met_wwjets",met_n_bin,0,met_h_bin);
	TH1F   *jets_wwjets = new TH1F("jets_wwjets_","Jets_wwjets",jet_n_bin,0,jet_h_bin);

	TH1F *mt2_h_ttwjets = new TH1F("mt2_ttwjets_","MT2_ttwjets",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_ttwjets = new  TH1F("mt_ttwjets_", "MT_ttwjets",mt_n_bin,0,mt_h_bin);
	TH1F   *met_ttwjets = new TH1F("met_ttwjets_","met_ttwjets",met_n_bin,0,met_h_bin);

	TH1F *mt2_h_ttzjets = new TH1F("mt2_ttzjets_","MT2_ttzjets",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_ttzjets = new  TH1F("mt_ttzjets_", "MT_ttzjets",mt_n_bin,0,mt_h_bin);
	TH1F   *met_ttzjets = new TH1F("met_ttzjets_","met_ttzjets",met_n_bin,0,met_h_bin);

	TH1F *mt2_h_dy = new TH1F("mt2_dy_","MT2_dy",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_dy = new  TH1F("mt_dy_", "MT_dy",mt_n_bin,0,mt_h_bin);
	TH1F   *met_dy = new TH1F("met_dy_","met_dy",met_n_bin,0,met_h_bin);
*/





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
int met50counter=0;
int vtx10counter=0; 	    
int file_count = 1;

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {
int mt2tailCounter = 0;
int totalEvents = 0;

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
/*
     //jet looper
      for (unsigned int k = 0; k < jets_p4().size(); k++){
        if (jets_p4().at(k).pt()*jets_p4Correction().at(k) < 30) continue;
        if (fabs(jets_p4().at(k).eta()) > 2.5) continue;
        if (ROOT::Math::VectorUtil::DeltaR(jets_p4().at(k), lr_p4()) < 0.4) continue;
        nJets++;
        if (btagDiscriminant().at(k) < 0.244) continue;
        nBtags++;
      }
*/
	if(lr_p4().pt() < 30) continue;

//////////// SELECTION ///////
if(file_count == 1){
met_data->Fill(met());
jets_data->Fill(njets());
	if (abs(lr_id()) == 13){
met_data_mu->Fill(met());
jets_data_mu->Fill(njets());
	}
}

else if(file_count == 2){
met_ttbar->Fill(met(),scale_1fb()*5.2);
jets_ttbar->Fill(njets(),scale_1fb()*5.2);
	if (abs(lr_id()) == 13){
met_ttbar_mu->Fill(met(),scale_1fb()*5.2);
jets_ttbar_mu->Fill(njets(),scale_1fb()*5.2);
	}
}

else if(file_count == 3){
met_wjets->Fill(met(),scale_1fb()*5.2);
jets_wjets->Fill(njets(),scale_1fb()*5.2);
	if (abs(lr_id()) == 13){
met_wjets_mu->Fill(met(),scale_1fb()*5.2);
jets_wjets_mu->Fill(njets(),scale_1fb()*5.2);
	}
}




/*
else if(file_count == 4){
met_wwjets->Fill(met(),scale_1fb()*5.2);
jets_wwjets->Fill(njets);
}
else if(file_count == 5){
met_dy->Fill(met(),scale_1fb()*5.2);
}
else if(file_count == 6){
met_ttwjets->Fill(met(),scale_1fb()*5.2);
}
else if(file_count == 7){
met_ttzjets->Fill(met(),scale_1fb()*5.2);
}
*/

if(met() < 40 ) met50counter++;
if(met() < 40 ) continue;
/*
*/
//////////// SELECTION ///////
/*
/////////////    ADD FAKE    ////////////

//magnitude of the vector
	float r = 40, x = 0, y = 0, phi =0 , theta = 0;
	float a[2];
//Create Lepton fake
	LorentzVector fake_lep;
	do{
//random numbers x,y
	TH1F   *jets_data= new TH1F("jets_data_","Jets_data",jet_n_bin,0,jet_h_bin);
	TH1F   *jets_data= new TH1F("jets_data_","Jets_data",jet_n_bin,0,jet_h_bin);
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
*/
////////  MT //////
float dphi = lr_p4().Phi() - metPhi();
if(dphi > 3.1415) dphi = 6.28 - dphi;
float MT = sqrt(2 * lr_p4().pt() * met() *(1 - cos(dphi))); 
////////  MT //////

////////  MT2 //////
double mt2_event = MT2(new_met.pt(),new_met.Phi(), lr_p4(), fake_lep);
if(mt2_event > 100) mt2tailCounter++;
////////  MT2 //////


////////  FILL //////
if(file_count == 1){
totalEvents++;
mt2_h_data->Fill(mt2_event);
mt_h_data->Fill(MT);

	if (abs(lr_id()) == 13){
mt2_h_data_mu->Fill(mt2_event);
mt_h_data_mu->Fill(MT);
	}

}
else if(file_count == 2){
totalEvents++;
mt2_h_ttbar->Fill(mt2_event,scale_1fb()*5.2);
mt_h_ttbar->Fill(MT,scale_1fb()*5.2);

	if (abs(lr_id()) == 13){
mt2_h_ttbar_mu->Fill(mt2_event,scale_1fb()*5.2);
mt_h_ttbar_mu->Fill(MT,scale_1fb()*5.2);
	}
}
else if(file_count == 3){
totalEvents++;
mt2_h_wjets->Fill(mt2_event,scale_1fb()*5.2);
mt_h_wjets->Fill(MT,scale_1fb()*5.2);

	if (abs(lr_id()) == 13){
mt2_h_wjets_mu->Fill(mt2_event,scale_1fb()*5.2);
mt_h_wjets_mu->Fill(MT,scale_1fb()*5.2);
	}
}

/*
else if(file_count == 4){
mt2_h_wwjets->Fill(mt2_event,scale_1fb()*5.2);
mt_h_wwjets->Fill(MT,scale_1fb()*5.2);
}
else if(file_count == 5){
mt2_h_ttzjets->Fill(mt2_event,scale_1fb()*5.2);
mt_h_ttzjets->Fill(MT,scale_1fb()*5.2);
}
else if(file_count == 6){
mt2_h_ttwjets->Fill(mt2_event,scale_1fb()*5.2);
mt_h_ttwjets->Fill(MT,scale_1fb()*5.2);
}
else if(file_count == 7){
mt2_h_dy->Fill(mt2_event,scale_1fb()*5.2);
mt_h_dy->Fill(MT,scale_1fb()*5.2);
}
*/

	}//event 

////////  DRAW  //////
/*
	TCanvas *c73 = new TCanvas("c73","c73");
	mt_h_data->Draw();
	TCanvas *c74 = new TCanvas("c74","c74");
	mt2_h_data->Draw();
	TCanvas *c75 = new TCanvas("c75","c75");
	met_data->Draw();
*/
if(file_count == 1){
     TFile* fout = new TFile("/home/users/sanil/single/may15hists/data_bv_met40.root","RECREATE");
	mt2_h_data->Write();
	mt_h_data->Write();
	met_data->Write();
	jets_data->Write();

	mt2_h_data_mu->Write();
	mt_h_data_mu->Write();
	met_data_mu->Write();
	jets_data_mu->Write();
	cout << "DATA mt2 tail:    " << mt2tailCounter << endl;
	cout << endl;
	cout << "DATA totalEvents:    " << totalEvents << endl;
	cout << endl;
  fout->Close();
}
else if(file_count == 2){
    TFile* fout = new TFile("/home/users/sanil/single/may15hists/ttbar_bv_met40.root","RECREATE");
	mt2_h_ttbar->Write();
	 mt_h_ttbar->Write();
	  met_ttbar->Write();
	jets_ttbar->Write();

	mt2_h_ttbar_mu->Write();
	 mt_h_ttbar_mu->Write();
	  met_ttbar_mu->Write();
	jets_ttbar_mu->Write();
	cout << "ttbar mt2 tail: " << mt2tailCounter*scale_1fb()*5.2 << endl;
	cout << endl;
	cout << "ttbar totalEvents:    " << totalEvents*scale_1fb()*5.2 << endl;
	cout << endl;
  fout->Close();
}
else if(file_count == 3){
    TFile* fout = new TFile("/home/users/sanil/single/may15hists/wjets_bv_met40.root","RECREATE");
	mt2_h_wjets->Write();
	 mt_h_wjets->Write();
	  met_wjets->Write();
	jets_wjets->Write();

	mt2_h_wjets_mu->Write();
	 mt_h_wjets_mu->Write();
	  met_wjets_mu->Write();
	jets_wjets_mu->Write();
	cout << "wjets mt2 tail: " << mt2tailCounter*scale_1fb()*5.2 << endl;
	cout << endl;
	cout << "wjets totalEvents:    " << totalEvents*scale_1fb()*5.2 << endl;
	cout << endl;
  fout->Close();
}
	/*
else if(file_count == 4){
   TFile* fout = new TFile("/home/users/sanil/single/wwjets_hists_met50.root","RECREATE");
	mt2_h_wwjets->Write();
	 mt_h_wwjets->Write();
	  met_wwjets->Write();
	jets_wwjets->Write();
	cout << "wwjets mt2 tail " << mt2tailCounter*scale_1fb()*5.2 << endl;
	cout << endl;
  fout->Close();
}
else if(file_count == 4){
       TFile* fout = new TFile("/home/users/sanil/single/dy_hists_met50.root","RECREATE");
	mt2_h_dy->Write();
	 mt_h_dy->Write();
	  met_dy->Write();
	cout << "dy mt2 tail " << mt2tailCounter*scale_1fb()*5.2 << endl;
	cout << endl;
  fout->Close();

else if(file_count == 2){
  TFile* fout = new TFile("/home/users/sanil/single/ttzjets_hists_met50.root","RECREATE");
	mt2_h_ttzjets->Write();
	 mt_h_ttzjets->Write();
	  met_ttzjets->Write();
	cout << "ttzjets mt2 tail " << mt2tailCounter*scale_1fb()*5.2 << endl;
	cout << endl;
  fout->Close();
}
else if(file_count == 6){
  TFile* fout = new TFile("/home/users/sanil/single/ttwjets_hists_met50.root","RECREATE");
	mt2_h_ttwjets->Write();
	 mt_h_ttwjets->Write();
	  met_ttwjets->Write();
	cout << "ttwjets mt2 tail " << mt2tailCounter*scale_1fb()*5.2 << endl;
	cout << endl;
  fout->Close();
}
*/


    // Clean Up
    delete tree;
    file->Close();
    delete file;

	file_count++;
  }	//file_loop
//met_data->Draw("P+E");
//met_wjets->Draw("same");
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
