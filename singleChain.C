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
#include "/home/users/sanil/single/NSBABY.h"
#include "/home/users/sanil/CORE/MT2/MT2.h"
#include "/home/users/sanil/CORE/MT2/MT2Utility.h"

#include "/home/users/sanil/CORE/lostlepton/ETHLeptonSelections.h"
#include "trackSelections.h"
#include "muonSelections.h"
#include "susySelections.h"
#include "electronSelections.h"
#include "eventSelections.h"


#include "../myheaders/myheader.h"

using namespace ETHSelection;
using namespace std;
using namespace nsbaby;
using namespace ROOT::Math;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

////////////////////////////////////////////////////////////////////////////////////////////

	TH1F *mt2_h = new TH1F("mt2_"," MT2",26,1,131);
	TH1F *random_phi_h = new TH1F("random phi","random phi_",80,-3.2,3.2);
	TH1F *random_theta_h = new TH1F("random_theta","random_theta _",80,-1,1);
	TH1F *fake_lepton_phi_h = new TH1F("fake_lepton_phi","fake_lepton_phi _",80,-3.2,3.2);
	TH1F *fake_lepton_px_h = new TH1F("fake_lepton_px","fake_lepton_px _",100,0,50);
	TH1F *fake_lepton_py_h = new TH1F("fake_lepton_py","fake_lepton_py _",100,0,50);
	TH1F *fake_lepton_pz_h = new TH1F("fake_lepton_pz","fake_lepton_pz _",100,0,50);
	TH1F *fake_met_phi_h = new TH1F("fake_met_phi","fake_met_phi _",80,-3.2,3.2);
	TH1F *fake_met_mag_h = new TH1F("fake_met_mag","fake_met_mag _",200,0,200);
	TH1F *new_met_phi_h = new TH1F("new_met_phi","new_met_phi _",80,-3.2,3.2);
	TH1F *new_met_mag_h = new TH1F("new_met_mag","new_met_mag _",200,0,200);
	TH1F *real_met_mag_h = new TH1F("real_met_mag","real_met_mag _",200,0,200);
TRandom  *d = new TRandom();

TRandom3 rx;
TRandom3 ry;

	real_met_mag_h	->SetLineColor(kRed);
/////////////////////////////////////////////////////////////////////////////////////////////

  // Loop over events to Analyze
  	unsigned int nEventsTotal = 0;
  	unsigned int nEventsChain = chain->GetEntries();
  	if( nEvents >= 0 ) nEventsChain = nEvents;
  	TObjArray *listOfFiles = chain->GetListOfFiles();
  	TIter fileIter(listOfFiles);
  	TFile *currentFile = 0;
//META
float test_new_met = 0;
float test_new_lep = 0;
int goodEvents = 0;
float scale = 0;
  // File Loop
	typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
  	while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    	TFile *file = new TFile( currentFile->GetTitle() );
    	TTree *tree = (TTree*)file->Get("tree");
    	if(fast) TTreeCache::SetLearnEntries(10);
    	if(fast) tree->SetCacheSize(128*1024*1024);
    	baby.Init(tree);
    
    // Loop over Events in current file
	if( nEventsTotal >= nEventsChain ) continue;
	unsigned int nEventsTree = tree->GetEntriesFast();
  	for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      baby.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      CMS2::progress( nEventsTotal, nEventsChain );

	//event based variables
	int mu_index = -1;
	int el_index = -1;
	float mu_max = 0.0;
	float el_max = 0.0;
	float lep_fake_phi = 0;
	float xc = 0 , yc = 0, zc = 0;
	//COUNTERS
	int el_count = 0;
	int mu_count = 0;
	int nJets = 0;
	int nBtags = 0;
int good_index = -1;
//scale = evt_scale1fb();
/*
// MUON SELECTION
      for (unsigned int index = 0; index< mus_p4().size(); index++){

  if (cms2.mus_p4()[index].pt() < 20.0)  continue;
  if (fabs(cms2.mus_p4()[index].eta()) > 2.4)  continue;

  const bool is_global  = ((cms2.mus_type().at(index) & (1<<1)) != 0);
  const bool is_pfmu    = ((cms2.mus_type().at(index) & (1<<5)) != 0);

  const int vtxidx = firstGoodVertex();

  if (!is_global)  continue;
  if (!is_pfmu)  continue;
  if (cms2.mus_gfit_validSTAHits().at(index) < 1)  continue;
  if (cms2.mus_numberOfMatchedStations().at(index) < 2)  continue;

  const int ctfidx = cms2.mus_trkidx().at(index);
  if (ctfidx < 0 || vtxidx < 0)  continue;
  const std::pair<double, double> cord0 = trks_d0_pv(ctfidx, vtxidx);
	el_count++;
  const std::pair<double, double> cordz = trks_dz_pv(ctfidx, vtxidx);
  if (fabs(cord0.first) > 0.2)  continue;
  if (fabs(cordz.first) > 0.5)  continue;
  if (cms2.trks_valid_pixelhits().at(ctfidx) < 1)  continue;
  if (cms2.trks_nlayers().at(ctfidx) < 6)  continue;

  if ( cms2.mus_gfit_chi2()[index] / cms2.mus_gfit_ndof()[index] >= 10 )  continue;
  if (muonIsoValuePF2012_deltaBeta(index) > 0.2)  continue;
	mu_count++;
	good_index = index;
	        }

// ELECTRON SELECTION
      for (unsigned int index = 0; index< els_p4().size(); index++){

if (cms2.els_p4()[index].pt() < 20.0)  continue;
  if (fabs(cms2.els_p4()[index].eta()) > 2.4)  continue;
  if (fabs(cms2.els_p4()[index].eta()) > 1.442 && fabs(cms2.els_p4()[index].eta()) < 1.566)  continue;
  if (fabs(electron_d0PV(index)) > 0.04)  continue;
  if (fabs(electron_dzPV_smurfV3(index)) > 0.2 )  continue;
  //if (electronIsoValuePF(index, 0) > 0.15)  continue;
  if(electronIsoValuePF2012_FastJetEffArea_v3(index) > 0.15)  continue;

  if (fabs(cms2.els_p4()[index].eta()) < 1.479) { //barrel
    if (cms2.els_dEtaIn()[index] > 0.007)  continue;
    if (cms2.els_dPhiIn()[index] > 0.8)  continue;
    if (cms2.els_sigmaIEtaIEta()[index] > 0.01)  continue;
    if (cms2.els_hOverE()[index] > 0.15)  continue;
  }
  else { //endcap
    if (cms2.els_dEtaIn()[index] > 0.01)  continue;
    if (cms2.els_dPhiIn()[index] > 0.7)  continue;
    if (cms2.els_sigmaIEtaIEta()[index] > 0.03)  continue;
	        }

	el_count++;
	good_index = index;
	}
	if(good_index == -1) continue;
if ( (mu_count + el_count) != 1) continue;
goodEvents++;

if(el_count == 1) lepton = els_p4().at(good_index);
if(mu_count == 1) lepton = mus_p4().at(good_index);
*/
LorentzVector lepton;
lepton = lr_p4();


	         //jet looper
      for (unsigned int k = 0; k < jets_p4().size(); k++){
        if (jets_p4().at(k).pt()*jets_p4Correction().at(k) < 30) continue;
        if (fabs(jets_p4().at(k).eta()) > 2.5) continue;
        if (ROOT::Math::VectorUtil::DeltaR(jets_p4().at(k), lepton) < 0.4) continue;
        nJets++;
        if (btagDiscriminant().at(k) < 0.244) continue;
        nBtags++;
      }


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
//figure out the phi angle of this fake lepton 
//phi is -pi to pi

lep_fake_phi = fake_lep.Phi();

//add neutrino by changing metPhi and magnitude of met
//calculate the phi of the fake neutrino from the fake lepton phi
//phi is -pi to pi
float fake_metPhi = -lep_fake_phi;
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
///////////////////////////////////////////////////
//FILL//

	random_phi_h	->Fill(phi);
	random_theta_h	->Fill(cos(theta));       
fake_lepton_phi_h	->Fill(fake_lep.Phi());
	fake_lepton_px_h->Fill(fake_lep.Px());
	fake_lepton_py_h->Fill(fake_lep.Py());
	fake_lepton_pz_h->Fill(fake_lep.Pz());

	fake_met_phi_h 	->Fill(fake_met.Phi());
	fake_met_mag_h	->Fill(fake_met.pt());
	new_met_phi_h 	->Fill(new_met.Phi());
	new_met_mag_h	->Fill(new_met.pt());


double mt2_event = MT2(new_met.pt(),new_met.Phi(), lepton, fake_lep);

mt2_h->Fill(mt2_event);
real_met_mag_h->Fill(met());
	} //loop over events in the current file

    
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  } //file loop
  
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
		}
	TCanvas *c73 = new TCanvas("c73","c73");
	mt2_h->Draw();
	TCanvas *c77 = new TCanvas("c77","c77");
	fake_lepton_px_h->Draw();
	TCanvas *c78 = new TCanvas("c78","c78");
	fake_lepton_py_h->Draw();
	TCanvas *c8 = new TCanvas("c8","c8");
	fake_lepton_pz_h->Draw();
	TCanvas *c74 = new TCanvas("c74","c74");
	random_phi_h	->Draw();
	TCanvas *c75 = new TCanvas("c75","c75");
	random_theta_h	->Draw();
	TCanvas *c76 = new TCanvas("c76","c76");
	fake_lepton_phi_h->Draw();
	TCanvas *c79 = new TCanvas("c79","c79");
	fake_lepton_pz_h->Draw();
	TCanvas *c72 = new TCanvas("c72","c72");
	fake_met_phi_h 	->Draw();
	TCanvas *c71 = new TCanvas("c71","c71");
	fake_met_mag_h	->Draw();
	TCanvas *c7 = new TCanvas("c7","c7");
	new_met_phi_h 	->Draw();
	TCanvas *c3 = new TCanvas("c3","c3");
	new_met_mag_h	->Draw();
	real_met_mag_h	->Draw("same");

  TFile* fout = new TFile("/home/users/sanil/single/single_hists_wjets.root","RECREATE");
	mt2_h		->Write();
	random_phi_h	->Write();
	random_theta_h	->Write();    
fake_lepton_phi_h	->Write();
	fake_lepton_px_h->Write();
	fake_lepton_py_h->Write();
	fake_lepton_pz_h->Write();
	fake_met_phi_h 	->Write();
	fake_met_mag_h	->Write();
	new_met_phi_h 	->Write();
	real_met_mag_h  ->Write();
	new_met_mag_h	->Write();
  fout->Close();


/*

	TCanvas *c3 = new TCanvas("c3","c3");
	new_met_h->Draw();
	//met_h->Draw("same");
	fake_px->Draw();
	TCanvas *c72 = new TCanvas("c72","c72");
	fake_py->Draw();
	fake_pz->Draw();
	TCanvas *c71 = new TCanvas("c71","c71");
	sphere_x->Draw();
	TCanvas *c55 = new TCanvas("c55","c55");
	sphere_y->Draw();
	TCanvas *c70 = new TCanvas("c70","c70");
	sphere_z->Draw();
	TCanvas *c71 = new TCanvas("c71","c71");
fake_lep_phi_h->Draw();
	TCanvas *c70 = new TCanvas("c70","c70");
fake_met_phi_h->Draw();
*/
	// return
//cout << goodEvents*scale << endl;
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << "Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:   " << Form( "%.0001f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:   " << Form( "%.0001f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return 0;
 }

