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
#include "TRandom.h"

// CMS2
#include "../CORE/CMS2.h"
#include "../CORE/MT2/MT2.h"
#include "../CORE/MT2/MT2Utility.h"
#include "muonSelections.h"
#include "susySelections.h"
#include "electronSelections.h"
#include "../myheaders/myheader.h"

using namespace std;
using namespace tas;
using namespace ROOT::Math;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

////////////////////////////////////////////////////////////////////////////////////////////
	TH1F *numberjets = new TH1F("numberjes ","numberjets",10,0,10);
	TH1F *lepton_pt = new TH1F("lepton_pt","lpton pt",90,20,200);
	TH1F *fake_px = new TH1F("fake_px","fake_px",40,0,40);
	TH1F *fake_py = new TH1F("fake_py","fake_py",40,0,40);
	TH1F *fake_pz = new TH1F("fake_pz","fake_pz",40,0,40);
	TH1F *el_number = new TH1F("# of mu","# of mu",10,0,10);
	TH1F *mu_number = new TH1F("# of el","# of el",10,0,10);
	TH1F *mt2_h 	= new TH1F("mt2 ","MT2",100,0,200);
	TH1F *fake_phi_h 	= new TH1F("fakePhi ","fake Phi  ",70,0,4);

TRandom  *d = new TRandom();

/////////////////////////////////////////////////////////////////////////////////////////////

  // Loop over events to Analyze
  	unsigned int nEventsTotal = 0;
  	unsigned int nEventsChain = chain->GetEntries();
  	if( nEvents >= 0 ) nEventsChain = nEvents;
  	TObjArray *listOfFiles = chain->GetListOfFiles();
  	TIter fileIter(listOfFiles);
  	TFile *currentFile = 0;
//META COUNTERS
float test_new_met = 0;
float test_new_lep = 0;
  // File Loop
	typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
  	while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    	TFile *file = new TFile( currentFile->GetTitle() );
    	TTree *tree = (TTree*)file->Get("Events");
    	if(fast) TTreeCache::SetLearnEntries(10);
    	if(fast) tree->SetCacheSize(128*1024*1024);
    	cms2.Init(tree);
    
    // Loop over Events in current file
	if( nEventsTotal >= nEventsChain ) continue;
	unsigned int nEventsTree = tree->GetEntriesFast();
  	for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      cms2.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      CMS2::progress( nEventsTotal, nEventsChain );

	//event based variables
	int mu_index = -1;
	int el_index = -1;
	float mu_max = 0.0;
	float el_max = 0.0;
	float lep_fake_phi = 0;
	//COUNTERS
	int el_count = 0;
	int mu_count = 0;
	int nJets = 0;
	int nBtags = 0;
// MUON SELECTION
      for (unsigned int i = 0; i< mus_p4().size(); i++){
	mu_count++;
	  //pT
          if (mus_p4().at(i).pt() < 20) continue;
	  //eta
          if (fabs(mus_p4().at(i).eta()) > 2.4) continue;

          // select the highest pt hypothesis
          float curPt = mus_p4().at(i).pt();
          if (curPt > mu_max){
              mu_max = curPt;
              mu_index = i;
          }
      }

// ELECTRON SELECTION
      for (unsigned int i = 0; i< els_p4().size(); i++){
	el_count++;
	  //pT
          if (els_p4().at(i).pt() < 20) continue;
	  //eta
          if (fabs(els_p4().at(i).eta()) > 2.4) continue;

          // select the highest pt hypothesis
          float currPt = els_p4().at(i).pt();
          if (currPt > el_max){
              el_max = currPt;
              el_index = i;
          }
      }

el_number->Fill(el_count);
mu_number->Fill(mu_count);

//skip if no good lepton
if((el_index == -1) && (mu_index == -1)) continue;
//skip if there are two good leptons
if((mu_index != -1) && (el_index != -1)) continue;
LorentzVector lepton;
if(el_index > -1) lepton = els_p4().at(el_index);
else if(mu_index > -1) lepton = mus_p4().at(mu_index);


	         //jet looper
      for (unsigned int k = 0; k < pfjets_p4().size(); k++){
        if (pfjets_p4().at(k).pt()*pfjets_corL1FastL2L3().at(k) < 30) continue;
        if (fabs(pfjets_p4().at(k).eta()) > 2.5) continue;
        if (ROOT::Math::VectorUtil::DeltaR(pfjets_p4().at(k), lepton) < 0.4) continue;
        nJets++;
        if (pfjets_combinedSecondaryVertexBJetTag().at(k) < 0.244) continue;
        nBtags++;
      }
numberjets->Fill(nJets);
lepton_pt->Fill(lepton.pt());



//SPHERE PICKING Sphere() from TRandom.h
double x, y, z;
//radius 40
d->Sphere(x,y,z,40);

//Create Lepton fake
LorentzVector fake_lep;
//set px py pz and E=40
fake_lep.SetPxPyPzE(x,y,z,40); 
//figure out the phi angle of this fake lepton
lep_fake_phi = fake_lep.Phi();


//add neutrino by changing metPhi and magnitude of met
float metPhi = evt_pfmetPhi_type1cor();
float met = evt_pfmet_type1cor();
//calculate the phi of the fake neutrino from the fake lepton phi
float fake_metPhi = wrapAnglePi(3.1415-lep_fake_phi);
//now vector addition using fake_metPhi, metPhi, met, and 40(momentum of the new neutrino)
float metx = cos(metPhi)*met;
float mety = sin(metPhi)*met;

float fake_met = 40;

float fake_metx = cos(fake_metPhi)*fake_met;
float fake_mety = sin(fake_metPhi)*fake_met;

float new_metx = metx + fake_metx;
float new_mety = mety + fake_mety;
//calculate the updated met
float new_met = 0;
new_met = sqrt(new_metx * new_metx + new_mety*new_mety);
float new_metPhi = wrapAnglePi(fake_metPhi + metPhi);
///////////////////////////////////////////////////
fake_phi_h->Fill(new_metPhi);
/*
fake_px->Fill(fake_lep.Px());
fake_py->Fill(fake_lep.Py());
fake_pz->Fill(fake_lep.Pz());
*/

test_new_met = new_metPhi;
test_new_lep = fake_lep.Phi();
	} //loop over events in the current file

    
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  } //file loop
  
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
		}
/*
	TCanvas *c73 = new TCanvas("c73","c73");
	fake_px->Draw();
	TCanvas *c72 = new TCanvas("c72","c72");
	fake_py->Draw();
	TCanvas *c71 = new TCanvas("c71","c71");
	fake_pz->Draw();
*/
	TCanvas *c70 = new TCanvas("c70","c70");
	fake_phi_h->Draw();
	cout << "new met phi "<< test_new_met << endl;
	cout << "new lep phi "<< test_new_lep << endl;
	// return
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

