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
#include "TRandom1.h"
#include "TRandom.h"
#include "TLorentzVector.h"

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


#include "/home/users/sanil/myheaders/myheader.h"
#include "/home/users/sanil/myheaders/triggerWeights.h"

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

	TH1F   *fake_wpt_h = new TH1F("fake_w_pt_","fake_W_pt",50,0,100);
	TH1F   *landau_h = new TH1F("landau_h_","landau_H",60,0,mt_h_bin);
	TH1F   *gammaM_h = new TH1F("gammaM_","gamma",120,1,6);

	TH1F *mt2_h_data = new TH1F("mt2_data_","MT2_data",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_data = new  TH1F("mt_data_", "MT_data",mt_n_bin,0,mt_h_bin);
	TH1F   *met_data = new TH1F("met_data_","met_data",met_n_bin,0,met_h_bin);
	TH1F   *jets_data= new TH1F("jets_data_","Jets_data",jet_n_bin,0,jet_h_bin);
	TH1F   *data_w_pt = new     TH1F("data_w_pt_","data_W_pt",60,0,mt_h_bin);
	TH1F   *data_lep_pt = new     TH1F("data_lep_pt_","data_Lep_pt",60,0,mt_h_bin);

	TH1F *mt2_h_ttbar = new TH1F("mt2_ttbar_","MT2_ttbar",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_ttbar = new  TH1F("mt_ttbar_", "MT_ttbar",mt_n_bin,0,mt_h_bin);
	TH1F   *met_ttbar = new TH1F("met_ttbar_","met_ttbar",met_n_bin,0,met_h_bin);
	TH1F   *jets_ttbar  = new TH1F("jets_ttbar_","Jets_ttbar",jet_n_bin,0,jet_h_bin);
	TH1F   *ttbar_w_pt = new     TH1F("ttbar_w_pt_","ttbar_W_pt",60,0,mt_h_bin);
	TH1F   *ttbar_lep_pt = new     TH1F("ttbar_lep_pt_","ttbar_Lep_pt",60,0,mt_h_bin);

	TH1F *mt2_h_wjets = new TH1F("mt2_wjets_","MT2_wjets",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_wjets = new  TH1F("mt_wjets_", "MT_wjets",mt_n_bin,0,mt_h_bin);
	TH1F   *met_wjets = new TH1F("met_wjets_","met_wjets",met_n_bin,0,met_h_bin);
	TH1F   *jets_wjets  = new TH1F("jets_wjets_","Jets_wjets",jet_n_bin,0,jet_h_bin);
	TH1F   *wjets_w_pt = new     TH1F("wjets_w_pt_","wjets_W_pt",60,0,mt_h_bin);
	TH1F   *wjets_lep_pt = new     TH1F("wjets_lep_pt_","wjets_Lep_pt",60,0,mt_h_bin);

	wjets_lep_pt->SetLineColor(kRed);
	data_lep_pt->SetLineColor(kYellow);
	ttbar_lep_pt->SetLineColor(kCyan);

	wjets_lep_pt->SetLineWidth(4);
	 data_lep_pt->SetLineWidth(4);
	ttbar_lep_pt->SetLineWidth(4);

////MUONS////
	TH1F *mt2_h_data_mu = new TH1F("mt2_data_mu_","MT2_data_mu",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_data_mu = new  TH1F("mt_data_mu_", "MT_data_mu",mt_n_bin,0,mt_h_bin);
	TH1F   *met_data_mu = new TH1F("met_data_mu_","met_data_mu",met_n_bin,0,met_h_bin);
	TH1F   *jets_data_mu= new TH1F("jets_data_mu_","Jets_data_mu",jet_n_bin,0,jet_h_bin);
	TH1F   *data_w_pt_mu = new     TH1F("data_w_pt_mu_","data_W_pt_mu",60,0,mt_h_bin);
	TH1F   *data_lep_pt_mu = new     TH1F("data_lep_pt_mu_","data_Lep_pt_mu",60,0,mt_h_bin);

	TH1F *mt2_h_ttbar_mu = new TH1F("mt2_ttbar_mu_","MT2_ttbar_mu",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_ttbar_mu = new  TH1F("mt_ttbar_mu_", "MT_ttbar_mu",mt_n_bin,0,mt_h_bin);
	TH1F   *met_ttbar_mu = new TH1F("met_ttbar_mu_","met_ttbar_mu",met_n_bin,0,met_h_bin);
	TH1F   *jets_ttbar_mu  = new TH1F("jets_ttbar_mu_","Jets_ttbar_mu",jet_n_bin,0,jet_h_bin);
	TH1F   *ttbar_w_pt_mu = new     TH1F("ttbar_w_pt_mu_","ttbar_W_pt_mu",60,0,mt_h_bin);
	TH1F   *ttbar_lep_pt_mu = new     TH1F("ttbar_lep_pt_mu_","ttbar_Lep_pt_mu",60,0,mt_h_bin);

	TH1F *mt2_h_wjets_mu = new TH1F("mt2_wjets_mu_","MT2_wjets_mu",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_wjets_mu = new  TH1F("mt_wjets_mu_", "MT_wjets_mu",mt_n_bin,0,mt_h_bin);
	TH1F   *met_wjets_mu = new TH1F("met_wjets_mu_","met_wjets_mu",met_n_bin,0,met_h_bin);
	TH1F   *jets_wjets_mu  = new TH1F("jets_wjets_mu_","Jets_wjets_mu",jet_n_bin,0,jet_h_bin);
	TH1F   *wjets_w_pt_mu = new     TH1F("wjets_w_pt_mu_","wjets_W_pt_mu",60,0,mt_h_bin);
	TH1F   *wjets_lep_pt_mu = new     TH1F("wjets_lep_pt_mu_","wjets_Lep_pt_mu",60,0,mt_h_bin);

////ELECTRONS////
	TH1F *mt2_h_data_el = new TH1F("mt2_data_el_","MT2_data_el",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_data_el = new  TH1F("mt_data_el_", "MT_data_el",mt_n_bin,0,mt_h_bin);
	TH1F   *met_data_el = new TH1F("met_data_el_","met_data_el",met_n_bin,0,met_h_bin);
	TH1F   *jets_data_el= new TH1F("jets_data_el_","Jets_data_el",jet_n_bin,0,jet_h_bin);
	TH1F   *data_w_pt_el = new     TH1F("data_w_pt_el_","data_W_pt_el",60,0,mt_h_bin);
	TH1F   *data_lep_pt_el = new     TH1F("data_lep_pt_el_","data_Lep_pt_el",60,0,mt_h_bin);

	TH1F *mt2_h_ttbar_el = new TH1F("mt2_ttbar_el_","MT2_ttbar_el",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_ttbar_el = new  TH1F("mt_ttbar_el_", "MT_ttbar_el",mt_n_bin,0,mt_h_bin);
	TH1F   *met_ttbar_el = new TH1F("met_ttbar_el_","met_ttbar_el",met_n_bin,0,met_h_bin);
	TH1F   *jets_ttbar_el  = new TH1F("jets_ttbar_el_","Jets_ttbar_el",jet_n_bin,0,jet_h_bin);
	TH1F   *ttbar_w_pt_el = new     TH1F("ttbar_w_pt_el_","ttbar_W_pt_el",60,0,mt_h_bin);
	TH1F   *ttbar_lep_pt_el = new     TH1F("ttbar_lep_pt_el_","ttbar_Lep_pt_el",60,0,mt_h_bin);

	TH1F *mt2_h_wjets_el = new TH1F("mt2_wjets_el_","MT2_wjets_el",mt2_n_bin,0,mt2_h_bin);
	TH1F  *mt_h_wjets_el = new  TH1F("mt_wjets_el_", "MT_wjets_el",mt_n_bin,0,mt_h_bin);
	TH1F   *met_wjets_el = new TH1F("met_wjets_el_","met_wjets_el",met_n_bin,0,met_h_bin);
	TH1F   *jets_wjets_el  = new TH1F("jets_wjets_el_","Jets_wjets_el",jet_n_bin,0,jet_h_bin);
	TH1F   *wjets_w_pt_el = new     TH1F("wjets_w_pt_el_","wjets_W_pt_el",60,0,mt_h_bin);
	TH1F   *wjets_lep_pt_el = new     TH1F("wjets_lep_pt_el_","wjets_Lep_pt_el",60,0,mt_h_bin);
/*
  TFile* f1 = new TFile("/home/users/sanil/single/may15hists/data_w_pt.root");
  TH1F* wp_dist   = (TH1F*) f1->Get("w_pt_"); // 
  f1->Close();
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


////////////// RANDOM NUMBERS GET W momentum


TRandom1  d;  
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
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVec;
//GLOBAL COUNTERS//
int file_count = 1;

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

float mt2_0_counter = 0, //0 bin counter
     mt2_10_counter = 0, //0-10 geV counter
     mt2_20_counter = 0, //10-20 geV counter etc.
     mt2_30_counter = 0,
     mt2_40_counter = 0,
     mt2_50_counter = 0,
     mt2_60_counter = 0,
     mt2_70_counter = 0,
     mt2_80_counter = 0,
     mt2_90_counter = 0,
    mt2_100_counter = 0,
    mt2_110_counter = 0,
    mt2_120_counter = 0,
    mt2_130_counter = 0,
    mt2_140_counter = 0,
    mt2_150_counter = 0,
    mt2_160_counter = 0,
    mt2_170_counter = 0,
    mt2_180_counter = 0,
    mt2_190_counter = 0,
    mt2_200_counter = 0,
    mt2_G200_counter = 0,
    totalEvents = 0;

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

	float lepton_weight = 0;
	
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      single.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      SL::progress( nEventsTotal, nEventsChain );

      // Analysis Code
//event based variables
	//COUNTERS
	int el_count = 0, mu_count = 0, nJets = 0, nBtags = 0, good_index = -1, mu_index = -1, el_index = -1;

      float metx = met() * cos( metPhi() );
      float mety = met() * sin( metPhi() );

///metPhi correction//

      float cmetx = met() * cos( metPhi() );
      float cmety = met() * sin( metPhi() );
      float shiftx = 0.;
      float shifty = 0.;

      shiftx = (! isRealData()) ? (0.1166 + 0.0200*nvtxs()) : (0.2661 + 0.3217*nvtxs());
      shifty = (! isRealData()) ? (0.2764 - 0.1280*nvtxs()) : (-0.2251 - 0.1747*nvtxs());

      cmetx -= shiftx;
      cmety -= shifty;

	//problem with this metPhi modulation correction right now
	//keep using cmet, too lazy to change all the variables
	float cmet = met();
     	float cmetphi = metPhi();
	 //float cmetphi = atan2( mety , metx );
      //float cmet = sqrt( cmetx*cmetx + cmety*cmety ); // cmet == corrected met
LorentzVector met_vec;
met_vec.SetPxPyPzE(metx,mety,0,0);

float lepy = lr_p4().py();
float lepx = lr_p4().px();

float w_x = sqrt ( metx*metx + lepx*lepx );
float w_y = sqrt ( mety*mety + lepy*lepy );

float w_pt_ = sqrt (w_x*w_x + w_y*w_y); 


//jet looper
      for (unsigned int k = 0; k < jets_p4().size(); k++){
        if (jets_p4().at(k).pt()*jets_p4Correction().at(k) < 30) continue;
        if (fabs(jets_p4().at(k).eta()) > 2.5) continue;
        if (ROOT::Math::VectorUtil::DeltaR(jets_p4().at(k), lr_p4()) < 0.4) continue;
        nJets++;
        if (btagDiscriminant().at(k) < 0.244) continue;
        nBtags++;
      }

       //metPhi Correction//
//////////// SELECTION ///////
//Event Requirements
//muon eta > 2.1
//electron eta > 2.4
if(lr_p4().pt() < 30) continue;
if(nJets < 2) continue;
if(nBtags != 0) continue;
if(cmet < 40 ) continue;
//Event Requirements

//TRIGGER WEIGH//
if(abs(lr_id()) == 13){
lepton_weight = electronTriggerWeight(lr_p4().pt(),lr_p4().eta());
}
else if(abs(lr_id()) == 11){ 
lepton_weight = muonTriggerWeight(lr_p4().pt(),lr_p4().eta());
}

//float w_pt_ = (lr_p4() + met_vec).pt();

//met FILL
if(file_count == 1){
met_data->Fill(cmet);
	if (abs(lr_id()) == 13){
met_data_mu->Fill(cmet);
	}
	if (abs(lr_id()) == 11){
met_data_el->Fill(cmet);
	}
}

else if(file_count == 2){
met_ttbar->Fill(cmet,scale_1fb()*lepton_weight*5.2);
	if (abs(lr_id()) == 13){
met_ttbar_mu->Fill(cmet,scale_1fb()*lepton_weight*5.2);
	}
	if (abs(lr_id()) == 11){
met_ttbar_el->Fill(cmet,scale_1fb()*lepton_weight*5.2);
	}
}

else if(file_count == 3){
met_wjets->Fill(cmet,scale_1fb()*lepton_weight*5.2);
	if (abs(lr_id()) == 13){
met_wjets_mu->Fill(cmet,scale_1fb()*lepton_weight*5.2);
	}
	if (abs(lr_id()) == 11){
met_wjets_el->Fill(cmet,scale_1fb()*lepton_weight*5.2);
	}
}
/*
		\\\~*~*~*~\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\|||||/////////////////////////////////////~*~*~*~///
		///~*~*~*~///////////////////////////////////////|||||\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\~*~*~*~\\\
*/
/////////////    ADD FAKE    ////////////
//W FAKE
//W vector variables
float xw = 0, yw = 0, phiw = 0 , thetaw = 0;
float w_p = -1, w_px = -1, w_py = -1, w_pz = -1;
float b[2];

//lepton vector variables
float xc = 0 , yc = 0, zc = 0;
float r = 40, x = 0, y = 0, phi =0 , theta = 0;
float a[2];
float gammaM;
	//random w momentum from the histogram
	//w_p = d.Landau(95.0291,10.3062);
float rand_pt = d.Landau(40.5,17);

	w_p = rand_pt;

	ry.RndmArray(2,b);
	xw = b[0];
	yw = b[1];

//angles for W vector intial direction
	thetaw = acos(2*xw-1);	   
 	phiw = 2*3.14156265359*(yw-0.5);       

//Change of coordinate to XYZ and scale by magnitude of p
      w_px = w_p* sin(thetaw)* cos(phiw); 
      w_py = w_p* sin(thetaw)* sin(phiw);
      w_pz = w_p* cos(thetaw);
//LEPTON
//Create Lepton fake
	LorentzVector fake_lep, fake_met;
//use TLorentVector to boost
	TLorentzVector lep_boost, met_boost;

//do while eta > 2.4 for the lepton
	do{
//random numbers x,y
//rx TRandom class
	rx.RndmArray(2,a);
	x = a[0];
	y = a[1];

//calculate phi and theta
	phi = (x-0.5)*2*3.14156265359;
	theta = acos(2*y - 1);

//calculate cartesian coordinates of the momentum
	 xc = r*sin(theta)*cos(phi);
	 yc = r*sin(theta)*sin(phi);
	 zc = r*cos(theta);

//lepton momentum and met momentum should be oppositely signed
	lep_boost.SetPxPyPzE(xc,yc,zc,40); 
	met_boost.SetPxPyPzE(-xc,-yc,-zc,40);

//BOOST
//gammaM is gamma * mass, reduces to sqrt(momentum^2 + mass^2)
	float gammaM_ = sqrt(w_p*w_p + 80.2*80.2);
	gammaM = gammaM_;

//Boost lepton
	lep_boost.Boost(w_px / (gammaM), w_py / (gammaM), w_pz / (gammaM) );

		} while  (lep_boost.Eta() > 2.4 ); //make sure lepton is within the eta requirement

//boost the neutrino
	met_boost.Boost(w_px / (gammaM), w_py / (gammaM), w_pz / (gammaM) );

//switch from boost-able TLorentzVector into LorentzVector, 
//some things work different in TLorentzVector
	fake_lep.SetPxPyPzE(lep_boost.Px(),lep_boost.Py(),lep_boost.Pz(),0);
	fake_met.SetPxPyPzE(met_boost.Px(),met_boost.Py(),met_boost.Pz(),0);

//now vector addition using fake_metPhi, metPhi, met, and 40(momentum of the new neutrino)
	float met_x = cos(cmetphi)*cmet;
	float met_y = sin(cmetphi)*cmet;

//real met vector
	LorentzVector real_met;
	real_met.SetPxPyPzE(met_x,met_y,0,0); 

//new met vector
	LorentzVector new_met;
	new_met = real_met + fake_met;
///////////////  END FAKE   ////////////////////////////////////
float fake_W_pt = (fake_lep+fake_met).pt();
if(file_count == 1) fake_wpt_h->Fill(fake_W_pt);
/*
		\\\~*~*~*~\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\|||||/////////////////////////////////////~*~*~*~///
		///~*~*~*~///////////////////////////////////////|||||\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\~*~*~*~\\\
*/

////////  MT //////
float dphi = lr_p4().Phi() - cmetphi;
if(dphi > 3.1415) dphi = 6.28 - dphi;
float MT = sqrt(2 * lr_p4().pt() * cmet *(1 - cos(dphi))); 
////////  MT //////

double mt2_event = MT2(new_met.Pt(),new_met.Phi(), lr_p4(), fake_lep);
////////  MT2 //////
if(file_count == 1){
if(mt2_event == 0)				  mt2_0_counter+=1;
else if(0 < mt2_event && mt2_event <= 10)	 mt2_10_counter+=1; 
else if(10 < mt2_event && mt2_event <= 20)	 mt2_20_counter+=1; 
else if(20 < mt2_event && mt2_event <= 30)	 mt2_30_counter+=1;
else if(30 < mt2_event && mt2_event <= 40)	 mt2_40_counter+=1; 
else if(40 < mt2_event && mt2_event <= 50)	 mt2_50_counter+=1;
else if(50 < mt2_event && mt2_event <= 60)	 mt2_60_counter+=1;
else if(60 < mt2_event && mt2_event <= 70)	 mt2_70_counter+=1;
else if(70 < mt2_event && mt2_event <= 80)	 mt2_80_counter+=1;
else if(80 < mt2_event && mt2_event <= 90)	 mt2_90_counter+=1;
else if(90 < mt2_event && mt2_event <= 100)	mt2_100_counter+=1;
else if(100 < mt2_event && mt2_event <= 110)	mt2_110_counter+=1;
else if(110 < mt2_event && mt2_event <= 120)	mt2_120_counter+=1;
else if(120 < mt2_event && mt2_event <= 130)	mt2_130_counter+=1;
else if(130 < mt2_event && mt2_event <= 140)	mt2_140_counter+=1;
else if(140 < mt2_event && mt2_event <= 150)	mt2_150_counter+=1;
else if(150 < mt2_event && mt2_event <= 160)	mt2_160_counter+=1;
else if(160 < mt2_event && mt2_event <= 170)	mt2_170_counter+=1;
else if(170 < mt2_event && mt2_event <= 180)	mt2_180_counter+=1;
else if(180 < mt2_event && mt2_event <= 190)	mt2_190_counter+=1;
else if(190 < mt2_event && mt2_event <= 200)	mt2_200_counter+=1;
else if( mt2_event > 200)		       mt2_G200_counter+=1;
	} //data

if(file_count == 2 || file_count == 3){
if(mt2_event == 0)			   	  mt2_0_counter += 1*lepton_weight;
else if(0 < mt2_event && mt2_event <= 10)	 mt2_10_counter += 1*lepton_weight; 
else if(10 < mt2_event && mt2_event <= 20)	 mt2_20_counter += 1*lepton_weight; 
else if(20 < mt2_event && mt2_event <= 30)	 mt2_30_counter += 1*lepton_weight;
else if(30 < mt2_event && mt2_event <= 40)	 mt2_40_counter += 1*lepton_weight; 
else if(40 < mt2_event && mt2_event <= 50)	 mt2_50_counter += 1*lepton_weight;
else if(50 < mt2_event && mt2_event <= 60)	 mt2_60_counter += 1*lepton_weight;
else if(60 < mt2_event && mt2_event <= 70)	 mt2_70_counter += 1*lepton_weight;
else if(70 < mt2_event && mt2_event <= 80)	 mt2_80_counter += 1*lepton_weight;
else if(80 < mt2_event && mt2_event <= 90)	 mt2_90_counter += 1*lepton_weight;
else if(90 < mt2_event && mt2_event <= 100)	mt2_100_counter += 1*lepton_weight;
else if(100 < mt2_event && mt2_event <= 110)	mt2_110_counter += 1*lepton_weight;
else if(110 < mt2_event && mt2_event <= 120)	mt2_120_counter += 1*lepton_weight;
else if(120 < mt2_event && mt2_event <= 130)	mt2_130_counter += 1*lepton_weight;
else if(130 < mt2_event && mt2_event <= 140)	mt2_140_counter += 1*lepton_weight;
else if(140 < mt2_event && mt2_event <= 150)	mt2_150_counter += 1*lepton_weight;
else if(150 < mt2_event && mt2_event <= 160)	mt2_160_counter += 1*lepton_weight;
else if(160 < mt2_event && mt2_event <= 170)	mt2_170_counter += 1*lepton_weight;
else if(170 < mt2_event && mt2_event <= 180)	mt2_180_counter += 1*lepton_weight;
else if(180 < mt2_event && mt2_event <= 190)	mt2_190_counter += 1*lepton_weight;
else if(190 < mt2_event && mt2_event <= 200)	mt2_200_counter += 1*lepton_weight;
else if( mt2_event > 200)		       mt2_G200_counter += 1*lepton_weight;
	} //MC 
////////  MT2 //////
float lep_pt = lr_p4().pt();

////////  FILL //////
// JETS, MT, MT2, W_pT
//DATA
if(file_count == 1){
totalEvents+=1;
	data_w_pt->Fill(w_pt_);
	data_lep_pt->Fill(lep_pt);
	jets_data->Fill(nJets);
	mt2_h_data->Fill(mt2_event);
	mt_h_data->Fill(MT);
if (abs(lr_id()) == 13){
	data_lep_pt_mu->Fill(lep_pt);
	mt2_h_data_mu->Fill(mt2_event);
	mt_h_data_mu->Fill(MT);
	jets_data_mu->Fill(nJets);
	}
if (abs(lr_id()) == 11){
	data_lep_pt_el->Fill(lep_pt);
	mt2_h_data_el->Fill(mt2_event);
	mt_h_data_el->Fill(MT);
	jets_data_el->Fill(nJets);
	}
}
//TTBAR
else if(file_count == 2){
totalEvents+=1*lepton_weight*scale_1fb()*5.2;
	ttbar_w_pt->Fill(w_pt_,scale_1fb()*lepton_weight*5.2);
	ttbar_lep_pt->Fill(lep_pt,scale_1fb()*lepton_weight*5.2);
	mt2_h_ttbar->Fill(mt2_event,scale_1fb()*lepton_weight*5.2);
	mt_h_ttbar->Fill(MT,scale_1fb()*lepton_weight*5.2);
	jets_ttbar->Fill(nJets,scale_1fb()*lepton_weight*5.2);
if (abs(lr_id()) == 13){
	ttbar_lep_pt_mu->Fill(lep_pt,scale_1fb()*lepton_weight*5.2);
	mt2_h_ttbar_mu->Fill(mt2_event,scale_1fb()*lepton_weight*5.2);
	mt_h_ttbar_mu->Fill(MT,scale_1fb()*lepton_weight*5.2);
	jets_ttbar_mu->Fill(nJets,scale_1fb()*lepton_weight*5.2);
	}
if (abs(lr_id()) == 11){
	ttbar_lep_pt_el->Fill(lep_pt,scale_1fb()*lepton_weight*5.2);
	mt2_h_ttbar_el->Fill(mt2_event,scale_1fb()*lepton_weight*5.2);
	mt_h_ttbar_el->Fill(MT,scale_1fb()*lepton_weight*5.2);
	jets_ttbar_el->Fill(nJets,scale_1fb()*lepton_weight*5.2);
	}
}
//WJETS
else if(file_count == 3){
totalEvents+=1*lepton_weight*scale_1fb()*5.2;
	wjets_w_pt ->Fill(w_pt_,scale_1fb()*lepton_weight*5.2);	
	wjets_lep_pt ->Fill(lep_pt,scale_1fb()*lepton_weight*5.2);	
	mt2_h_wjets->Fill(mt2_event,scale_1fb()*lepton_weight*5.2);
	mt_h_wjets ->Fill(MT,scale_1fb()*lepton_weight*5.2);
	jets_wjets ->Fill(nJets,scale_1fb()*lepton_weight*5.2);
if (abs(lr_id()) == 13){
	wjets_lep_pt_mu->Fill(lep_pt,scale_1fb()*lepton_weight*5.2);
	mt2_h_wjets_mu->Fill(mt2_event,scale_1fb()*lepton_weight*5.2);
	mt_h_wjets_mu->Fill(MT,scale_1fb()*lepton_weight*5.2);
	jets_wjets_mu->Fill(nJets,scale_1fb()*lepton_weight*5.2);
	}
if (abs(lr_id()) == 11){
	wjets_lep_pt_el->Fill(lep_pt,scale_1fb()*lepton_weight*5.2);
	mt2_h_wjets_el->Fill(mt2_event,scale_1fb()*lepton_weight*5.2);
	mt_h_wjets_el->Fill(MT,scale_1fb()*lepton_weight*5.2);
	jets_wjets_el->Fill(nJets,scale_1fb()*lepton_weight*5.2);
	}
}
	/////FILL END //////
	}//event 
///////////////// write histograms ////////////
char* date = "may30";
if(file_count == 1){
     TFile* fout = new TFile(Form("/home/users/sanil/single/%shists/data_sl5_bv_lep30.root",date),"RECREATE");
	mt2_h_data->Write();
	mt_h_data->Write();
	met_data->Write();
	jets_data->Write();
	data_w_pt->Write();
	data_lep_pt->Write();

	data_lep_pt_mu->Write();
	mt2_h_data_mu->Write();
	mt_h_data_mu->Write();
	met_data_mu->Write();
	data_w_pt_mu->Write();
	jets_data_mu->Write();

	data_lep_pt_el->Write();
	mt2_h_data_el->Write();
	mt_h_data_el->Write();
	met_data_el->Write();
	data_w_pt_el->Write();
	jets_data_el->Write();

	//cout << "DATA mt2 tail:    " << mt2tailCounter << endl;
	cout << endl;
	cout << "DATA totalEvents:    " << totalEvents << endl;
	cout << endl;
  fout->Close();
}
else if(file_count == 2){
    TFile* fout = new TFile(Form("/home/users/sanil/single/%shists/ttbar_sl5_bv_lep30.root",date),"RECREATE");
	mt2_h_ttbar->Write();
	mt_h_ttbar->Write();
	met_ttbar->Write();
	jets_ttbar->Write();
	ttbar_w_pt->Write();
	ttbar_lep_pt->Write();


	ttbar_lep_pt_mu->Write();
	mt2_h_ttbar_mu->Write();
	mt_h_ttbar_mu->Write();
	met_ttbar_mu->Write();
	ttbar_w_pt_mu->Write();
	jets_ttbar_mu->Write();

	ttbar_lep_pt_el->Write();
	mt2_h_ttbar_el->Write();
	mt_h_ttbar_el->Write();
	met_ttbar_el->Write();
	ttbar_w_pt_el->Write();
	jets_ttbar_el->Write();
	//cout << "ttbar mt2 tail: " << mt2tailCounter*scale_1fb()*lepton_weight*5.2 << endl;
	cout << endl;
	cout << "ttbar totalEvents:    " << totalEvents << endl;
	cout << endl;
  fout->Close();
}
else if(file_count == 3){
	//TCanvas *c75 = new TCanvas("c75","c75");
    TFile* fout = new TFile(Form("/home/users/sanil/single/%shists/wjets_sl5_bv_lep30.root",date),"RECREATE");
	mt2_h_wjets->Write();
	mt_h_wjets->Write();
	met_wjets->Write();
	jets_wjets->Write();
	wjets_w_pt->Write();
	wjets_lep_pt->Write();

	wjets_lep_pt_mu->Write();
	mt2_h_wjets_mu->Write();
	mt_h_wjets_mu->Write();
	met_wjets_mu->Write();
	wjets_w_pt_mu->Write();
	jets_wjets_mu->Write();

	wjets_lep_pt_el->Write();
	mt2_h_wjets_el->Write();
	mt_h_wjets_el->Write();
	met_wjets_el->Write();
	wjets_w_pt_el->Write();
	jets_wjets_el->Write();
	//cout << "wjets mt2 tail: " << mt2tailCounter*scale_1fb()*lepton_weight*5.2 << endl;
	cout << endl;
	cout << "wjets totalEvents:    " << totalEvents << endl;
	cout << endl;
  fout->Close();
}
///////////////// write histograms END ////////////
		// +++++++++ //

       	//MT2 count writing//

	///txt output///
if(file_count == 1){
ofstream file_d(Form("/home/users/sanil/single/%shists/lep30data_mt2_bin.txt",date));
	if(!file_d.is_open()){return 0;}
	if( file_d.is_open()){
file_d << "total DATA events: " << totalEvents << endl; 
file_d << "starts with 0 < MT2 <= 10, and goes on increments of 10 geV; " << endl;
file_d << mt2_0_counter   << endl;
file_d << mt2_10_counter  << endl;
file_d << mt2_20_counter  << endl;
file_d << mt2_30_counter  << endl;
file_d << mt2_40_counter  << endl;
file_d << mt2_50_counter  << endl;
file_d << mt2_60_counter  << endl;
file_d << mt2_70_counter  << endl;
file_d << mt2_80_counter  << endl;
file_d << mt2_90_counter  << endl;
file_d << mt2_100_counter << endl;
file_d << mt2_110_counter << endl;
file_d << mt2_120_counter << endl;
file_d << mt2_130_counter << endl;
file_d << mt2_140_counter << endl;
file_d << mt2_150_counter << endl;
file_d << mt2_160_counter << endl;
file_d << mt2_170_counter << endl;
file_d << mt2_180_counter << endl;
file_d << mt2_190_counter << endl;
file_d << mt2_200_counter << endl;
file_d << mt2_G200_counter << endl;
file_d << "^^^ MT2 > 200: " << endl;
	}
}
else if (file_count == 3){
ofstream file_w(Form("/home/users/sanil/single/%shists/lep30wjets_mt2_bin.txt",date));
	if(!file_w.is_open()){return 0;}
	if( file_w.is_open()){
file_w << "total W+Jets events: " << totalEvents << endl; 
file_w << "starts with 0 < MT2 <= 10, and goes on increments of 10 geV; " << endl;
file_w << mt2_0_counter   *scale_1fb()*5.2 << endl;
file_w << mt2_10_counter  *scale_1fb()*5.2 << endl;
file_w << mt2_20_counter  *scale_1fb()*5.2 << endl;
file_w << mt2_30_counter  *scale_1fb()*5.2 << endl;
file_w << mt2_40_counter  *scale_1fb()*5.2 << endl;
file_w << mt2_50_counter  *scale_1fb()*5.2 << endl;
file_w << mt2_60_counter  *scale_1fb()*5.2 << endl;
file_w << mt2_70_counter  *scale_1fb()*5.2 << endl;
file_w << mt2_80_counter  *scale_1fb()*5.2 << endl;
file_w << mt2_90_counter  *scale_1fb()*5.2 << endl;
file_w << mt2_100_counter *scale_1fb()*5.2 << endl;
file_w << mt2_110_counter *scale_1fb()*5.2 << endl;
file_w << mt2_120_counter *scale_1fb()*5.2 << endl;
file_w << mt2_130_counter *scale_1fb()*5.2 << endl;
file_w << mt2_140_counter *scale_1fb()*5.2 << endl;
file_w << mt2_150_counter *scale_1fb()*5.2 << endl;
file_w << mt2_160_counter *scale_1fb()*5.2 << endl;
file_w << mt2_170_counter *scale_1fb()*5.2 << endl;
file_w << mt2_180_counter *scale_1fb()*5.2 << endl;
file_w << mt2_190_counter *scale_1fb()*5.2 << endl;
file_w << mt2_200_counter *scale_1fb()*5.2 << endl;
file_w << mt2_G200_counter*scale_1fb()*5.2 << endl;
file_w << "^^^ MT2 > 200 ^^^" << endl;
	}
}
else if (file_count == 2){
ofstream file_t(Form("/home/users/sanil/single/%shists/lep30ttbar_mt2_bin.txt",date));
	if(!file_t.is_open()){return 0;}
	if( file_t.is_open()){
file_d << "total tt~ events: " << totalEvents << endl; 
file_t << "starts with 0 < MT2 <= 10, and goes on increments of 10 geV; " << endl;
file_t << mt2_0_counter  *scale_1fb()*5.2 << endl;
file_t << mt2_10_counter *scale_1fb()*5.2 << endl;
file_t << mt2_20_counter *scale_1fb()*5.2 << endl;
file_t << mt2_30_counter *scale_1fb()*5.2 << endl;
file_t << mt2_40_counter *scale_1fb()*5.2 << endl;
file_t << mt2_50_counter *scale_1fb()*5.2 << endl;
file_t << mt2_60_counter *scale_1fb()*5.2 << endl;
file_t << mt2_70_counter *scale_1fb()*5.2 << endl;
file_t << mt2_80_counter *scale_1fb()*5.2 << endl;
file_t << mt2_90_counter *scale_1fb()*5.2 << endl;
file_t << mt2_100_counter*scale_1fb()*5.2 << endl;
file_t << mt2_110_counter*scale_1fb()*5.2 << endl;
file_t << mt2_120_counter*scale_1fb()*5.2 << endl;
file_t << mt2_130_counter*scale_1fb()*5.2 << endl;
file_t << mt2_140_counter*scale_1fb()*5.2 << endl;
file_t << mt2_150_counter*scale_1fb()*5.2 << endl;
file_t << mt2_160_counter*scale_1fb()*5.2 << endl;
file_t << mt2_170_counter*scale_1fb()*5.2 << endl;
file_t << mt2_180_counter*scale_1fb()*5.2 << endl;
file_t << mt2_190_counter*scale_1fb()*5.2 << endl;
file_t << mt2_200_counter*scale_1fb()*5.2 << endl;
file_t << mt2_G200_counter*scale_1fb()*5.2 << endl;
file_t << "^^^ MT2 > 200 ^^^" << endl;
	}
	///txt output end///
}
       	//MT2 count writing END//

    // Clean Up
    delete tree;
    file->Close();
    delete file;

	file_count++;
  }	//file_loop
fake_wpt_h->Draw();

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
