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
#include "TRandom1.h"
#include "TPad.h"
#include "TH2F.h"
#include "TH1.h"
#include "THStack.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"

// CMS2
#include "CMS2.h"
#include "muonSelections.h"
#include "susySelections.h"
#include "ssSelections.h"
#include "eventSelections.h"
#include "MT2/MT2.h"
#include "MT2/MT2Utility.h"


// Good run list

#include "/home/users/jgran/CMSSW_5_3_2_patch4_V05-03-23/src/CMS2/NtupleMacros/Tools/goodrun.cc"

// My includes
#include "myheader.h"
#include "NSBABY.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVec;

using namespace std;
//using namespace tas;
using namespace ROOT::Math;
using namespace nsbaby;


int ScanChain( TChain* chain, char* suffix = "", bool ismc = true, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  const int FIRTBIN = 0;
  const int LASTBIN = 150;
  const int BINNUM = 60;
  const int INVMLASTBIN = 300;
  const int METLASTBIN = 300;
  const int MTLASTBIN = 200;

  TFile* f1 = new TFile("./hists/mt2_hists_odata_nmcrb_3.root");
  // TH1F* RWPT   = (TH1F*) f1->Get("wpt"); 
  f1->Close();

////////////////////////////////////////////////////////////////////////////////////////////
  TH1F* MT2_hist  = new TH1F("MT2","MT2 distribution", BINNUM, FIRTBIN, LASTBIN+1);
  TH1F* MT2_el    = new TH1F("MT2e","MT2 distribution for e  events", BINNUM, FIRTBIN, LASTBIN+1);
  TH1F* MT2_mu    = new TH1F("MT2m","MT2 distribution for mu events", BINNUM, FIRTBIN, LASTBIN+1);
  TH1F* h2_nozero = new TH1F("hnz","MT2 without 0 bin", BINNUM, 1, LASTBIN+1);
  TH1F* MT2_Cut50 = new TH1F("h50","MT2 without 0 bin with met cut 50", BINNUM, 1, LASTBIN+1);
        
  TH1F* MT_hist  = new TH1F("MT",  "MT distribution", BINNUM, FIRTBIN, MTLASTBIN+1);
  TH1F* MTc_hist = new TH1F("MTc", "MT distribution", BINNUM, FIRTBIN, MTLASTBIN+1);
  TH1F* MT_el    = new TH1F("MTe", "MT distribution for e  events", BINNUM, FIRTBIN, MTLASTBIN+1);
  TH1F* MT_mu    = new TH1F("MTm", "MT distribution for mu events", BINNUM, FIRTBIN, MTLASTBIN+1);
        
  TH1F* JetMult_a = new TH1F("jeta", "Jet mutiplicity all", 7, 0, 7);
  TH1F* Jeta_el   = new TH1F("jae",  "Jet mutiplicity for e  events", 7, 0, 7);
  TH1F* Jeta_mu   = new TH1F("jam",  "Jet mutiplicity for mu events", 7, 0, 7);
  TH1F* JetMult_b = new TH1F("jetb", "b-Jet mutiplicity", 7, 0, 7);
  TH1F* Jetb_el   = new TH1F("jbe",  "b-Jet mutiplicity for e  events", 7, 0, 7);
  TH1F* Jetb_mu   = new TH1F("jbm",  "b-Jet mutiplicity for mu events", 7, 0, 7);
        
  TH1F* Met_all = new TH1F("cmet",  "MET for all events", 70, 0, METLASTBIN+1);
  TH1F* Met_el  = new TH1F("mete",  "MET for e  events",  70, 0, METLASTBIN+1);
  TH1F* Met_mu  = new TH1F("metm",  "MET for mu events",  70, 0, METLASTBIN+1);
  TH1F* FMet    = new TH1F("fmet",  "Real Met + Fake \nu",70, 0, METLASTBIN+1);
        
  TH1F* MetPhi    = new TH1F("metphi",  "MetPhi ",		   90, -3.3, 3.3);
  TH1F* MetPhiCor = new TH1F("cmetphi", "Corrected MetPhi ",       90, -3.3, 3.3);
  TH1F* FMetPhi  = new TH1F("fmetphi", "Metphi + fake nu",         90, -3.3, 3.3);
  TH1F* Metx    = new TH1F("metx", "Met x component",           120, -240, 240);
  TH1F* Mety    = new TH1F("mety", "Met y component",           120, -240, 240);
  TH1F* MetCorx = new TH1F("cmetx","Corrected Met x component", 120, -240, 240);
  TH1F* MetCory = new TH1F("cmety","Corrected Met y component", 120, -240, 240);
        
  TH1F* InvM_lep  = new TH1F("invm","Dileptons InvMass all",   70, 0, INVMLASTBIN+1);
  TH1F* W_PT  = new TH1F("wpt","Pt value of lep + met ",   70, 0, 351);
  TH1F* FAKEW_MT = new TH1F("fakeWmt","Pt value of lep + met ",   70, 0, 351);
  TH1F* TEST  = new TH1F("test","Testing histo for an interesting variable",   70, 0, 351);
  TH1F* TEST2  = new TH1F("test2","Testing2 histo for an interesting variable",   70, 0, 251);

  MT_hist   ->Sumw2(); 
  MTc_hist  ->Sumw2(); 
  MT2_hist  ->Sumw2(); 
  h2_nozero ->Sumw2();
  MT2_Cut50 ->Sumw2();

  InvM_lep  ->Sumw2();
  JetMult_a ->Sumw2();
  JetMult_b ->Sumw2();
  Met_all   ->Sumw2();
  FMet      ->Sumw2();


///////////////////////////////////////////////////////////////////////////////////////////

  int file_count = 0;

  int less_jets = 0;
  int nGoodEvents = 0;

  float bTagDiscriminant = 0.244; 

  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  // Set Good Run List
  // set_goodrun_file("/home/users/jgran/analysis/sswh/fakes/json/final_19p49fb_cms2.txt");

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    
    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("tree");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    // cms2.Init(tree);
    baby.Init(tree);  
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      // cms2.GetEntry(event);
      baby.GetEntry(event);  
      ++nEventsTotal;
    
      // Progress	
      CMS2::progress( nEventsTotal, nEventsChain );

      // Select Good Runs

      Int_t n_jets = 0;
      Int_t n_bTag = 0;
	
      // Analysis Code for MT2
      if( lr_p4().pt() < 30 ) continue;


      //// tt~ event selection ////
      // if(jets_p4().size() < 2) {less_jets++; continue;} // pre-selection check

      // for(unsigned int c = 0; c < jets_p4().size(); c++) {
      // 	float _jetPt = jets_p4().at(c).pt() * jets_p4Correction().at(c);
      // 	// jet pt times correction to jets
      // 	float dr_lr = VectorUtil::DeltaR(jets_p4().at(c), lt_p4() );
      // 	// delta R is the distance in the eta-phi plane
      // 	if(_jetPt < 30) continue;       // disgard those with small pt

      // 	if(dr_lr < 0.4) continue;      // disgard small delta R 

      // 	if(fabs(jets_p4().at(c).eta()) > 2.5) continue; // May here cause less tt~ events as dumping higher eta jets

      // 	float _bTag = btagDiscriminant().at(c);
      // 	if ( _bTag > bTagDiscriminant)  n_bTag++;	 // mark as the has qualified b-jet, _bTag should be 0~1

      // 	// float _jetdp = abs(pfjets_p4().at(c).pt() - _jetPt);  
      // 	// should be diff in correction, tested equal to pt*abs(1-correc)
      // 	n_jets++;					 // as qualified jet multiplicity

      // }//pfjets_p4().size()

      // if(n_jets < 2 || n_bTag == 0) {less_jets++; continue;}


      // if( met() < 80 ) continue;
      //      if( njets() < 4 ) continue;

      // float ht = 0;
      // for(unsigned int c = 0; c < jets_p4().size(); c++) {
      // 	float jetPt = jets_p4().at(c).pt() * jets_p4Correction().at(c);
      // 	ht += jetPt;
      // }      
      // if( ht < 200 ) continue;

      //// Jet requirements

    //// To be done: ////
      //// Later ///    GenRandomW:   //////			
      TRandom1 rand; 

      float M_W = 80.2;

      ///////////////// fake w initialization  ///////////////////
      float fakeW_eta = -4;	
      float fakeW_theta = -4;
      float fakeW_phi = -4;
      float fakeW_p  = -1;
      Float_t fakeW_px = 0, fakeW_py = 0 , fakeW_pz = 0;

      ///////////////// fake lepton initialization  ///////////////////
      Float_t fakeLep_px = 0, fakeLep_py = 0 , fakeLep_pz = 0;
      float fakeLep_p  = M_W/2;		// Later can apply a width here
      float fakeLep_eta = -4;
      float fakeLep_theta = -4;
      float fakeLep_phi = -4;

      LorentzVec fakeLep_p4;
      
      float ww = 0; 

      // do{
      // 	ww = rand.Gaus(60,20);
      // }while(ww < 10);

      fakeW_p = 30;
      // do{
      	float w1 = rand.Rndm();
      	float w2 = rand.Rndm();
      	fakeW_theta = acos(2*w2-1);
      	fakeW_eta = - log(0.5*fakeW_theta);
      	fakeW_phi = 2*3.14156265359*(w1-0.5);
      // } while(fabs(fakeW_eta) > 2.4); 	

      fakeW_px = fakeW_p* sin(fakeW_theta)* cos(fakeW_phi);
      fakeW_py = fakeW_p* sin(fakeW_theta)* sin(fakeW_phi);
      fakeW_pz = fakeW_p* cos(fakeW_theta);

      // fakeW_px = 10;
      // fakeW_py = 30;
      // fakeW_pz = 10;

      TLorentzVector fakeLepP4, fakeNeuP4;

      // float rdmWpt = RWPT->GetRandom();
      // TEST2->Fill(rdmWpt);
      
     // Generate rest decay 
      do{
	float r1 = rand.Rndm();
	float r2 = rand.Rndm();
	fakeLep_theta = acos(2*r2-1);
	fakeLep_eta = - log(0.5*fakeLep_theta);
	fakeLep_phi = 2*3.14156265359*(r1-0.5);
	//} while(fabs(fakeLep_eta) > 2.4); 	// Want those would pass selection

	// fakeW_p = sqrt(1100);
	/*
	  fakeLep_p4.SetE(fakeLep_p);
	  fakeLep_p4.SetPhi(fakeLep_phi);
	  fakeLep_p4.SetEta(fakeLep_eta);
	*/
	fakeLep_px = fakeLep_p* sin(fakeLep_theta)* cos(fakeLep_phi);
	fakeLep_py = fakeLep_p* sin(fakeLep_theta)* sin(fakeLep_phi);
	fakeLep_pz = fakeLep_p* cos(fakeLep_theta);

	// float v2 = 1/((M_W/(fakeW_p))*(M_W/(fakeW_p))+1);
	// float gamma = 1/(sqrt(1-v2));
	float gammaM = sqrt(fakeW_p*fakeW_p + M_W*M_W); // gammaM == gamma times Mass of W

	// Initialize lep and nu at W rest frame
	fakeLepP4.SetPxPyPzE( fakeLep_px,  fakeLep_py,  fakeLep_pz, fakeLep_p);
	fakeNeuP4.SetPxPyPzE(-fakeLep_px, -fakeLep_py, -fakeLep_pz, fakeLep_p);

	// Boost from W rest frame to Lab frame as W_p / gamma*Mass == v 
	fakeLepP4.Boost(fakeW_px / gammaM, fakeW_py / (gammaM), fakeW_pz / (gammaM) );
	fakeNeuP4.Boost(fakeW_px / (gammaM), fakeW_py / (gammaM), fakeW_pz / (gammaM) );

      } while(fakeLepP4.Eta() < 2.4);


      TEST->Fill( (fakeLepP4 + fakeNeuP4).Pt(), scale_1fb());

      fakeLep_p4.SetPxPyPzE(fakeLep_px, fakeLep_py, fakeLep_pz, fakeLep_p);

      float nux = - fakeLep_px;   
      float nuy = - fakeLep_py;   


      // Make MET correction and add nu_pt to met 

      // for metphi correction
      float metx = met() * cos( metPhi() );
      float mety = met() * sin( metPhi() );
      float shiftx = 0.;
      float shifty = 0.;

      Metx->Fill(metx, scale_1fb());
      Mety->Fill(mety, scale_1fb());
      
      shiftx = (! isRealData()) ? (0.1166 + 0.0200*nvtxs()) : (0.2661 + 0.3217*nvtxs());
      shifty = (! isRealData()) ? (0.2764 - 0.1280*nvtxs()) : (-0.2251 - 0.1747*nvtxs());

      metx -= shiftx;
      mety -= shifty;

      MetCorx->Fill(metx, scale_1fb());
      MetCory->Fill(mety, scale_1fb());

      float cmet = sqrt( metx*metx + mety*mety ); // cmet == corrected met
      float cmetphi = atan2( mety , metx );

      // end of metphi correction

      MetPhi->Fill(metPhi(), scale_1fb());
      MetPhiCor->Fill(cmetphi, scale_1fb());

      float fmetx = metx + nux;	                 // fmet == met after adding fake neutrinos
      float fmety = mety + nuy;

      float fmet = sqrt( fmetx*fmetx + fmety*fmety );
      float fmetphi = atan2( fmety , fmetx );

      if(fmet > (float)METLASTBIN)  FMet->Fill((float)METLASTBIN, scale_1fb());
      else       FMet->Fill(fmet, scale_1fb());
      if(cmet > (float)METLASTBIN) {
	Met_all->Fill((float)METLASTBIN, scale_1fb());
	if(abs(lr_id()) == 11 ) Met_el->Fill((float)METLASTBIN, scale_1fb());
	if(abs(lr_id()) == 13 ) Met_mu->Fill((float)METLASTBIN, scale_1fb());
      }
      else {
	Met_all->Fill(cmet, scale_1fb());
	if(abs(lr_id()) == 11 ) Met_el->Fill(cmet, scale_1fb());
	if(abs(lr_id()) == 13 ) Met_mu->Fill(cmet, scale_1fb());
      }
      
      FMetPhi->Fill(fmetphi, scale_1fb());


      // Plot PT of real W
      float wr_px = metx + lr_p4().Px();
      float wr_py = mety + lr_p4().Py();
      float wr_p = sqrt(wr_px* wr_px  + wr_py*wr_py);
      

      Fill1F(W_PT, wr_p, scale_1fb());


      // Fill the jet mutiplicity to the histogram
      if(njets() > 7)  JetMult_a->Fill((float) 7, scale_1fb());
      else {
	JetMult_a->Fill(njets(), scale_1fb());
	if(abs(lr_id()) == 11 ) Jeta_el->Fill(njets(), scale_1fb());
	if(abs(lr_id()) == 13 ) Jeta_mu->Fill(njets(), scale_1fb());
      }

      if(nbTag() > 7)  JetMult_b->Fill((float) 7, scale_1fb());
      else {
	JetMult_b->Fill(nbTag(), scale_1fb());
	if(abs(lr_id()) == 11 ) Jetb_el->Fill(nbTag(), scale_1fb());
	if(abs(lr_id()) == 13 ) Jetb_mu->Fill(nbTag(), scale_1fb());
      }
      
      float mt2 = MT2(fmet, fmetphi, lr_p4() , fakeLep_p4);

      if(mt2 > (float)LASTBIN){
	MT2_hist->Fill((float)LASTBIN, scale_1fb());
	if(abs(lr_id()) == 11 ) MT2_el->Fill((float)LASTBIN, scale_1fb());
	if(abs(lr_id()) == 13 ) MT2_mu->Fill((float)LASTBIN, scale_1fb());
	h2_nozero->Fill((float)LASTBIN, scale_1fb());
      }
      else{
	MT2_hist->Fill(mt2, scale_1fb());
	if(abs(lr_id()) == 11 ) MT2_el->Fill(mt2, scale_1fb());
	if(abs(lr_id()) == 13 ) MT2_mu->Fill(mt2, scale_1fb());
	h2_nozero->Fill(mt2, scale_1fb());
      }

      // if(M_ll > (float)INVMLASTBIN)  InvM_lep->Fill((float)INVMLASTBIN, scale_1fb());
      // else InvM_lep->Fill(M_ll, scale_1fb());

      float mt  = sqrt(2*lr_p4().Pt()* met()*(1 - cos(lr_p4().Phi() - metPhi())));
      float mtc = sqrt(2*lr_p4().Pt()* cmet *(1 - cos(lr_p4().Phi() - cmetphi )));
      //cout << lr_p4().Pt() <<"\t" << met() <<"\t" << 1- cos(lr_p4().Phi() - cmetphi ) <<"\t" << mt << endl;
      
      if(mt  > (float) MTLASTBIN)  MT_hist ->Fill((float) MTLASTBIN, scale_1fb());
      else MT_hist ->Fill(mt, scale_1fb());

      if(mtc > (float) MTLASTBIN)  {
	MTc_hist ->Fill((float) MTLASTBIN, scale_1fb());
	if(abs(lr_id()) == 11 ) MT_el->Fill((float) MTLASTBIN, scale_1fb());
	if(abs(lr_id()) == 13 ) MT_mu->Fill((float) MTLASTBIN, scale_1fb());
      }
      else {
	MTc_hist ->Fill(mtc, scale_1fb());
	if(abs(lr_id()) == 11 ) MT_el->Fill(mtc, scale_1fb());
	if(abs(lr_id()) == 13 ) MT_mu->Fill(mtc, scale_1fb());
      }


      // GenMet->Fill(gen_met());
      // RecMet->Fill(metcor);

    } //loop over events in the current file

    file_count++;

    // Clean Up
    delete tree;
    file->Close();
    delete file;
  } //file loop
  
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }

  cout << "\nNumber of total Events: " << nEventsTotal 
       //<< " at 20GeV cut, " << nttbarEvents-metCut30 
       <<endl <<endl;
  cout << "For the data samples: "
       //<< "\n# notGoodRun: " << notGoodRun 
       //<< "\n# goodRun: " << goodRun
       //<< "\n# noGoodVtx: " << noGoodVtx
       << endl;


  // nttbarScaled *= nEvtFdr/nEventsChain; // no need for running on baby

  // char* suffix = "_1_mcdy";
  TFile* fout = new TFile(Form("./hists/mt2_hists%s.root",suffix),"RECREATE");

  MT2_hist ->Write();
  MT2_el   ->Write();
  MT2_mu   ->Write();
  MT_hist  ->Write();
  MTc_hist ->Write();
  MT_el    ->Write();
  MT_mu    ->Write();
  MT2_Cut50->Write();
  h2_nozero->Write();
  Met_all  ->Write();
  Met_el   ->Write();
  Met_mu   ->Write();
  FMet     ->Write();
  InvM_lep ->Write();
  JetMult_a->Write();
  JetMult_b->Write();
  Jeta_el  ->Write();
  Jeta_mu  ->Write();
  Jetb_el  ->Write();
  Jetb_mu  ->Write();
  W_PT     ->Write();


  fout->Close();

  //TCanvas* c1 = new TCanvas;
  // TLegend* l1 = new TLegend(0.4,0.1,1,0.4,"Legend");
  // l1->AddEntry(MT2_hist,"MT2","f");

  // MT2_ttbar->Draw();
  // // c1->BuildLegend(0.2, 0.1, 0.3, 0.2);
  // c1->SaveAs(Form("./hists/mt2%s.png", suffix));
  // h2_nozero->Draw();
  // c1->SaveAs(Form("./hists/mt2_nonzero%s.png", suffix));

  TCanvas* c2 = new TCanvas;		   
  MetCorx->Draw();
  Metx->Draw("Same");
  MetCorx->SetLineColor(kRed+1);
  //c2->SaveAs(Form("./hists/metcx%s.png", suffix));
  MetCory->Draw();
  Mety->Draw("Same");
  MetCory->SetLineColor(kRed+1);
  //c2->SaveAs(Form("./hists/metcy%s.png", suffix));
  MetPhi->Draw();
  MetPhiCor->Draw("same");
  MetPhiCor->SetLineColor(kRed+1);
  //c2->SaveAs(Form("./hists/metPhi_cor%s.png", suffix));

  FMet->Draw();
  Met_all->Draw("Same");
  Met_all->SetLineColor(kRed+1);
  //c2->SaveAs(Form("./hists/fmet%s.png", suffix));

  cout << "Entries FMet vs Met: " << FMet->GetEntries() << "\t" << Met_all->GetEntries() << endl;

  // FMetPhi->Draw();
  // MetPhiCor->Draw("same");
  //FMetPhi->SetLineColor(kRed+1);
  //c2->SaveAs(Form("./hists/fmetphi_cor%s.png", suffix));

  // MT_hist ->Draw("hist");			   
  // MTc_hist->Draw("samehist");			   
  // MTc_hist->SetLineColor(kRed+1);

  // TEST->Draw();

  TCanvas* c3 = new TCanvas;		   
  MT2_hist->SetLineStyle(1);			   
  MT2_hist->Draw("hist");			   
  // TEST2->Draw();
  TCanvas* c4 = new TCanvas;		   
  //RWPT->Draw();

  //MT2_Cut50->Draw("same");
  //MT2_Cut50->SetLineColor(kRed+1);
  // c3->SaveAs(Form("./hists/mt2_emuonly%s.png", suffix));
					   
  //  c1->SetLogy();			   
  // c1->Close();			   
  					   

  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << "Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:   " << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:   " << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;

  return 0;
 }

