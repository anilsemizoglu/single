///// A lot of includes /////
#include "TLorentzVector.h"	

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVec;

//using namespace tas;
using namespace ROOT::Math;

int ScanChain( TChain* chain, char* suffix = "", bool ismc = true, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {


  /////////// Get Histogram of W_p distribution
  TFile* f1 = new TFile("./hists/mt2_hists_odata_nmcrb_3.root");
  TH1F* Wp_dist   = (TH1F*) f1->Get("wp"); // 
  f1->Close();

  TH1F* MT2_hist  = new TH1F("MT2","MT2 distribution", BINNUM, FIRTBIN, LASTBIN+1);
  MT_hist   ->Sumw2(); 

  // Loop over events to Analyze
  // File loop codes are omitted as unimportant here
      
  /**************************************
   ***** Lepton and Jet Selections ******
   **************************************/

   // Fake l nu generation start here

      TRandom1 rand; 		          // Random number generator

      float M_W = 80.2;			  // Later can give a width here

      // Variable Initialization

      float fakeW_eta = -4;	
      float fakeW_theta = -4;
      float fakeW_phi = -4;
      float fakeW_p  = -1;
      Float_t fakeW_px = 0, fakeW_py = 0 , fakeW_pz = 0;

      Float_t fakeLep_px = 0, fakeLep_py = 0 , fakeLep_pz = 0;
      float fakeLep_p  = M_W/2;		  // lep_p ==  nu_p  == W_mass / 2
      float fakeLep_eta = -4;
      float fakeLep_theta = -4;
      float fakeLep_phi = -4;

      LorentzVec fakeLep_p4;
      
      float ww = 0; 

      ///////////////// Fake W p4 generation ///////////////////
      fakeW_p = Wp_dist->GetRandom();	  // Get the magnitude of p from histogram

      float w1 = rand.Rndm();		  // Get some random number in (0,1) for flat phi
      float w2 = rand.Rndm();		  // Get some random number in (0,1) for flat cosTheta
      fakeW_theta = acos(2*w2-1);	   
      fakeW_phi = 2*3.14156265359*(w1-0.5);       

      // Change of coordinate to XYZ and scale by magnitude of p
      fakeW_px = fakeW_p* sin(fakeW_theta)* cos(fakeW_phi); 
      fakeW_py = fakeW_p* sin(fakeW_theta)* sin(fakeW_phi);
      fakeW_pz = fakeW_p* cos(fakeW_theta);


      TLorentzVector fakeLepP4, fakeNeuP4;	  // Declare 4 vector for lepton and neutrino to be boosted
      
      // Generate rest decay 
      do{
	//////// W rest frame lepton neutrino generation  ////////////

	float r1 = rand.Rndm();			  // Get some number in (0,1) for flat phi
	float r2 = rand.Rndm();			  // Get some number in (0,1) for flat cosTheta
	fakeLep_theta = acos(2*r2-1);		  
	fakeLep_phi = 2*3.14156265359*(r1-0.5);

	// Change of coordinate ot cartesian and scale to 40 GeV
	fakeLep_px = fakeLep_p* sin(fakeLep_theta)* cos(fakeLep_phi);
	fakeLep_py = fakeLep_p* sin(fakeLep_theta)* sin(fakeLep_phi);
	fakeLep_pz = fakeLep_p* cos(fakeLep_theta);

	// Initialize lep and nu at W rest frame
	fakeLepP4.SetPxPyPzE( fakeLep_px,  fakeLep_py,  fakeLep_pz, fakeLep_p);
	fakeNeuP4.SetPxPyPzE(-fakeLep_px, -fakeLep_py, -fakeLep_pz, fakeLep_p);

	float gammaM = sqrt(fakeW_p*fakeW_p + M_W*M_W);   // gammaM == gamma_of_W * Mass_of_W

	// Boost from W rest frame to Lab frame using   v = W_p / gamma*Mass 
	fakeLepP4.Boost(fakeW_px / (gammaM), fakeW_py / (gammaM), fakeW_pz / (gammaM) );
	fakeNeuP4.Boost(fakeW_px / (gammaM), fakeW_py / (gammaM), fakeW_pz / (gammaM) );

      } while(fakeLepP4.Eta() < 2.4);    	  // Keep the same selection as "real" lepton

      // Sanity Check
      if( fabs( (fakeLepP4 + fakeNeuP4).P() - fakeW_p ) > 1 )  cerr << "Boosting fail!! \n";
      if( fabs( (fakeLepP4 + fakeNeuP4).M() - M_W ) > 1 )      cerr << "Boosting fail!! \n";

      // Cast the vector from TLorentzVector to LorentzVec
      fakeLep_p4.SetPxPyPzE(fakeLep_px, fakeLep_py, fakeLep_pz, fakeLep_p);

      // Only need the transverse component to add into the met
      float nux = fakeNeu.Px();   
      float nuy = fakeNeu.Py();   

      // Make MET correction and add nu_pt to met 

      //////// for metphi correction, credit to Verena /////
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

      float cmet = sqrt( metx*metx + mety*mety ); // cmet == corrected met
      float cmetphi = atan2( mety , metx );

      /////// end of metphi correction  ///////

      float fmetx = metx + nux;	                 // fmet == met after adding fake neutrinos
      float fmety = mety + nuy;

      float fmet = sqrt( fmetx*fmetx + fmety*fmety );
      float fmetphi = atan2( fmety , fmetx );

      // Calculate mt2
      float mt2 = MT2(fmet, fmetphi, lr_p4() , fakeLep_p4); // lr stands for "real" lepton 

      // Calculate the original MT
      float mt  = sqrt(2*lr_p4().Pt()* met()*(1 - cos(lr_p4().Phi() - metPhi())));

      //// The rest are only histograms filling 

    } //loop over events in the current file


  TFile* fout = new TFile(Form("./hists/mt2_hists%s.root",suffix),"RECREATE");
  MT2_hist ->Write();
  fout->Close();

  //// The rest are only histograms filling 

  return 0;
 }

