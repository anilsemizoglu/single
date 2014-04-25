#include "TStyle.h"

int mergeDataMC_1()
{
  TFile* wjets = new TFile("/home/users/sanil/single/single_hists.root");
  TFile* data = new TFile("/home/users/sanil/single/data_hists.root");

//get files
  TH1F* data_hee   = (TH1F*) data->Get("hee"); 
//scale
  mctt_hem  ->Scale(lumi);      mcdy_hem  ->Scale(lumi);
//color
  mctt_hem  ->SetFillColor(kRed+1);       mcdy_hem  ->SetFillColor(kOrange+1);
//stack
  THStack* mcs_hem   = new THStack("mchem","Stacked mc for tt~ and DY");
//add
  mcs_hem  ->Add(mctt_hem  );	  mcs_hem  ->Add(mcdy_hem  );

  TCanvas* c1 = new TCanvas;
//draw
  data_hem  ->Draw("P+E");	mcs_hem  ->Draw("same");      data_hem  ->Draw("sameP+E");
  c1->SaveAs("./hists/combinedDataMC_hem_1.png");	                                    

  return 0;
}
