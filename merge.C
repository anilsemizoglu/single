
#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <iomanip>

#include "TStyle.h"
#include "TFile.h"
#include "dataMCplotMaker.h"
#include "/home/users/sanil/myheaders/myheader.h"

int merge()
{
  TFile* data =       new TFile("/home/users/sanil/single/data_hists_met60.root");
  TFile* wjets =     new TFile("/home/users/sanil/single/wjets_hists_met60.root");
  TFile* dy =           new TFile("/home/users/sanil/single/dy_hists_met60.root");
  TFile* ttbar =     new TFile("/home/users/sanil/single/ttbar_hists_met60.root");
  TFile* wwjets =   new TFile("/home/users/sanil/single/wwjets_hists_met60.root");
  TFile* ttwjets = new TFile("/home/users/sanil/single/ttwjets_hists_met60.root");
  TFile* ttzjets = new TFile("/home/users/sanil/single/ttzjets_hists_met60.root");

//get files
  TH1F* mt2_h_data   = (TH1F*) data->Get("mt2_data_"); 
  TH1F* mt_h_data    = (TH1F*) data->Get("mt_data_"); 
  TH1F* met_data    = (TH1F*)  data->Get("met_data_"); 

  TH1F* mt2_h_wjets   = (TH1F*) wjets->Get("mt2_wjets_"); 
  TH1F* mt_h_wjets    = (TH1F*)wjets-> Get("mt_wjets_"); 
  TH1F* met_wjets    = (TH1F*)wjets->Get("met_wjets_"); 

  TH1F* mt2_h_wwjets   = (TH1F*) wwjets->Get("mt2_wwjets_"); 
  TH1F* mt_h_wwjets    = (TH1F*)wwjets-> Get("mt_wwjets_"); 
  TH1F* met_wwjets    = (TH1F*)wwjets->Get("met_wwjets_"); 

  TH1F* mt2_h_ttbar   = (TH1F*) ttbar->Get("mt2_ttbar_"); 
  TH1F* mt_h_ttbar    = (TH1F*)ttbar-> Get("mt_ttbar_"); 
  TH1F* met_ttbar    = (TH1F*)ttbar->Get("met_ttbar_"); 

  TH1F* mt2_h_dy   = (TH1F*) dy->Get("mt2_dy_"); 
  TH1F* mt_h_dy    = (TH1F*)dy-> Get("mt_dy_"); 
  TH1F* met_dy    = (TH1F*)dy->Get("met_dy_"); 

  TH1F* mt2_h_ttwjets   = (TH1F*) ttwjets->Get("mt2_ttwjets_"); 
  TH1F* mt_h_ttwjets    = (TH1F*)ttwjets-> Get("mt_ttwjets_"); 
  TH1F* met_ttwjets    = (TH1F*)ttwjets->Get("met_ttwjets_"); 

  TH1F* mt2_h_ttzjets   = (TH1F*) ttzjets->Get("mt2_ttzjets_"); 
  TH1F* mt_h_ttzjets   = (TH1F*) ttzjets->Get("mt_ttzjets_"); 
  TH1F* met_ttzjets   = (TH1F*) ttzjets->Get("met_ttzjets_"); 
/*
Normalize(mt2_h_data);
//6 stacked
Normalize6(mt2_h_wjets);
Normalize6(mt2_h_dy);
Normalize6(mt2_h_ttbar);
Normalize6(mt2_h_wwjets);
Normalize6(mt2_h_ttzjets);

Normalize(mt_h_data);
//6 stacked
Normalize6(mt_h_wjets);
Normalize6(mt_h_ttbar);
Normalize6(mt_h_wwjets);
Normalize6(mt_h_ttzjets);
Normalize6(mt_h_dy);
met_data->Scale(0.5);
*/
//mt2_h_data->Scale(0.5);
//mt_h_data->Scale(0.5);

vector<TH1F*>	 mt2_background;
vector<TH1F*>	 mt_background;
vector<TH1F*>	 met_background;
vector<char*>	 titles;
vector<Color_t > colors;

mt2_background.push_back(mt2_h_wjets);
mt2_background.push_back(mt2_h_ttbar);
mt2_background.push_back(mt2_h_wwjets);
//mt2_background.push_back(mt2_h_ttwjets);
//mt2_background.push_back(mt2_h_dy);
//mt2_background.push_back(mt2_h_ttzjets);

mt_background.push_back(mt_h_wjets);
mt_background.push_back(mt_h_ttbar);
mt_background.push_back(mt_h_wwjets);
//mt_background.push_back(mt_h_ttwjets);
//mt_background.push_back(mt_h_dy);
//mt_background.push_back(mt_h_ttzjets);

met_background.push_back(met_wjets);
met_background.push_back(met_ttbar);
met_background.push_back(met_wwjets);
//met_background.push_back(met_ttwjets);
//met_background.push_back(met_dy);
//met_background.push_back(met_ttzjets);

titles.push_back("wjets");
titles.push_back("ttbar");
titles.push_back("wwjets");
//titles.push_back("ttwjets");
//titles.push_back("dy");
//titles.push_back("ttzjets");

colors.push_back(kRed-3);
colors.push_back(kCyan-3);
colors.push_back(kViolet-3);
colors.push_back(kBlue+1);
colors.push_back(kYellow-3);
colors.push_back(kGreen-2);

//dataMCplotMaker(mt2_h_data,mt2_background, titles,"MT2, met > 75","","--lumi 5.2, -- outputName MT2_met75");
//dataMCplotMaker(mt_h_data, mt_background, titles,"MT, met > 60","","--lumi 5.2, -- outputName MT_met60");
dataMCplotMaker(met_data, met_background, titles,"Met","","--lumi 5.2 -- outputName MET");

  return 0;
}
