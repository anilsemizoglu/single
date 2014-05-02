{
  gSystem->AddIncludePath(Form("-I%s/CORE", gSystem->Getenv("HOME")));
  gSystem->Load(Form("%s/CORE/libCMS2NtupleMacrosCORE.so", gSystem->Getenv("HOME")));  
gROOT->ProcessLine(".L /home/users/sanil/myheaders/SL.cc+");
  gROOT->ProcessLine(".L ScanChain.C+");
 
  TChain *ch = new TChain("tree"); 
  //ch->Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_HT-200To250_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/V05-03-28/merged_ntuple_1.root");

ch->Add("/home/users/sicheng/makebaby/babies/data_sl_3.root");
ch->Add("/home/users/sicheng/makebaby/babies/ttzjets_sl.root");
//ch->Add("/hadoop/cms/store/user/sicheng/wjets_baby.root");  
//ch->Add("/home/users/sicheng/makebaby/babies/wjets_sl.root");
ch->Add("/home/users/sicheng/makebaby/babies/wjets_sl_3.root");
ch->Add("/home/users/sicheng/makebaby/babies/dy_sl.root");
ch->Add("/home/users/sicheng/makebaby/babies/ttbar_sl_3.root");
ch->Add("/home/users/sicheng/makebaby/babies/ttwjets_sl.root");
ch->Add("/home/users/sicheng/makebaby/babies/wwjets_sl.root");
	ScanChain(ch);
}
