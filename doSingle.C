{
  gSystem->AddIncludePath(Form("-I%s/CORE", gSystem->Getenv("HOME")));
  gSystem->Load(Form("%s/CORE/libCMS2NtupleMacrosCORE.so", gSystem->Getenv("HOME")));  
  gROOT->ProcessLine(".L singleChain.C+");
 
  TChain *ch = new TChain("Events"); 
  ch->Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_HT-200To250_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/V05-03-28/merged_ntuple_1.root");
  
	ScanChain(ch);
}
