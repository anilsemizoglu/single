{
gSystem->Load("~/CORE/libCMS2NtupleMacrosCORE.so");
  gROOT->ProcessLine(".L dataMCplotMaker.cc+");
gROOT->ProcessLine(".L merge.C+");
gROOT->ProcessLine("merge()");
}
