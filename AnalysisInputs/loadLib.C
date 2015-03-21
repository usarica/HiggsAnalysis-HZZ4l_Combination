{
	gSystem->AddIncludePath("-I$ROOFITSYS/include/");
	gSystem->Load("libHiggsAnalysisCombinedLimit.so");
	gROOT->ProcessLine(".L ../CreateDatacards/include/tdrstyle.cc");
}
