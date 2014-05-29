/**************REMINDER OF HYPOTHESES*************************/
/*
const float gi_phi2_phi4[kNumHypo][13]={ // Coefficients as in US_MEweights branch of ZZAnalusis/AnalysisStep/test/Macros/HZZ4l.h
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0, 0, 0, 0, 0},									 // Pure SM 0
	{	0,	1.0,	0,	0,		0,	0,	1.0,	0,	0, 0, 0, 0, 0},								 // fa2=1 1
	{	0,	0,	0,	1.0,		0,	0,	0,	1.0,	0, 0, 0, 0, 0},								 // fa3=1 2
	{	0,	0,	0,	0,		0,	0,	0,	0,	1.0, 0, 0, 0, 0},									 // fL1=1 3
	{	1.0,	1.638,	0,	0,		0,	0,	0.5,	0,	0, 0, 0, 0, 0},							 // fa2=0.5 4
	{	1.0,	0,	0,	2.521,		0,	0,	0,	0.5,	0, 0, 0, 0, 0},							 // fa3=0.5 5
	{	0,	0.650,	0,	1.0,		0,	0,	0.5,	0.5,	0, 0, 0, 0, 0},						 // fa2=fa3=0.5 6
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	12046.01, 0, 0, 0, 0},							 // fLambda1=-0.5, for T3 templates 7
	{	1.0,	1.638,	0,	2.521,		0,	0,	1.0/3.0,	1.0/3.0,	0, 0, 0, 0, 0},			 // fa2=fa3=1/3 8
	{	1.0,	0.546,	0,	0,		0,	0,	0.1,	0,	0, 0, 0, 0, 0},							 // fa2=0.1 9
	{	1.0,	0,	0,	0.840,		0,	0,	0,	0.1,	0, 0, 0, 0, 0},							 // fa3=0.1 10
	{	1.0,	0.579,	0,	0.891,		0,	0,	0.1,	0.1,	0, 0, 0, 0, 0},					 // fa2=fa3=0.1 11
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	-12046.01, 0, 0, 0, 0},							 // fLambda1=0.5 12
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	-7885.965, 0, 0, 0, 0},	 						 // flambda1=0.3 13  <<<<<<<<<<-------------- gapZZHypo has no dedicated file
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	-4015.337, 0, 0, 0, 0},							 // flambda1=0.1 14
	{	1.0,	1.638,	0,	0,		PI_VAL/2.0,	0,	0.5,	0,	0, 0, 0, 0, 0},					 // fa2=0.5, phia2=90 15
	{	1.0,	0,	0,	2.521,		0,	PI_VAL/2.0,	0,	0.5,	0, 0, 0, 0, 0},					 // fa3=0.5, phia3=90 16
	{	0,	0.650,	0,	1.0,		0,	PI_VAL/2.0,	0.5,	0.5,	0, 0, 0, 0, 0},				 // fa2=fa3=0.5, phia3=90 17
	{	1.0,	1.638,	0,	2.521,		0,	PI_VAL/2.0,	1.0/3.0,	1.0/3.0,	0, 0, 0, 0, 0},	 // fa2=fa3=1/3, phia3=90 18
	{	1.0,	0.546,	0,	0,		PI_VAL/2.0,	0,	0.1,	0,	0, 0, 0, 0, 0},					 // fa2=0.1, phia2=90 19
	{	1.0,	0,	0,	0.840,		0,	PI_VAL/2.0,	0,	0.1,	0, 0, 0, 0, 0},					 // fa3=0.1, phia3=90 20
	{	1.0,	0.579,	0,	0.891,		0,	PI_VAL/2.0,	0.1,	0.1,	0, 0, 0, 0, 0},			 // fa2=0.1, fa3=0.1, phia3=90 21
	{	1.0,	1.638,	0,	0,		PI_VAL,	0,	0.5,	0,	0, 0, 0, 0, 0},						 // fa2=-0.5 22
	{	1.0,	0,	0,	2.521,		0,	PI_VAL,	0,	0.5,	0, 0, 0, 0, 0},						 // fa3=-0.5 23
	{	0,	0.650,	0,	1.0,		PI_VAL,	0,	0.5,	0.5,	0, 0, 0, 0, 0},					 // fa2=-0.5, fa3=0.5 24    <<<<<<<---------- kNumFiles+1 is the last dedicated file
	{	0,	1.638,	0,	0,		0,	0,	0,	0,	12046.01, 0, 0, 0, 0},							 // fa2=0.5, fLambda1=-0.5, for templates 25
	{	0,	1.638,	0,	0,		PI_VAL*3.0/2.0,	0,	0,	0,	12046.01, 0, 0, 0, 0},				 // fa2=-0.5i, fLambda1=-0.5, for templates 26
	{	0,	1.638,	0,	0,		0,	0,	0,	0,	-12046.01, 0, 0, 0, 0},							 // fa2=0.5, fLambda1=0.5 27
	{	1.0,	1.638,	0,	0,		PI_VAL,	0,	0,	0,	-12046.01, 0, 0, 0, 0},					 // fa2=-0.33, fLambda1=0.33 28
	{	0,	0,	0,	2.521,		0,	0,	0,	0,	12046.01, 0, 0, 0, 0},							 // fa3=0.5, fLambda1=-0.5, for templates 29
	{	0,	0,	0,	2.521,		PI_VAL*3.0/2.0,	0,	0,	0,	12046.01, 0, 0, 0, 0},				 // fa3=-0.5i, fLambda1=-0.5, for templates 30
	{	0,	0,	0,	2.521,		0,	0,	0,	0,	-12046.01, 0, 0, 0, 0},							 // fa3=0.5, fLambda1=0.5 31
	{	1.0,	0,	0,	2.521,		PI_VAL,	0,	0,	0,	-12046.01, 0, 0, 0, 0},					 // fa3=-0.33, fLambda1=0.33 32
	{	0,	0,	0,	0,		0,	0,	0,	0,	0, 1, 0, 0, 0},										 // Pure PM-ZG 33    <<<<<<<<<<<<<<---------- firstNonZZHypo
	{	0,	0,	0,	0,		0,	0,	0,	0,	0, 0, 1, 0, 0},										 // Pure PM-GG 34
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0, 0.0473, 0, 0, 0},							 // fPM-ZG=0.5 35
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0, 0, -0.0531, 0, 0},							 // fPM-GG=0.5 36
	{	0,	0,	0,	0,		0,	0,	0,	0,	0, 0, 0, 1, 0},										 // Pure M-ZG 37
	{	0,	0,	0,	0,		0,	0,	0,	0,	0, 0, 0, 0, 1},										 // Pure M-GG 38
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0, 0, 0, 0.052161, 0},							 // fM-ZG=0.5 39
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0, 0, 0, 0, -0.053666},							 // fM-GG=0.5 40
	{	0,	0,	0,	0,		0,	0,	0,	0,	0, 0.0473, 0, 0.052161, 0},							 // fZG=fM-ZG=0.5 41
	{	0,	0,	0,	0,		0,	0,	0,	0,	0, 0, -0.0531, 0, -0.053666},						 // fGG=fM-GG=0.5 42
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0, 0.00175, -0.002, 0, 0}						 // SM ZZ, ZG. GG 43
};

********** NOTES **********
*** Be careful about this ordering in your datacards below and name your channels appropriately when combining datacards:

	TString chan[6]={"ch1","ch2","ch3","ch4","ch5","ch6"};
	int order[6]={3,4,5,0,1,2};
	RooDataSet *data[6];
	RooDataSet *toyf03[6];
	TString channameO[6]={"7TeV/4e","7TeV/4mu","7TeV/2mu2e","8TeV/4mu","8TeV/2mu2e","8TeV/4e"};
	float normO[6]=		{	0.7194,	1.2859,	1.7264,	6.0802,	7.9085,	3.1441	}; 
	float normzxO[6]=	{	0.6226,	0.2230,	1.0628,	1.1878,	4.2929,	2.7676	};
	float normqqO[6]=	{	0.8386,	1.7971,	2.2456,	7.6478,	8.8585,	2.9364	};
	float normggO[6]=	{	0.0341,	0.0625,	0.0741,	0.4131,	0.5005,	0.2041	};

*** Depending on the datacard maker implementation, you might need to change these strings:

	TString toyKDName[] = {
		"CMS_zz4l_KD1",
		"CMS_zz4l_KD2",
		"CMS_zz4l_smd"
	};
	TString toyKDNameF[] = {
		"CMS_zz4l_KD1/F",
		"CMS_zz4l_KD2/F",
		"CMS_zz4l_smd/F"
	};

*** Supported KDs are

	TString KDlist[] = {
		"D_g1Q2_phi0", // 0
		"D_g1_vs_g2_phi0", // 1
		"D_g1_vs_g4_phi0", // 2
		"D_ZG", // 3
		"D_GG", // 4
		"D_ZG_PS", // 5
		"D_GG_PS", // 6
		"D_g1Q2int_phi0", // 7
		"D_g2int_phi0", // 8
		"D_g4int_phi0", // 9
		"D_g2int_phi90", // 10
		"D_g4int_phi90", // 11
		"D_ZGint", // 12
		"D_GGint", // 13
		"D_ZG_PSint", // 14
		"D_GG_PSint" // 15
	};
	int firstIntKD=7; <<<<<------- Increment this number if you add a pure hypothesis KD into the array above

*/
using namespace RooFit;

void writeShuffledAsimovToys(int BSMSample=0, int iKD1=0, int iKD2=1, int iKD3=-1, double signalStrength=1, int writeToys=0){
	const int kNumHypo=44;
	const int kNumFiles=24;
	const int gapZZHypo=13;
	const int firstNonZZHypo=33;

/***** RUN CONDITIONS *****/
	if(signalStrength<0) assert(0);
	if(writeToys<0) assert(0);
	if(iKD1<0) assert(0);
	if(iKD3>=0 && iKD2<0) assert(0);
	if(BSMSample>=kNumHypo) assert(0);
/*** END RUN CONDITIONS ***/
	int maxNtoysPossible=1000000;

	TString KDlist[] = {
		"D_g1Q2_phi0", // 0
		"D_g1_vs_g2_phi0", // 1
		"D_g1_vs_g4_phi0", // 2
		"D_ZG", // 3
		"D_GG", // 4
		"D_ZG_PS", // 5
		"D_GG_PS", // 6
		"D_g1Q2int_phi0", // 7
		"D_g2int_phi0", // 8
		"D_g4int_phi0", // 9
		"D_g2int_phi90", // 10
		"D_g4int_phi90", // 11
		"D_ZGint", // 12
		"D_GGint", // 13
		"D_ZG_PSint", // 14
		"D_GG_PSint" // 15
	};
	int firstIntKD=7;

	TString toyKDName[] = {
		"CMS_zz4l_KD1",
		"CMS_zz4l_KD2",
		"CMS_zz4l_smd",
		"CMS_zz4l_KD3"
	};
	TString toyKDNameF[] = {
		"CMS_zz4l_KD1/F",
		"CMS_zz4l_KD2/F",
		"CMS_zz4l_smd/F",
		"CMS_zz4l_KD3/F"
	};
	char* cMCwgt[2] = {
		"MC_CV_weight",
		"MC_CV_4GeVcut_weight"
	};

	TString chan[6]={"ch1","ch2","ch3","ch4","ch5","ch6"};
//	int order[6]={1,2,0,5,3,4};
	int order[6]={2,1,0,4,3,5};
	RooDataSet *data[6];
	RooDataSet *toyf03[6];
	TString channameO[6]={"7TeV/4e","7TeV/4mu","7TeV/2mu2e","8TeV/4mu","8TeV/2mu2e","8TeV/4e"};
	double normO[6]=	{	0.7194,	1.2859,	1.7264,	6.0802,	7.9085,	3.1441	}; 
	double normzxO[6]=	{	0.6226,	0.2230,	1.0628,	1.1878,	4.2929,	2.7676	};
	double normqqO[6]=	{	0.8386,	1.7971,	2.2456,	7.6478,	8.8585,	2.9364	};
	double normggO[6]=	{	0.0341,	0.0625,	0.0741,	0.4131,	0.5005,	0.2041	};

	double norm[6] = { 0 };
	double normgg[6] = { 0 };
	double normzx[6] = { 0 };
	double normqq[6] = { 0 };
	TString channame[6];
	TTree* tsig[6];
	TTree* tggzz[6];
	TTree* tqqzz[6];
	TTree* tzx[6];

	RooCategory cat("CMS_channel","CMS_channel");
	for(int j=0;j<6;j++){
		cat.defineType(chan[j],j);
		cat.setLabel(chan[j]);
		int idx=order[j];
		channame[j]=channameO[idx];
		norm[j]=normO[idx];
		normgg[j]=normggO[idx];
		normqq[j]=normqqO[idx];
		normzx[j]=normzxO[idx];

		TString ctoys_sig="ToyEvents_Sig_";
		ctoys_sig.Append(chan[j]);
		TString ctoys_ggzz="ToyEvents_ggZZ_";
		ctoys_ggzz.Append(chan[j]);
		TString ctoys_qqzz="ToyEvents_qqZZ_";
		ctoys_qqzz.Append(chan[j]);
		TString ctoys_zx="ToyEvents_ZX_";
		ctoys_zx.Append(chan[j]);

		tsig[j] = new TTree(ctoys_sig,"");
		tggzz[j] = new TTree(ctoys_ggzz,"");
		tqqzz[j] = new TTree(ctoys_qqzz,"");
		tzx[j] = new TTree(ctoys_zx,"");
	}
	 
	RooRealVar* rweightFit = new RooRealVar("_weight_","_weight_",0.,1.e20);
	RooRealVar* rCMS_zz4l_KD1 = new RooRealVar(toyKDName[0],toyKDName[0],0.,1.);
	RooRealVar* rCMS_zz4l_smd = new RooRealVar(toyKDName[2],toyKDName[2],0.,1.);
	RooRealVar* rCMS_zz4l_mass = new RooRealVar("CMS_zz4l_mass","CMS_zz4l_mass",100.,1000.);
	RooRealVar* rCMS_zz4l_KD2;
	RooRealVar* rCMS_zz4l_KD3;
	RooArgSet rContainer(*rweightFit,*rCMS_zz4l_KD1,*rCMS_zz4l_smd,*rCMS_zz4l_mass);

	if (iKD2 != iKD1 && iKD2 >= 0){
		if(iKD2>=firstIntKD) rCMS_zz4l_KD2 = new RooRealVar(toyKDName[1], toyKDName[1], -1., 1.);
		else rCMS_zz4l_KD2 = new RooRealVar(toyKDName[1], toyKDName[1], 0., 1.);
		rContainer.add(*rCMS_zz4l_KD2);
	};
	if (iKD3 != iKD1 && iKD3 != iKD2 && iKD3 >= 0){
		if(iKD3>=firstIntKD) rCMS_zz4l_KD3 = new RooRealVar(toyKDName[3], toyKDName[3], -1., 1.);
		else rCMS_zz4l_KD3 = new RooRealVar(toyKDName[3], toyKDName[3], 0., 1.);
		rContainer.add(*rCMS_zz4l_KD3);
	};

	for (int j=0; j<6; j++){
//	for (int j=0; j<1; j++){
		cout<<chan[j]<<endl;
		TChain *sig = new TChain("SelectedTree");
		TChain *gg=new TChain("SelectedTree_ggZZ");
		TChain *qq;
		if(writeToys==0) qq = new TChain("SelectedTree_qqZZ");
		else qq = new TChain("SelectedTree_qqZZ_Dedicated");
		TChain *zx=new TChain("SelectedTree_ZX");
		 sig->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_"+channame[j]+"/HZZ4lTree_H125p6_ShuffledSignalBkg.root");
		 gg->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_"+channame[j]+"/HZZ4lTree_H125p6_ShuffledSignalBkg.root");
		 qq->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_"+channame[j]+"/HZZ4lTree_H125p6_ShuffledSignalBkg.root");

		 zx->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/2mu2e/HZZ4lTree_H125p6_ShuffledSignalBkg.root");
		 zx->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/4e/HZZ4lTree_H125p6_ShuffledSignalBkg.root");
		 zx->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/4mu/HZZ4lTree_H125p6_ShuffledSignalBkg.root");
		 zx->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_7TeV/2mu2e/HZZ4lTree_H125p6_ShuffledSignalBkg.root");
		 zx->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_7TeV/4e/HZZ4lTree_H125p6_ShuffledSignalBkg.root");
		 zx->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_7TeV/4mu/HZZ4lTree_H125p6_ShuffledSignalBkg.root");

	int BSMFile=BSMSample;
	if (BSMSample<0 || BSMSample>=kNumHypo){
		cerr << "No hypothesis" << endl;
		assert(0);
	}
	else if (BSMSample>gapZZHypo && BSMSample < (kNumFiles + 1)){
		cout << "Passed gap hypothesis " << gapZZHypo << " while file range is still usable." << endl;
		BSMFile--;
	}
	else if (BSMSample == gapZZHypo || BSMSample >= (kNumFiles + 1)){
		BSMFile = 0;
		cout << "No dedicated sample is available, defaulting to " << BSMFile << endl;
	};
	cout << "BSM Hypothesis is " << BSMSample << ", and  dedicated file is " << BSMFile << endl;

	int chooseMCwgt=0;
	if(BSMSample>=firstNonZZHypo) chooseMCwgt=1;
	char cMCwgtchosen[30];
	sprintf(cMCwgtchosen,"%s",cMCwgt[chooseMCwgt]);
	char cMCwgtchosen_indexed[30];
	sprintf(cMCwgtchosen_indexed,"%s[%i]",cMCwgtchosen,BSMSample);
	cout << "Signal MC weight string: " << cMCwgtchosen_indexed << endl;

	char cBSMSampleAll_cut[100];
	sprintf(cBSMSampleAll_cut,"(ZZMass<140.6&&ZZMass>=105.6)*%s",cMCwgtchosen_indexed);
	cout << "BSM cut for all events: " << cBSMSampleAll_cut << endl;

	char cBSMSampleDedicated_cut[100];
	sprintf(cBSMSampleDedicated_cut,"(ZZMass<140.6&&ZZMass>=105.6 && EventSample==%i)*%s",BSMFile,cMCwgtchosen_indexed);
	cout << "BSM cut for dedicated events: " << cBSMSampleDedicated_cut << endl;

	TH1D *htempqqZZ=new TH1D("htqqZZ","",1,0.0,1000.);
	TH1D *htempggZZ=new TH1D("htggZZ","",1,0.0,1000.);
	TH1D *htempZX=new TH1D("htZX","",1,0.0,1000.);
	double bkgSum=0;

	TH1D *htempSM=new TH1D("htSM","",1,0.0,1000.);
	sig->Draw("ZZMass>>htSM","(ZZMass<140.6&&ZZMass>=105.6)*MC_CV_weight[0]");
	double sumSM = htempSM->GetSumOfWeights();
	cout << "Rsignal for SM hypothesis all files (without 4 GeV cut): " << sumSM << endl;

	TH1D *htempSMCut=new TH1D("htSMCut","",1,0.0,1000.);
	sig->Draw("ZZMass>>htSMCut","(ZZMass<140.6&&ZZMass>=105.6)*MC_CV_4GeVcut_weight[0]");
	double sumSMCut = htempSMCut->GetSumOfWeights();
	cout << "Rsignal for SM hypothesis all files (with 4 GeV cut): " << sumSMCut << endl;


	TH1D *htempBSM = new TH1D("htBSM","",1,0.0,1000.);
	sig->Draw("ZZMass>>htBSM",cBSMSampleAll_cut);
	double sumBSM = htempBSM->GetSumOfWeights();
	cout << "Rsignal for BSM hypothesis all files: " << sumBSM << endl;
	if(chooseMCwgt==1) sumBSM *= (sumSM/sumSMCut);

	TH1D *htempBSMDedicated = new TH1D("htBSMDedicated","",1,0.0,1000.);
	sig->Draw("ZZMass>>htBSMDedicated",cBSMSampleDedicated_cut);
	double sumBSMDedicated = htempBSMDedicated->GetSumOfWeights();//(EventSample==BSMFile)
	cout << "Rsignal for BSM hypothesis dedicated file: " << sumBSMDedicated << endl;

	double signalNormScale = (norm[j]*sumBSM/sumSM)/sumBSMDedicated;
	cout << "Nsignal for BSM hypothesis: " << signalNormScale*sumBSMDedicated << endl;
	double NBSMexpected = signalNormScale*sumBSMDedicated;
	if (signalStrength != 1){
		signalNormScale *= signalStrength;
		cout << "Nsignal for BSM hypothesis after signal strength modification: " << NBSMexpected << endl;
	};

	//double D_bkg, weightp3L;
	float D_bkg;
	float D_KD1, D_KD3, D_KD2, ZZMass;
	int sample;
	double weightqqZZL=0, weightggZZL=0, weightZX=0, weightFit=0, weightFit2=0;
	float MC_CV_weight[kNumHypo];
	float MC_weight_QQBGGProper[2];
	float MC_weight_noxsec;
	float MC_weight_Kfactor,ZXfake_weightProper;
	int EventSample;
	int ToyNumber=0;
	int eventIterator=0;
	double sumIndToyWgt=0;
	int store_NToyEvents=0;
	int counter_IndToyEvent=0;

	tsig[j]->Branch(toyKDName[0],&D_KD1,toyKDNameF[0]);
	if (iKD2 != iKD1 && iKD2 >= 0) tsig[j]->Branch(toyKDName[1], &D_KD2, toyKDNameF[1]);
	if (iKD3 != iKD1 && iKD3 != iKD2 && iKD3 >= 0) tsig[j]->Branch(toyKDName[3], &D_KD3, toyKDNameF[3]);
	tsig[j]->Branch(toyKDName[2],&D_bkg,toyKDNameF[2]);
	tsig[j]->Branch("CMS_zz4l_mass",&ZZMass,"CMS_zz4l_mass/F");
	tsig[j]->Branch("_weight_",&weightFit2,"_weight_/D");
	tsig[j]->Branch("ToyNumber",&ToyNumber,"ToyNumber/I");
	tsig[j]->Branch("NToyEvents",&store_NToyEvents,"NToyEvents/I");
	tggzz[j]->Branch(toyKDName[0],&D_KD1,toyKDNameF[0]);
	if (iKD2 != iKD1 && iKD2 >= 0) tggzz[j]->Branch(toyKDName[1], &D_KD2, toyKDNameF[1]);
	if (iKD3 != iKD1 && iKD3 != iKD2 && iKD3 >= 0) tggzz[j]->Branch(toyKDName[3], &D_KD3, toyKDNameF[3]);
	tggzz[j]->Branch(toyKDName[2],&D_bkg,toyKDNameF[2]);
	tggzz[j]->Branch("CMS_zz4l_mass",&ZZMass,"CMS_zz4l_mass/F");
	tggzz[j]->Branch("_weight_",&weightFit2,"_weight_/D");
	tggzz[j]->Branch("ToyNumber",&ToyNumber,"ToyNumber/I");
	tggzz[j]->Branch("NToyEvents",&store_NToyEvents,"NToyEvents/I");
	tqqzz[j]->Branch(toyKDName[0],&D_KD1,toyKDNameF[0]);
	if (iKD2 != iKD1 && iKD2 >= 0) tqqzz[j]->Branch(toyKDName[1], &D_KD2, toyKDNameF[1]);
	if (iKD3 != iKD1 && iKD3 != iKD2 && iKD3 >= 0) tqqzz[j]->Branch(toyKDName[3], &D_KD3, toyKDNameF[3]);
	tqqzz[j]->Branch(toyKDName[2],&D_bkg,toyKDNameF[2]);
	tqqzz[j]->Branch("CMS_zz4l_mass",&ZZMass,"CMS_zz4l_mass/F");
	tqqzz[j]->Branch("_weight_",&weightFit2,"_weight_/D");
	tqqzz[j]->Branch("ToyNumber",&ToyNumber,"ToyNumber/I");
	tqqzz[j]->Branch("NToyEvents",&store_NToyEvents,"NToyEvents/I");
	tzx[j]->Branch(toyKDName[0],&D_KD1,toyKDNameF[0]);
	if (iKD2 != iKD1 && iKD2 >= 0) tzx[j]->Branch(toyKDName[1], &D_KD2, toyKDNameF[1]);
	if (iKD3 != iKD1 && iKD3 != iKD2 && iKD3 >= 0) tzx[j]->Branch(toyKDName[3], &D_KD3, toyKDNameF[3]);
	tzx[j]->Branch(toyKDName[2],&D_bkg,toyKDNameF[2]);
	tzx[j]->Branch("CMS_zz4l_mass",&ZZMass,"CMS_zz4l_mass/F");
	tzx[j]->Branch("_weight_",&weightFit2,"_weight_/D");
	tzx[j]->Branch("ToyNumber",&ToyNumber,"ToyNumber/I");
	tzx[j]->Branch("NToyEvents",&store_NToyEvents,"NToyEvents/I");
	TRandom3 rand_sig(j+100);
	TRandom3 rand_ggzz(j+200);
	TRandom3 rand_qqzz(j+300);
	TRandom3 rand_zx(j+400);
	std::vector<int> count_sig;
	std::vector<int> count_ggzz;
	std::vector<int> count_qqzz;
	std::vector<int> count_zx;
	std::vector<double> fill_sig;
	std::vector<double> fill_ggzz;
	std::vector<double> fill_qqzz;
	std::vector<double> fill_zx;



			  //sig->SetBranchAddress("D_g1_vs_g4_phi0",&D_KD1);
			  //sig->SetBranchAddress("D_g4int_phi0",&D_KD2);
			  //sig->SetBranchAddress("D_bkg",&D_bkg);
			  sig->SetBranchAddress(KDlist[iKD1],&D_KD1);
	  		  cout << "sig KD1 name: " << KDlist[iKD1] << endl;
			  if (iKD2 != iKD1 && iKD2 >= 0){
				  sig->SetBranchAddress(KDlist[iKD2], &D_KD2);
				  cout << "sig KD2 name: " << KDlist[iKD2] << endl;
			  };
			  if (iKD3 != iKD1 && iKD3 != iKD2 && iKD3 >= 0) {
				  sig->SetBranchAddress(KDlist[iKD3], &D_KD3);
				  cout << "sig KD3 name: " << KDlist[iKD3] << endl;
			  };
			  sig->SetBranchAddress("D_bkg",&D_bkg);
			  sig->SetBranchAddress("ZZMass",&ZZMass);
			  sig->SetBranchAddress(cMCwgtchosen,MC_CV_weight);
		//      sig->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec);
		//      sig->SetBranchAddress("sample",&sample);
			  sig->SetBranchAddress("EventSample",&EventSample);
			  //sig->SetBranchAddress("weight1L",&weightp3L);

	TTree *outTree= new TTree("SelectedTree","SelectedTree");
	outTree->Branch(toyKDName[0],&D_KD1,toyKDNameF[0]);
	cout << "outTree KD1 name: " << toyKDName[0] << endl;
	if (iKD2 != iKD1 && iKD2 >= 0){
		outTree->Branch(toyKDName[1], &D_KD2, toyKDNameF[1]);
		cout << "outTree KD2 name: " << toyKDName[1] << endl;
	};
	if (iKD3 != iKD1 && iKD3 != iKD2 && iKD3 >= 0){
		outTree->Branch(toyKDName[3], &D_KD3, toyKDNameF[2]);
		cout << "outTree KD3 name: " << toyKDName[3] << endl;
	};
	outTree->Branch(toyKDName[2],&D_bkg,toyKDNameF[2]);
	cout << "outTree D_bkg name: " << toyKDName[2] << endl;
	outTree->Branch("CMS_zz4l_mass",&ZZMass,"CMS_zz4l_mass/F");
	outTree->Branch("_weight_",&weightFit,"_weight_/D");


			double sumSigValidate=0;
			int nentry = sig->GetEntries();
			int nFilled=0;
			for(int i=0;i<nentry;i++){
				sig->GetEntry(i);
				if(ZZMass<140.6 &&ZZMass>=105.6){
					if(EventSample==BSMFile && MC_CV_weight[BSMSample]>0){
						weightFit = signalNormScale*MC_CV_weight[BSMSample];
						sumSigValidate += weightFit;
						outTree->Fill();
						if(nFilled==0) cout << "Typical weight is " << weightFit << endl;
						nFilled++;

						if (counter_IndToyEvent == store_NToyEvents){
							store_NToyEvents = rand_sig.Poisson(NBSMexpected);
							counter_IndToyEvent = 0;
							sumIndToyWgt=0;
							if (store_NToyEvents == 0){
								count_sig.push_back(store_NToyEvents);
								fill_sig.push_back(sumIndToyWgt);
							};
						};
						if (counter_IndToyEvent < store_NToyEvents){
							sumIndToyWgt += MC_CV_weight[BSMSample];
							if (counter_IndToyEvent == store_NToyEvents - 1){
								count_sig.push_back(store_NToyEvents);
								fill_sig.push_back(sumIndToyWgt);
							};
							counter_IndToyEvent++;
						};
					}
				}
			}
			cout << "Nsignal filled for BSM: " << sumSigValidate << " with number of events filled " << nFilled << endl;
			cout << "Size of fill_sig: " << fill_sig.size() << endl;
			cout << "Size of count_sig: " << count_sig.size() << endl;
			int sum_of_countsig = 0;
			for (int ms = 0; ms < count_sig.size(); ms++){
				sum_of_countsig += count_sig[ms];
				if (count_sig[ms] == 0) sum_of_countsig++;
			};
			cout << "Sum of count_sig: " << sum_of_countsig << endl;
			double sum_of_fillsig = 0;
			for(int ms=0;ms<fill_sig.size();ms++) sum_of_fillsig+=fill_sig[ms];
			cout << "Sum of fill_sig: " << sum_of_fillsig << endl;

			sumIndToyWgt=0;
			store_NToyEvents=0;
			counter_IndToyEvent=0;
			for(int i=0;i<nentry;i++){
				sig->GetEntry(i);
				if(ZZMass<140.6 &&ZZMass>=105.6){
					if(EventSample==BSMFile && MC_CV_weight[BSMSample]>0){
						if(count_sig.size()<=ToyNumber) continue;
						store_NToyEvents = count_sig[ToyNumber];
						sumIndToyWgt = fill_sig[ToyNumber];
						double nexp=store_NToyEvents;
						double sumtoywgt=sumIndToyWgt;
						if(nexp==0) sumtoywgt=1;
						double scale = nexp/sumtoywgt;

						weightFit2 = scale*MC_CV_weight[BSMSample];
						tsig[j]->Fill();
						counter_IndToyEvent++;
						if (counter_IndToyEvent >= store_NToyEvents){
							counter_IndToyEvent = 0;
							ToyNumber++;
							store_NToyEvents=0;
							sumIndToyWgt=0;
						};
					}
				}
			}
			cout << "Ntoy possible from signal: " << ToyNumber << endl;
			if(maxNtoysPossible>ToyNumber) maxNtoysPossible=ToyNumber;
			ToyNumber=0;
			store_NToyEvents=0;
			sumIndToyWgt=0;
			counter_IndToyEvent=0;
			nFilled=0;

	/*******bkg**********/
		  qq->SetBranchAddress(KDlist[iKD1],&D_KD1);
  	  	  cout << "qqZZ KD1 name: " << KDlist[iKD1] << endl;
		  if (iKD2 != iKD1 && iKD2 >= 0){
			  qq->SetBranchAddress(KDlist[iKD2], &D_KD2);
	  	  	  cout << "qqZZ KD2 name: " << KDlist[iKD2] << endl;
		  };
		  if (iKD3 != iKD1 && iKD3 != iKD2 && iKD3 >= 0){
			  qq->SetBranchAddress(KDlist[iKD3], &D_KD3);
	  	  	  cout << "qqZZ KD3 name: " << KDlist[iKD3] << endl;
		  };
		  qq->SetBranchAddress("D_bkg",&D_bkg);
		  qq->SetBranchAddress("ZZMass",&ZZMass);
	   //   qq->SetBranchAddress("sample",&sample);
		  qq->SetBranchAddress("MC_weight_QQBGGProper",MC_weight_QQBGGProper);
		  qq->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec);
		  qq->SetBranchAddress("MC_weight_Kfactor",&MC_weight_Kfactor);
	  //    qq->SetBranchAddress("weightqqZZL",&weightqqZZL);
			qq->Draw("ZZMass>>htqqZZ","MC_weight_QQBGGProper[1]*MC_weight_noxsec*(ZZMass<140.6&&ZZMass>=105.6)");
			//qq->Draw("ZZMass>>ht","weightqqZZL*(ZZMass<140.6&&ZZMass>=105.6&&sample<1)");
			bkgSum = htempqqZZ->GetSumOfWeights();
			cout << "qqZZ dedicated sum: " << bkgSum << endl;
			nentry = qq->GetEntries();
			double sumQQZZValidate=0;
			for(int i=0;i<nentry;i++){
				qq->GetEntry(i);
				if(ZZMass<140.6 &&ZZMass>=105.6){
					//if(sample<1){
					//weightFit= normqq[j]/sumweight*weightqqZZL;
					weightFit = normqq[j]/bkgSum*MC_weight_QQBGGProper[1]*MC_weight_noxsec;
					outTree->Fill();
					sumQQZZValidate += weightFit;
					if(nFilled==0) cout << "Typical weight is " << weightFit << endl;
					nFilled++;

					if (counter_IndToyEvent == store_NToyEvents){
						store_NToyEvents = rand_qqzz.Poisson(normqq[j]);
						counter_IndToyEvent = 0;
						sumIndToyWgt=0;
						if (store_NToyEvents == 0){
							count_qqzz.push_back(store_NToyEvents);
							fill_qqzz.push_back(sumIndToyWgt);
						};
					};
					if (counter_IndToyEvent < store_NToyEvents){
						sumIndToyWgt += MC_weight_QQBGGProper[1]*MC_weight_noxsec;
						if (counter_IndToyEvent == store_NToyEvents - 1){
							count_qqzz.push_back(store_NToyEvents);
							fill_qqzz.push_back(sumIndToyWgt);
						};
						counter_IndToyEvent++;
					};
				}
			}
			cout << "qqZZ filled: " << sumQQZZValidate << " with number of events filled " << nFilled << endl;
			cout << "Size of fill_qqzz: " << fill_qqzz.size() << endl;
			cout << "Size of count_qqzz: " << count_qqzz.size() << endl;
			int sum_of_countqqzz = 0;
			for (int ms = 0; ms < count_qqzz.size(); ms++){
				sum_of_countqqzz += count_qqzz[ms];
				if (count_qqzz[ms] == 0) sum_of_countqqzz++;
			};
			cout << "Sum of count_qqzz: " << sum_of_countqqzz << endl;
			double sum_of_fillqqzz = 0;
			for(int ms=0;ms<fill_qqzz.size();ms++) sum_of_fillqqzz+=fill_qqzz[ms];
			cout << "Sum of fill_qqzz: " << sum_of_fillqqzz << endl;

			sumIndToyWgt=0;
			store_NToyEvents=0;
			counter_IndToyEvent=0;
			for(int i=0;i<nentry;i++){
				qq->GetEntry(i);
				if(ZZMass<140.6 &&ZZMass>=105.6){
					if(count_qqzz.size()<=ToyNumber) continue;
					store_NToyEvents = count_qqzz[ToyNumber];
					sumIndToyWgt = fill_qqzz[ToyNumber];
					double nexp=store_NToyEvents;
					double sumtoywgt=sumIndToyWgt;
					if(nexp==0) sumtoywgt=1;
					double scale = nexp/sumtoywgt;

					weightFit2 = scale*MC_weight_QQBGGProper[1]*MC_weight_noxsec;
					tqqzz[j]->Fill();
					counter_IndToyEvent++;
					if (counter_IndToyEvent >= store_NToyEvents){
						counter_IndToyEvent = 0;
						ToyNumber++;
						store_NToyEvents=0;
						sumIndToyWgt=0;
					}
				}
			}
			cout << "Ntoy possible from qqZZ: " << ToyNumber << endl;
			if(maxNtoysPossible>ToyNumber) maxNtoysPossible=ToyNumber;
			ToyNumber=0;
			store_NToyEvents=0;
			sumIndToyWgt=0;
			counter_IndToyEvent=0;
			nFilled=0;

		  gg->SetBranchAddress(KDlist[iKD1],&D_KD1);
  	  	  cout << "ggZZ KD1 name: " << KDlist[iKD1] << endl;
		  if (iKD2 != iKD1 && iKD2 >= 0){
			  gg->SetBranchAddress(KDlist[iKD2], &D_KD2);
			  cout << "ggZZ KD2 name: " << KDlist[iKD2] << endl;
		  };
		  if (iKD3 != iKD1 && iKD3 != iKD2 && iKD3 >= 0){
			  gg->SetBranchAddress(KDlist[iKD3], &D_KD3);
	  	  	  cout << "ggZZ KD3 name: " << KDlist[iKD3] << endl;
		  };
		  gg->SetBranchAddress("D_bkg",&D_bkg);
		  gg->SetBranchAddress("ZZMass",&ZZMass);
		  gg->SetBranchAddress("MC_weight_QQBGGProper",MC_weight_QQBGGProper);
		  gg->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec);
		  gg->SetBranchAddress("MC_weight_Kfactor",&MC_weight_Kfactor);
		  //gg->SetBranchAddress("weightggZZL",&weightggZZL);
	//      gg->SetBranchAddress("weightggZZL",&weightggZZL);
	//      gg->SetBranchAddress("sample",&sample);
		//gg->Draw("ZZMass>>ht","weightggZZL*(ZZMass<140.6&&ZZMass>=105.6&&sample>1&&sample<5)");
			gg->Draw("ZZMass>>htggZZ","(MC_weight_QQBGGProper[0]*MC_weight_noxsec*MC_weight_Kfactor)*(ZZMass<140.6&&ZZMass>=105.6)");
			bkgSum= htempggZZ->GetSumOfWeights();
			cout << "ggZZ dedicated sum: " << bkgSum << endl;
			nentry = gg->GetEntries();
			double sumGGZZValidate=0;
			for(int i=0;i<nentry;i++){
				gg->GetEntry(i);
				if(ZZMass<140.6 &&ZZMass>=105.6){
					//if(sample>1&&sample<5){
					//weightFit= normgg[j]/sumweight*weightggZZL;
					weightFit = normgg[j]/bkgSum*MC_weight_QQBGGProper[0]*MC_weight_noxsec*MC_weight_Kfactor;
					outTree->Fill();
					sumGGZZValidate += weightFit;
					if(nFilled==0) cout << "Typical weight is " << weightFit << endl;
					nFilled++;

					if (counter_IndToyEvent == store_NToyEvents){
						store_NToyEvents = rand_ggzz.Poisson(normgg[j]);
						counter_IndToyEvent = 0;
						sumIndToyWgt=0;
						if (store_NToyEvents == 0){
							count_ggzz.push_back(store_NToyEvents);
							fill_ggzz.push_back(sumIndToyWgt);
						};
					};
					if (counter_IndToyEvent < store_NToyEvents){
						sumIndToyWgt += MC_weight_QQBGGProper[0]*MC_weight_noxsec*MC_weight_Kfactor;
						if (counter_IndToyEvent == store_NToyEvents - 1){
							count_ggzz.push_back(store_NToyEvents);
							fill_ggzz.push_back(sumIndToyWgt);
						};
						counter_IndToyEvent++;
					};
				}
			}
			cout << "ggZZ filled: " << sumGGZZValidate << " with number of events filled " << nFilled << endl;
			cout << "Size of fill_ggzz: " << fill_ggzz.size() << endl;
			cout << "Size of count_ggzz: " << count_ggzz.size() << endl;
			int sum_of_countggzz = 0;
			for (int ms = 0; ms < count_ggzz.size(); ms++){
				sum_of_countggzz += count_ggzz[ms];
				if (count_ggzz[ms] == 0) sum_of_countggzz++;
			};
			cout << "Sum of count_ggzz: " << sum_of_countggzz << endl;
			double sum_of_fillggzz = 0;
			for(int ms=0;ms<fill_ggzz.size();ms++) sum_of_fillggzz+=fill_ggzz[ms];
			cout << "Sum of fill_ggzz: " << sum_of_fillggzz << endl;

			sumIndToyWgt=0;
			store_NToyEvents=0;
			counter_IndToyEvent=0;
			for(int i=0;i<nentry;i++){
				gg->GetEntry(i);
				if(ZZMass<140.6 &&ZZMass>=105.6){
					if(count_ggzz.size()<=ToyNumber) continue;
					store_NToyEvents = count_ggzz[ToyNumber];
					sumIndToyWgt = fill_ggzz[ToyNumber];
					double nexp=store_NToyEvents;
					double sumtoywgt=sumIndToyWgt;
					if(nexp==0) sumtoywgt=1;
					double scale = nexp/sumtoywgt;

					weightFit2 = scale*MC_weight_QQBGGProper[0]*MC_weight_noxsec*MC_weight_Kfactor;
					tggzz[j]->Fill();
					counter_IndToyEvent++;
					if (counter_IndToyEvent >= store_NToyEvents){
						counter_IndToyEvent = 0;
						ToyNumber++;
						store_NToyEvents=0;
						sumIndToyWgt=0;
					}
				}
			}
			cout << "Ntoy possible from ggZZ: " << ToyNumber << endl;
			if(maxNtoysPossible>ToyNumber) maxNtoysPossible=ToyNumber;
			ToyNumber=0;
			store_NToyEvents=0;
			sumIndToyWgt=0;
			counter_IndToyEvent=0;
			nFilled=0;

		  zx->SetBranchAddress(KDlist[iKD1],&D_KD1);
		  cout << "ZX KD1 name: " << KDlist[iKD1] << endl;
		  if (iKD2 != iKD1 && iKD2 >= 0){
			  zx->SetBranchAddress(KDlist[iKD2], &D_KD2);
			  cout << "ZX KD2 name: " << KDlist[iKD2] << endl;
		  };
		  if (iKD3 != iKD1 && iKD3 != iKD2 && iKD3 >= 0){
			  zx->SetBranchAddress(KDlist[iKD3], &D_KD3);
	  	  	  cout << "ZX KD3 name: " << KDlist[iKD3] << endl;
		  };
		  zx->SetBranchAddress("D_bkg",&D_bkg);
		  zx->SetBranchAddress("ZZMass",&ZZMass);
		  //zx->SetBranchAddress("FakeRate_w",&FakeRate_w);
	//      zx->SetBranchAddress("weightZX",&weightZX);
	//      zx->SetBranchAddress("sample",&sample);
		  zx->SetBranchAddress("ZXfake_weightProper",&ZXfake_weightProper);
	//zx->Draw("ZZMass>>ht","(ZZMass<140.6&&ZZMass>=105.6)*weightZX");
			zx->Draw("ZZMass>>htZX","(ZZMass<140.6&&ZZMass>=105.6 && ZXfake_weightProper>0)*ZXfake_weightProper");
			bkgSum= htempZX->GetSumOfWeights();
			cout << "ZX dedicated sum: " << bkgSum << endl;
			nentry = zx->GetEntries();
			double sumZXValidate=0;
			for(int i=0;i<nentry;i++){
				zx->GetEntry(i);
				if(ZZMass<140.6 && ZZMass>=105.6 && ZXfake_weightProper>0){
					//weightFit= weightZX*normzx[j]/sumweight;
					weightFit= ZXfake_weightProper*normzx[j]/bkgSum;
					outTree->Fill();
					sumZXValidate += weightFit;
					if(nFilled==0) cout << "Typical weight is " << weightFit << endl;
					nFilled++;

					if (counter_IndToyEvent == store_NToyEvents){
						store_NToyEvents = rand_zx.Poisson(normzx[j]);
						counter_IndToyEvent = 0;
						sumIndToyWgt=0;
						if (store_NToyEvents == 0){
							count_zx.push_back(store_NToyEvents);
							fill_zx.push_back(sumIndToyWgt);
						};
					};
					if (counter_IndToyEvent < store_NToyEvents){
						sumIndToyWgt += ZXfake_weightProper;
						if (counter_IndToyEvent == store_NToyEvents - 1){
							count_zx.push_back(store_NToyEvents);
							fill_zx.push_back(sumIndToyWgt);
						};
						counter_IndToyEvent++;
					};
				}
			}
			cout << "ZX filled: " << sumZXValidate << " with number of events filled " << nFilled << endl;
			cout << "Size of fill_zx: " << fill_zx.size() << endl;
			cout << "Size of count_zx: " << count_zx.size() << endl;
			int sum_of_countzx = 0;
			for (int ms = 0; ms < count_zx.size(); ms++){
				sum_of_countzx += count_zx[ms];
				if (count_zx[ms] == 0) sum_of_countzx++;
			};
			cout << "Sum of count_zx: " << sum_of_countzx << endl;
			double sum_of_fillzx = 0;
			for(int ms=0;ms<fill_zx.size();ms++) sum_of_fillzx+=fill_zx[ms];
			cout << "Sum of fill_zx: " << sum_of_fillzx << endl;

			sumIndToyWgt=0;
			store_NToyEvents=0;
			counter_IndToyEvent=0;
			for(int i=0;i<nentry;i++){
				zx->GetEntry(i);
				if(ZZMass<140.6 &&ZZMass>=105.6 && ZXfake_weightProper>0){
					if(count_zx.size()<=ToyNumber) continue;
					store_NToyEvents = count_zx[ToyNumber];
					sumIndToyWgt = fill_zx[ToyNumber];
					double nexp=store_NToyEvents;
					double sumtoywgt=sumIndToyWgt;
					if(nexp==0) sumtoywgt=1;
					double scale = nexp/sumtoywgt;

					weightFit2 = scale*ZXfake_weightProper;
					tzx[j]->Fill();
					counter_IndToyEvent++;
					if (counter_IndToyEvent >= store_NToyEvents){
						counter_IndToyEvent = 0;
						ToyNumber++;
						store_NToyEvents=0;
						sumIndToyWgt=0;
					}
				}
			}
			cout << "Ntoy possible from ZX: " << ToyNumber << endl;
			if(maxNtoysPossible>ToyNumber) maxNtoysPossible=ToyNumber;
			ToyNumber=0;
			store_NToyEvents=0;
			sumIndToyWgt=0;
			counter_IndToyEvent=0;
			nFilled=0;


		data[j]=new RooDataSet(Form("data%d",j),Form("data%d",j),outTree,rContainer,"","_weight_");
		data[j]->weightError(RooAbsData::None);

		RooArgSet toyArgs(rContainer);
		toyArgs.add(cat);
		toyf03[j]= new RooDataSet (Form("toy_asimov%d",j),Form("toy_asimov%d",j), toyArgs, Index(cat),WeightVar("_weight_"), Import(chan[j],*data[j]) );
		if(toyf03[j]->isWeighted()) cout << "Asimov dataset is weighted." << endl;
	}
	char cDir[]="./AsimovToys";
	char coutput[50];
/*
	char coutput_main[50];
	char cSignalStrength[30];
	int signalStrength_I = signalStrength;
	int signalStrength_K = 100*signalStrength;
	int signalStrength_Mod = signalStrength_K%100;
	if(signalStrength_Mod!=0) sprintf(cSignalStrength,"SigStr_%ip%i_",signalStrength_I,signalStrength_Mod);
	else sprintf(cSignalStrength,"%s_%i","SigStr",signalStrength_I);
	cout << cSignalStrength << endl;
	if(signalStrength!=1) sprintf(coutput_main,"%s/Asimov_Hypo%i_%s",cDir,BSMSample,cSignalStrength);
	else sprintf(coutput_main,"%s/Asimov_Hypo%i",cDir,BSMSample);
	cout << coutput_main << endl;
	if(iKD2>=0){
		sprintf(coutput,"%s_KD1_%i_vs_KD2_%i.root",coutput_main,iKD1,iKD2);
	}
	else{
		sprintf(coutput,"%s_KD1_%i.root",coutput_main,iKD1);
	};
*/


	if(iKD3>=0){
		if(writeToys==0) sprintf(coutput,"%s/Asimov_Hypo%i_KD1_%i_vs_KD2_%i_vs_KD3_%i.root",cDir,BSMSample,iKD1,iKD2,iKD3);
		else sprintf(coutput,"%s/Toys_Hypo%i_KD1_%i_vs_KD2_%i_vs_KD3_%i.root",cDir,BSMSample,iKD1,iKD2,iKD3);
	}
	else if(iKD2>=0){
		if(writeToys==0) sprintf(coutput,"%s/Asimov_Hypo%i_KD1_%i_vs_KD2_%i.root",cDir,BSMSample,iKD1,iKD2);
		else sprintf(coutput,"%s/Toys_Hypo%i_KD1_%i_vs_KD2_%i.root",cDir,BSMSample,iKD1,iKD2);
	}
	else{
		if(writeToys==0) sprintf(coutput,"%s/Asimov_Hypo%i_KD1_%i.root",cDir,BSMSample,iKD1);
		else sprintf(coutput,"%s/Toys_Hypo%i_KD1_%i.root",cDir,BSMSample,iKD1);
	};
	cout << coutput << endl;
	TFile *toyfile = new TFile(coutput,"recreate");
	cout << coutput << endl;
	TDirectory *cdtof = toyfile->mkdir("toys");

	RooDataSet *toyf03_total = toyf03[0]; 
	for(int j=0;j<6;j++){
		if(j!=0) toyf03_total->append(*(toyf03[j]));
		toyf03[j]->Print("v");		
		RooCategory *cata = dynamic_cast<RooCategory *>(toyf03_total->get()->find("CMS_channel"));
		int na = cata->numBins((const char *)0);
		cout<<na<<endl;

		TString cValPlot = "Validate_";
		cValPlot.Append(chan[j]);
		TString cCanValPlot = "canvasValidate_";
		cCanValPlot.Append(chan[j]);
		TCanvas* can = new TCanvas(cCanValPlot,chan[j],800,800);
		RooPlot* valPlot = new RooPlot(cValPlot,chan[j],*rCMS_zz4l_mass,105.6,140.6,80);
		toyf03[j]->plotOn(valPlot);
		can->cd();
		valPlot->Draw();
		toyfile->WriteTObject(can);
		delete valPlot;
		can->Close();
		cdtof->cd();

		cout << "Sum of weights in dataset is " << toyf03[j]->sumEntries("CMS_zz4l_mass<140.6 && CMS_zz4l_mass>=105.6") << endl;
	}
	cout << "Final toy name: " << toyf03_total->GetName() << endl;
	toyf03_total->SetName("toy_asimov");
	toyf03_total->SetTitle("toy_asimov");
	if(toyf03_total->isWeighted()) cout << "Total asimov dataset is weighted." << endl;
	cout << "Sum of weights in final dataset is " << toyf03_total->sumEntries("CMS_zz4l_mass<140.6 && CMS_zz4l_mass>=105.6") << endl;
	/****** pure bkg toys ********/
	 //TFile workspace("cards_dcp_8TeV/HCG/125.6/bkgonly.root");
	//TFile workspace("cards_NoRebinAdapSmoothReweightMirror_8tev/HCG/125.6/higgsCombinebkgonly.GenerateOnly.mH125.6.123456.root");
	//TFile workspace("cards_NoRebinAdapSmoothReweightMirror_8tev/HCG/125.6/higgsCombinebkgOnly.GenerateOnly.mH125.6.123456.root");


	//TFile workspace("cards_test6_8TeV/HCG/125.6/higgsCombinebkgOnly.GenerateOnly.mH125.6.123456.root");
	// RooDataSet * toyDataset = (RooDataSet*) workspace.Get("toys/toy_asimov");
	// int d1_entries = toyDataset->numEntries();
	//RooCategory *catb = dynamic_cast<RooCategory *>(toyDataset->get()->find("CMS_channel"));
	//    int nb = catb->numBins((const char *)0);
	//cout<<nb<<endl; 

		RooArgSet* set;
	//for(int i=0;i< d2_entries;i++){
	//  set = toyf03->get(i);
	//toyDataset->add(*set,toyf03->weight());
	////cout<<toyf03->weight()<<endl;
	//}
	//toyDataset->Print("v");

	RooCategory *cata = dynamic_cast<RooCategory *>(toyf03_total->get()->find("CMS_channel"));
	int na = cata->numBins((const char *)0);
	cout<<na<<endl; 
	/*** add pure bkg toys***/
	//toyf03_total->append(*toyDataset);

	cout<<toyf03_total->numEntries()<<endl;
	cdtof->cd();
//	cout << "cdtof" << endl;
	if(writeToys==0) toyf03_total->Write("toy_asimov");
//	cout << "write" << endl;

	toyfile->cd();
	for(int j=0;j<6;j++){
		toyfile->WriteTObject(tsig[j]);
		toyfile->WriteTObject(tggzz[j]);
		toyfile->WriteTObject(tqqzz[j]);
		toyfile->WriteTObject(tzx[j]);
		delete tsig[j];
		delete tggzz[j];
		delete tqqzz[j];
		delete tzx[j];
	};

	cout << "Successfully closed the toys file!" << endl;
	cout << "Maximum possible number of toys: " << maxNtoysPossible << endl;

	toyfile->Close();
}
