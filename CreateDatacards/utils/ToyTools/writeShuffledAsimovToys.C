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

void writeShuffledAsimovToys(int BSMSample=6, int iKD1=0, int iKD2=1, double signalStrength=1){
	const int kNumHypo=44;
	const int kNumFiles=24;
	const int gapZZHypo=13;
	const int firstNonZZHypo=33;

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
		"CMS_zz4l_smd"
	};
	TString toyKDNameF[] = {
		"CMS_zz4l_KD1/F",
		"CMS_zz4l_KD2/F",
		"CMS_zz4l_smd/F"
	};
	char* cMCwgt[2] = {
		"MC_CV_weight",
		"MC_CV_4GeVcut_weight"
	};

	TString chan[6]={"ch1","ch2","ch3","ch4","ch5","ch6"};
	int order[6]={2,1,0,4,3,5};
	RooDataSet *data[6];
	RooDataSet *toyf03[6];
	TString channameO[6]={"7TeV/4e","7TeV/4mu","7TeV/2mu2e","8TeV/4mu","8TeV/2mu2e","8TeV/4e"};
	double normO[6]=		{	0.7194,	1.2859,	1.7264,	6.0802,	7.9085,	3.1441	}; 
	double normzxO[6]=	{	0.6226,	0.2230,	1.0628,	1.1878,	4.2929,	2.7676	};
	double normqqO[6]=	{	0.8386,	1.7971,	2.2456,	7.6478,	8.8585,	2.9364	};
	double normggO[6]=	{	0.0341,	0.0625,	0.0741,	0.4131,	0.5005,	0.2041	};

	double norm[6] = { 0 };
	double normgg[6] = { 0 };
	double normzx[6] = { 0 };
	double normqq[6] = { 0 };
	TString channame[6];

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
	}
	for (int j=0; j<6; j++){
//	for (int j=0; j<1; j++){
		cout<<chan[j]<<endl;
		TChain *sig = new TChain("SelectedTree");
		TChain *gg=new TChain("SelectedTree_ggZZ");
		TChain *qq=new TChain("SelectedTree_qqZZ");
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


	// sig->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_5_3_9/src/ZZMatrixElement/MELA/test/"+channame[j]+"/signal.root");
	//		zx->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_5_3_9/src/ZZMatrixElement/MELA/test/2mu2e/bkg.root");
	//		zx->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_5_3_9/src/ZZMatrixElement/MELA/test/4e/bkg.root");
	//		zx->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_5_3_9/src/ZZMatrixElement/MELA/test/4mu/bkg.root");
	//		gg->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_5_3_9/src/ZZMatrixElement/MELA/test/"+channame[j]+"/bkg.root");
	//		qq->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_5_3_9/src/ZZMatrixElement/MELA/test/"+channame[j]+"/bkg.root");

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
	if (signalStrength != 1){
		signalNormScale *= signalStrength;
		cout << "Nsignal for BSM hypothesis after signal strength modification: " << signalNormScale*sumBSMDedicated << endl;
	};


	//double D_bkg, weightp3L;
	float D_bkg;
	float D_g1g4, weightFit, D_cp, ZZMass;
	int sample;
	double weightqqZZL, weightggZZL, weightZX;
	float MC_CV_weight[kNumHypo];
	float MC_weight_QQBGGProper[2];
	float MC_weight_noxsec;
	float MC_weight_Kfactor,ZXfake_weightProper;
	int EventSample;

			  //sig->SetBranchAddress("D_g1_vs_g4_phi0",&D_g1g4);
			  //sig->SetBranchAddress("D_g4int_phi0",&D_cp);
			  //sig->SetBranchAddress("D_bkg",&D_bkg);
			  sig->SetBranchAddress(KDlist[iKD1],&D_g1g4);
	  		  cout << "sig KD1 name: " << KDlist[iKD1] << endl;
			  if (iKD2 != iKD1 && iKD2 >= 0){
				  sig->SetBranchAddress(KDlist[iKD2], &D_cp);
				  cout << "sig KD2 name: " << KDlist[iKD2] << endl;
			  };
			  sig->SetBranchAddress("D_bkg",&D_bkg);
			  sig->SetBranchAddress("ZZMass",&ZZMass);
			  sig->SetBranchAddress(cMCwgtchosen,MC_CV_weight);
		//      sig->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec);
		//      sig->SetBranchAddress("sample",&sample);
			  sig->SetBranchAddress("EventSample",&EventSample);
			  //sig->SetBranchAddress("weight1L",&weightp3L);

	  RooRealVar rweightFit("_weight_","_weight_",0.,1.e20);
	  RooRealVar rCMS_zz4l_pseudoKD(toyKDName[0],toyKDName[0],0.,1.);
	  if (iKD2 != iKD1 && iKD2 >= 0){
		  if(iKD2>=firstIntKD) RooRealVar rCMS_zz4l_dcp(toyKDName[1], toyKDName[1], -1., 1.);
		  else RooRealVar rCMS_zz4l_dcp(toyKDName[1], toyKDName[1], 0., 1.);
	  };
	  RooRealVar rCMS_zz4l_smd(toyKDName[2],toyKDName[2],0.,1.);
	//  RooRealVar rCMS_zz4l_smd("CMS_zz4l_smd","CMS_zz4l_smd",0.,1.);
	  RooRealVar rCMS_zz4l_mass("CMS_zz4l_mass","CMS_zz4l_mass",100.,1000.);
	TTree *outTree= new TTree("SelectedTree","SelectedTree");
	outTree->Branch(toyKDName[0],&D_g1g4,toyKDNameF[0]);
	cout << "outTree KD1 name: " << toyKDName[0] << endl;
	if (iKD2 != iKD1 && iKD2 >= 0){
		outTree->Branch(toyKDName[1], &D_cp, toyKDNameF[1]);
		cout << "outTree KD2 name: " << toyKDName[1] << endl;
	};
	outTree->Branch(toyKDName[2],&D_bkg,toyKDNameF[2]);
	cout << "outTree D_bkg name: " << toyKDName[2] << endl;
	outTree->Branch("CMS_zz4l_mass",&ZZMass,"CMS_zz4l_mass/F");
	outTree->Branch("_weight_",&weightFit,"_weight_/F");


			double sumSigValidate=0;
			int nentry = sig->GetEntries();
			for(int i=0;i<nentry;i++){
				sig->GetEntry(i);
				if(ZZMass<140.6 &&ZZMass>=105.6){
					if(EventSample==BSMFile && MC_CV_weight[BSMSample]>0){
						weightFit = signalNormScale*MC_CV_weight[BSMSample];
						sumSigValidate += weightFit;
						outTree->Fill();
					}
				}
			}
			cout << "Nsignal filled for BSM: " << sumSigValidate << endl;


	/*******bkg**********/
		  qq->SetBranchAddress(KDlist[iKD1],&D_g1g4);
  	  	  cout << "qqZZ KD1 name: " << KDlist[iKD1] << endl;
		  if (iKD2 != iKD1 && iKD2 >= 0){
			  qq->SetBranchAddress(KDlist[iKD2], &D_cp);
	  	  	  cout << "qqZZ KD2 name: " << KDlist[iKD2] << endl;
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
					//}
				}
			}
			cout << "qqZZ filled: " << sumQQZZValidate << endl;

		  gg->SetBranchAddress(KDlist[iKD1],&D_g1g4);
  	  	  cout << "ggZZ KD1 name: " << KDlist[iKD1] << endl;
		  if (iKD2 != iKD1 && iKD2 >= 0){
			  gg->SetBranchAddress(KDlist[iKD2], &D_cp);
			  cout << "ggZZ KD2 name: " << KDlist[iKD2] << endl;
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
				//}
				}
			}
			cout << "ggZZ filled: " << sumGGZZValidate << endl;

		  zx->SetBranchAddress(KDlist[iKD1],&D_g1g4);
		  cout << "ZX KD1 name: " << KDlist[iKD1] << endl;
		  if (iKD2 != iKD1 && iKD2 >= 0){
			  zx->SetBranchAddress(KDlist[iKD2], &D_cp);
			  cout << "ZX KD2 name: " << KDlist[iKD2] << endl;
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
			//if(sample==5){
				if(ZZMass<140.6 && ZZMass>=105.6 && ZXfake_weightProper>0){
					//weightFit= weightZX*normzx[j]/sumweight;
					weightFit= ZXfake_weightProper*normzx[j]/bkgSum;
					outTree->Fill();
					sumZXValidate += weightFit;
					//}
				}
			}
			cout << "ZX filled: " << sumZXValidate << endl;

		data[j]=new RooDataSet(Form("data%d",j),Form("data%d",j),outTree,RooArgSet(rCMS_zz4l_pseudoKD,rCMS_zz4l_dcp,rCMS_zz4l_smd,rCMS_zz4l_mass,rweightFit),"","_weight_");

		toyf03[j]= new RooDataSet (Form("toy_asimov%d",j),Form("toy_asimov%d",j), RooArgSet(rCMS_zz4l_pseudoKD,rCMS_zz4l_dcp,rCMS_zz4l_smd,rCMS_zz4l_mass,rweightFit,cat), Index(cat),WeightVar("_weight_"), Import(chan[j],*data[j]) );
		//toyf03[j]->Print("v");
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
	if(iKD2>=0){
		sprintf(coutput,"%s/Asimov_Hypo%i_KD1_%i_vs_KD2_%i.root",cDir,BSMSample,iKD1,iKD2);
	}
	else{
		sprintf(coutput,"%s/Asimov_Hypo%i_KD1_%i.root",cDir,BSMSample,iKD1);
	};
	cout << coutput << endl;
	TFile *toyfile = new TFile(coutput,"recreate");
	cout << coutput << endl;
	TDirectory *cdtof = toyfile->mkdir("toys");
	RooDataSet *toyf03_total = toyf03[0]; 
	toyf03[0]->Print("v");
	for(int j=1;j<6;j++){
		toyf03_total->append(*(toyf03[j]));
		toyf03[j]->Print("v");
		RooCategory *cata = dynamic_cast<RooCategory *>(toyf03_total->get()->find("CMS_channel"));
		int na = cata->numBins((const char *)0);
		cout<<na<<endl;
	}
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
	toyf03_total->Write("toy_asimov");
//	cout << "write" << endl;
	toyfile->Close();
	cout << "Successfully closed the toys file!" << endl;
}
