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


void writetoys_ShuffledIndividual(int BSMSample=0, int iKD1=0, int iKD2=1, int iKD3=-1, const int maxNToysPossible=1500, int weightedEvents=1, int writeToys=0){
	const int kNumHypo=44;
	const int kNumFiles=24;
	const int gapZZHypo=13;
	const int firstNonZZHypo=33;

/***** RUN CONDITIONS *****/
	if(!(weightedEvents==0 || weightedEvents==1)) assert(0);
	if(writeToys<0) assert(0);
	if(iKD1<0) assert(0);
	if(iKD3>=0 && iKD2>=0) assert(0);
	if(BSMSample>=kNumHypo) assert(0);
/*** END RUN CONDITIONS ***/
	char cDir[]="./AsimovToys";
	char cinput[100];
	char coutput[100];

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
	TString toyKDName_Alt[] = {
		"CMS_zz4l_KD1",
		"CMS_zz4l_KD2",
		"CMS_zz4l_smd",
		"CMS_zz4l_KD3"
	};
	TString toyKDNameF_Alt[] = {
		"CMS_zz4l_KD1/F",
		"CMS_zz4l_KD2/F",
		"CMS_zz4l_smd/F",
		"CMS_zz4l_KD3/F"
	};


	TString chan[6]={"ch1","ch2","ch3","ch4","ch5","ch6"};
	int order[6]={0,1,2,3,4,5};
	RooDataSet *data[6];
	RooDataSet *toyf03[6];

	RooCategory cat("CMS_channel","CMS_channel");

	for(int j=0;j<6;j++){
		cat.defineType(chan[j],j);
		cat.setLabel(chan[j]);
	}

	TChain *sig[6]; 
	TChain *gg[6];  
	TChain *qq[6];  
	TChain *zx[6];  

	cout << "Starting..." << endl;

	if(iKD3>=0){
		if(writeToys==0) sprintf(cinput,"%s/Asimov_Hypo%i_KD1_%i_vs_KD2_%i_vs_KD3_%i.root",cDir,BSMSample,iKD1,-iKD2,iKD3);
		else sprintf(cinput,"%s/Toys_Hypo%i_KD1_%i_vs_KD2_%i_vs_KD3_%i.root",cDir,BSMSample,iKD1,-iKD2,iKD3);
	}
	else if(iKD2>=0){
		if(writeToys==0) sprintf(cinput,"%s/Asimov_Hypo%i_KD1_%i_vs_KD2_%i_vs_KD3_%i.root",cDir,BSMSample,iKD1,iKD2,-iKD3);
		else sprintf(cinput,"%s/Toys_Hypo%i_KD1_%i_vs_KD2_%i_vs_KD3_%i.root",cDir,BSMSample,iKD1,iKD2,-iKD3);
	}
	else{
		if(writeToys==0) sprintf(cinput,"%s/Asimov_Hypo%i_KD1_%i.root",cDir,BSMSample,iKD1);
		else sprintf(cinput,"%s/Toys_Hypo%i_KD1_%i.root",cDir,BSMSample,iKD1);
	};
	cout << cinput << endl;

	if (weightedEvents == 1){
		if (iKD3 >= 0){
			if (writeToys == 0) sprintf(coutput, "%s/WeightedSynchToys_Hypo%i_KD1_%i_vs_KD2_%i.root", cDir, BSMSample, iKD1, iKD3);
			else sprintf(coutput, "%s/WeightedSynchToys_NoAsimov_Hypo%i_KD1_%i_vs_KD2_%i.root", cDir, BSMSample, iKD1, iKD3);
		}
		else if (iKD2 >= 0){
			if (writeToys == 0) sprintf(coutput, "%s/WeightedSynchToys_Hypo%i_KD1_%i_vs_KD2_%i.root", cDir, BSMSample, iKD1, iKD2);
			else sprintf(coutput, "%s/WeightedSynchToys_NoAsimov_Hypo%i_KD1_%i_vs_KD2_%i.root", cDir, BSMSample, iKD1, iKD2);
		}
		else{
			if (writeToys == 0) sprintf(coutput, "%s/WeightedSynchToys_Hypo%i_KD1_%i.root", cDir, BSMSample, iKD1);
			else sprintf(coutput, "%s/WeightedSynchToys_NoAsimov_Hypo%i_KD1_%i.root", cDir, BSMSample, iKD1);
		};
	}
	else{
		if(iKD3>=0){
			if(writeToys==0) sprintf(coutput,"%s/UnweightedSynchToys_Hypo%i_KD1_%i_vs_KD2_%i.root",cDir,BSMSample,iKD1,iKD3);
			else sprintf(coutput,"%s/UnweightedSynchToys_NoAsimov_Hypo%i_KD1_%i_vs_KD2_%i.root",cDir,BSMSample,iKD1,iKD3);
		}
		else if(iKD2>=0){
			if(writeToys==0) sprintf(coutput,"%s/UnweightedSynchToys_Hypo%i_KD1_%i_vs_KD2_%i.root",cDir,BSMSample,iKD1,iKD2);
			else sprintf(coutput,"%s/UnweightedSynchToys_NoAsimov_Hypo%i_KD1_%i_vs_KD2_%i.root",cDir,BSMSample,iKD1,iKD2);
		}
		else{
			if(writeToys==0) sprintf(coutput,"%s/UnweightedSynchToys_Hypo%i_KD1_%i.root",cDir,BSMSample,iKD1);
			else sprintf(coutput,"%s/UnweightedSynchToys_NoAsimov_Hypo%i_KD1_%i.root",cDir,BSMSample,iKD1);
		};
	};
	cout << coutput << endl;

	TFile *toyfile = new TFile(coutput,"recreate");
	TDirectory *cdtof = toyfile->mkdir("toys");

	for (int j=0; j<6; j++){
		TString ctoys_sig="ToyEvents_Sig_";
		ctoys_sig.Append(chan[j]);
		TString ctoys_ggzz="ToyEvents_ggZZ_";
		ctoys_ggzz.Append(chan[j]);
		TString ctoys_qqzz="ToyEvents_qqZZ_";
		ctoys_qqzz.Append(chan[j]);
		TString ctoys_zx="ToyEvents_ZX_";
		ctoys_zx.Append(chan[j]);

		sig[j] = new TChain(ctoys_sig);
		gg[j] = new TChain(ctoys_ggzz);
		qq[j] = new TChain(ctoys_qqzz);
		zx[j] = new TChain(ctoys_zx);

		sig[j]->Add(cinput);
		gg[j]->Add(cinput);
		qq[j]->Add(cinput);
		zx[j]->Add(cinput);
	}

	RooRealVar* rweightFit = new RooRealVar("_weight_","_weight_",0.,1.e20);
	RooRealVar* rCMS_zz4l_KD1 = new RooRealVar(toyKDName_Alt[0],toyKDName_Alt[0],0.,1.);
	RooRealVar* rCMS_zz4l_smd = new RooRealVar(toyKDName_Alt[2],toyKDName_Alt[2],0.,1.);
	RooRealVar* rCMS_zz4l_mass = new RooRealVar("CMS_zz4l_mass","CMS_zz4l_mass",100.,1000.);
	RooRealVar* rCMS_zz4l_KD2;
	RooArgSet rContainer(*rCMS_zz4l_KD1,*rCMS_zz4l_smd,*rCMS_zz4l_mass);
	if(weightedEvents==1) rContainer.add(*rweightFit);
	if ( (iKD2 != iKD1 && iKD2 >= 0) ){
		if(iKD2>=firstIntKD) rCMS_zz4l_KD2 = new RooRealVar(toyKDName_Alt[1], toyKDName_Alt[1], -1., 1.);
		else rCMS_zz4l_KD2 = new RooRealVar(toyKDName_Alt[1], toyKDName_Alt[1], 0., 1.);
		rContainer.add(*rCMS_zz4l_KD2);
	};
	if ( (iKD3 != iKD1 && iKD3 >= 0) ){
		if(iKD3>=firstIntKD) rCMS_zz4l_KD2 = new RooRealVar(toyKDName_Alt[1], toyKDName_Alt[1], -1., 1.);
		else rCMS_zz4l_KD2 = new RooRealVar(toyKDName_Alt[1], toyKDName_Alt[1], 0., 1.);
		rContainer.add(*rCMS_zz4l_KD2);
	};

	float D_bkg;
	float D_KD1, D_KD2, ZZMass;
	int ToyNumber=0;
	double weightqqZZL = 0, weightggZZL = 0, weightZX = 0, weightFit = 0, weightFit2 = 0;

	TTree* outTree[6*maxNToysPossible];
	char cTreeName[]="SelectedTree";
	for (int iToy = 0; iToy < maxNToysPossible; iToy++){
		for (int j = 0; j < 6; j++){
			char cOutTree[50];
			sprintf(cOutTree, "%s_%i_%i", cTreeName, iToy, j);
			int iArray = 6 * iToy + j;
			outTree[iArray] = new TTree(cOutTree, cOutTree);
			outTree[iArray]->Branch(toyKDName_Alt[0], &D_KD1, toyKDNameF_Alt[0]);
			outTree[iArray]->Branch(toyKDName_Alt[1], &D_KD2, toyKDNameF_Alt[1]);
			outTree[iArray]->Branch(toyKDName_Alt[2], &D_bkg, toyKDNameF_Alt[2]);
			outTree[iArray]->Branch("CMS_zz4l_mass", &ZZMass, "CMS_zz4l_mass/F");
			if (weightedEvents == 1) outTree[iArray]->Branch("_weight_", &weightFit, "_weight_/D");
		};
	};


/*** SMARTER LOOP ***/

		for (int j = 0; j < 6; j++){
			sig[j]->SetBranchAddress(toyKDName[0], &D_KD1);
			if (iKD2 >= 0) sig[j]->SetBranchAddress("CMS_zz4l_KD2", &D_KD2);
			if (iKD3 >= 0) sig[j]->SetBranchAddress("CMS_zz4l_KD3", &D_KD2);
			sig[j]->SetBranchAddress("CMS_zz4l_smd", &D_bkg);
			sig[j]->SetBranchAddress("CMS_zz4l_mass", &ZZMass);
			sig[j]->SetBranchAddress("_weight_", &weightFit);
			sig[j]->SetBranchAddress("ToyNumber", &ToyNumber);
			int nentry = sig[j]->GetEntries();
			for (int i = 0; i < nentry; i++){
				sig[j]->GetEntry(i);
				if (ZZMass<140.6 && ZZMass >= 105.6 && weightFit>0 && ToyNumber < maxNToysPossible){
/*					cout << "Tree: " << 6 * ToyNumber + j << endl;
					cout << D_KD1 << endl;
					cout << D_KD2 << endl;
					cout << weightFit << endl;
*/					outTree[6 * ToyNumber + j]->Fill();
				};
			}
			/*******bkg**********/

			qq[j]->SetBranchAddress(toyKDName[0], &D_KD1);
			if (iKD2 >= 0) qq[j]->SetBranchAddress(toyKDName[1], &D_KD2);
			if (iKD3 >= 0) qq[j]->SetBranchAddress(toyKDName[3], &D_KD2);
			qq[j]->SetBranchAddress(toyKDName[2], &D_bkg);
			qq[j]->SetBranchAddress("CMS_zz4l_mass", &ZZMass);
			qq[j]->SetBranchAddress("_weight_", &weightFit);
			qq[j]->SetBranchAddress("ToyNumber", &ToyNumber);
			nentry = qq[j]->GetEntries();
			for (int i = 0; i < nentry; i++){
				qq[j]->GetEntry(i);
				if (ZZMass<140.6 && ZZMass >= 105.6 && weightFit>0 && ToyNumber < maxNToysPossible) outTree[6 * ToyNumber + j]->Fill();
			}
			gg[j]->SetBranchAddress(toyKDName[0], &D_KD1);
			if (iKD2 >= 0) gg[j]->SetBranchAddress(toyKDName[1], &D_KD2);
			if (iKD3 >= 0) gg[j]->SetBranchAddress(toyKDName[3], &D_KD2);
			gg[j]->SetBranchAddress(toyKDName[2], &D_bkg);
			gg[j]->SetBranchAddress("CMS_zz4l_mass", &ZZMass);
			gg[j]->SetBranchAddress("_weight_", &weightFit);
			gg[j]->SetBranchAddress("ToyNumber", &ToyNumber);
			nentry = gg[j]->GetEntries();
			for (int i = 0; i < nentry; i++){
				gg[j]->GetEntry(i);
				if (ZZMass<140.6 && ZZMass >= 105.6 && weightFit>0 && ToyNumber < maxNToysPossible) outTree[6 * ToyNumber + j]->Fill();
			}
			zx[j]->SetBranchAddress(toyKDName[0], &D_KD1);
			if (iKD2 >= 0) zx[j]->SetBranchAddress(toyKDName[1], &D_KD2);
			if (iKD3 >= 0) zx[j]->SetBranchAddress(toyKDName[3], &D_KD2);
			zx[j]->SetBranchAddress(toyKDName[2], &D_bkg);
			zx[j]->SetBranchAddress("CMS_zz4l_mass", &ZZMass);
			zx[j]->SetBranchAddress("_weight_", &weightFit);
			zx[j]->SetBranchAddress("ToyNumber", &ToyNumber);
			nentry = zx[j]->GetEntries();
			for (int i = 0; i < nentry; i++){
				zx[j]->GetEntry(i);
				if (ZZMass<140.6 && ZZMass >= 105.6 && weightFit>0 && ToyNumber < maxNToysPossible) outTree[6 * ToyNumber + j]->Fill();
			}
		}

	cdtof->cd();
	for (int iToy = 0; iToy<maxNToysPossible; iToy++){
		if (iToy%100==0) cout<<iToy<<endl;
		for (int j=0; j<6; j++){
			int iArray = 6 * iToy + j;
			TString cOutTreeStore = outTree[iArray]->GetName();
			outTree[iArray]->SetName(cTreeName);

			if (weightedEvents == 1){
				data[j] = new RooDataSet(Form("data%d", j), Form("data%d", j), outTree[iArray], rContainer, "", "_weight_");
				data[j]->weightError(RooAbsData::None);
			}
			else{
				data[j] = new RooDataSet(Form("data%d", j), Form("data%d", j), outTree[iArray], rContainer);
			};

			RooArgSet toyArgs(rContainer);
			toyArgs.add(cat);
			if(weightedEvents==1) toyf03[j] = new RooDataSet (Form("toy_asimov%d",j),Form("toy_asimov%d",j), toyArgs, Index(cat),WeightVar("_weight_"), Import(chan[j],*data[j]) );
			else toyf03[j] = new RooDataSet (Form("toy_asimov%d",j),Form("toy_asimov%d",j), toyArgs, Index(cat), Import(chan[j],*data[j]) );
			if(toyf03[j]->isWeighted()) cout << "Asimov dataset is weighted." << endl;

			outTree[iArray]->SetName(cOutTreeStore);
		}

		RooDataSet *toyf03_total = toyf03[0]; 
		for(int j=0;j<6;j++){
			if(j!=0) toyf03_total->append(*(toyf03[j]));
//			toyf03[j]->Print("v");
			RooCategory *cata = dynamic_cast<RooCategory *>(toyf03_total->get()->find("CMS_channel"));
			int na = cata->numBins((const char *)0);
		};
		toyf03_total->SetName(Form("toy_%d",iToy));

	    RooArgSet* set;

		RooCategory *cata = dynamic_cast<RooCategory *>(toyf03_total->get()->find("CMS_channel"));
		int na = cata->numBins((const char *)0);
		cdtof->cd();
		toyf03_total->Write(Form("toy_%d",iToy));
		toyf03_total->Print("v");

		if (iToy == 0){
			toyfile->cd();
			TString cValPlot = "Validate_";
			cValPlot.Append("ZZMass");
			TString cCanValPlot = "canvasValidate_";
			cCanValPlot.Append("ZZMass");
			TCanvas* can = new TCanvas(cCanValPlot,"",800,800);
			RooPlot* valPlot = new RooPlot(cValPlot,"",*rCMS_zz4l_mass,105.6,140.6,80);
			toyf03_total->plotOn(valPlot);
			can->cd();
			valPlot->Draw();
			toyfile->WriteTObject(can);
			delete valPlot;
			can->Close();

			cValPlot = "Validate_";
			cValPlot.Append("KD1");
			cCanValPlot = "canvasValidate_";
			cCanValPlot.Append("KD1");
			can = new TCanvas(cCanValPlot,"",800,800);
			valPlot = new RooPlot(cValPlot,"",*rCMS_zz4l_KD1,0,1,40);
			toyf03_total->plotOn(valPlot);
			can->cd();
			valPlot->Draw();
			toyfile->WriteTObject(can);
			delete valPlot;
			can->Close();

			cValPlot = "Validate_";
			cValPlot.Append("KD2");
			cCanValPlot = "canvasValidate_";
			cCanValPlot.Append("KD2");
			can = new TCanvas(cCanValPlot,"",800,800);
			valPlot = new RooPlot(cValPlot,"",*rCMS_zz4l_KD2,rCMS_zz4l_KD2->getMin(),rCMS_zz4l_KD2->getMax(),40);
			toyf03_total->plotOn(valPlot);
			can->cd();
			valPlot->Draw();
			toyfile->WriteTObject(can);
			delete valPlot;
			can->Close();

			if (weightedEvents==1){
				cValPlot = "Validate_";
				cValPlot.Append("weight");
				cCanValPlot = "canvasValidate_";
				cCanValPlot.Append("weight");
				can = new TCanvas(cCanValPlot,"",800,800);
				valPlot = new RooPlot(cValPlot,"",*rweightFit,0,2.5,50);
				toyf03_total->plotOn(valPlot);
				can->cd();
				valPlot->Draw();
				toyfile->WriteTObject(can);
				delete valPlot;
				can->Close();
			};

			cdtof->cd();
		};
	};
	for (int iToy = 0; iToy < maxNToysPossible; iToy++){
		for (int j = 0; j < 6; j++){
			int iArray = 6 * iToy + j;
			delete outTree[iArray];
		};
	};

	toyfile->Close();
}
