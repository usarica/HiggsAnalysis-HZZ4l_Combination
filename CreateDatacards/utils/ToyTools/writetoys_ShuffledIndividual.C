/**************REMINDER OF HYPOTHESES*************************/
/*
const double gi_phi2_phi4[kNumHypo][p_NumArray_MEweight_SpinZero] = { // Coefficients as in US_MEweights branch of ZZAnalusis/AnalysisStep/test/Macros/HZZ4l.h
		//	mH	GH	Re-g1	Im-g1	Re-g2	Im-g2	Re-g3	Im-g3	Re-g4	Im-g4	Re-L1	Im-L1	Re-LQ	Im-LQ	Re-a2ZG	Im-a2ZG	Re-a2GG	Im-a2GG	Re-a3ZG	Im-a3ZG	Re-a3GG	Im-a3GG	Re-ZGL1	Im-ZGL1
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	Pure	SM	0		
		{	125.6,	0.00415,	0,	0,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=1	1			
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa3=1	2			
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	0,	0,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fL1=1	3			
		{	125.6,	0.00415,	1.0,	0,	1.638,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=0.5	4			
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	2.521,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa3=0.5	5			
		{	125.6,	0.00415,	0,	0,	0.650,	0,	0,	0,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=fa3=0.5	6			
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	12046.01,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fLambda1=-0.5,	for	T3	templates	7
		{	125.6,	0.00415,	1.0,	0,	1.638,	0,	0,	0,	2.521,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=fa3=1/3	8			
		{	125.6,	0.00415,	1.0,	0,	0.546,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=0.1	9			
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0.840,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa3=0.1	10			
		{	125.6,	0.00415,	1.0,	0,	0.579,	0,	0,	0,	0.891,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=fa3=0.1	11			
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	-12046.01,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fLambda1=0.5	12			
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	-7885.965,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	flambda1=0.3	13	--> No dedicated sample		
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	-4015.337,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	flambda1=0.1	14			
		{	125.6,	0.00415,	1.0,	0,	0,	1.638,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=0.5,	phia2=90	15		
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	2.521,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa3=0.5,	phia3=90	16		
		{	125.6,	0.00415,	0,	0,	0.650,	0,	0,	0,	0,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=fa3=0.5,	phia3=90	17		
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	-12046.01,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fLambda1=0.5,	phiL1=90	18	--> No dedicated sample					
		{	125.6,	0.00415,	1.0,	0,	1.638,	0,	0,	0,	0,	2.521,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=fa3=1/3,	phia3=90	19		
		{	125.6,	0.00415,	1.0,	0,	0,	0.546,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=0.1,	phia2=90	20		
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0.840,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa3=0.1,	phia3=90	21		
		{	125.6,	0.00415,	1.0,	0,	0.579,	0,	0,	0,	0,	0.891,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=0.1,	fa3=0.1,	phia3=90	22	
		{	125.6,	0.00415,	1.0,	0,	-1.638,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=-0.5	23			
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	-2.521,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa3=-0.5	24			
		{	125.6,	0.00415,	0,	0,	-0.650,	0,	0,	0,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=-0.5,	fa3=0.5	25		

		{	125.6,	0.00415,	0,	0,	1.638,	0,	0,	0,	0,	0,	12046.01,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=0.5,	fLambda1=-0.5,	for	templates	26
		{	125.6,	0.00415,	0,	0,	0,	-1.638,	0,	0,	0,	0,	12046.01,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=-0.5i,	fLambda1=-0.5,	for	templates	27
		{	125.6,	0.00415,	0,	0,	1.638,	0,	0,	0,	0,	0,	-12046.01,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=0.5,	fLambda1=0.5	28		
		{	125.6,	0.00415,	1.0,	0,	-1.638,	0,	0,	0,	0,	0,	-12046.01,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa2=-0.33,	fLambda1=0.33	29		
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	2.521,	0,	12046.01,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa3=0.5,	fLambda1=-0.5,	for	templates	30
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	0,	-2.521,	12046.01,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa3=-0.5i,	fLambda1=-0.5,	for	templates	31
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	2.521,	0,	-12046.01,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa3=0.5,	fLambda1=0.5	32		
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	-2.521,	0,	-12046.01,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fa3=-0.33,	fLambda1=0.33	33		
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-6338.9995537,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fLambdaQ=1	34			
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-6338.9995537,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fLambdaQ=1,	phiLQ=90	35			
		{	125.6,	0.00415,	0,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fg1=1,	phig1=90	36			
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-6338.9995537,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fg1=fLambdaQ=0.5,	phiLQ=0		37		
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-6338.9995537,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fg1=fLambdaQ=0.5,	phiLQ=90	38			
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	Pure	PM-ZG	39		
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1.0,	0,	0,	0,	0,	0,	0,	0	},	//	Pure	PM-GG	40		
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0.0473,	0,	0,	0,	0,	0,	0,	0,	0,	0	},	//	fPM-ZG=0.5	41			
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-0.0531,	0,	0,	0,	0,	0,	0,	0	},	//	fPM-GG=0.5	42			
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1.0,	0,	0,	0,	0,	0	},	//	Pure	M-ZG	43		
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1.0,	0,	0,	0	},	//	Pure	M-GG	44		
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0.052161,	0,	0,	0,	0,	0	},	//	fM-ZG=0.5	45			
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-0.053666,	0,	0,	0	},	//	fM-GG=0.5	46			
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0.0473,	0,	0,	0,	0.052161,	0,	0,	0,	0,	0	},	//	fZG=fM-ZG=0.5	47			
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-0.0531,	0,	0,	0,	-0.053666,	0,	0,	0	},	//	fGG=fM-GG=0.5	48			
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0.00175,	0,	-0.002,	0,	0,	0,	0,	0,	0,	0	},	//	SM	ZZ,	ZG.	GG	49
		{	125.6,	0.00415,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1.0,	0	},	//	Pure	ZG-L1	50		
		{	125.6,	0.00415,	1.0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-7591.914,	0	}	//	fZG-L1=0.5	51			
};

********** NOTES **********
*** Be careful about this ordering in your datacards below; tree orderings are

	TString chan[6]={"ch1","ch2","ch3","ch4","ch5","ch6"};
	int order[6]={2,1,0,4,3,5};
	TString channameO[6]={"7TeV/4e","7TeV/4mu","7TeV/2mu2e","8TeV/4mu","8TeV/2mu2e","8TeV/4e"};
	double normO[6]=	{	0.6951,	1.2439,	1.6662,	5.9471,	7.6807,	3.0898	}; 
	double normzxO[6]=	{	0.6226,	0.2230,	1.0628,	1.1878,	4.2929,	2.7676	};
	double normqqO[6]=	{	0.8386,	1.7971,	2.2456,	7.6478,	8.8585,	2.9364	};
	double normggO[6]=	{	0.0341,	0.0625,	0.0741,	0.4131,	0.5005,	0.2041	};

*** Change ordering via

	TString chan[nCat_Tree]={"ch1","ch2","ch3","ch4","ch5","ch6"}; // Input categories: Should match input file
	TString chan_Write[nCat_Write]={"ch1","ch2","ch3","ch4","ch5","ch6"}; // Written categories: Whatever names and count you wish



*** Variables are recorded in trees as
	TString KDlist_Trees[N_KDs] = {
		"CMS_zz4l_smd",
		"CMS_zz4l_D_g1Q2_phi0",
		"CMS_zz4l_D_g1_vs_g2_phi0",
		"CMS_zz4l_D_g1_vs_g4_phi0",
		"CMS_zz4l_D_ZG",
		"CMS_zz4l_D_GG",
		"CMS_zz4l_D_ZG_PS",
		"CMS_zz4l_D_GG_PS",
		"CMS_zz4l_D_ZG_L1",
		"CMS_zz4l_D_g1Q2int_phi0",
		"CMS_zz4l_D_g2int_phi0",
		"CMS_zz4l_D_g4int_phi0",
		"CMS_zz4l_D_g2int_phi90",
		"CMS_zz4l_D_g4int_phi90",
		"CMS_zz4l_D_ZGint",
		"CMS_zz4l_D_GGint",
		"CMS_zz4l_D_ZG_PSint",
		"CMS_zz4l_D_GG_PSint",
		"CMS_zz4l_D_ZG_L1int_phi0",
		"CMS_zz4l_D_ZG_L1int_phi90"
	};
	TString Kinematicslist_Trees[N_Kinematics] = {
		"CMS_zz4l_mass",
		"CMS_zz4l_Z1Mass",
		"CMS_zz4l_Z2Mass",
		"CMS_zz4l_costhetastar",
		"CMS_zz4l_helcosthetaZ1",
		"CMS_zz4l_helcosthetaZ2",
		"CMS_zz4l_helphi",
		"CMS_zz4l_phistarZ1"
	};

*** You may change the following as you wish. Please change both sets consistently:

	const int nMaxSupportedKDs=3;
	TString KDlist_Write[nMaxSupportedKDs] = { // CHANGE THESE TO MATCH DATACARDS
		"CMS_zz4l_smd",
		"CMS_zz4l_KD1",
		"CMS_zz4l_KD2"
	};
	TString Kinematicslist_Write[N_Kinematics] = { // CHANGE THESE TO MATCH DATACARDS
		"CMS_zz4l_mass",
		"CMS_zz4l_mZ1",
		"CMS_zz4l_mZ2",
		"CMS_zz4l_hs",
		"CMS_zz4l_h1",
		"CMS_zz4l_h2",
		"CMS_zz4l_phi",
		"CMS_zz4l_phi1"
	};
	TString KDlist_WriteD[nMaxSupportedKDs] = { // CHANGE THESE TO MATCH DATACARDS
		"CMS_zz4l_smd/D",
		"CMS_zz4l_KD1/D",
		"CMS_zz4l_KD2/D"
	};
	TString Kinematicslist_WriteD[N_Kinematics] = { // CHANGE THESE TO MATCH DATACARDS
		"CMS_zz4l_mass/D",
		"CMS_zz4l_mZ1/D",
		"CMS_zz4l_mZ2/D",
		"CMS_zz4l_hs/D",
		"CMS_zz4l_h1/D",
		"CMS_zz4l_h2/D",
		"CMS_zz4l_phi/D",
		"CMS_zz4l_phi1/D"
	};



*/
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include "RooFitResult.h"
#include "TMath.h"
#include "TList.h"
#include "TChain.h"
#include "TString.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TRandom3.h"


using namespace RooFit;
using namespace std;


void writetoys_ShuffledIndividual(int BSMSample=0, int iKD1=1, int iKD2=2, int maxNToysPossible=1550, int weightedEvents=1, int useAsimovV2=0){ // <<<--- Choose filename
	const int kNumHypo=52;
	const int kNumFiles=24;
	const int gapZZHypo=13;
	const int gapZZHypo_2=18;
	const int firstNonZZHypo=39;
	double ZZMassRange[2] = {	105.6	,	140.6	}; // SHOULD CHANGE FOR 8D METHOD

	bool useKinematics=false; // SET TRUE WITH NEGATIVE iKD1 AND iKD2

	const int N_Kinematics=8;
	const int N_KDs=20;
	const int firstIntKD=9;

/***** RUN CONDITIONS *****/
	if(!(weightedEvents==0 || weightedEvents==1)) assert(0);
	if(useAsimovV2<0) assert(0);
	if (iKD1<=0){
		if (iKD2>=0) assert(0);
		else useKinematics = true;
	}
	else if(iKD1>=N_KDs) assert(0);
	else if(iKD2==iKD1) iKD2=-iKD2; // Turns off KD2 if identical to KD1
	if(BSMSample>=kNumHypo) assert(0);
/*** END RUN CONDITIONS ***/
	char cDir[]="./AsimovToys";
	char cinput[100];
	char coutput[100];

	TString KDlist_Trees[N_KDs] = {
		"CMS_zz4l_smd",
		"CMS_zz4l_D_g1Q2_phi0",
		"CMS_zz4l_D_g1_vs_g2_phi0",
		"CMS_zz4l_D_g1_vs_g4_phi0",
		"CMS_zz4l_D_ZG",
		"CMS_zz4l_D_GG",
		"CMS_zz4l_D_ZG_PS",
		"CMS_zz4l_D_GG_PS",
		"CMS_zz4l_D_ZG_L1",
		"CMS_zz4l_D_g1Q2int_phi0",
		"CMS_zz4l_D_g2int_phi0",
		"CMS_zz4l_D_g4int_phi0",
		"CMS_zz4l_D_g2int_phi90",
		"CMS_zz4l_D_g4int_phi90",
		"CMS_zz4l_D_ZGint",
		"CMS_zz4l_D_GGint",
		"CMS_zz4l_D_ZG_PSint",
		"CMS_zz4l_D_GG_PSint",
		"CMS_zz4l_D_ZG_L1int_phi0",
		"CMS_zz4l_D_ZG_L1int_phi90"
	};
	TString Kinematicslist_Trees[N_Kinematics] = {
		"CMS_zz4l_mass",
		"CMS_zz4l_Z1Mass",
		"CMS_zz4l_Z2Mass",
		"CMS_zz4l_costhetastar",
		"CMS_zz4l_helcosthetaZ1",
		"CMS_zz4l_helcosthetaZ2",
		"CMS_zz4l_helphi",
		"CMS_zz4l_phistarZ1"
	};

	const int nMaxSupportedKDs=3;
	TString KDlist_Write[nMaxSupportedKDs] = { // CHANGE THESE TO MATCH DATACARDS
		"CMS_zz4l_smd",
		"CMS_zz4l_KD1",
		"CMS_zz4l_KD2"
	};
	TString Kinematicslist_Write[N_Kinematics] = { // CHANGE THESE TO MATCH DATACARDS
		"CMS_zz4l_mass",
		"CMS_zz4l_mZ1",
		"CMS_zz4l_mZ2",
		"CMS_zz4l_hs",
		"CMS_zz4l_h1",
		"CMS_zz4l_h2",
		"CMS_zz4l_phi",
		"CMS_zz4l_phi1"
	};
	TString KDlist_WriteD[nMaxSupportedKDs] = { // CHANGE THESE TO MATCH DATACARDS
		"CMS_zz4l_smd/D",
		"CMS_zz4l_KD1/D",
		"CMS_zz4l_KD2/D"
	};
	TString Kinematicslist_WriteD[N_Kinematics] = { // CHANGE THESE TO MATCH DATACARDS
		"CMS_zz4l_mass/D",
		"CMS_zz4l_mZ1/D",
		"CMS_zz4l_mZ2/D",
		"CMS_zz4l_hs/D",
		"CMS_zz4l_h1/D",
		"CMS_zz4l_h2/D",
		"CMS_zz4l_phi/D",
		"CMS_zz4l_phi1/D"
	};

	const int nCat_Tree=6; // Input categories: Should match input file
	const int nCat_Write=6; // Written categories: Whatever names and count you wish
	TString chan[nCat_Tree]={"ch1","ch2","ch3","ch4","ch5","ch6"}; // Input categories: Should match input file
	TString chan_Write[nCat_Write]={"ch1","ch2","ch3","ch4","ch5","ch6"}; // Written categories: Whatever names and count you wish
	bool split2e2mu = (nCat_Tree==(nCat_Write-2) && useKinematics); // Split 2e2mu when using kinematics? Add extra 2 channels to chan_Write for chan::ch1 and ch4 to fill with 2e2mu

/*** INTERMEDIATE DATASETS TO DO THE DIRTY WORK ****/
	RooDataSet *data[nCat_Write];
	RooDataSet *toyf03[nCat_Write];
	RooDataSet *data_asimov[nCat_Write];
	RooDataSet *toyf03_asimov[nCat_Write];

/*** THE BOSS ***/
	RooCategory cat("CMS_channel","CMS_channel");
	for(int j=0;j<nCat_Write;j++){ // Define written categories
		cat.defineType(chan_Write[j],j);
		cat.setLabel(chan_Write[j]);
	};

	cout << "Starting..." << endl;

	if(useAsimovV2==0) sprintf(cinput,"%s/Asimov_v1_Hypo%i.root",cDir,BSMSample);
	else sprintf(cinput,"%s/Asimov_v2_Hypo%i.root",cDir,BSMSample);
	cout << "Input file: " << cinput << endl;

	if (weightedEvents == 1){
		if (useKinematics){
			if (useAsimovV2 == 0) sprintf(coutput, "%s/WeightedSynchToys_v1_Hypo%i_Kinematics.root", cDir, BSMSample);
			else sprintf(coutput, "%s/WeightedSynchToys_v2_Hypo%i_Kinematics.root", cDir, BSMSample);
		}
		else if (iKD2 >= 0){
			if (useAsimovV2 == 0) sprintf(coutput, "%s/WeightedSynchToys_v1_Hypo%i_KD1_%i_vs_KD2_%i.root", cDir, BSMSample, iKD1, iKD2);
			else sprintf(coutput, "%s/WeightedSynchToys_v2_Hypo%i_KD1_%i_vs_KD2_%i.root", cDir, BSMSample, iKD1, iKD2);
		}
		else{
			if (useAsimovV2 == 0) sprintf(coutput, "%s/WeightedSynchToys_v1_Hypo%i_KD1_%i.root", cDir, BSMSample, iKD1);
			else sprintf(coutput, "%s/WeightedSynchToys_v2_Hypo%i_KD1_%i.root", cDir, BSMSample, iKD1);
		};
	}
	else{
		if(useKinematics){
			if(useAsimovV2==0) sprintf(coutput,"%s/UnweightedSynchToys_v1_Hypo%i_Kinematics.root",cDir,BSMSample);
			else sprintf(coutput,"%s/UnweightedSynchToys_v2_Hypo%i_Kinematics.root",cDir,BSMSample);
		}
		else if(iKD2>=0){
			if(useAsimovV2==0) sprintf(coutput,"%s/UnweightedSynchToys_v1_Hypo%i_KD1_%i_vs_KD2_%i.root",cDir,BSMSample,iKD1,iKD2);
			else sprintf(coutput,"%s/UnweightedSynchToys_v2_Hypo%i_KD1_%i_vs_KD2_%i.root",cDir,BSMSample,iKD1,iKD2);
		}
		else{
			if(useAsimovV2==0) sprintf(coutput,"%s/UnweightedSynchToys_v1_Hypo%i_KD1_%i.root",cDir,BSMSample,iKD1);
			else sprintf(coutput,"%s/UnweightedSynchToys_v2_Hypo%i_KD1_%i.root",cDir,BSMSample,iKD1);
		};
	};
	cout << "Output file: " << coutput << endl;

	TFile *toyfile = new TFile(coutput,"recreate"); // This is where the asimov will reside
	TDirectory *cdtof = toyfile->mkdir("toys"); // This is where the toys will reside

	toyfile->cd();
	TChain *asimov[nCat_Tree];
	TChain *sig[nCat_Tree];
	TChain *gg[nCat_Tree];
	TChain *qq[nCat_Tree];
	TChain *zx[nCat_Tree];
	for (int j=0; j<nCat_Tree; j++){ // Grab the trees
		TString cAsimov=chan[j];
		TString ctoys_sig="ToyEvents_Sig_";
		ctoys_sig.Append(chan[j]);
		TString ctoys_ggzz="ToyEvents_ggZZ_";
		ctoys_ggzz.Append(chan[j]);
		TString ctoys_qqzz="ToyEvents_qqZZ_";
		ctoys_qqzz.Append(chan[j]);
		TString ctoys_zx="ToyEvents_ZX_";
		ctoys_zx.Append(chan[j]);

		asimov[j] = new TChain(cAsimov);
		sig[j] = new TChain(ctoys_sig);
		gg[j] = new TChain(ctoys_ggzz);
		qq[j] = new TChain(ctoys_qqzz);
		zx[j] = new TChain(ctoys_zx);

		asimov[j]->Add(cinput);
		sig[j]->Add(cinput);
		gg[j]->Add(cinput);
		qq[j]->Add(cinput);
		zx[j]->Add(cinput);
	}

	RooRealVar* rweightFit = new RooRealVar("_weight_","_weight_",0.,1.e20);
	RooRealVar* rCMS_zz4l_mass = new RooRealVar(Kinematicslist_Write[0],Kinematicslist_Write[0],100.,1000.); // Always write ZZMass
	RooRealVar* rCMS_zz4l_VarContainer[nMaxSupportedKDs+N_Kinematics-1] = { 0 }; // (-1) since mass is always added separately above, notice that no Z1 or Z2 id is added here
	RooArgSet rContainer(*rCMS_zz4l_mass,*rweightFit);
	if (!useKinematics){ // Add 2 or 3 KDs
		rCMS_zz4l_VarContainer[0] = new RooRealVar(KDlist_Write[0],KDlist_Write[0], 0., 1.);
		rContainer.add(*rCMS_zz4l_VarContainer[0]); // Always add superMELA

		if(iKD1 < firstIntKD) rCMS_zz4l_VarContainer[1] = new RooRealVar(KDlist_Write[1],KDlist_Write[1], 0., 1.);
		else rCMS_zz4l_VarContainer[1] = new RooRealVar(KDlist_Write[1],KDlist_Write[1], -1., 1.);
		rContainer.add(*rCMS_zz4l_VarContainer[1]);

		if(iKD2 < firstIntKD && iKD2>0) rCMS_zz4l_VarContainer[2] = new RooRealVar(KDlist_Write[2],KDlist_Write[2], 0., 1.);
		else if(iKD2>0) rCMS_zz4l_VarContainer[2] = new RooRealVar(KDlist_Write[2],KDlist_Write[2], -1., 1.);
		rContainer.add(*rCMS_zz4l_VarContainer[2]);
	}
	else{ // Add the kinematics
		rCMS_zz4l_VarContainer[0] = new RooRealVar(Kinematicslist_Write[1],Kinematicslist_Write[1], 40., 120.); // Z1Mass
		rCMS_zz4l_VarContainer[1] = new RooRealVar(Kinematicslist_Write[2],Kinematicslist_Write[2], 12., 120.); // Z2Mass
		rCMS_zz4l_VarContainer[2] = new RooRealVar(Kinematicslist_Write[3],Kinematicslist_Write[3], -1., 1.); // costhetastar
		rCMS_zz4l_VarContainer[3] = new RooRealVar(Kinematicslist_Write[4],Kinematicslist_Write[4], -1., 1.); // costheta1
		rCMS_zz4l_VarContainer[4] = new RooRealVar(Kinematicslist_Write[5],Kinematicslist_Write[5], -1., 1.); // costheta2
		rCMS_zz4l_VarContainer[5] = new RooRealVar(Kinematicslist_Write[6],Kinematicslist_Write[6], -TMath::Pi(), TMath::Pi()); // phi
		rCMS_zz4l_VarContainer[6] = new RooRealVar(Kinematicslist_Write[7],Kinematicslist_Write[7], -TMath::Pi(), TMath::Pi()); // phi1

		for(int k=1;k<N_Kinematics;k++) rContainer.add(*rCMS_zz4l_VarContainer[k-1]);
	};

	double KD_pass[N_KDs];
	double Kinematics_pass[N_Kinematics];
	int ToyNumber=0;
	double weightFit = 0, weightFit2 = 0;
	int Z1ids,Z2ids;

/*** CHECK TO SEE IF THERE ARE AS MANY TOYS AS REQUESTED ***/

	for (int j = 0; j < nCat_Tree; j++){
		for (int k = 0; k < N_KDs; k++) asimov[j]->SetBranchAddress(KDlist_Trees[k], (KD_pass+k));
		for (int k = 0; k < N_Kinematics; k++) asimov[j]->SetBranchAddress(Kinematicslist_Trees[k], (Kinematics_pass+k));
		asimov[j]->SetBranchAddress("_weight_", &weightFit);
		asimov[j]->SetBranchAddress("Z1ids", &Z1ids);
		asimov[j]->SetBranchAddress("Z2ids", &Z2ids);


		for (int k = 0; k < N_KDs; k++) sig[j]->SetBranchAddress(KDlist_Trees[k], (KD_pass+k));
		for (int k = 0; k < N_Kinematics; k++) sig[j]->SetBranchAddress(Kinematicslist_Trees[k], (Kinematics_pass+k));
		sig[j]->SetBranchAddress("_weight_", &weightFit2);
		sig[j]->SetBranchAddress("ToyNumber", &ToyNumber);
		sig[j]->SetBranchAddress("Z1ids", &Z1ids);
		sig[j]->SetBranchAddress("Z2ids", &Z2ids);

		for (int k = 0; k < N_KDs; k++) qq[j]->SetBranchAddress(KDlist_Trees[k], (KD_pass+k));
		for (int k = 0; k < N_Kinematics; k++) qq[j]->SetBranchAddress(Kinematicslist_Trees[k], (Kinematics_pass+k));
		qq[j]->SetBranchAddress("_weight_", &weightFit2);
		qq[j]->SetBranchAddress("ToyNumber", &ToyNumber);
		qq[j]->SetBranchAddress("Z1ids", &Z1ids);
		qq[j]->SetBranchAddress("Z2ids", &Z2ids);

		for (int k = 0; k < N_KDs; k++) gg[j]->SetBranchAddress(KDlist_Trees[k], (KD_pass+k));
		for (int k = 0; k < N_Kinematics; k++) gg[j]->SetBranchAddress(Kinematicslist_Trees[k], (Kinematics_pass+k));
		gg[j]->SetBranchAddress("_weight_", &weightFit2);
		gg[j]->SetBranchAddress("ToyNumber", &ToyNumber);
		gg[j]->SetBranchAddress("Z1ids", &Z1ids);
		gg[j]->SetBranchAddress("Z2ids", &Z2ids);

		for (int k = 0; k < N_KDs; k++) zx[j]->SetBranchAddress(KDlist_Trees[k], (KD_pass+k));
		for (int k = 0; k < N_Kinematics; k++) zx[j]->SetBranchAddress(Kinematicslist_Trees[k], (Kinematics_pass+k));
		zx[j]->SetBranchAddress("_weight_", &weightFit2);
		zx[j]->SetBranchAddress("ToyNumber", &ToyNumber);
		zx[j]->SetBranchAddress("Z1ids", &Z1ids);
		zx[j]->SetBranchAddress("Z2ids", &Z2ids);

		sig[j]->GetEntry(sig[j]->GetEntries()-1);
		cout << "Signal " << j << ": " << sig[j]->GetEntries() << " entries" << endl;
		if(ToyNumber<maxNToysPossible-1) maxNToysPossible = ToyNumber+1;
		qq[j]->GetEntry(qq[j]->GetEntries()-1);
		cout << "qqZZ " << j << ": " << qq[j]->GetEntries() << " entries" << endl;
		if(ToyNumber<maxNToysPossible-1) maxNToysPossible = ToyNumber+1;
		gg[j]->GetEntry(gg[j]->GetEntries()-1);
		cout << "ggZZ " << j << ": " << gg[j]->GetEntries() << " entries" << endl;
		if(ToyNumber<maxNToysPossible-1) maxNToysPossible = ToyNumber+1;
		zx[j]->GetEntry(zx[j]->GetEntries()-1);
		cout << "Z+X " << j << ": " << zx[j]->GetEntries() << " entries" << endl;
		if(ToyNumber<maxNToysPossible-1) maxNToysPossible = ToyNumber+1;

		cout << "Asimov " << j << ": " << asimov[j]->GetEntries() << " entries" << endl;
	};
	cout << "Maximum number of toys possible: " << maxNToysPossible << endl;

	TTree* outAsimovTree[nCat_Write];
	TTree** outTree = new TTree* [nCat_Write*maxNToysPossible];
	char cTreeName[]="SelectedTree";
	char cAsimovTreeName[]="SelectedAsimovTree";
	for (int j = 0; j < nCat_Write; j++){
		cdtof->cd();
		for (int iToy = 0; iToy < maxNToysPossible; iToy++){
			char cOutTree[50];
			sprintf(cOutTree, "%s_%i_%i", cTreeName, iToy, j);
			int iArray = nCat_Write * iToy + j;
			outTree[iArray] = new TTree(cOutTree, cOutTree);
			outTree[iArray]->Branch(KDlist_Write[0], (KD_pass+0), KDlist_WriteD[0]);
			if(iKD1>0) outTree[iArray]->Branch(KDlist_Write[1], (KD_pass+iKD1), KDlist_WriteD[1]);
			if(iKD2>0) outTree[iArray]->Branch(KDlist_Write[2], (KD_pass+iKD2), KDlist_WriteD[2]);
			for(int k=0;k<N_Kinematics;k++) outTree[iArray]->Branch(Kinematicslist_Write[k], (Kinematics_pass+k), Kinematicslist_WriteD[k]);
			outTree[iArray]->Branch("_weight_", &weightFit2, "_weight_/D");
//			outTree[iArray]->Branch("Z1ids", &Z1ids);
//			outTree[iArray]->Branch("Z2ids", &Z2ids);
		};
		toyfile->cd();
		char cOutAsimovTree[50];
		sprintf(cOutAsimovTree, "%s_%i", cAsimovTreeName, j);
		outAsimovTree[j] = new TTree(cOutAsimovTree, cOutAsimovTree);
		outAsimovTree[j]->Branch(KDlist_Write[0], (KD_pass+0), KDlist_WriteD[0]);
		if(iKD1>0) outAsimovTree[j]->Branch(KDlist_Write[1], (KD_pass+iKD1), KDlist_WriteD[1]);
		if(iKD2>0) outAsimovTree[j]->Branch(KDlist_Write[2], (KD_pass+iKD2), KDlist_WriteD[2]);
		for(int k=0;k<N_Kinematics;k++) outAsimovTree[j]->Branch(Kinematicslist_Write[k], (Kinematics_pass+k), Kinematicslist_WriteD[k]);
		outAsimovTree[j]->Branch("_weight_", &weightFit, "_weight_/D");
//		outAsimovTree[j]->Branch("Z1ids", &Z1ids);
//		outAsimovTree[j]->Branch("Z2ids", &Z2ids);
	};


/*** SMARTER LOOP ***/

		for (int j = 0; j < nCat_Tree; j++){
			int nentry = asimov[j]->GetEntries();
			for (int i = 0; i < nentry; i++){
				asimov[j]->GetEntry(i);
				if (Kinematics_pass[0]<ZZMassRange[1] && Kinematics_pass[0] >= ZZMassRange[0] && weightFit>0){
					if(!(split2e2mu && (j==0 || j==3) ) ) outAsimovTree[j]->Fill();
					else{
						if(Z1ids==-169) outAsimovTree[j]->Fill();
						else outAsimovTree[nCat_Tree + (j % 2) ]->Fill();
					};
				}
			}


			nentry = sig[j]->GetEntries();
			for (int i = 0; i < nentry; i++){
				sig[j]->GetEntry(i);
				if (Kinematics_pass[0]<ZZMassRange[1] && Kinematics_pass[0] >= ZZMassRange[0] && weightFit2>0 && ToyNumber < maxNToysPossible){
					if(weightedEvents==0) weightFit2=1;
					if(!(split2e2mu && (j==0 || j==3) ) ) outTree[ nCat_Write * ToyNumber + j ]->Fill();
					else{
						if(Z1ids==-169) outTree[ nCat_Write * ToyNumber + j ]->Fill();
						else outTree[ nCat_Write * ToyNumber + nCat_Tree + (j % 2) ]->Fill();
					};
				};
			}
			/*******bkg**********/

			nentry = qq[j]->GetEntries();
			for (int i = 0; i < nentry; i++){
				qq[j]->GetEntry(i);
				if (Kinematics_pass[0]<ZZMassRange[1] && Kinematics_pass[0] >= ZZMassRange[0] && weightFit2>0 && ToyNumber < maxNToysPossible){
					if(weightedEvents==0) weightFit2=1;
					if(!(split2e2mu && (j==0 || j==3) ) ) outTree[ nCat_Write * ToyNumber + j ]->Fill();
					else{
						if(Z1ids==-169) outTree[ nCat_Write * ToyNumber + j ]->Fill();
						else outTree[ nCat_Write * ToyNumber + nCat_Tree + (j % 2) ]->Fill();
					};
				};
			}

			nentry = gg[j]->GetEntries();
			for (int i = 0; i < nentry; i++){
				gg[j]->GetEntry(i);
				if (Kinematics_pass[0]<ZZMassRange[1] && Kinematics_pass[0] >= ZZMassRange[0] && weightFit2>0 && ToyNumber < maxNToysPossible){
					if(weightedEvents==0) weightFit2=1;
					if(!(split2e2mu && (j==0 || j==3) ) ) outTree[ nCat_Write * ToyNumber + j ]->Fill();
					else{
						if(Z1ids==-169) outTree[ nCat_Write * ToyNumber + j ]->Fill();
						else outTree[ nCat_Write * ToyNumber + nCat_Tree + (j % 2) ]->Fill();
					};
				};
			}

			nentry = zx[j]->GetEntries();
			for (int i = 0; i < nentry; i++){
				zx[j]->GetEntry(i);
				if (Kinematics_pass[0]<ZZMassRange[1] && Kinematics_pass[0] >= ZZMassRange[0] && weightFit2>0 && ToyNumber < maxNToysPossible){
					if(weightedEvents==0) weightFit2=1;
					if(!(split2e2mu && (j==0 || j==3) ) ) outTree[ nCat_Write * ToyNumber + j ]->Fill();
					else{
						if(Z1ids==-169) outTree[ nCat_Write * ToyNumber + j ]->Fill();
						else outTree[ nCat_Write * ToyNumber + nCat_Tree + (j % 2) ]->Fill();
					};
				};
			}
		}

	toyfile->cd();
	RooDataSet* toyf03_total_asimov; 
	for (int j=0; j<nCat_Write; j++){
		TString cOutTreeStore = outAsimovTree[j]->GetName();
		outAsimovTree[j]->SetName(cTreeName);
		cout << outAsimovTree[j]->GetEntries() << endl;
		data_asimov[j] = new RooDataSet(Form("data%d", j), Form("data%d", j), outAsimovTree[j], rContainer, "", "_weight_");
		data_asimov[j]->weightError(RooAbsData::None);

		RooArgSet toyArgs(rContainer);
		toyArgs.add(cat);
		toyf03_asimov[j] = new RooDataSet (Form("toy_asimov%d",j),Form("toy_asimov%d",j), toyArgs, Index(cat),WeightVar("_weight_"), Import(chan[j],*data_asimov[j]) );
		if(toyf03_asimov[j]->isWeighted()) cout << "Asimov dataset is weighted." << endl;
		outAsimovTree[j]->SetName(cOutTreeStore);

		toyf03_asimov[j]->Print("v");
		if(j==0) toyf03_total_asimov = toyf03_asimov[j];
		else toyf03_total_asimov->append(*(toyf03_asimov[j]));
	};
	toyf03_total_asimov->SetName("toy_asimov");

	RooCategory* cata = dynamic_cast<RooCategory *>(toyf03_total_asimov->get()->find("CMS_channel"));
	int na = cata->numBins((const char *)0);
	toyfile->cd();
	toyf03_total_asimov->Write("toy_asimov");
	toyf03_total_asimov->Print("v");

	TString cValPlot = "ValidateAsimov_";
	cValPlot.Append("ZZMass");
	TString cCanValPlot = "canvasValidateAsimov_";
	cCanValPlot.Append("ZZMass");
	TCanvas* can = new TCanvas(cCanValPlot,"",800,800);
	RooPlot* valPlot = new RooPlot(cValPlot,"",*rCMS_zz4l_mass,ZZMassRange[0],ZZMassRange[1],80);
	toyf03_total_asimov->plotOn(valPlot);
	can->cd();
	valPlot->Draw();
	toyfile->WriteTObject(can);
	delete valPlot;
	can->Close();

	if (!useKinematics){
		cValPlot = "ValidateAsimov_";
		cValPlot.Append(rCMS_zz4l_VarContainer[1]->GetName());
		cCanValPlot = "canvasValidateAsimov_";
		cCanValPlot.Append(rCMS_zz4l_VarContainer[1]->GetName());
		can = new TCanvas(cCanValPlot, "", 800, 800);
		valPlot = new RooPlot(cValPlot, "", *rCMS_zz4l_VarContainer[1], rCMS_zz4l_VarContainer[1]->getMin(), rCMS_zz4l_VarContainer[1]->getMax(), 40);
		toyf03_total_asimov->plotOn(valPlot);
		can->cd();
		valPlot->Draw();
		toyfile->WriteTObject(can);
		delete valPlot;
		can->Close();

		cValPlot = "ValidateAsimov_";
		cValPlot.Append(rCMS_zz4l_VarContainer[2]->GetName());
		cCanValPlot = "canvasValidateAsimov_";
		cCanValPlot.Append(rCMS_zz4l_VarContainer[2]->GetName());
		can = new TCanvas(cCanValPlot, "", 800, 800);
		valPlot = new RooPlot(cValPlot, "", *rCMS_zz4l_VarContainer[2], rCMS_zz4l_VarContainer[2]->getMin(), rCMS_zz4l_VarContainer[2]->getMax(), 40);
		toyf03_total_asimov->plotOn(valPlot);
		can->cd();
		valPlot->Draw();
		toyfile->WriteTObject(can);
		delete valPlot;
		can->Close();
	};	

	if(weightedEvents==0) rContainer.remove(*rweightFit);
	cdtof->cd();
	for (int iToy = 0; iToy<maxNToysPossible; iToy++){
		for (int j=0; j<nCat_Write; j++){
			int iArray = nCat_Write * iToy + j;
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
//			if(toyf03[j]->isWeighted()) cout << "Toy dataset is weighted." << endl;

			outTree[iArray]->SetName(cOutTreeStore);
		}

		RooDataSet *toyf03_total = toyf03[0]; 
		for(int j=0;j<nCat_Write;j++){
			if(j!=0) toyf03_total->append(*(toyf03[j]));
			RooCategory *cata = dynamic_cast<RooCategory *>(toyf03_total->get()->find("CMS_channel"));
			int na = cata->numBins((const char *)0);
		};
		toyf03_total->SetName(Form("toy_%d",iToy));

		cata = dynamic_cast<RooCategory *>(toyf03_total->get()->find("CMS_channel"));
		na = cata->numBins((const char *)0);
		cdtof->cd();
		toyf03_total->Write(Form("toy_%d",iToy));
//		toyf03_total->Print("v");

		if (iToy == 0){
			toyfile->cd();
			cValPlot = "ValidateToy0_";
			cValPlot.Append("ZZMass");
			cCanValPlot = "canvasValidateToy0_";
			cCanValPlot.Append("ZZMass");
			can = new TCanvas(cCanValPlot,"",800,800);
			valPlot = new RooPlot(cValPlot,"",*rCMS_zz4l_mass,ZZMassRange[0],ZZMassRange[1],80);
			toyf03_total->plotOn(valPlot);
			can->cd();
			valPlot->Draw();
			toyfile->WriteTObject(can);
			delete valPlot;
			can->Close();

			if (!useKinematics){
				cValPlot = "ValidateToy0_";
				cValPlot.Append(rCMS_zz4l_VarContainer[1]->GetName());
				cCanValPlot = "canvasValidateToy0_";
				cCanValPlot.Append(rCMS_zz4l_VarContainer[1]->GetName());
				can = new TCanvas(cCanValPlot, "", 800, 800);
				valPlot = new RooPlot(cValPlot, "", *rCMS_zz4l_VarContainer[1], rCMS_zz4l_VarContainer[1]->getMin(), rCMS_zz4l_VarContainer[1]->getMax(), 40);
				toyf03_total->plotOn(valPlot);
				can->cd();
				valPlot->Draw();
				toyfile->WriteTObject(can);
				delete valPlot;
				can->Close();

				cValPlot = "ValidateToy0_";
				cValPlot.Append(rCMS_zz4l_VarContainer[2]->GetName());
				cCanValPlot = "canvasValidateToy0_";
				cCanValPlot.Append(rCMS_zz4l_VarContainer[2]->GetName());
				can = new TCanvas(cCanValPlot, "", 800, 800);
				valPlot = new RooPlot(cValPlot, "", *rCMS_zz4l_VarContainer[2], rCMS_zz4l_VarContainer[2]->getMin(), rCMS_zz4l_VarContainer[2]->getMax(), 40);
				toyf03_total->plotOn(valPlot);
				can->cd();
				valPlot->Draw();
				toyfile->WriteTObject(can);
				delete valPlot;
				can->Close();
			};

			if (weightedEvents==1){
				cValPlot = "ValidateToy0_";
				cValPlot.Append("weight");
				cCanValPlot = "canvasValidateToy0_";
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


	for (int j = 0; j < nCat_Write; j++){
		delete outAsimovTree[j];
		for (int iToy = 0; iToy < maxNToysPossible; iToy++){
			int iArray = nCat_Write * iToy + j;
			delete outTree[iArray];
		};
	};
	delete[] outTree;

	toyfile->Close();
}
