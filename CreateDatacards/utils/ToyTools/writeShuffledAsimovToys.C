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
*** Be careful about this ordering in your datacards below and name your channels appropriately when combining datacards:

	TString chan[6]={"ch1","ch2","ch3","ch4","ch5","ch6"};
	int order[6]={2,1,0,4,3,5};
	TString channameO[6]={"7TeV/4e","7TeV/4mu","7TeV/2mu2e","8TeV/4mu","8TeV/2mu2e","8TeV/4e"};
	double normO[6]=	{	0.6951,	1.2439,	1.6662,	5.9471,	7.6807,	3.0898	}; 
	double normzxO[6]=	{	0.6226,	0.2230,	1.0628,	1.1878,	4.2929,	2.7676	};
	double normqqO[6]=	{	0.8386,	1.7971,	2.2456,	7.6478,	8.8585,	2.9364	};
	double normggO[6]=	{	0.0341,	0.0625,	0.0741,	0.4131,	0.5005,	0.2041	};

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

	const int N_Kinematics=8;
	const int N_KDs = 20;
	const int firstIntKD=9; <<<<<------- Increment this number if you add a pure hypothesis KD into the array above
	TString KDlist[N_KDs] = {
		"D_bkg", // 0
		"D_g1Q2_phi0", // 1
		"D_g1_vs_g2_phi0", // 2
		"D_g1_vs_g4_phi0", // 3
		"D_ZG", // 4
		"D_GG", // 5
		"D_ZG_PS", // 6
		"D_GG_PS", // 7
		"D_ZG_L1", // 8
		"D_g1Q2int_phi0", // 9
		"D_g2int_phi0", // 10
		"D_g4int_phi0", // 11
		"D_g2int_phi90", // 12
		"D_g4int_phi90", // 13
		"D_ZGint", // 14
		"D_GGint", // 15
		"D_ZG_PSint", // 16
		"D_GG_PSint", // 17
		"D_ZG_L1int_phi0", // 18
		"D_ZG_L1int_phi90" // 19
	};
	TString Kinematicslist[N_Kinematics] = {
		"ZZMass", // 0
		"Z1Mass", // 1
		"Z2Mass", // 2
		"costhetastar", // 3
		"helcosthetaZ1", // 4
		"helcosthetaZ2", // 5
		"helphi", // 6
		"phistarZ1" // 7
	};

*** KDs and kinematics are written with names

	TString KDlist_CMS[N_KDs] = {
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
	TString Kinematicslist_CMS[N_Kinematics] = {
		"CMS_zz4l_mass",
		"CMS_zz4l_Z1Mass",
		"CMS_zz4l_Z2Mass",
		"CMS_zz4l_costhetastar",
		"CMS_zz4l_helcosthetaZ1",
		"CMS_zz4l_helcosthetaZ2",
		"CMS_zz4l_helphi",
		"CMS_zz4l_phistarZ1"
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

void writeShuffledAsimovToys(int BSMSample=0, double signalStrength=1, int use_qqZZ_Dedicated=0){
	RooFit::Verbose(false);
	RooFit::PrintLevel(-1000);
	RooMsgService::instance().setStreamStatus(1,false);
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const int kNumHypo=52;
	const int kNumFiles=24;
	const int gapZZHypo=13;
	const int gapZZHypo_2=18;
	const int firstNonZZHypo=39;

/***** RUN CONDITIONS *****/
	if(signalStrength<0) assert(0);
	if(use_qqZZ_Dedicated<0) assert(0);
	if(BSMSample>=kNumHypo) assert(0);
/*** END RUN CONDITIONS ***/
	int maxNtoysPossible=1000000;

	const int N_Kinematics=8;
	const int N_KDs = 20;
	const int firstIntKD=9;
	TString KDlist[N_KDs] = {
		"D_bkg", // 0
		"D_g1Q2_phi0", // 1
		"D_g1_vs_g2_phi0", // 2
		"D_g1_vs_g4_phi0", // 3
		"D_ZG", // 4
		"D_GG", // 5
		"D_ZG_PS", // 6
		"D_GG_PS", // 7
		"D_ZG_L1", // 8
		"D_g1Q2int_phi0", // 9
		"D_g2int_phi0", // 10
		"D_g4int_phi0", // 11
		"D_g2int_phi90", // 12
		"D_g4int_phi90", // 13
		"D_ZGint", // 14
		"D_GGint", // 15
		"D_ZG_PSint", // 16
		"D_GG_PSint", // 17
		"D_ZG_L1int_phi0", // 18
		"D_ZG_L1int_phi90" // 19
	};
	TString Kinematicslist[N_Kinematics] = {
		"ZZMass", // 0
		"Z1Mass", // 1
		"Z2Mass", // 2
		"costhetastar", // 3
		"helcosthetaZ1", // 4
		"helcosthetaZ2", // 5
		"helphi", // 6
		"phistarZ1" // 7
	};

	TString KDlist_CMS[N_KDs] = {
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
	TString Kinematicslist_CMS[N_Kinematics] = {
		"CMS_zz4l_mass",
		"CMS_zz4l_Z1Mass",
		"CMS_zz4l_Z2Mass",
		"CMS_zz4l_costhetastar",
		"CMS_zz4l_helcosthetaZ1",
		"CMS_zz4l_helcosthetaZ2",
		"CMS_zz4l_helphi",
		"CMS_zz4l_phistarZ1"
	};

	char* cMCwgt[2] = {
		"MC_CV_weight",
		"MC_CV_4GeVcut_weight"
	};


	TString chan[6]={"ch1","ch2","ch3","ch4","ch5","ch6"};
//	int order[6]={1,2,0,5,3,4};
	int order[6]={2,1,0,4,3,5};
	TString channameO[6]={"7TeV/4e","7TeV/4mu","7TeV/2mu2e","8TeV/4mu","8TeV/2mu2e","8TeV/4e"};
//	double normO[6]=	{	0.7194,	1.2859,	1.7264,	6.0802,	7.9085,	3.1441	}; 
	double normO[6]=	{	0.6951,	1.2439,	1.6662,	5.9471,	7.6807,	3.0898	}; 
	double normzxO[6]=	{	0.6226,	0.2230,	1.0628,	1.1878,	4.2929,	2.7676	};
	double normqqO[6]=	{	0.8386,	1.7971,	2.2456,	7.6478,	8.8585,	2.9364	};
	double normggO[6]=	{	0.0341,	0.0625,	0.0741,	0.4131,	0.5005,	0.2041	};

	double norm[6] = { 0 };
	double normgg[6] = { 0 };
	double normzx[6] = { 0 };
	double normqq[6] = { 0 };
	TString channame[6];
	TTree* outTree[6];
	TTree* tsig[6];
	TTree* tggzz[6];
	TTree* tqqzz[6];
	TTree* tzx[6];


/**** PRODUCE INTERMEDIATE TREES ****/

	float KD_vars[N_KDs] = { 0 };
	double KD_pass[N_KDs] = { 0 };
	float Kinematics_vars[N_Kinematics] = { 0 };
	double Kinematics_pass[N_Kinematics] = { 0 };

	short Z1ids,Z2ids;	
	int rZ1ids,rZ2ids;	

	double weightFit=0, weightFit2=0;
	float MC_CV_weight[kNumHypo];
	float MC_weight_QQBGGProper[2];
	float MC_weight_noxsec;
	float MC_weight_Kfactor,ZXfake_weightProper;
	int EventSample;
	int ToyNumber=0;
	double sumIndToyWgt=0;
	int store_NToyEvents=0;
	int counter_IndToyEvent=0;

	for(int j=0;j<6;j++){
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
		outTree[j] = new TTree(chan[j],channame[j]);

		outTree[j]->Branch("Z1ids", &rZ1ids);
		outTree[j]->Branch("Z2ids", &rZ2ids);
		outTree[j]->Branch("_weight_",&weightFit,"_weight_/D");

		tsig[j]->Branch("Z1ids", &rZ1ids);
		tsig[j]->Branch("Z2ids", &rZ2ids);
//		tsig[j]->Branch("_Asimov_weight_",&weightFit,"_Asimov_weight_/D");
		tsig[j]->Branch("_weight_",&weightFit2,"_weight_/D");
		tsig[j]->Branch("ToyNumber",&ToyNumber,"ToyNumber/I");
		tsig[j]->Branch("NToyEvents",&store_NToyEvents,"NToyEvents/I");

		tggzz[j]->Branch("Z1ids", &rZ1ids);
		tggzz[j]->Branch("Z2ids", &rZ2ids);
//		tggzz[j]->Branch("_Asimov_weight_",&weightFit,"_Asimov_weight_/D");
		tggzz[j]->Branch("_weight_",&weightFit2,"_weight_/D");
		tggzz[j]->Branch("ToyNumber",&ToyNumber,"ToyNumber/I");
		tggzz[j]->Branch("NToyEvents",&store_NToyEvents,"NToyEvents/I");

		tqqzz[j]->Branch("Z1ids", &rZ1ids);
		tqqzz[j]->Branch("Z2ids", &rZ2ids);
//		tqqzz[j]->Branch("_Asimov_weight_",&weightFit,"_Asimov_weight_/D");
		tqqzz[j]->Branch("_weight_",&weightFit2,"_weight_/D");
		tqqzz[j]->Branch("ToyNumber",&ToyNumber,"ToyNumber/I");
		tqqzz[j]->Branch("NToyEvents",&store_NToyEvents,"NToyEvents/I");

		tzx[j]->Branch("Z1ids", &rZ1ids);
		tzx[j]->Branch("Z2ids", &rZ2ids);
//		tzx[j]->Branch("_Asimov_weight_",&weightFit,"_Asimov_weight_/D");
		tzx[j]->Branch("_weight_",&weightFit2,"_weight_/D");
		tzx[j]->Branch("ToyNumber",&ToyNumber,"ToyNumber/I");
		tzx[j]->Branch("NToyEvents",&store_NToyEvents,"NToyEvents/I");

		for (int k = 0; k < N_KDs; k++){
			tsig[j]->Branch(KDlist_CMS[k], &KD_pass[k]);
			tggzz[j]->Branch(KDlist_CMS[k], &KD_pass[k]);
			tqqzz[j]->Branch(KDlist_CMS[k], &KD_pass[k]);
			tzx[j]->Branch(KDlist_CMS[k], &KD_pass[k]);
			outTree[j]->Branch(KDlist_CMS[k], &KD_pass[k]);
		};
		for(int k=0;k<N_Kinematics;k++){
			tsig[j]->Branch(Kinematicslist_CMS[k],&Kinematics_pass[k]);
			tggzz[j]->Branch(Kinematicslist_CMS[k],&Kinematics_pass[k]);
			tqqzz[j]->Branch(Kinematicslist_CMS[k],&Kinematics_pass[k]);
			tzx[j]->Branch(Kinematicslist_CMS[k],&Kinematics_pass[k]);
			outTree[j]->Branch(Kinematicslist_CMS[k],&Kinematics_pass[k]);
		};
	}

/***** NEED TO SHUFFLE ZX ****/

	TString files_ZX[6]{
		"/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/2mu2e/HZZ4lTree_H125p6_ShuffledSignalBkg.root",
		"/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/4e/HZZ4lTree_H125p6_ShuffledSignalBkg.root",
		"/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/4mu/HZZ4lTree_H125p6_ShuffledSignalBkg.root",
		"/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_7TeV/2mu2e/HZZ4lTree_H125p6_ShuffledSignalBkg.root",
		"/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_7TeV/4e/HZZ4lTree_H125p6_ShuffledSignalBkg.root",
		"/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_7TeV/4mu/HZZ4lTree_H125p6_ShuffledSignalBkg.root"
	};
	TChain* zx_Container[6];
	int nevents_ZX[6]={0};
	int naccevents_ZX[6]={0};
	TRandom3 shuffler_ZX(123456);
	TRandom3 assigner_ZX(891011);
	for (int c = 0; c < 6; c++){
		zx_Container[c] = new TChain("SelectedTree_ZX");
		zx_Container[c]->Add(files_ZX[c]);
		zx_Container[c]->SetBranchAddress("Z1ids",&Z1ids);
		zx_Container[c]->SetBranchAddress("Z2ids",&Z2ids);
		zx_Container[c]->SetBranchAddress("ZXfake_weightProper",&ZXfake_weightProper);
		for (int k = 0; k < N_KDs; k++){
			zx_Container[c]->SetBranchAddress(KDlist[k], &KD_vars[k]);
		};
		for (int k = 0; k < N_Kinematics; k++){
			zx_Container[c]->SetBranchAddress(Kinematicslist[k], &Kinematics_vars[k]);
		}
		nevents_ZX[c] = zx_Container[c]->GetEntries();
		naccevents_ZX[c] = nevents_ZX[c];
//		cout << naccevents_ZX[c] << endl;
	};
	for (int c = 1; c < 6; c++) naccevents_ZX[c] += naccevents_ZX[c-1];
	TTree* zx = new TTree("SelectedTree_ZX_Combined","");
	zx->Branch("Z1ids",&Z1ids);
	zx->Branch("Z2ids",&Z2ids);
	zx->Branch("ZXfake_weightProper",&ZXfake_weightProper);
	for (int k = 0; k < N_KDs; k++){
		zx->Branch(KDlist[k], &KD_vars[k]);
	};
	for (int k = 0; k < N_Kinematics; k++){
		zx->Branch(Kinematicslist[k], &Kinematics_vars[k]);
	}
	int nUsed_ZX[6]={0};
	int ctr=0;
	int ntotal_ZX=naccevents_ZX[5];
	cout << "Total ZX before: " << ntotal_ZX << endl;
	while (ctr < ntotal_ZX){
		int coin = shuffler_ZX.Integer(naccevents_ZX[5]);
		int lucky=-1;
		for(int smp=0;smp<6;smp++){
			if(coin<naccevents_ZX[smp] && (nevents_ZX[smp]-nUsed_ZX[smp])!=0){lucky=smp;break;};
		};
		if(lucky<0) continue;
		zx_Container[lucky]->GetEntry(nUsed_ZX[lucky]);
		nUsed_ZX[lucky] += 1;
		zx->Fill();
		for(int smp=lucky;smp<6;smp++) naccevents_ZX[smp] -= 1;
		ctr++;
	};
	for (int c = 0; c < 6; c++) delete zx_Container[c];
	cout << "Total ZX after: " << zx->GetEntries() << endl;

/**** BEGIN *****/

	for (int j=0; j<6; j++){
		cout<<chan[j]<<endl;
		TChain *sig = new TChain("SelectedTree");
		TChain *gg=new TChain("SelectedTree_ggZZ");
		TChain *qq;
		if(use_qqZZ_Dedicated==0) qq = new TChain("SelectedTree_qqZZ");
		else qq = new TChain("SelectedTree_qqZZ_Dedicated");
		sig->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_"+channame[j]+"/HZZ4lTree_H125p6_ShuffledSignalBkg.root");
		gg->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_"+channame[j]+"/HZZ4lTree_H125p6_ShuffledSignalBkg.root");
		qq->Add("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_"+channame[j]+"/HZZ4lTree_H125p6_ShuffledSignalBkg.root");

/**** SET SIGNAL YIELDS ****/

		int BSMFile=BSMSample;
		if (BSMSample<0 || BSMSample>=kNumHypo){
			cerr << "No hypothesis" << endl;
			assert(0);
		}
		else if (BSMSample>gapZZHypo && BSMSample<gapZZHypo_2){
			cout << "Passed gap hypothesis no. 1, " << gapZZHypo << ", while file range is still usable." << endl;
			BSMFile-=1;
		}
		else if (BSMSample>gapZZHypo_2 && BSMSample < (kNumFiles + 1)){
			cout << "Passed gap hypothesis no. 2, " << gapZZHypo << ", while file range is still usable." << endl;
			BSMFile-=2;
		}
		else if (BSMSample == gapZZHypo || BSMSample == gapZZHypo_2 || BSMSample >= (kNumFiles + 1)){
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

/**** SET RANDOM NUMBER GENERATORS ****/

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

		sig->SetBranchAddress("Z1ids",&Z1ids);
		sig->SetBranchAddress("Z2ids",&Z2ids);
		qq->SetBranchAddress("Z1ids",&Z1ids);
		qq->SetBranchAddress("Z2ids",&Z2ids);
		gg->SetBranchAddress("Z1ids",&Z1ids);
		gg->SetBranchAddress("Z2ids",&Z2ids);

		sig->SetBranchAddress(cMCwgtchosen,MC_CV_weight);
		sig->SetBranchAddress("EventSample",&EventSample);
		qq->SetBranchAddress("MC_weight_QQBGGProper",MC_weight_QQBGGProper);
		qq->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec);
		qq->SetBranchAddress("MC_weight_Kfactor",&MC_weight_Kfactor);
		gg->SetBranchAddress("MC_weight_QQBGGProper",MC_weight_QQBGGProper);
		gg->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec);
		gg->SetBranchAddress("MC_weight_Kfactor",&MC_weight_Kfactor);
		
		for (int k = 0; k < N_KDs; k++){
			sig->SetBranchAddress(KDlist[k], &KD_vars[k]);
			gg->SetBranchAddress(KDlist[k], &KD_vars[k]);
			qq->SetBranchAddress(KDlist[k], &KD_vars[k]);
		};
		for (int k = 0; k < N_Kinematics; k++){
			sig->SetBranchAddress(Kinematicslist[k], &Kinematics_vars[k]);
			gg->SetBranchAddress(Kinematicslist[k], &Kinematics_vars[k]);
			qq->SetBranchAddress(Kinematicslist[k], &Kinematics_vars[k]);
		}

		int nZerosThrown=0;

		double sumSigValidate=0;
			int nentry = sig->GetEntries();
			int nFilled=0;
			for(int i=0;i<nentry;i++){
				sig->GetEntry(i);
				if(Kinematics_vars[0]<140.6 && Kinematics_vars[0]>=105.6){
					if(EventSample==BSMFile && MC_CV_weight[BSMSample]>0){
						weightFit = signalNormScale*MC_CV_weight[BSMSample];
						rZ1ids=Z1ids;
						rZ2ids=Z2ids;
						for (int k = 0; k < N_KDs; k++) KD_pass[k]=KD_vars[k];
						for (int k = 0; k < N_Kinematics; k++) Kinematics_pass[k]=Kinematics_vars[k];
						outTree[j]->Fill();

						sumSigValidate += weightFit;
						if(nFilled==0) cout << "Typical weight is " << weightFit << endl;
						nFilled++;

						if (counter_IndToyEvent == store_NToyEvents){
							store_NToyEvents = (int) rand_sig.PoissonD(NBSMexpected)+0.5;
							counter_IndToyEvent = 0;
							sumIndToyWgt=0;
							if (store_NToyEvents == 0){
								count_sig.push_back(store_NToyEvents);
								fill_sig.push_back(sumIndToyWgt);
								nZerosThrown++;
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
			cout << "Number of 0s thrown: " << nZerosThrown << endl;
			nZerosThrown=0;

			sumIndToyWgt=0;
			store_NToyEvents=0;
			counter_IndToyEvent=0;
			for(int i=0;i<nentry;i++){
				sig->GetEntry(i);
				if(Kinematics_vars[0]<140.6 && Kinematics_vars[0]>=105.6){
					if(EventSample==BSMFile && MC_CV_weight[BSMSample]>0){
						if(count_sig.size()<=ToyNumber) continue;
						store_NToyEvents = count_sig[ToyNumber];
						sumIndToyWgt = fill_sig[ToyNumber];
						double nexp=store_NToyEvents;
						double sumtoywgt=sumIndToyWgt;
						if(nexp==0) sumtoywgt=1;
						double scale = nexp/sumtoywgt;

						weightFit2 = scale*MC_CV_weight[BSMSample];
						rZ1ids=Z1ids;
						rZ2ids=Z2ids;
						for (int k = 0; k < N_KDs; k++) KD_pass[k]=KD_vars[k];
						for (int k = 0; k < N_Kinematics; k++) Kinematics_pass[k]=Kinematics_vars[k];
						if(weightFit2>0) tsig[j]->Fill();
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

		    qq->Draw("ZZMass>>htqqZZ","MC_weight_QQBGGProper[1]*MC_weight_noxsec*(ZZMass<140.6&&ZZMass>=105.6)");
			bkgSum = htempqqZZ->GetSumOfWeights();
			cout << "qqZZ dedicated sum: " << bkgSum << endl;
			nentry = qq->GetEntries();
			double sumQQZZValidate=0;
			for(int i=0;i<nentry;i++){
				qq->GetEntry(i);
				if(Kinematics_vars[0]<140.6 && Kinematics_vars[0]>=105.6){
					weightFit = normqq[j]/bkgSum*MC_weight_QQBGGProper[1]*MC_weight_noxsec;
					rZ1ids=Z1ids;
					rZ2ids=Z2ids;
					for (int k = 0; k < N_KDs; k++) KD_pass[k]=KD_vars[k];
					for (int k = 0; k < N_Kinematics; k++) Kinematics_pass[k]=Kinematics_vars[k];
					outTree[j]->Fill();

					sumQQZZValidate += weightFit;
					if(nFilled==0) cout << "Typical weight is " << weightFit << endl;
					nFilled++;

					if (counter_IndToyEvent == store_NToyEvents){
						store_NToyEvents = (int) rand_qqzz.PoissonD(normqq[j])+0.5;
						counter_IndToyEvent = 0;
						sumIndToyWgt=0;
						if (store_NToyEvents == 0){
							count_qqzz.push_back(store_NToyEvents);
							fill_qqzz.push_back(sumIndToyWgt);
							nZerosThrown++;
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
			cout << "Number of 0s thrown: " << nZerosThrown << endl;
			nZerosThrown=0;

			sumIndToyWgt=0;
			store_NToyEvents=0;
			counter_IndToyEvent=0;
			for(int i=0;i<nentry;i++){
				qq->GetEntry(i);
				if(Kinematics_vars[0]<140.6 && Kinematics_vars[0]>=105.6){
					if(count_qqzz.size()<=ToyNumber) continue;
					store_NToyEvents = count_qqzz[ToyNumber];
					sumIndToyWgt = fill_qqzz[ToyNumber];
					double nexp=store_NToyEvents;
					double sumtoywgt=sumIndToyWgt;
					if(nexp==0) sumtoywgt=1;
					double scale = nexp/sumtoywgt;

					weightFit2 = scale*MC_weight_QQBGGProper[1]*MC_weight_noxsec;
					rZ1ids=Z1ids;
					rZ2ids=Z2ids;
					for (int k = 0; k < N_KDs; k++) KD_pass[k]=KD_vars[k];
					for (int k = 0; k < N_Kinematics; k++) Kinematics_pass[k]=Kinematics_vars[k];
					if(weightFit2>0) tqqzz[j]->Fill();
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


			gg->Draw("ZZMass>>htggZZ","(MC_weight_QQBGGProper[0]*MC_weight_noxsec*MC_weight_Kfactor)*(ZZMass<140.6&&ZZMass>=105.6)");
			bkgSum= htempggZZ->GetSumOfWeights();
			cout << "ggZZ dedicated sum: " << bkgSum << endl;
			nentry = gg->GetEntries();
			double sumGGZZValidate=0;
			for(int i=0;i<nentry;i++){
				gg->GetEntry(i);
				if(Kinematics_vars[0]<140.6 && Kinematics_vars[0]>=105.6){
					weightFit = normgg[j]/bkgSum*MC_weight_QQBGGProper[0]*MC_weight_noxsec*MC_weight_Kfactor;
					rZ1ids=Z1ids;
					rZ2ids=Z2ids;
					for (int k = 0; k < N_KDs; k++) KD_pass[k]=KD_vars[k];
					for (int k = 0; k < N_Kinematics; k++) Kinematics_pass[k]=Kinematics_vars[k];
					outTree[j]->Fill();

					sumGGZZValidate += weightFit;
					if(nFilled==0) cout << "Typical weight is " << weightFit << endl;
					nFilled++;

					if (counter_IndToyEvent == store_NToyEvents){
						store_NToyEvents = (int) rand_ggzz.PoissonD(normgg[j])+0.5;
						counter_IndToyEvent = 0;
						sumIndToyWgt=0;
						if (store_NToyEvents == 0){
							count_ggzz.push_back(store_NToyEvents);
							fill_ggzz.push_back(sumIndToyWgt);
							nZerosThrown++;
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
			cout << "Number of 0s thrown: " << nZerosThrown << endl;
			nZerosThrown=0;

			sumIndToyWgt=0;
			store_NToyEvents=0;
			counter_IndToyEvent=0;
			for(int i=0;i<nentry;i++){
				gg->GetEntry(i);
				if(Kinematics_vars[0]<140.6 && Kinematics_vars[0]>=105.6){
					if(count_ggzz.size()<=ToyNumber) continue;
					store_NToyEvents = count_ggzz[ToyNumber];
					sumIndToyWgt = fill_ggzz[ToyNumber];
					double nexp=store_NToyEvents;
					double sumtoywgt=sumIndToyWgt;
					if(nexp==0) sumtoywgt=1;
					double scale = nexp/sumtoywgt;

					weightFit2 = scale*MC_weight_QQBGGProper[0]*MC_weight_noxsec*MC_weight_Kfactor;
					rZ1ids=Z1ids;
					rZ2ids=Z2ids;
					for (int k = 0; k < N_KDs; k++) KD_pass[k]=KD_vars[k];
					for (int k = 0; k < N_Kinematics; k++) Kinematics_pass[k]=Kinematics_vars[k];
					if(weightFit2>0) tggzz[j]->Fill();
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

			zx->Draw("ZZMass>>htZX","(ZZMass<140.6&&ZZMass>=105.6 && ZXfake_weightProper>0)*ZXfake_weightProper");
			bkgSum= htempZX->GetSumOfWeights();
			cout << "ZX dedicated sum: " << bkgSum << endl;
			nentry = zx->GetEntries();
			double sumZXValidate=0;
			for(int i=0;i<nentry;i++){
				zx->GetEntry(i);
				if(Kinematics_vars[0]<140.6 && Kinematics_vars[0]>=105.6 && ZXfake_weightProper>0){
					if (j == 0 || j == 3){
						assigner_ZX.SetSeed(i);
						int crossAssign = assigner_ZX.Integer(2);
						if (crossAssign == 0){
							Z1ids = -121;
							Z2ids = -169;
						}
						else{
							Z2ids = -121;
							Z1ids = -169;
						};
					}
					else{
						if (j == 1 || j == 4){
							Z1ids = -169;
							Z2ids = -169;
						}
						else if (j == 2 || j == 5){
							Z1ids = -121;
							Z2ids = -121;
						}
					};

					weightFit= ZXfake_weightProper*normzx[j]/bkgSum;
					rZ1ids=Z1ids;
					rZ2ids=Z2ids;
					for (int k = 0; k < N_KDs; k++) KD_pass[k]=KD_vars[k];
					for (int k = 0; k < N_Kinematics; k++) Kinematics_pass[k]=Kinematics_vars[k];
					outTree[j]->Fill();

					sumZXValidate += weightFit;
					if(nFilled==0) cout << "Typical weight is " << weightFit << endl;
					nFilled++;

					if (counter_IndToyEvent == store_NToyEvents){
						store_NToyEvents = (int) rand_zx.PoissonD(normzx[j])+0.5;
						counter_IndToyEvent = 0;
						sumIndToyWgt=0;
						if (store_NToyEvents == 0){
							count_zx.push_back(store_NToyEvents);
							fill_zx.push_back(sumIndToyWgt);
							nZerosThrown++;
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
			cout << "Number of 0s thrown: " << nZerosThrown << endl;
			nZerosThrown=0;

			sumIndToyWgt=0;
			store_NToyEvents=0;
			counter_IndToyEvent=0;
			for(int i=0;i<nentry;i++){
				zx->GetEntry(i);
				if(Kinematics_vars[0]<140.6 && Kinematics_vars[0]>=105.6 && ZXfake_weightProper>0){
					if(count_zx.size()<=ToyNumber) continue;
					store_NToyEvents = count_zx[ToyNumber];
					sumIndToyWgt = fill_zx[ToyNumber];
					double nexp=store_NToyEvents;
					double sumtoywgt=sumIndToyWgt;
					if(nexp==0) sumtoywgt=1;
					double scale = nexp/sumtoywgt;

					if (j == 0 || j == 3){
						assigner_ZX.SetSeed(i);
						int crossAssign = assigner_ZX.Integer(2);
						if (crossAssign == 0){
							Z1ids = -121;
							Z2ids = -169;
						}
						else{
							Z2ids = -121;
							Z1ids = -169;
						};
					}
					else{
						if (j == 1 || j == 4){
							Z1ids = -169;
							Z2ids = -169;
						}
						else if (j == 2 || j == 5){
							Z1ids = -121;
							Z2ids = -121;
						}
					};

					weightFit2 = scale*ZXfake_weightProper;
					rZ1ids=Z1ids;
					rZ2ids=Z2ids;
					for (int k = 0; k < N_KDs; k++) KD_pass[k]=KD_vars[k];
					for (int k = 0; k < N_Kinematics; k++) Kinematics_pass[k]=Kinematics_vars[k];
					if(weightFit2>0) tzx[j]->Fill();
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


		delete htempqqZZ;
		delete htempggZZ;
		delete htempZX;
		delete htempSM;
		delete htempSMCut;
		delete htempBSM;
		delete htempBSMDedicated;
		delete sig;
		delete gg;
		delete qq;
	}
	delete zx;

	char cDir[]="./AsimovToys";
	char coutput[50];
	if (use_qqZZ_Dedicated == 0) sprintf(coutput, "%s/Asimov_v1_Hypo%i.root", cDir, BSMSample);
	else sprintf(coutput, "%s/Asimov_v2_Hypo%i.root", cDir, BSMSample);

	cout << coutput << endl;
	TFile *toyfile = new TFile(coutput,"recreate");
	cout << coutput << endl;

	toyfile->cd();
	for(int j=0;j<6;j++){
		toyfile->WriteTObject(tsig[j]);
		toyfile->WriteTObject(tggzz[j]);
		toyfile->WriteTObject(tqqzz[j]);
		toyfile->WriteTObject(tzx[j]);
		toyfile->WriteTObject(outTree[j]);
		delete tsig[j];
		delete tggzz[j];
		delete tqqzz[j];
		delete tzx[j];
		delete outTree[j];
	};

	cout << "Maximum possible number of toys: " << maxNtoysPossible << endl;
	toyfile->Close();
}
