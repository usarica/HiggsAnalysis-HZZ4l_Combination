/* 
 * Compute ZZ background rates and write them in a card fragment.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b ZZbackgroundRate.C+
 * This runs on 7 and 8 TeV and writes the output in a file (see stdout).
 *
 */


#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TAxis.h"
#include "TFile.h"
#include "TLegend.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include "TTree.h"
#include "TText.h"
#include "TROOT.h"
#include "TStyle.h"

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

//VBFtag convention: 0->0/1 Jets; 1->2+ jets; 2->No tagging

const int SelectionScheme=4;
double massRange[2] ={ 105.6, 140.6 };

using namespace std;
void compute(TString name, int sqrts, double lumi, int VBFtag);
double sumWeights(TString name, double lumi, int selVBF, int isQQBZZ=0);


// Run both sqrts in one go
void ZZbackgroundRate_onshellMassRange() {

  compute(filePath7TeV, 7, lumi7TeV, 1);
  compute(filePath8TeV, 8, lumi8TeV, 1);

  compute(filePath7TeV, 7, lumi7TeV, 0);
  compute(filePath8TeV, 8, lumi8TeV, 0);

  compute(filePath7TeV, 7, lumi7TeV, 2);
  compute(filePath8TeV, 8, lumi8TeV, 2);

}

bool testSelection(
	float Lep1SIP, float Lep2SIP, float Lep3SIP, float Lep4SIP,
	float Lep1_Z1SIP, float Lep2_Z1SIP, float Lep3_Z1SIP, float Lep4_Z1SIP,
	float KalmanCandVtx_chi2,
	int scheme=4
	){
	bool pass=true;
	switch (scheme){
	case 0:
		if (Lep1SIP>=4 || Lep2SIP>=4 || Lep3SIP>=4 || Lep4SIP>=4) pass = false;
		break;
	case 1:
		break;
	case 2:
		if (Lep1_Z1SIP>=4 || Lep2_Z1SIP>=4 || Lep3_Z1SIP>=5 || Lep4_Z1SIP>=5) pass = false;
		break;
	case 3:
		if (KalmanCandVtx_chi2>=30) pass = false;
		break;
	case 4:
		if (Lep1_Z1SIP>=4 || Lep2_Z1SIP>=4 || Lep3_Z1SIP>=5 || Lep4_Z1SIP>=5 || KalmanCandVtx_chi2>=30) pass = false;
		break;
	default:
		break;
	}
	return pass;
}

// The actual job
void compute(TString filePath, int sqrts, double lumi, int VBFtag){

  double qqZZ[3];
  double ggZZ[3];
  
	setprecision(9);
  cout<<"qqZZ 4mu"<<endl;  
  double qqZZ_4mu_ev4mu = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo4mu_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_4mu_ev4e = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo4e_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_4mu_ev2e2mu = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo2e2mu_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_4mu_ev4tau = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo4tau_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_4mu_ev2mu2tau = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo2mu2tau_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_4mu_ev2e2tau = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo2e2tau_Reprocessed.root",lumi,VBFtag,1);
  qqZZ[0]=qqZZ_4mu_ev4mu+qqZZ_4mu_ev4e+qqZZ_4mu_ev2e2mu+qqZZ_4mu_ev4tau+qqZZ_4mu_ev2mu2tau+qqZZ_4mu_ev2e2tau;
  cout<<qqZZ_4mu_ev4mu<<" "<<qqZZ_4mu_ev4e<<" "<<qqZZ_4mu_ev2e2mu<<" "<<qqZZ_4mu_ev4tau<<" "<<qqZZ_4mu_ev2mu2tau<<" "<<qqZZ_4mu_ev2e2tau
      <<" -> "<<qqZZ_4mu_ev4mu<<" + "
      <<qqZZ_4mu_ev4e+qqZZ_4mu_ev2e2mu+qqZZ_4mu_ev4tau+qqZZ_4mu_ev2mu2tau+qqZZ_4mu_ev2e2tau
      <<" -> "<< qqZZ[0]<< endl;
  

  cout<<"qqZZ 4e"<<endl;
  double qqZZ_4e_ev4e = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo4e_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_4e_ev4mu = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo4mu_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_4e_ev2e2mu = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo2e2mu_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_4e_ev4tau = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo4tau_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_4e_ev2mu2tau = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo2mu2tau_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_4e_ev2e2tau = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo2e2tau_Reprocessed.root",lumi,VBFtag,1);
  qqZZ[1]=qqZZ_4e_ev4e+qqZZ_4e_ev4mu+qqZZ_4e_ev2e2mu+qqZZ_4e_ev4tau+qqZZ_4e_ev2mu2tau+qqZZ_4e_ev2e2tau;
  cout<<qqZZ_4e_ev4e<<" "<<qqZZ_4e_ev4mu<<" "<<qqZZ_4e_ev2e2mu<<" "<<qqZZ_4e_ev4tau<<" "<<qqZZ_4e_ev2mu2tau<<" "<<qqZZ_4e_ev2e2tau
      <<" -> "<<qqZZ_4e_ev4e<<" + "
      <<qqZZ_4e_ev4mu+qqZZ_4e_ev2e2mu+qqZZ_4e_ev4tau+qqZZ_4e_ev2mu2tau+qqZZ_4e_ev2e2tau
      <<" -> "<< qqZZ[1] <<endl;

  cout<<"qqZZ 2e2mu"<<endl;
  double qqZZ_2e2mu_ev2e2mu = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo2e2mu_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_2e2mu_ev4mu = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo4mu_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_2e2mu_ev4e = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo4e_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_2e2mu_ev4tau = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo4tau_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_2e2mu_ev2mu2tau = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo2mu2tau_Reprocessed.root",lumi,VBFtag,1);
  double qqZZ_2e2mu_ev2e2tau = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo2e2tau_Reprocessed.root",lumi,VBFtag,1);
  qqZZ[2] = qqZZ_2e2mu_ev2e2mu+qqZZ_2e2mu_ev4e+qqZZ_2e2mu_ev4mu+qqZZ_2e2mu_ev4tau+qqZZ_2e2mu_ev2mu2tau+qqZZ_2e2mu_ev2e2tau; 
  cout<<qqZZ_2e2mu_ev2e2mu<<" "<<qqZZ_2e2mu_ev4e<<" "<<qqZZ_2e2mu_ev4mu<<" "<<qqZZ_2e2mu_ev4tau<<" "<<qqZZ_2e2mu_ev2mu2tau<<" "<<qqZZ_2e2mu_ev2e2tau
      <<" -> "<<qqZZ_2e2mu_ev2e2mu<<" + "
      <<qqZZ_2e2mu_ev4e+qqZZ_2e2mu_ev4mu+qqZZ_2e2mu_ev4tau+qqZZ_2e2mu_ev2mu2tau+qqZZ_2e2mu_ev2e2tau
      <<" -> "<<qqZZ[2]<<endl;
  
  cout<<"ggZZ 4mu"<<endl;
  double ggZZ_4mu_ev2l2l = sumWeights(filePath + "/4mu/HZZ4lTree_ggZZ2l2l_Reprocessed.root",lumi,VBFtag);
  double ggZZ_4mu_ev4l = sumWeights(filePath + "/4mu/HZZ4lTree_ggZZ4l_Reprocessed.root",lumi,VBFtag);
  ggZZ[0] = ggZZ_4mu_ev4l + ggZZ_4mu_ev2l2l;
  cout<<ggZZ_4mu_ev4l<<" + "<<ggZZ_4mu_ev2l2l<<" -> "<< ggZZ[0] <<endl;

  cout<<"ggZZ 4e"<<endl;
  double ggZZ_4e_ev2l2l = sumWeights(filePath + "/4e/HZZ4lTree_ggZZ2l2l_Reprocessed.root",lumi,VBFtag);
  double ggZZ_4e_ev4l = sumWeights(filePath + "/4e/HZZ4lTree_ggZZ4l_Reprocessed.root",lumi,VBFtag);
  ggZZ[1] = ggZZ_4e_ev4l + ggZZ_4e_ev2l2l;
  cout<<ggZZ_4e_ev4l<<" + "<<ggZZ_4e_ev2l2l<<" -> "<< ggZZ[1] <<endl;

  cout<<"ggZZ 2e2mu"<<endl;
  double ggZZ_2e2mu_ev2l2l = sumWeights(filePath + "/2mu2e/HZZ4lTree_ggZZ2l2l_Reprocessed.root",lumi,VBFtag);
  double ggZZ_2e2mu_ev4l = sumWeights(filePath + "/2mu2e/HZZ4lTree_ggZZ4l_Reprocessed.root",lumi,VBFtag);
  ggZZ[2] = ggZZ_2e2mu_ev4l + ggZZ_2e2mu_ev2l2l;
  cout<<ggZZ_2e2mu_ev4l<<" + "<<ggZZ_2e2mu_ev2l2l<<" -> "<< ggZZ[2] <<endl;

  TString schannel[3] = {"4mu","4e","2e2mu"};
  TString ssqrts = (long) sqrts + TString("TeV");
  for (int i=0; i<3; ++i) {
    TString outfile;
		TString crange = Form("M4l_%.1f-%.1f", massRange[0], massRange[1]);
		if (VBFtag<2) outfile = "CardFragments/ZZRates_" + crange + "_" + ssqrts + "_" + schannel[i] + "_" + Form("%d", int(VBFtag)) + ".txt";
		if (VBFtag==2) outfile = "CardFragments/ZZRates_" + crange + "_" + ssqrts + "_" + schannel[i] + ".txt";
    ofstream of(outfile,ios_base::out);
    of << "## rates --- format = chan N lumi ##" << endl
       << "## if lumi is blank, lumi for cards used ##" << endl;
    of << "rate qqZZ  " << qqZZ[i] << endl;
    of << "rate ggZZ  " << ggZZ[i] << endl;
    of.close();
    cout << "Output written to: " << outfile << endl;
  }
}


double sumWeights(TString name, double lumi, int selVBF, int isQQBZZ){
  TChain* tree = new TChain("SelectedTree");
  tree->Add((TString)(name));

	Float_t Lep1SIP=0, Lep2SIP=0, Lep3SIP=0, Lep4SIP=0;
	Float_t Lep1_Z1SIP=0, Lep2_Z1SIP=0, Lep3_Z1SIP=0, Lep4_Z1SIP=0;
	Float_t KalmanCandVtx_chi2=0;
	float ZZMass;
	float MC_weight;
	float MC_weight_extra;
	Short_t NJets;
  double totEvents=0;

	tree->SetBranchAddress("ZZMass", &ZZMass);
	tree->SetBranchAddress("MC_weight", &MC_weight);
	if (isQQBZZ==1) tree->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_extra);
	else tree->SetBranchAddress("MC_weight_Kfactor", &MC_weight_extra);
	tree->SetBranchAddress("NJets30", &NJets);
	if (SelectionScheme!=0){
		tree->SetBranchAddress("Lep1SIP", &Lep1SIP);
		tree->SetBranchAddress("Lep2SIP", &Lep2SIP);
		tree->SetBranchAddress("Lep3SIP", &Lep3SIP);
		tree->SetBranchAddress("Lep4SIP", &Lep4SIP);
		tree->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
		tree->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
		tree->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
		tree->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
		tree->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
	}
	
	for (int iEvt=0; iEvt<tree->GetEntries(); iEvt++){
    tree->GetEntry(iEvt);
		if (ZZMass<massRange[0] || ZZMass>=massRange[1]) continue;
		if ((selVBF == 1 && NJets > 1) || (selVBF == 0 && NJets < 2) || (selVBF==2)) totEvents += MC_weight*MC_weight_extra;
  }

  return totEvents*lumi;
}



