/* 
 * Fit ggZZ background shapes and write parameters in a card fragment.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b backgroundFits_ggzz_1Dw.C
 * This runs on the 3 final states for 7 and 8 TeV and writes the output in a file (see stdout).
 *
 */

/*
 #ifndef __CINT__
 #include "RooGlobalFunc.h"
 #endif
 #include "RooRealVar.h"
 #include "RooDataSet.h"
 #include "RooGaussian.h"
 #include "RooConstVar.h"
 #include "RooChebychev.h"
 #include "RooAddPdf.h"
 #include "RooWorkspace.h"
 #include "RooPlot.h"
 #include "TCanvas.h"
 #include "TAxis.h"
 #include "TFile.h"
 #include "TH1.h"
*/

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"

#include "../../CombinedLimit/interface/HZZ4LRooPdfs.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;
using namespace RooFit;



//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

const int SelectionScheme=4;

void backgroundFits_ggzz_1Dw(int channel, int sqrts, int VBFtag);

// Run all final states and sqrts in one go
void backgroundFits_ggzz_1Dw_useMCFM() {
  gSystem->Exec("mkdir -p bkgFigs7TeV");
  gSystem->Exec("mkdir -p bkgFigs8TeV");

//  backgroundFits_ggzz_1Dw(1,7,1);
//  backgroundFits_ggzz_1Dw(2,7,1);
//  backgroundFits_ggzz_1Dw(3,7,1);  
  backgroundFits_ggzz_1Dw(1,8,1);
  backgroundFits_ggzz_1Dw(2,8,1);
  backgroundFits_ggzz_1Dw(3,8,1);

//  backgroundFits_ggzz_1Dw(1,7,0);
//  backgroundFits_ggzz_1Dw(2,7,0);
//  backgroundFits_ggzz_1Dw(3,7,0);  
  backgroundFits_ggzz_1Dw(1,8,0);
  backgroundFits_ggzz_1Dw(2,8,0);
  backgroundFits_ggzz_1Dw(3,8,0);

//  backgroundFits_ggzz_1Dw(1,7,2);
//  backgroundFits_ggzz_1Dw(2,7,2);
//  backgroundFits_ggzz_1Dw(3,7,2);  
  backgroundFits_ggzz_1Dw(1,8,2);
  backgroundFits_ggzz_1Dw(2,8,2);
  backgroundFits_ggzz_1Dw(3,8,2);
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
void backgroundFits_ggzz_1Dw(int channel, int sqrts, int VBFtag)
{
  TString schannel;
  if      (channel == 1) schannel = "4mu";
  else if (channel == 2) schannel = "4e";
  else if (channel == 3) schannel = "2e2mu";
  else cout << "Not a valid channel: " << schannel << endl;

  TString ssqrts = (long) sqrts + TString("TeV");

  cout << "schannel = " << schannel << "  sqrts = " << sqrts << " VBFtag = "<< VBFtag << endl;

  TString outfile;
  if(VBFtag<2) outfile = "CardFragments/ggzzMCFMBackgroundFit_" + ssqrts + "_" + schannel + "_" + Form("%d",int(VBFtag)) + ".txt";
  if(VBFtag==2) outfile = "CardFragments/ggzzMCFMBackgroundFit_" + ssqrts + "_" + schannel + ".txt";
  ofstream of(outfile,ios_base::out);


	//  setTDRStyle(false);
	gStyle->SetPadLeftMargin(0.16);

	TString filepath;
	if (sqrts==7) {
		filepath = "/scratch0/hep/usarical/HiggsLifetime/No_SIP/LHC_7TeV/";
		//		filepath = filePath7TeV;
	}
	else if (sqrts==8) {
		filepath = "/scratch0/hep/usarical/HiggsLifetime/No_SIP/LHC_8TeV/";
		//		filepath = filePath8TeV;
	}

  TChain* tree = new TChain("SelectedTree");
	tree->Add(filepath+ "/" + (schannel=="2e2mu" ? "2mu2e" : schannel) + "/HZZ4lTree_ggTo4mu_Contin-MCFM67_Reprocessed.root");
	tree->Add(filepath+ "/" + (schannel=="2e2mu" ? "2mu2e" : schannel) + "/HZZ4lTree_ggTo4e_Contin-MCFM67_Reprocessed.root");
	tree->Add(filepath+ "/" + (schannel=="2e2mu" ? "2mu2e" : schannel) + "/HZZ4lTree_ggTo2e2mu_Contin-MCFM67_Reprocessed.root");


	RooRealVar* MC_weight = new RooRealVar("MC_weight", "MC_weight", 0., 100.);
	RooRealVar* ZZMass = new RooRealVar("ZZMass", "ZZMass", 100., 1600.);
	RooRealVar* NJets30 = new RooRealVar("NJets30", "NJets30", 0., 100.);
	RooArgSet ntupleVarSet(*ZZMass, *NJets30, *MC_weight);
  RooDataSet *set = new RooDataSet("set","set",ntupleVarSet,WeightVar("MC_weight"));
	
	Float_t myMC, myMass, myMC_extra;
  Short_t myNJets;
  int nentries = tree->GetEntries();

	tree->SetBranchAddress("ZZMass", &myMass);
	tree->SetBranchAddress("MC_weight", &myMC);
	if (tree->GetBranchStatus("MC_weight_Kfactor")) tree->SetBranchAddress("MC_weight_Kfactor", &myMC_extra);
	tree->SetBranchAddress("NJets30", &myNJets);

	Float_t Lep1SIP=0, Lep2SIP=0, Lep3SIP=0, Lep4SIP=0;
	Float_t Lep1_Z1SIP=0, Lep2_Z1SIP=0, Lep3_Z1SIP=0, Lep4_Z1SIP=0;
	Float_t KalmanCandVtx_chi2=0;

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

  for(int i =0;i<nentries;i++) {
    tree->GetEntry(i);
    if(VBFtag==1 && myNJets<2)continue;
    if(VBFtag==0 && myNJets>1)continue;
		if (myMass<100 || myMass>=1600) continue;
		bool passSelection = testSelection(
			Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP,
			Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP,
			KalmanCandVtx_chi2,
			SelectionScheme
			);
		if (!passSelection) continue;

		ntupleVarSet.setRealValue("ZZMass", myMass);
		ntupleVarSet.setRealValue("MC_weight", myMC*myMC_extra);
		ntupleVarSet.setRealValue("NJets30", (double)myNJets);

		set->add(ntupleVarSet, myMC*myMC_extra);
	}

  //RooRealVar* ZZLD = new RooRealVar("ZZLD","ZZLD",0.,1.);
  //char cut[10];
  //sprintf(cut,"ZZLD>0.5");
  //RooDataSet* set = new RooDataSet("set","set",tree,RooArgSet(*ZZMass,*MC_weight,*ZZLD),cut,"MC_weight");

  double totalweight = 0.;
  for (int i=0 ; i<set->numEntries() ; i++) { 
    set->get(i) ; 
    totalweight += set->weight();
    //cout << CMS_zz4l_mass->getVal() << " = " << set->weight() << endl ; 
  } 
  cout << "nEntries: " << set->numEntries() << ", totalweight: " << totalweight << endl;
		
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
	
  //// ---------------------------------------
  //Background
  RooRealVar CMS_qqzzbkg_a0("CMS_qqzzbkg_a0","CMS_qqzzbkg_a0",115.3,0.,200.);
  RooRealVar CMS_qqzzbkg_a1("CMS_qqzzbkg_a1","CMS_qqzzbkg_a1",21.96,0.,200.);
  RooRealVar CMS_qqzzbkg_a2("CMS_qqzzbkg_a2","CMS_qqzzbkg_a2",122.8,0.,200.);
  RooRealVar CMS_qqzzbkg_a3("CMS_qqzzbkg_a3","CMS_qqzzbkg_a3",0.03479,0.,1.);
  RooRealVar CMS_qqzzbkg_a4("CMS_qqzzbkg_a4","CMS_qqzzbkg_a4",185.5,0.,200.);
  RooRealVar CMS_qqzzbkg_a5("CMS_qqzzbkg_a5","CMS_qqzzbkg_a5",12.67,0.,200.);
  RooRealVar CMS_qqzzbkg_a6("CMS_qqzzbkg_a6","CMS_qqzzbkg_a6",34.81,0.,100.);
  RooRealVar CMS_qqzzbkg_a7("CMS_qqzzbkg_a7","CMS_qqzzbkg_a7",0.1393,0.,1.);
  RooRealVar CMS_qqzzbkg_a8("CMS_qqzzbkg_a8","CMS_qqzzbkg_a8",66.,0.,200.);
  RooRealVar CMS_qqzzbkg_a9("CMS_qqzzbkg_a9","CMS_qqzzbkg_a9",0.07191,0.,1.);
	
  RooggZZPdf_v2* bkg_ggzz = new RooggZZPdf_v2("bkg_ggzz","bkg_ggzz",*ZZMass,
					      CMS_qqzzbkg_a0,CMS_qqzzbkg_a1,CMS_qqzzbkg_a2,CMS_qqzzbkg_a3,CMS_qqzzbkg_a4,
					      CMS_qqzzbkg_a5,CMS_qqzzbkg_a6,CMS_qqzzbkg_a7,CMS_qqzzbkg_a8,CMS_qqzzbkg_a9);
	
  //// ---------------------------------------
	
  RooFitResult *r1 = bkg_ggzz->fitTo( *set, Save(kTRUE), SumW2Error(kTRUE) );//, Save(kTRUE), SumW2Error(kTRUE)) ;

  cout << endl;
  cout << "------- Parameters for " << schannel << " sqrts=" << sqrts << endl;
  cout << "  a0_bkgd = " << CMS_qqzzbkg_a0.getVal() << endl;
  cout << "  a1_bkgd = " << CMS_qqzzbkg_a1.getVal() << endl;
  cout << "  a2_bkgd = " << CMS_qqzzbkg_a2.getVal() << endl;
  cout << "  a3_bkgd = " << CMS_qqzzbkg_a3.getVal() << endl;
  cout << "  a4_bkgd = " << CMS_qqzzbkg_a4.getVal() << endl;
  cout << "  a5_bkgd = " << CMS_qqzzbkg_a5.getVal() << endl;
  cout << "  a6_bkgd = " << CMS_qqzzbkg_a6.getVal() << endl;
  cout << "  a7_bkgd = " << CMS_qqzzbkg_a7.getVal() << endl;
  cout << "  a8_bkgd = " << CMS_qqzzbkg_a8.getVal() << endl;
  cout << "  a9_bkgd = " << CMS_qqzzbkg_a9.getVal() << endl;
  cout << "---------------------------" << endl << endl;  

  of << "ggZZshape a0_bkgd  " << CMS_qqzzbkg_a0.getVal() << endl;
  of << "ggZZshape a1_bkgd  " << CMS_qqzzbkg_a1.getVal() << endl;
  of << "ggZZshape a2_bkgd  " << CMS_qqzzbkg_a2.getVal() << endl;
  of << "ggZZshape a3_bkgd  " << CMS_qqzzbkg_a3.getVal() << endl;
  of << "ggZZshape a4_bkgd  " << CMS_qqzzbkg_a4.getVal() << endl;
  of << "ggZZshape a5_bkgd  " << CMS_qqzzbkg_a5.getVal() << endl;
  of << "ggZZshape a6_bkgd  " << CMS_qqzzbkg_a6.getVal() << endl;
  of << "ggZZshape a7_bkgd  " << CMS_qqzzbkg_a7.getVal() << endl;
  of << "ggZZshape a8_bkgd  " << CMS_qqzzbkg_a8.getVal() << endl;
  of << "ggZZshape a9_bkgd  " << CMS_qqzzbkg_a9.getVal() << endl;
  of << endl;
  of.close();

  cout << endl << "Output written to: " << outfile << endl;

  int iLineColor = 1;
  string lab = "blah";
  if (channel == 1) { iLineColor = 2; lab = "4#mu"; }
  if (channel == 3) { iLineColor = 4; lab = "2e2#mu"; }
  if (channel == 2) { iLineColor = 6; lab = "4e"; }
  char lname[192];
  sprintf(lname,"gg #rightarrow ZZ #rightarrow %s", lab.c_str() );
  char lname2[192];
  sprintf(lname2,"Shape Model, %s", lab.c_str() );
  // dummy!                                                                                                                                               
  TF1* dummyF = new TF1("dummyF","1",0.,1.);
  TH1F* dummyH = new TH1F("dummyH","",1, 0.,1.);
  dummyF->SetLineColor( iLineColor );
  dummyF->SetLineWidth( 2 );
  
  TLegend * box2 = new TLegend(0.5,0.70,0.90,0.90);
  box2->SetFillColor(0);
  box2->SetBorderSize(0);
  box2->AddEntry(dummyH,"Simulation (GG2ZZ)  ","pe");
  box2->AddEntry(dummyH,lname,"");
  box2->AddEntry(dummyH,"","");
  box2->AddEntry(dummyF,lname2,"l");

  TPaveText *pt = new TPaveText(0.15,0.955,0.4,0.99,"NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->AddText("CMS Preliminary 2012");
  TPaveText *pt2 = new TPaveText(0.84,0.955,0.99,0.99,"NDC");
  pt2->SetFillColor(0);
  pt2->SetBorderSize(0);

  // Plot m4l and 
	ZZMass->setBins((1600-100)/8);
	RooPlot* frameM4l = ZZMass->frame(Title("M4L"), Bins((1600-100)/8));
  set->plotOn(frameM4l, MarkerStyle(24)) ;
	TH1F* bkg_ggzz_histo = (TH1F*)bkg_ggzz->createHistogram(Form("%s_histo", bkg_ggzz->GetName()), *ZZMass);
	bkg_ggzz_histo->SetLineColor(iLineColor);
	bkg_ggzz_histo->Scale(totalweight/bkg_ggzz_histo->Integral());
	set->plotOn(frameM4l) ;

  //comaprison with different shape, if needed (uncommenting also the code above)
  //bkg_ggzz_bkgd->plotOn(frameM4l,LineColor(1),NormRange("largerange")) ;

  frameM4l->GetXaxis()->SetTitle("m_{4l} [GeV]");
  frameM4l->GetYaxis()->SetTitle("a.u.");
/*  frameM4l->GetYaxis()->SetRangeUser(0,0.03);
  if(channel == 3)frameM4l->GetYaxis()->SetRangeUser(0,0.05);
  if(VBFtag<2){
    if(channel == 3)frameM4l->GetYaxis()->SetRangeUser(0,0.01);
    else frameM4l->GetYaxis()->SetRangeUser(0,0.005);
  }
*/  frameM4l->GetXaxis()->SetRangeUser(100,600);
  TCanvas *c = new TCanvas("c","c",800,600);
  c->cd();
  frameM4l->Draw();
	bkg_ggzz_histo->Draw("same");
  box2->Draw();
  pt->Draw();
  pt2->Draw();

  TString outputPath = "bkgFigs";
  outputPath = outputPath+ (long) sqrts + "TeV/";
  TString outputName;
  if(VBFtag<2) outputName =  outputPath + "bkgMCFMggzz_" + schannel + "_" + Form("%d",int(VBFtag));
  if(VBFtag==2) outputName =  outputPath + "bkgMCFMggzz_" + schannel;
  c->SaveAs(outputName + ".eps");
  c->SaveAs(outputName + ".png");
  c->SaveAs(outputName + ".root");
  delete c;
} 

