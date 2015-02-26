/* 
 * Prepare root files containing data events.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b prepareData.C+
 * This creates all files for 3 final states, 7 and 8 TeV and stores them in the final destination directory
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
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "TTree.h"
#include "TText.h"
#include "TStyle.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#endif

using namespace std;

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

//--- Flags to control re-computation of KD
bool recompute_ = false;       // Recompute KD (instead of taking it from the tree); when true, the following flags apply:
bool usePowhegTemplate=false;  // false use analytic bg
bool withPt_ = false;          // Include pT in KD
bool withY_  = false;          //    "    Y  "  "
int sqrts    = 8;              // sqrts, used only for withPt_/withY_

bool onlyICHEPStat = false;
bool dumpEvents = true;

int convertTreeForDatacards(TString inFile, TString outfile, bool useJET, bool VBFtag, int CoM);

// Run all final states and sqrts in one go
void prepareData_Lifetime_BS() {

  string appendDataName = "_Lifetime_BS";
  gSystem->Exec("mkdir -p "+ DataRootFilePath + appendDataName);
  Int_t nOut[18];

  nOut[0]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath + appendDataName+"/hzz4mu_"  +lumistr7TeV+"_1.root",true, true, 7);
  nOut[1]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath + appendDataName+"/hzz4e_"   +lumistr7TeV+"_1.root", true, true, 7);
  nOut[2]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleOr.root", DataRootFilePath + appendDataName+"/hzz2e2mu_"+lumistr7TeV+"_1.root", true, true, 7);
  nOut[3]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleMu.root", DataRootFilePath + appendDataName+"/hzz4mu_"  +lumistr8TeV+"_1.root", true, true, 8);
  nOut[4]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath + appendDataName+"/hzz4e_"   +lumistr8TeV+"_1.root", true, true, 8);
  nOut[5]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleOr.root", DataRootFilePath + appendDataName+"/hzz2e2mu_"+lumistr8TeV+"_1.root", true, true, 8);

  nOut[6]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleMu.root", DataRootFilePath + appendDataName+"/hzz4mu_"  +lumistr7TeV+"_0.root", true, false, 7);
  nOut[7]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath + appendDataName+"/hzz4e_"   +lumistr7TeV+"_0.root", true, false, 7);
  nOut[8]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleOr.root", DataRootFilePath + appendDataName+"/hzz2e2mu_"+lumistr7TeV+"_0.root", true, false, 7);
  nOut[9]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleMu.root", DataRootFilePath + appendDataName+"/hzz4mu_"  +lumistr8TeV+"_0.root", true, false, 8);
  nOut[10]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath + appendDataName+"/hzz4e_"   +lumistr8TeV+"_0.root", true, false, 8);
  nOut[11]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleOr.root", DataRootFilePath + appendDataName+"/hzz2e2mu_"+lumistr8TeV+"_0.root", true, false, 8);

  nOut[12]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleMu.root", DataRootFilePath + appendDataName+"/hzz4mu_"  +lumistr7TeV+".root", false, false, 7);
  nOut[13]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath + appendDataName+"/hzz4e_"   +lumistr7TeV+".root", false, false, 7);
  nOut[14]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleOr.root", DataRootFilePath + appendDataName+"/hzz2e2mu_"+lumistr7TeV+".root", false, false, 7);
  nOut[15]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleMu.root", DataRootFilePath + appendDataName+"/hzz4mu_"  +lumistr8TeV+".root", false, false, 8);
  nOut[16]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath + appendDataName+"/hzz4e_"   +lumistr8TeV+".root", false, false, 8);
  nOut[17]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleOr.root", DataRootFilePath + appendDataName+"/hzz2e2mu_"+lumistr8TeV+".root", false, false, 8);

  TString jetsT[3]={"DiJet     ",
		    "Untagged  ",
		    "Inclusive "};
  TString energies[2]={"7TeV", "8TeV "};
  TString channels[3]={"4mu   ",
		       "4e    ",
		       "2e2mu "};
  cout<<"Summary  "<<endl;
  for(int ijet=0;ijet<3;ijet++)
    for(int ien = 0;ien<2;ien++)
      for(int ich=0;ich<3;ich++)
	cout<<"N Events out for "<<jetsT[ijet]<<energies[ien].Data()<<channels[ich].Data()<<": "<<nOut[ijet*6+ien*3+ich]<<endl;
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
    if (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5) pass = false;
    break;
  case 3:
    if (KalmanCandVtx_chi2>=30) pass = false;
    break;
  case 4:
    if (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2>=30) pass = false;
    break;
  default:
    break;
  }
  return pass;
}


// The actual job
int convertTreeForDatacards(TString inFile, TString outfile, bool useJET, bool VBFtag, int CoM){
  TString comstring;
  comstring.Form("%iTeV", CoM);
  int EnergyIndex=0; if (CoM==8) EnergyIndex=1;

  TChain* treedata ;
  treedata= new TChain("SelectedTree");
  treedata->Add(inFile);

  int neventOut=0;
  Int_t run;
  Int_t ls=0;
  Long64_t event;
  float mzz, mzzErr;
  float m1, m2;
  float pt4l, fisher;
  int NJets_30;

  float ZZMass, ZZPt, ZZEta, ZZPhi;
  float Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP;
  float Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP;

  float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
  float OfflinePrimaryVtx_x, OfflinePrimaryVtx_y, OfflinePrimaryVtx_z;
  float CandVtx_x, CandVtx_y, CandVtx_z;

  float p0plus_VAJHU;
  float bkg_VAMCFM;
  float p0plus_m4l;
  float bkg_m4l;

  double Dxy=-99;
  double Txy=-99;
  double KD=-99;
  double D_bkg=-99;

  // BeamSpot info
  int RunNumber_Ref, LumiNumber_Ref;
  float BeamPosX, BeamPosY, BeamPosZ, BeamPosXErr, BeamPosYErr, BeamPosZErr;
  TString strBeamSpot[2]={
    "data_GR_R_44_V15C",
    "data_FT_53_V21_AN4"
  };
  TChain* tBeam = new TChain("BeamSpotRecord");
  TString cinput_BS = "BeamSpotRecord_" + strBeamSpot[EnergyIndex] + "_" + comstring + ".root";
  tBeam->Add(cinput_BS);
  tBeam->SetBranchAddress("RunNumber", &RunNumber_Ref);
  tBeam->SetBranchAddress("LumiNumber", &LumiNumber_Ref);
  tBeam->SetBranchAddress("BeamPosX", &BeamPosX);
  tBeam->SetBranchAddress("BeamPosY", &BeamPosY);
  tBeam->SetBranchAddress("BeamPosZ", &BeamPosZ);
  tBeam->SetBranchAddress("BeamPosXErr", &BeamPosXErr);
  tBeam->SetBranchAddress("BeamPosYErr", &BeamPosYErr);
  tBeam->SetBranchAddress("BeamPosZErr", &BeamPosZErr);

  treedata->SetBranchAddress("RunNumber",&run);
  treedata->SetBranchAddress("LumiNumber",&ls);
  treedata->SetBranchAddress("EventNumber",&event);
  treedata->SetBranchAddress("ZZMass", &ZZMass);
  treedata->SetBranchAddress("ZZPt", &ZZPt);
  treedata->SetBranchAddress("ZZEta", &ZZEta);
  treedata->SetBranchAddress("ZZPhi", &ZZPhi);
  treedata->SetBranchAddress("ZZMassErrCorr", &mzzErr);

  treedata->SetBranchAddress("ZZFisher",&fisher);
  treedata->SetBranchAddress("NJets30", &NJets_30);

  //Spin Disc variables
			treedata->SetBranchAddress("p0plus_VAJHU",&p0plus_VAJHU);
			treedata->SetBranchAddress("bkg_VAMCFM",&bkg_VAMCFM);
			treedata->SetBranchAddress("p0plus_m4l",&p0plus_m4l);
			treedata->SetBranchAddress("bkg_m4l",&bkg_m4l);

      treedata->SetBranchAddress("Lep1SIP", &Lep1SIP);
      treedata->SetBranchAddress("Lep2SIP", &Lep2SIP);
      treedata->SetBranchAddress("Lep3SIP", &Lep3SIP);
      treedata->SetBranchAddress("Lep4SIP", &Lep4SIP);
      treedata->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
      treedata->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
      treedata->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
      treedata->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
      treedata->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
      treedata->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
      treedata->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
      treedata->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
      treedata->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
      treedata->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
      treedata->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);

   
  TFile* newFile  = new TFile(outfile, "RECREATE");
  newFile->cd();
  TTree* newTree = new TTree("data_obs","data_obs"); 
  Double_t CMS_zz4l_mass, melaLD, CMS_zz4l_massErr, CMS_zz4l_massRelErr;
  Double_t pt = -99, Fisher = -99;

  const double cm_to_microns = 10000.;
  const double KD_scale = 1250.;

  newTree->Branch("CMS_zz4l_mass",&CMS_zz4l_mass,"CMS_zz4l_mass/D");
  newTree->Branch("CMS_zz4l_massErr",&CMS_zz4l_massErr,"CMS_zz4l_massErr/D");
  newTree->Branch("CMS_zz4l_massRelErr",&CMS_zz4l_massRelErr,"CMS_zz4l_massRelErr/D");
  newTree->Branch("CMS_zz4l_Fisher",&Fisher,"CMS_zz4l_Fisher/D");
  newTree->Branch("CMS_zz4l_Pt",&pt,"CMS_zz4l_Pt/D");

  //Spin Analysis Branches
			newTree->Branch("CMS_zz4l_KD",&KD);
			newTree->Branch("CMS_zz4l_smd",&D_bkg);


  cout << inFile << " entries: " << treedata->GetEntries() << endl;
  for(int iEvt=0; iEvt<treedata->GetEntries(); iEvt++){
    run = 1;
    ls = 1;

    treedata->GetEntry(iEvt);

    bool passSelection = testSelection(
      Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP,
      Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP,
      KalmanCandVtx_chi2,
      4
//      scheme
      );
    if (!passSelection) continue;

    //cout << run << endl;
    
    if (onlyICHEPStat && run>=198049) continue;

    //if( NJets_30 > 0) cout << "NJets_30: " << NJets_30 << endl;

    if ((useJET && VBFtag && NJets_30 < 2) || (useJET && !VBFtag && NJets_30 >= 2)) continue;

    mzz = ZZMass;
    pt4l = ZZPt;

    CMS_zz4l_mass = double(mzz);
    CMS_zz4l_massErr = double(mzzErr);
    CMS_zz4l_massRelErr = double(mzzErr/mzz);

    if(useJET && !VBFtag) 
      {
	if(pt4l > 200.)
	  {
	    //if pt out of range set to middle of higest bin
	    pt4l = 198.;
	  }
	pt = double(pt4l);
      }
    if(useJET && VBFtag) 
      {
	if(fisher > 2.)
	  {
	    //if fisher out of range set to middle of higest bin
	    fisher = 1.98;
	  }
	Fisher = double(fisher);
      }

			D_bkg = p0plus_VAJHU*p0plus_m4l / ( p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l );
/*
      CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
      CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
      CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;
*/
      bool matchFound = false;
      for (int evBS=0; evBS < tBeam->GetEntries(); evBS++){
        tBeam->GetEntry(evBS);
        if (run == RunNumber_Ref && ls == LumiNumber_Ref){
          matchFound = true;
          break;
        }
        else continue;
      }
      if (!matchFound) cerr << "No beamspot matching possible!" << endl;

      CandVtx_x = KalmanCandVtx_x - BeamPosX;
      CandVtx_y = KalmanCandVtx_y - BeamPosY;
      CandVtx_z = KalmanCandVtx_z - BeamPosZ;

      Dxy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi)) * cm_to_microns;
      Txy = Dxy*ZZMass / ZZPt;
      KD = Txy;
      KD = tanh(KD / KD_scale); // Scaled to due to tanh function
      if (KD >= 1.) KD = 1. - 1.0e-4;
      if (KD <= -1.) KD = -1. + 1.0e-4;

    ++neventOut;
    newTree->Fill();
    if (dumpEvents) {
      if (!useJET && !VBFtag) {
	cout.setf(ios::fixed);
	int ifs=-1;
	if (outfile.Contains("4e")) ifs = 0;
	else if (outfile.Contains("4mu")) ifs = 1;
	else if (outfile.Contains("2e2mu")) ifs = 2;
	cout << run << ":" << ls << ":" << event << " " << ifs << " " << setprecision(2) << CMS_zz4l_mass << " " << m1 << " " << m2 << " " << CMS_zz4l_massErr << " "  <<  setprecision(3) << melaLD << " "  <<  setprecision(2)  << ((pt4l>0)?pt4l:-1) << " " << setprecision(3) << ((NJets_30 >= 2)?fisher:-1) << endl;
      }
    }
  }

  delete tBeam;

  newTree->Write("data_obs"); 
  newFile->Close();
  
  cout << "written: " << outfile << " entries: " << neventOut << endl << endl;
  return  neventOut;
  
}


