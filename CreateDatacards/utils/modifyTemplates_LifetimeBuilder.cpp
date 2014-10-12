#include <iostream>
#include <cmath>
#include <string>
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

const int nCTau=41;
const int maxCTau=1000;


void modifySigTemplates(TString dir, TString sqrts = "7TeV", TString channame = "2e2mu"){
	TString cinput_Dbkg = dir + "/" + sqrts + "/" + channame + "_templates_Modified_Nominal_ScaleResUpDown.root";
	TString cinput_Txymain = dir + "/" + sqrts + "/" + channame + "_templates_TxyUpDown_CTau";
	TString coutput_Dbkg = dir + "/" + sqrts + "/" + channame + "_templates_SignalScaleResSyst.root";

	TFile* finput_Dbkg = new TFile(cinput_Dbkg, "read");
	if(finput_Dbkg->IsZombie()) return;
	TFile* foutput = new TFile(coutput_Dbkg, "recreate");

	const int kNumScaleResSysts=3;
	const int kNumTxySysts=3;
	char* cScaleResSyst[kNumScaleResSysts] = {
		"_ScaleResNominal", "_ScaleResUp", "_ScaleResDown"
	};
	TH1F* hDbkg[3];
	TH3F* hDbkg_3D[3] = {
		(TH3F*)finput_Dbkg->Get("T_3D_1"),
		(TH3F*)finput_Dbkg->Get("T_3D_1_ScaleResUp"),
		(TH3F*)finput_Dbkg->Get("T_3D_1_ScaleResDown")
	};
	foutput->cd();
	for (int ss = 0; ss < kNumScaleResSysts; ss++){
		hDbkg[ss] = (TH1F*)hDbkg_3D[ss]->ProjectionZ();
		hDbkg[ss]->SetName(Form("T_2D%s", cScaleResSyst[ss]));
		hDbkg[ss]->Scale(1./hDbkg[ss]->Integral());
		foutput->WriteTObject(hDbkg[ss]);
	}
	foutput->Close();
	finput_Dbkg->Close();

	for (int ct = 0; ct < nCTau; ct++){
		int genctau = ct*maxCTau / (nCTau - 1);
		TString cinput = cinput_Txymain;
		TString coutput = cinput_Txymain;
		cinput.Append(Form("%i%s", genctau, ".root"));
		coutput.Append(Form("%i%s", genctau, "_Modified.root"));
		TFile* finput = new TFile(cinput,"read");
		if(finput->IsZombie()) return;
		foutput = new TFile(coutput,"recreate");

		TH1F* h_1D[3];
		TH2F* h_2D[3] = {
			(TH2F*)finput->Get("T_2D"),
			(TH2F*)finput->Get("T_2D_TxyUp"),
			(TH2F*)finput->Get("T_2D_TxyDown")
		};
		h_2D[0]->SetName("T_2D_TxyNominal");
		for (int ss = 0; ss < kNumTxySysts; ss++){
			TString templateName = h_2D[ss]->GetName();
			h_2D[ss]->SetName(Form("%s_Original",h_2D[ss]->GetName()));
			h_1D[ss] = (TH1F*)h_2D[ss]->ProjectionX();
			h_1D[ss]->SetName(templateName);
			if(ss>0) h_1D[ss]->Scale(h_1D[0]->Integral()/h_1D[ss]->Integral());
			foutput->WriteTObject(h_1D[ss]);
		}

		foutput->Close();
		finput->Close();
	}
}


void modifyBkgTemplates (TString dir, TString sqrts = "7TeV", TString channame = "2e2mu"){
	TString cinput_TxyNominal = dir + "/" + sqrts + "/" + channame + "_templates_Nominal_bkg.root";
	TString cinput_TxyUp = dir + "/" + sqrts + "/" + channame + "_templates_TxyUp_bkg.root";
	TString cinput_TxyDown = dir + "/" + sqrts + "/" + channame + "_templates_TxyDown_bkg.root";
	TString cinput_Dbkg = dir + "/" + sqrts + "/" + channame + "_templates_Modified_Nominal_bkg.root";
	TString coutput = dir + "/" + sqrts + "/" + channame + "_templates_Merged_bkg.root";
	
	TFile* finput[4] = {
		new TFile(cinput_Dbkg, "read"),
		new TFile(cinput_TxyNominal, "read"),
		new TFile(cinput_TxyUp, "read"),
		new TFile(cinput_TxyDown, "read")
	};
	if (finput[1]->IsZombie()){ finput[0]->Close(); return; }
	TFile* foutput = new TFile(coutput,"recreate");
	TH3F* ggzz_dbkg = (TH3F*) finput[0]->Get("template_ggZZ");
	TH3F* qqzz_dbkg = (TH3F*) finput[0]->Get("template_qqZZ");
	TH3F* zx_dbkg[3] = {
		(TH3F*)finput[0]->Get("template_ZX"),
		(TH3F*)finput[0]->Get("template_ZX_Up"),
		(TH3F*)finput[0]->Get("template_ZX_Down")
	};
	cout << ggzz_dbkg->GetName() << endl;
	cout << ggzz_dbkg->Integral() << endl;
	ggzz_dbkg->Scale(1./ggzz_dbkg->Integral());
	cout << qqzz_dbkg->GetName() << endl;
	cout << qqzz_dbkg->Integral() << endl;
	qqzz_dbkg->Scale(1./qqzz_dbkg->Integral());
	for (int ii = 0; ii < 3; ii++){
		cout << zx_dbkg[ii]->GetName() << endl;
		cout << zx_dbkg[ii]->Integral() << endl;
		zx_dbkg[ii]->Scale(1. / zx_dbkg[ii]->Integral());
	}
	TH1F* ggzz_1d_dbkg = (TH1F*) ggzz_dbkg->ProjectionZ();
	ggzz_1d_dbkg->SetName(Form("%s_Dbkg",ggzz_dbkg->GetName()));
	TH1F* qqzz_1d_dbkg = (TH1F*) qqzz_dbkg->ProjectionZ();
	qqzz_1d_dbkg->SetName(Form("%s_Dbkg",qqzz_dbkg->GetName()));
	TH1F* zx_1d_dbkg[3];
	zx_1d_dbkg[0] = (TH1F*) zx_dbkg[0]->ProjectionZ();
	zx_1d_dbkg[0]->SetName(Form("%s_Dbkg_Nominal",zx_dbkg[0]->GetName()));
	zx_1d_dbkg[1] = (TH1F*) zx_dbkg[1]->ProjectionZ();
	zx_1d_dbkg[1]->SetName(Form("%s_Dbkg_ScaleResUp",zx_dbkg[0]->GetName()));
	zx_1d_dbkg[2] = (TH1F*) zx_dbkg[2]->ProjectionZ();
	zx_1d_dbkg[2]->SetName(Form("%s_Dbkg_ScaleResDown",zx_dbkg[0]->GetName()));

	foutput->WriteTObject(ggzz_1d_dbkg);
	foutput->WriteTObject(qqzz_1d_dbkg);
	for (int ii = 0; ii < 3; ii++) foutput->WriteTObject(zx_1d_dbkg[ii]);

	TH2F* ggzz_txy[3];
	TH2F* qqzz_txy[3];
	TH2F* zx_txy[3];
	TH1F* ggzz_1d_txy[3];
	TH1F* qqzz_1d_txy[3];
	TH1F* zx_1d_txy[3];
	for (int ii = 0; ii < 3; ii++){
		ggzz_txy[ii] = (TH2F*)finput[ii + 1]->Get("template_ggZZ");
		qqzz_txy[ii] = (TH2F*)finput[ii + 1]->Get("template_qqZZ");
		zx_txy[ii] = (TH2F*)finput[ii + 1]->Get("template_ZX");

		TString templateName=ggzz_txy[ii]->GetName();
		if(ii==0) templateName.Append("_TxyNominal");
		if(ii==1) templateName.Append("_TxyUp");
		if(ii==2) templateName.Append("_TxyDown");
		ggzz_1d_txy[ii] = (TH1F*) ggzz_txy[ii]->ProjectionX();
		ggzz_1d_txy[ii]->SetName(templateName);

		templateName=qqzz_txy[ii]->GetName();
		if(ii==0) templateName.Append("_TxyNominal");
		if(ii==1) templateName.Append("_TxyUp");
		if(ii==2) templateName.Append("_TxyDown");
		qqzz_1d_txy[ii] = (TH1F*) qqzz_txy[ii]->ProjectionX();
		qqzz_1d_txy[ii]->SetName(templateName);

		templateName=zx_txy[ii]->GetName();
		if(ii==0) templateName.Append("_TxyNominal");
		if(ii==1) templateName.Append("_TxyUp");
		if(ii==2) templateName.Append("_TxyDown");
		zx_1d_txy[ii] = (TH1F*) zx_txy[ii]->ProjectionX();
		zx_1d_txy[ii]->SetName(templateName);
	}
	for (int ii = 1; ii < 3; ii++){
		ggzz_1d_txy[ii]->Scale(ggzz_1d_txy[0]->Integral() / ggzz_1d_txy[ii]->Integral());
		qqzz_1d_txy[ii]->Scale(qqzz_1d_txy[0]->Integral() / qqzz_1d_txy[ii]->Integral());
		zx_1d_txy[ii]->Scale(zx_1d_txy[0]->Integral() / zx_1d_txy[ii]->Integral());
	}
	for (int ii = 0; ii < 3; ii++){
		foutput->WriteTObject(ggzz_1d_txy[ii]);
		foutput->WriteTObject(qqzz_1d_txy[ii]);
		foutput->WriteTObject(zx_1d_txy[ii]);
	}

	foutput->Close();
	for(int ii=0;ii<4;ii++) finput[ii]->Close();
}

void modifyTemplates_LifetimeBuilder(TString dir, int processSig=1, int processBkg=1){
	TString channame[3] = { "4mu","4e","2e2mu" };
	TString sqrts[2] = { "7TeV","8TeV" };
	for (int s = 0; s < 2; s++){
		for (int f = 0; f < 3; f++){
			if(processBkg==1) modifyBkgTemplates(dir, sqrts[s], channame[f]);
			if(processSig==1) modifySigTemplates(dir, sqrts[s], channame[f]);
		};
	};
}