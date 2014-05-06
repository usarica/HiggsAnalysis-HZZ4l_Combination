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


void modifySigTemplates(TString dir, TString sqrts = "7TeV", TString channame = "2e2mu"){
	TString cinputNominal = dir + "/" + sqrts + "/" + channame + "_templates_Nominal.root";
	TString cinputResUp = dir + "/" + sqrts + "/" + channame + "_templates_ResUp.root";
	TString cinputScaleUp = dir + "/" + sqrts + "/" + channame + "_templates_ScaleUp.root";
	TString cinputResDown = dir + "/" + sqrts + "/" + channame + "_templates_ResDown.root";
	TString cinputScaleDown = dir + "/" + sqrts + "/" + channame + "_templates_ScaleDown.root";
	TString cinput[5] = {
		cinputNominal,
		cinputResUp,
		cinputResDown,
		cinputScaleUp,
		cinputScaleDown
	};

	TString coutput = dir + "/" + sqrts + "/" + channame + "_templates_Modified_Nominal_ScaleResUpDown.root";
	TFile* foutput = new TFile(coutput, "recreate");

	TFile* finput[5];
	TH3F* hsig[5 * 9];
	TH3F* hsig_scaleres[2 * 9];
	for (int f = 0; f < 5; f++){
		finput[f] = new TFile(cinput[f], "read");
		cout << "Reading file: " << finput[f]->GetName() << "..." << endl;

		TString ctemplate_main = "T_3D_";
		for (int t = 0; t < 9; t++){
			char tcode[2];
			sprintf(tcode, "%i", t+1);
			TString ctemplate = ctemplate_main + tcode;
			TString ctemplate_scaleres = ctemplate + "_ScaleRes";
			if ((f - 1) % 2 == 0) ctemplate_scaleres = ctemplate_scaleres + "Up";
			if ((f - 1) % 2 == 1) ctemplate_scaleres = ctemplate_scaleres + "Down";
			hsig[9 * f + t] = (TH3F*) finput[f]->Get(ctemplate);

			cout << "Receive template: " << hsig[9 * f + t]->GetName() << endl;

			if (f == 1 || f == 2){
				hsig_scaleres[9 * (f - 1) + t] = (TH3F*) hsig[9 * f + t]->Clone(ctemplate_scaleres);
				hsig_scaleres[9 * (f - 1) + t]->SetTitle(ctemplate_scaleres);
				cout << "Clone template " << hsig[9 * f + t]->GetName() << " as " << hsig_scaleres[9 * (f - 1) + t]->GetName() << endl;
			};
			if (f == 3 || f == 4){
				hsig_scaleres[9 * (f - 3) + t]->Add(hsig[9 * f + t], 1.0);
				cout << "Add template " << hsig[9 * f + t]->GetName() << " to " << hsig_scaleres[9 * (f - 3) + t]->GetName() << endl;
			};
		};
	};
	for (int t = 0; t < 9; t++){
		TString ccanvas = "validateSig_";
		char tcode[2];
		sprintf(tcode, "%i", t+1);
		ccanvas = ccanvas + tcode;

		foutput->WriteTObject(hsig[9 * 0 + t]);
		for (int f = 0; f < 2; f++){
			hsig_scaleres[9 * f + t]->Add(hsig[9 * 0 + t], -1.0);
			cout << "Add template (-1) x " << hsig[9 * 0 + t]->GetName() << " to " << hsig_scaleres[9 * f + t]->GetName() << endl;
			foutput->WriteTObject(hsig_scaleres[9 * f + t]);
		};


		TString ccanvas_xy = ccanvas + "_xy";
		TCanvas* c = new TCanvas(ccanvas_xy,"",1600,1200);
		c->Divide(2,2);
		c->cd(1);
		hsig[9 * 0 + t]->Project3D("xy")->Draw("colz");
		c->cd(2);
		hsig_scaleres[9 * 0 + t]->Project3D("xy")->Draw("colz");
		c->cd(3);
		hsig_scaleres[9 * 1 + t]->Project3D("xy")->Draw("colz");
		foutput->WriteTObject(c);
		delete c;
		TString ccanvas_xz = ccanvas + "_xz";
		c = new TCanvas(ccanvas_xz,"",1600,1200);
		c->Divide(2,2);
		c->cd(1);
		hsig[9 * 0 + t]->Project3D("xz")->Draw("colz");
		c->cd(2);
		hsig_scaleres[9 * 0 + t]->Project3D("xz")->Draw("colz");
		c->cd(3);
		hsig_scaleres[9 * 1 + t]->Project3D("xz")->Draw("colz");
		foutput->WriteTObject(c);
		delete c;
		TString ccanvas_yz = ccanvas + "_yz";
		c = new TCanvas(ccanvas_yz,"",1600,1200);
		c->Divide(2,2);
		c->cd(1);
		hsig[9 * 0 + t]->Project3D("yz")->Draw("colz");
		c->cd(2);
		hsig_scaleres[9 * 0 + t]->Project3D("yz")->Draw("colz");
		c->cd(3);
		hsig_scaleres[9 * 1 + t]->Project3D("yz")->Draw("colz");
		foutput->WriteTObject(c);
		delete c;
	};

	foutput->Close();
	for (int f = 0; f < 5; f++) finput[f]->Close();
}


void modifyBkgTemplates (TString dir, TString sqrts = "7TeV", TString channame = "2e2mu"){
	TString cinput = dir + "/" + sqrts + "/" + channame + "_templates_Nominal_bkg.root";
	TString coutput = dir + "/" + sqrts + "/" + channame + "_templates_Modified_Nominal_bkg.root";
	
	TFile* finput = new TFile(cinput,"read");
	TFile* foutput = new TFile(coutput,"recreate");
	TH3F* ggzz = (TH3F*) finput->Get("template_ggZZ");
	TH3F* qqzz = (TH3F*) finput->Get("template_qqZZ");
	TH3F* zx= (TH3F*) finput->Get("template_ZX");

	qqzz->Scale(1.0/qqzz->Integral());
	zx->Scale(1.0/zx->Integral());

	TH3F* zx_up = (TH3F*) qqzz->Clone("template_ZX_Up");
	TH3F* zx_down = (TH3F*) zx->Clone("template_ZX_Down");
	zx_down->Multiply(zx);
	zx_down->Divide(qqzz);
	zx_down->Scale(1.0/zx_down->Integral());

	foutput->WriteTObject(ggzz);
	foutput->WriteTObject(qqzz);
	foutput->WriteTObject(zx);
	foutput->WriteTObject(zx_up);
	foutput->WriteTObject(zx_down);

	TString ccanvas = "validateBkg";
	TString ccanvas_xy = ccanvas + "_xy";
	TCanvas* c = new TCanvas(ccanvas_xy,"",1600,1200);
	c->Divide(2,2);
	c->cd(2);
	gStyle->SetPalette(1);
	zx->SetLineColor(1);
	zx->Project3D("xy")->Draw("colz");
	c->cd(1);
	qqzz->SetLineColor(2);
	qqzz->Project3D("xy")->Draw("colz");
	c->cd(3);
	zx_up->SetLineColor(4);
	zx_up->Project3D("xy")->Draw("colz");
	c->cd(4);
	zx_down->SetLineColor(4);
	zx_down->Project3D("xy")->Draw("colz");
	foutput->WriteTObject(c);
	delete c;
	TString ccanvas_xz = ccanvas + "_xz";
	c = new TCanvas(ccanvas_xz,"",1600,1200);
	c->Divide(2,2);
	c->cd(2);
	gStyle->SetPalette(1);
	zx->SetLineColor(1);
	zx->Project3D("xz")->Draw("colz");
	c->cd(1);
	qqzz->SetLineColor(2);
	qqzz->Project3D("xz")->Draw("colz");
	c->cd(3);
	zx_up->SetLineColor(4);
	zx_up->Project3D("xz")->Draw("colz");
	c->cd(4);
	zx_down->SetLineColor(4);
	zx_down->Project3D("xz")->Draw("colz");
	foutput->WriteTObject(c);
	delete c;
	TString ccanvas_yz = ccanvas + "_yz";
	c = new TCanvas(ccanvas_yz,"",1600,1200);
	c->Divide(2,2);
	c->cd(2);
	gStyle->SetPalette(1);
	zx->SetLineColor(1);
	zx->Project3D("yz")->Draw("colz");
	c->cd(1);
	qqzz->SetLineColor(2);
	qqzz->Project3D("yz")->Draw("colz");
	c->cd(3);
	zx_up->SetLineColor(4);
	zx_up->Project3D("yz")->Draw("colz");
	c->cd(4);
	zx_down->SetLineColor(4);
	zx_down->Project3D("yz")->Draw("colz");
	foutput->WriteTObject(c);
	delete c;

	foutput->Close();
	finput->Close();
}

void modifyTemplates(TString dir){
	TString channame[3] = { "4mu","4e","2e2mu" };
	TString sqrts[2] = { "7TeV","8TeV" };
	for (int s = 0; s < 2; s++){
		for (int f = 0; f < 3; f++){
			modifyBkgTemplates(dir, sqrts[s], channame[f]);
			modifySigTemplates(dir, sqrts[s], channame[f]);
		};
	};
}