#include <iostream>
#include <cmath>
#include <string>
#include "TSystem.h"
#include "TROOT.h"
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

const int nCTau=101;
const int maxCTau=1000;


void extractSignalDbkg(TH1F* hsig[5]){
	const int kNumTemplates=1;
	const int kNumSysts=5;

	TH1F* hsig_scaleres[(kNumSysts-3) * kNumTemplates];
	TString ctemplate_main = "T_2D_";
	for (int f = 0; f < kNumSysts; f++){
		for (int t = 0; t < kNumTemplates; t++){
			char tcode[2];
			sprintf(tcode, "%i", t+1);
			TString ctemplate = ctemplate_main + tcode;
			TString ctemplate_scaleres = ctemplate + "_ScaleRes";
			if ((f - 1) % 2 == 0) ctemplate_scaleres = ctemplate_scaleres + "Up";
			if ((f - 1) % 2 == 1) ctemplate_scaleres = ctemplate_scaleres + "Down";

			if (f == 1){
				hsig_scaleres[kNumTemplates * (f - 1) + t] = (TH1F*)hsig[kNumTemplates * f + t]->Clone(ctemplate_scaleres);
				hsig_scaleres[kNumTemplates * (f - 1) + t]->SetTitle(ctemplate_scaleres);
			}
			if (f == 2){
				hsig_scaleres[kNumTemplates * (f - 1) + t] = (TH1F*)hsig[kNumTemplates * (f-1) + t]->Clone(ctemplate_scaleres);
				hsig_scaleres[kNumTemplates * (f - 1) + t]->SetTitle(ctemplate_scaleres);
			}
			if (f == 3){
				hsig_scaleres[kNumTemplates * (f - 3) + t]->Add(hsig[kNumTemplates * f + t], 1.0);
			}
			if (f == 4){
				hsig_scaleres[kNumTemplates * (f - 3) + t]->Add(hsig[kNumTemplates * (f-1) + t], 1.0);
			}
		}
	}
	for (int t = 0; t < kNumTemplates; t++){
		for (int f = 0; f < 1; f++){
			hsig_scaleres[kNumTemplates * f + t]->Add(hsig[kNumTemplates * 0 + t], -1.0);
		}
		for (int f = 1; f < (kNumSysts-3); f++){
			hsig_scaleres[kNumTemplates * f + t]->Scale(-1.0);
			hsig_scaleres[kNumTemplates * f + t]->Add(hsig[kNumTemplates * 0 + t], 3.0);
		}
		for (int f = 0; f < (kNumSysts-3); f++){
			for (int binx=1; binx<=hsig_scaleres[kNumTemplates * f + t]->GetNbinsX(); binx++){
				double bincontent = hsig_scaleres[kNumTemplates * f + t]->GetBinContent(binx);
				if (bincontent<=0) bincontent=1.0e-20;
				hsig[kNumTemplates * (f+1) + t]->SetBinContent(binx,bincontent);
			}
			hsig[kNumTemplates * (f+1) + t]->Scale(1./hsig[kNumTemplates * (f+1) + t]->Integral());
		}
	}
}


void regularizeCTauSlice(TH1F* hSlice);

TGraph* regularizeSignalNominal(const int nbins, double xx[], double yy[], int nIter, int iProd);

void regularizeSignalSystRatio(const int nbins, double xx[], double yy[], double yy_up[], double yy_dn[], int nIter);


void modifySigTemplates(TString dir, TString sqrts = "7TeV", TString channame = "2e2mu"){
	TString cinput_Dbkg = dir + "/" + sqrts + "/" + channame + "_templates_TxyUpDown_CTau0.root";
	TString coutput_Dbkg = dir + "/" + sqrts + "/" + channame + "_templates_SignalScaleResSyst.root";
  TString cinput_Txymain = dir + "/" + sqrts + "/" + channame + "_templates_";
  const int nProd = 5;
  TString cProduction[nProd]={ "", "VBFH", "WH", "ZH", "ttH" };
  double prodXSEC[nProd][2]={
    { 14.99, 19.09 },
    { 1.214, 1.572 },
    { 0.5688, 0.6931 },
    { 0.3299, 0.4091 },
    { 0.08508, 0.1274 }
  };

  int EnergyIndex = 0;
  if (sqrts.Contains("8TeV")) EnergyIndex=1;
  double sum_XSEC=0;
  for (int iProd=0; iProd<nProd; iProd++) sum_XSEC += prodXSEC[iProd][EnergyIndex];

	TFile* finput_Dbkg = new TFile(cinput_Dbkg, "read");
	if(finput_Dbkg->IsZombie()) return;
	TFile* foutput = new TFile(coutput_Dbkg, "recreate");

	cout << "File are opened" << endl;

	const int kNumScaleResSysts=3;
	const int kNumScaleResSystsInput=5;
	const int kNumTxySysts=3;
	char* cScaleResSystInput[kNumScaleResSystsInput] ={
		"_Nominal", "_ResUp", "_ResDown", "_ScaleUp", "_ScaleDown"
	};
	char* cScaleResSyst[kNumScaleResSysts] ={
		"_ScaleResNominal", "_ScaleResUp", "_ScaleResDown"
	};
	TH1F* hDbkg[kNumScaleResSystsInput];
	TH2F* hDbkg_2D[kNumScaleResSystsInput];
	for (int s=0; s<kNumScaleResSystsInput; s++){
		TString hname = "T_2D_Dbkg";
		if (s>0) hname.Append(cScaleResSystInput[s]);
		hDbkg_2D[s] = (TH2F*)finput_Dbkg->Get(hname);
		cout << hname << endl;
		hDbkg[s] = (TH1F*)hDbkg_2D[s]->ProjectionX();
	}
	cout << "Extracting signal bkg." << endl;
	extractSignalDbkg(hDbkg);
	cout << "Extracted!" << endl;
	foutput->cd();
	for (int ss = 0; ss < kNumScaleResSysts; ss++){
		hDbkg[ss]->SetName(Form("T_2D%s", cScaleResSyst[ss]));
		hDbkg[ss]->Scale(1./hDbkg[ss]->Integral());
		foutput->WriteTObject(hDbkg[ss]);
	}
	foutput->Close();
	finput_Dbkg->Close();
	cout << "Signal bkg. recorded." << endl;

  double xx[nCTau];
  double yy_all[nProd][100][nCTau];
  double yy_all_up[nProd][100][nCTau];
  double yy_all_dn[nProd][100][nCTau];
  int nbins_KD;

  TH1F* h_KDarray[nCTau][nProd*kNumTxySysts];
  double totalyield[nProd] ={ 0 };

  for (int iProd=0; iProd<nProd; iProd++){
    for (int ct = 0; ct < nCTau; ct++){
      int genctau = ct*maxCTau / (nCTau - 1);
      xx[ct] = (double)genctau;

      TString cinput = cinput_Txymain;
      if (iProd>0) cinput.Append(Form("%s_", cProduction[iProd].Data()));
      cinput.Append("TxyUpDown_CTau");
      cinput.Append(Form("%i%s", genctau, ".root"));
      TFile* finput = new TFile(cinput, "read");
      if (finput->IsZombie()) return;

      TH2F* h_2D[3] ={
        (TH2F*)finput->Get("T_2D"),
        (TH2F*)finput->Get("T_2D_TxyUp"),
        (TH2F*)finput->Get("T_2D_TxyDown")
      };
      h_2D[0]->SetName("T_2D_TxyNominal");
      TH1F* h_1D[kNumTxySysts];
      for (int ss = 0; ss < kNumTxySysts; ss++){
        TString templateName = h_2D[ss]->GetName();
        h_2D[ss]->SetName(Form("%s_Original", h_2D[ss]->GetName()));
        h_1D[ss] = (TH1F*)h_2D[ss]->ProjectionX();
        h_1D[ss]->SetName(templateName);
        if (ss>0) h_1D[ss]->Scale(h_1D[0]->Integral()/h_1D[ss]->Integral());
//        if(iProd>0) regularizeCTauSlice(h_1D[ss]);
        if (ss==0 && ct==0){
          if (iProd==0) nbins_KD = h_1D[ss]->GetNbinsX();
          totalyield[iProd] = h_1D[ss]->Integral();
        }
        for (int bin=1; bin<=nbins_KD; bin++){
          double bincontent = h_1D[ss]->GetBinContent(bin);
          if (ss==0) yy_all[iProd][bin-1][ct] = bincontent;
          else if (ss==1) yy_all_up[iProd][bin-1][ct] = bincontent;
          else if (ss==2) yy_all_dn[iProd][bin-1][ct] = bincontent;
        }
      }
      for (int ss = 0; ss < kNumTxySysts; ss++){
        gROOT->cd();
        h_KDarray[ct][kNumTxySysts*iProd + ss] = (TH1F*)h_1D[ss]->Clone(Form("%s_CTau%i", h_1D[ss]->GetName(), genctau));
        finput->cd();
        delete h_1D[ss];
        delete h_2D[ss];
      }
      finput->Close();
    }

    TString plotDir = dir + "/" + sqrts + "/";
    plotDir.Append("ValidationPlots/");
    TString mkdirCommand = "mkdir -p ";
    mkdirCommand.Append(plotDir);
    gSystem->Exec(mkdirCommand);

    gROOT->cd();
    gStyle->SetTitleFont(62, "t");
    gStyle->SetOptStat(0);
    gROOT->SetStyle(gStyle->GetName());
    gROOT->ForceStyle();

    for (int bin=0; bin<nbins_KD; bin++){
      double max_plot = 0;
      double yy[nCTau];
      double yy_up[nCTau];
      double yy_dn[nCTau];
      for (int ct=0; ct<nCTau; ct++){
        yy[ct] = yy_all[iProd][bin][ct];
        yy_up[ct] = yy_all_up[iProd][bin][ct];
        yy_dn[ct] = yy_all_dn[iProd][bin][ct];

        double max_yy = max(max(yy_up[ct], yy_dn[ct]), yy[ct]);
        if (max_plot<max_yy) max_plot = max_yy;
      }
      TString tg_name = Form("Bin%i_", bin+1);
      TString tg_up_name = Form("Bin%i_", bin+1);
      TString tg_dn_name = Form("Bin%i_", bin+1);
      if (iProd==0){
        tg_name.Append("ggH_");
        tg_up_name.Append("ggH_");
        tg_dn_name.Append("ggH_");
      }
      else{
        tg_name.Append(Form("%s_", cProduction[iProd].Data()));
        tg_up_name.Append(Form("%s_", cProduction[iProd].Data()));
        tg_dn_name.Append(Form("%s_", cProduction[iProd].Data()));
      }
      tg_name = tg_name + channame + "_" + sqrts + "_TxyNominal";
      tg_up_name = tg_up_name + channame + "_" + sqrts + "_TxyUp";
      tg_dn_name = tg_dn_name + channame + "_" + sqrts + "_TxyDown";

      gROOT->cd();
      TGraph* tg = new TGraph(nCTau, xx, yy);
      tg->SetName(tg_name);
      tg->SetLineColor(kBlack);
      TGraph* tg_up = new TGraph(nCTau, xx, yy_up);
      tg_up->SetName(tg_up_name);
      tg_up->SetLineColor(kRed);
      TGraph* tg_dn = new TGraph(nCTau, xx, yy_dn);
      tg_dn->SetName(tg_dn_name);
      tg_dn->SetLineColor(kBlue);
      TGraph* tg_all[3] ={ tg, tg_up, tg_dn };

      int nIter=20;
      TGraph* tg_fixed = regularizeSignalNominal(nCTau, xx, yy, nIter, iProd);
      tg_fixed->SetName(Form("%s_fixed", tg->GetName()));
      tg_fixed->SetLineColor(kBlack);
      tg_fixed->SetLineStyle(7);
      TGraph* tg_up_fixed = regularizeSignalNominal(nCTau, xx, yy_up, nIter, iProd);
      tg_up_fixed->SetName(Form("%s_fixed", tg_up->GetName()));
      tg_up_fixed->SetLineColor(kRed);
      tg_up_fixed->SetLineStyle(7);
      TGraph* tg_dn_fixed = regularizeSignalNominal(nCTau, xx, yy_dn, nIter, iProd);
      tg_dn_fixed->SetName(Form("%s_fixed", tg_dn->GetName()));
      tg_dn_fixed->SetLineColor(kBlue);
      tg_dn_fixed->SetLineStyle(7);
      TGraph* tg_fixed_all[3] ={ tg_fixed, tg_up_fixed, tg_dn_fixed };

      TGraph* tg_up_regratio;
      TGraph* tg_dn_regratio;
      regularizeSignalSystRatio(nCTau, xx, yy, yy_up, yy_dn, nIter/2);
      tg_up_regratio = new TGraph(nCTau, xx, yy_up);
      tg_up_regratio->SetName(Form("%s_regular", tg_up->GetName()));
      tg_up_regratio->SetLineColor(kRed);
      tg_up_regratio->SetLineStyle(2);
      tg_up_regratio->SetLineWidth(2);
      tg_dn_regratio = new TGraph(nCTau, xx, yy_dn);
      tg_dn_regratio->SetName(Form("%s_regular", tg_dn->GetName()));
      tg_dn_regratio->SetLineColor(kBlue);
      tg_dn_regratio->SetLineStyle(2);
      tg_dn_regratio->SetLineWidth(2);


      for (int ct=0; ct<nCTau; ct++){
        yy_all[iProd][bin][ct] = yy[ct];
        yy_all_up[iProd][bin][ct] = yy_up[ct];
        yy_all_dn[iProd][bin][ct] = yy_dn[ct];

        double max_yy = max(max(yy_up[ct], yy_dn[ct]), yy[ct]);
        if (max_plot<max_yy) max_plot = max_yy;
      }

      gROOT->cd();
      TString canvasname_2D = "cCanvas_";
      canvasname_2D.Append(Form("Bin%i_", bin+1));
      if (iProd==0) canvasname_2D.Append("ggH_");
      else canvasname_2D.Append(Form("%s_", cProduction[iProd].Data()));
      canvasname_2D = canvasname_2D + channame + "_" + sqrts + "_SignalTxyComparison";
      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
      c2D->cd();
      gStyle->SetOptStat(0);
      c2D->SetFillColor(0);
      c2D->SetBorderMode(0);
      c2D->SetBorderSize(2);
      c2D->SetTickx(1);
      c2D->SetTicky(1);
      c2D->SetLeftMargin(0.17);
      c2D->SetRightMargin(0.05);
      c2D->SetTopMargin(0.07);
      c2D->SetBottomMargin(0.13);
      c2D->SetFrameFillStyle(0);
      c2D->SetFrameBorderMode(0);
      c2D->SetFrameFillStyle(0);
      c2D->SetFrameBorderMode(0);

      TLegend *l2D = new TLegend(0.20, 0.57, 0.58, 0.90);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      for (int ss=0; ss<3; ss++){
        tg_all[ss]->SetLineWidth(2);
        tg_all[ss]->SetMarkerColor(tg_all[ss]->GetLineColor());
        tg_all[ss]->SetTitle("");
        tg_all[ss]->GetXaxis()->SetTitle("c#tau_{H} (#mum)");
        tg_all[ss]->GetYaxis()->SetTitle(Form("Events / bin %i", bin+1));

        tg_fixed_all[ss]->SetLineWidth(2);
        tg_fixed_all[ss]->SetMarkerColor(tg_fixed_all[ss]->GetLineColor());
        tg_fixed_all[ss]->SetTitle("");
        tg_fixed_all[ss]->GetXaxis()->SetTitle("c#tau_{H} (#mum)");
        tg_fixed_all[ss]->GetYaxis()->SetTitle(Form("Events / bin %i", bin+1));

        tg_all[ss]->GetXaxis()->SetNdivisions(505);
        tg_all[ss]->GetXaxis()->SetLabelFont(42);
        tg_all[ss]->GetXaxis()->SetLabelOffset(0.007);
        tg_all[ss]->GetXaxis()->SetLabelSize(0.04);
        tg_all[ss]->GetXaxis()->SetTitleSize(0.06);
        tg_all[ss]->GetXaxis()->SetTitleOffset(0.9);
        tg_all[ss]->GetXaxis()->SetTitleFont(42);
        tg_all[ss]->GetYaxis()->SetNdivisions(505);
        tg_all[ss]->GetYaxis()->SetLabelFont(42);
        tg_all[ss]->GetYaxis()->SetLabelOffset(0.007);
        tg_all[ss]->GetYaxis()->SetLabelSize(0.04);
        tg_all[ss]->GetYaxis()->SetTitleSize(0.06);
        tg_all[ss]->GetYaxis()->SetTitleOffset(1.1);
        tg_all[ss]->GetYaxis()->SetTitleFont(42);
        tg_all[ss]->GetYaxis()->SetRangeUser(0, max_plot);
        //      tg_all[ss]->GetXaxis()->SetRangeUser(0, 150);
      }

      l2D->AddEntry(tg, "Nominal", "l");
      l2D->AddEntry(tg_up, "Res. up", "l");
      l2D->AddEntry(tg_dn, "Res. down", "l");

      tg->Draw("al");
      tg_up->Draw("lsame");
      tg_dn->Draw("lsame");
      //    for (int ss=0; ss<3; ss++) tg_fixed_all[ss]->Draw("lsame");
      for (int ss=0; ss<1; ss++) tg_fixed_all[ss]->Draw("lsame");
      tg_up_regratio->Draw("lsame");
      tg_dn_regratio->Draw("lsame");

      canvasname_2D.Prepend(plotDir);
      TString canvasname_2D_pdf = canvasname_2D;
      TString canvasname_2D_eps = canvasname_2D;
      TString canvasname_2D_png = canvasname_2D;
      TString canvasname_2D_root = canvasname_2D;
      TString canvasname_2D_c = canvasname_2D;
      canvasname_2D_pdf.Append(".pdf");
      canvasname_2D_eps.Append(".eps");
      canvasname_2D_png.Append(".png");
      canvasname_2D_root.Append(".root");
      canvasname_2D_c.Append(".C");
      c2D->SaveAs(canvasname_2D_pdf);
      //    c2D->SaveAs(canvasname_2D_eps);
      c2D->SaveAs(canvasname_2D_png);
      c2D->SaveAs(canvasname_2D_root);
      //    c2D->SaveAs(canvasname_2D_c);

      delete l2D;
      c2D->Close();
      delete tg_up_regratio;
      delete tg_dn_regratio;
      for (int ss=0; ss<3; ss++) delete tg_fixed_all[ss];
      for (int ss=0; ss<3; ss++) delete tg_all[ss];
      gROOT->cd();
    }
  }

  for (int ct = 0; ct < nCTau; ct++){
    int genctau = ct*maxCTau / (nCTau - 1);
    xx[ct] = (double)genctau;

    TString coutput = cinput_Txymain;
    coutput.Append(Form("%i%s", genctau, "_Modified.root"));
    foutput = new TFile(coutput, "recreate");

    for (int iProd=0; iProd<nProd; iProd++){
      for (int ss = 0; ss < kNumTxySysts; ss++){
        for (int bin=1; bin<=nbins_KD; bin++){
          if (ss==0) h_KDarray[ct][kNumTxySysts*iProd + ss]->SetBinContent(bin, yy_all[iProd][bin-1][ct]);
          else if (ss==1) h_KDarray[ct][kNumTxySysts*iProd + ss]->SetBinContent(bin, yy_all_up[iProd][bin-1][ct]);
          else if (ss==2) h_KDarray[ct][kNumTxySysts*iProd + ss]->SetBinContent(bin, yy_all_dn[iProd][bin-1][ct]);
        }
        double templateyield = h_KDarray[ct][kNumTxySysts*iProd + ss]->Integral();
        double scale = totalyield[iProd] / templateyield;
        h_KDarray[ct][kNumTxySysts*iProd + ss]->Scale(scale);
      }
    }
    TString strTxySystAppend[kNumTxySysts]={ "TxyNominal", "TxyUp", "TxyDown" };
    double xsecScale = 1;
    for (int iProd=0; iProd<nProd; iProd++){
      for (int ss = 0; ss < kNumTxySysts; ss++){
        foutput->cd();

        TString hname = "T_2D_";
        hname.Append(strTxySystAppend[ss]);
        if (iProd>0){
          hname.Append(Form("_%s", cProduction[iProd].Data()));
          h_KDarray[ct][kNumTxySysts*iProd + ss]->SetName(hname);
          xsecScale = sum_XSEC / prodXSEC[iProd][EnergyIndex];
          h_KDarray[ct][kNumTxySysts*iProd + ss]->Scale(xsecScale);
          if (ss==0 && ct==0) cout << "Xsec scale: " << xsecScale << endl;
          if (ss==0 && ct==0) cout << cProduction[iProd] << " rate: " << h_KDarray[ct][kNumTxySysts*iProd + ss]->Integral() << endl;
          foutput->WriteTObject(h_KDarray[ct][kNumTxySysts*iProd + ss]);
          gROOT->cd();
          delete h_KDarray[ct][kNumTxySysts*iProd + ss];
        }
        else{
          h_KDarray[ct][kNumTxySysts*iProd + ss]->SetName(hname);

          hname.Append("_ggH");
          TH1F* hggH = (TH1F*)h_KDarray[ct][kNumTxySysts*iProd + ss]->Clone(hname);
          xsecScale = sum_XSEC / prodXSEC[iProd][EnergyIndex];
          hggH->Scale(xsecScale);
          if (ss==0 && ct==0) cout << "Xsec scale: " << xsecScale << endl;
          if (ss==0 && ct==0) cout << "ggH rate: " << hggH->Integral() << endl;
          foutput->WriteTObject(hggH);
          delete hggH;

          for (int jProd=1; jProd<nProd; jProd++) h_KDarray[ct][kNumTxySysts*iProd + ss]->Add(h_KDarray[ct][kNumTxySysts*jProd + ss], 1);
          foutput->WriteTObject(h_KDarray[ct][kNumTxySysts*iProd + ss]);
          if (ss==0 && ct==0) cout << "Total rate: " << h_KDarray[ct][kNumTxySysts*iProd + ss]->Integral() << endl;
          gROOT->cd();
          delete h_KDarray[ct][kNumTxySysts*iProd + ss];
        }
      }
    }
    foutput->cd();
    foutput->Close();
  }
}


void modifyBkgTemplates (TString dir, TString sqrts = "7TeV", TString channame = "2e2mu"){
	TString cinput_TxyNominal = dir + "/" + sqrts + "/" + channame + "_templates_Nominal_bkg.root";
	TString cinput_TxyUp = dir + "/" + sqrts + "/" + channame + "_templates_TxyUp_bkg.root";
	TString cinput_TxyDown = dir + "/" + sqrts + "/" + channame + "_templates_TxyDown_bkg.root";
//	TString cinput_TxyUp = dir + "/" + sqrts + "/" + channame + "_templates_Nominal_bkg.root";
//	TString cinput_TxyDown = dir + "/" + sqrts + "/" + channame + "_templates_Nominal_bkg.root";
	TString cinput_Dbkg = dir + "/" + sqrts + "/" + channame + "_templates_Nominal_bkg.root";
	TString coutput = dir + "/" + sqrts + "/" + channame + "_templates_Merged_bkg.root";
	
	const int TxySyst=3;
	TFile* finput[TxySyst+1] ={
		new TFile(cinput_Dbkg, "read"),
		new TFile(cinput_TxyNominal, "read"),
		new TFile(cinput_TxyUp, "read"),
		new TFile(cinput_TxyDown, "read")
	};
	if (finput[1]->IsZombie()){ finput[0]->Close(); return; }
	TFile* foutput = new TFile(coutput,"recreate");
	TH2F* ggzz_dbkg = (TH2F*) finput[0]->Get("template_ggZZ_Dbkg");
	TH2F* qqzz_dbkg = (TH2F*) finput[0]->Get("template_qqZZ_Dbkg");
	TH2F* zx_dbkg[3] = {
		(TH2F*)finput[0]->Get("template_ZX_Dbkg"),
		(TH2F*)qqzz_dbkg->Clone("ZX_Up"),
		(TH2F*)qqzz_dbkg->Clone("ZX_Down")
	};
	cout << ggzz_dbkg->GetName() << '\t';
	cout << ggzz_dbkg->Integral() << endl;
	ggzz_dbkg->Scale(1./ggzz_dbkg->Integral());
	cout << qqzz_dbkg->GetName() << '\t';
	cout << qqzz_dbkg->Integral() << endl;
	qqzz_dbkg->Scale(1./qqzz_dbkg->Integral());
	for (int ii = 0; ii < 3; ii++){
		cout << zx_dbkg[ii]->GetName() << '\t';
		cout << zx_dbkg[ii]->Integral() << endl;
		zx_dbkg[ii]->Scale(1. / zx_dbkg[ii]->Integral());
	}
	TH1F* ggzz_1d_dbkg = (TH1F*) ggzz_dbkg->ProjectionX();
	ggzz_1d_dbkg->SetName(ggzz_dbkg->GetName());
	TH1F* qqzz_1d_dbkg = (TH1F*) qqzz_dbkg->ProjectionX();
	qqzz_1d_dbkg->SetName(qqzz_dbkg->GetName());
	TH1F* zx_1d_dbkg[3];
	zx_1d_dbkg[0] = (TH1F*) zx_dbkg[0]->ProjectionX();
	zx_1d_dbkg[0]->SetName(Form("%s_Nominal",zx_dbkg[0]->GetName()));
	zx_1d_dbkg[1] = (TH1F*) zx_dbkg[1]->ProjectionX();
	zx_1d_dbkg[1]->SetName(Form("%s_ScaleResUp",zx_dbkg[0]->GetName()));
	zx_1d_dbkg[2] = (TH1F*) zx_dbkg[2]->ProjectionX();
	zx_1d_dbkg[2]->SetName(Form("%s_ScaleResDown",zx_dbkg[0]->GetName()));
	zx_1d_dbkg[2]->Scale(-1.);
	zx_1d_dbkg[2]->Add(zx_1d_dbkg[0], 2.);
	for (int binx=1; binx<=(zx_1d_dbkg[2]->GetNbinsX()); binx++){
		double bincontent = zx_1d_dbkg[2]->GetBinContent(binx);
		if (bincontent<=0) zx_1d_dbkg[2]->SetBinContent(binx,1.0e-20);
	}
	zx_1d_dbkg[2]->Scale(1. / zx_1d_dbkg[2]->Integral());


	foutput->WriteTObject(ggzz_1d_dbkg);
	foutput->WriteTObject(qqzz_1d_dbkg);
	for (int ii = 0; ii < 3; ii++) foutput->WriteTObject(zx_1d_dbkg[ii]);

	cout << "Finished bkg Dbkg templates"<<endl;
	foutput->cd();
	TH2F* ggzz_txy[TxySyst];
	TH2F* qqzz_txy[TxySyst];
	TH2F* zx_txy[TxySyst];
	TH1F* ggzz_1d_txy[TxySyst];
	TH1F* qqzz_1d_txy[TxySyst];
	TH1F* zx_1d_txy[TxySyst];
	for (int ii = 0; ii < TxySyst; ii++){
		ggzz_txy[ii] = (TH2F*)finput[ii + 1]->Get("template_ggZZ");
		qqzz_txy[ii] = (TH2F*)finput[ii + 1]->Get("template_qqZZ");
		zx_txy[ii] = (TH2F*)finput[ii + 1]->Get("template_ZX");

		TString templateName=ggzz_txy[ii]->GetName();
		cout << templateName << endl;
		if(ii==0) templateName.Append("_TxyNominal");
		if(ii==1) templateName.Append("_TxyUp");
		if(ii==2) templateName.Append("_TxyDown");
		ggzz_1d_txy[ii] = (TH1F*) ggzz_txy[ii]->ProjectionX();
		ggzz_1d_txy[ii]->SetName(templateName);

		templateName=qqzz_txy[ii]->GetName();
		cout << templateName << endl;
		if (ii==0) templateName.Append("_TxyNominal");
		if(ii==1) templateName.Append("_TxyUp");
		if(ii==2) templateName.Append("_TxyDown");
		qqzz_1d_txy[ii] = (TH1F*) qqzz_txy[ii]->ProjectionX();
		qqzz_1d_txy[ii]->SetName(templateName);

		templateName=zx_txy[ii]->GetName();
		cout << templateName << endl;
		if (ii==0) templateName.Append("_TxyNominal");
		if(ii==1) templateName.Append("_TxyUp");
		if(ii==2) templateName.Append("_TxyDown");
		zx_1d_txy[ii] = (TH1F*) zx_txy[ii]->ProjectionX();
		zx_1d_txy[ii]->SetName(templateName);
	}
	for (int ii = 1; ii < TxySyst; ii++){
		ggzz_1d_txy[ii]->Scale(ggzz_1d_txy[0]->Integral() / ggzz_1d_txy[ii]->Integral());
		qqzz_1d_txy[ii]->Scale(qqzz_1d_txy[0]->Integral() / qqzz_1d_txy[ii]->Integral());
		zx_1d_txy[ii]->Scale(zx_1d_txy[0]->Integral() / zx_1d_txy[ii]->Integral());
	}
  if (TxySyst>=2){
    zx_1d_txy[2]->Scale(2.);
    zx_1d_txy[2]->Add(zx_1d_txy[0],-1.);
    for (int binx=1; binx<=zx_1d_txy[2]->GetNbinsX(); binx++){
      double bincontent = zx_1d_txy[2]->GetBinContent(binx);
      if (bincontent<0){
        bincontent = 1e-15;
        zx_1d_txy[2]->SetBinContent(binx, bincontent);
      }
    }
  }
  for (int ii = 0; ii < TxySyst; ii++){
    double integral_zx = zx_1d_txy[0]->Integral();
    for (int binx=1; binx<=zx_1d_txy[ii]->GetNbinsX(); binx++){
      double bincontent = zx_1d_txy[ii]->GetBinContent(binx);
      if (bincontent<=0){
        bincontent = 1e-15;
        zx_1d_txy[ii]->SetBinContent(binx, bincontent);
      }
    }
    zx_1d_txy[ii]->Scale(integral_zx / zx_1d_txy[ii]->Integral());
  }
	for (int ii = 0; ii < TxySyst; ii++){
    cout << "Integral qq: " << qqzz_1d_txy[ii]->Integral() << endl;
    cout << "Integral gg: " << ggzz_1d_txy[ii]->Integral() << endl;
    cout << "Integral Z+X: " << zx_1d_txy[ii]->Integral() << endl;
    foutput->WriteTObject(ggzz_1d_txy[ii]);
		foutput->WriteTObject(qqzz_1d_txy[ii]);
		foutput->WriteTObject(zx_1d_txy[ii]);
	}

	foutput->Close();
	for (int ii=0; ii<TxySyst+1; ii++) finput[ii]->Close();
}

void modifyTemplates_LifetimeBuilder(TString dir, int processSig=1, int processBkg=1){
	TString channame[3] = { "4mu","4e","2e2mu" };
	TString sqrts[2] = { "7TeV","8TeV" };
//	for (int s = 1; s < 2; s++){
	for (int s = 0; s < 2; s++){
		for (int f = 0; f < 3; f++){
			if(processBkg==1) modifyBkgTemplates(dir, sqrts[s], channame[f]);
			if(processSig==1) modifySigTemplates(dir, sqrts[s], channame[f]);
		}
	}
}


void regularizeCTauSlice(TH1F* hSlice){
  int nbins_slice = hSlice->GetNbinsX();
  double integral_in = hSlice->Integral();
  bool isEven = (nbins_slice % 2 == 0);
  int bin_mid = 0, bin_midleft = 0, bin_midright = 0;
  if (isEven){
    bin_mid = nbins_slice / 2 + 1;
    bin_midright = bin_mid+1;
  }
  else{
    bin_mid = (nbins_slice+1) / 2;
    bin_midright = bin_mid+2;
  }
  bin_midleft = bin_mid-2;
  int bin_first = 2, bin_last = nbins_slice - 1;

  int nbins_second = nbins_slice - bin_mid + 1;
  double* xx_second = new double[nbins_second-1];
  double* yy_second = new double[nbins_second-1];

  double threshold = 0.05;
  const int nIter = 5;
  for (int it=0; it<nIter; it++){
    for (int binIt = bin_midright; binIt<bin_last; binIt++){
      int ctr = 0;
      for (int bin = bin_mid; bin<=nbins_slice; bin++){
        if (bin==binIt) continue;
        xx_second[ctr] = hSlice->GetXaxis()->GetBinCenter(bin);
        yy_second[ctr] = hSlice->GetBinContent(bin);
        ctr++;
      }

      TGraph* interpolator = new TGraph(nbins_second-1, xx_second, yy_second);
      double derivative_first = (yy_second[1]-yy_second[0])/(xx_second[1]-xx_second[0]);
      double derivative_last = (yy_second[nbins_second-2]-yy_second[nbins_second-3])/(xx_second[nbins_second-2]-xx_second[nbins_second-3]);
      TSpline3* spline = new TSpline3("spline", interpolator, "b1e1", derivative_first, derivative_last);

      double val = spline->Eval(hSlice->GetXaxis()->GetBinCenter(binIt));
      if (fabs(hSlice->GetBinContent(binIt)-val)>threshold*val) hSlice->SetBinContent(binIt, val);

      delete spline;
      delete interpolator;
    }
  }
  delete[] yy_second;
  delete[] xx_second;

  double integral_out = hSlice->Integral();
  double scale = integral_out / integral_in;
  hSlice->Scale(scale);
}

TGraph* regularizeSignalNominal(const int nbins, double xx[], double yy[], int nIter, int iProd){
  double minContent = 1e-6;
  if (iProd==1) minContent *= 5e-2;
  else if (iProd==2) minContent *= 2e-2;
  else if (iProd==3) minContent *= 1e-2;
  else if (iProd==4) minContent *= 2e-3;
  for (int binx=2; binx<nbins; binx++){
    if (yy[binx]!=yy[binx]) yy[binx] = yy[binx-1];
    if (yy[binx]<minContent){
      if (yy[binx-1]>yy[binx]){
        cout << binx << " is about 0. Setting " << yy[binx]  << " to " << yy[binx-1] << endl;
        yy[binx] = yy[binx-1];
      }
      else{
        cout << binx << " is about 0. Setting " << yy[binx]  << " to " << yy[binx-2] << endl;
        yy[binx] = yy[binx-2];
      }
    }
  }
  for (int it=0; it<1000; it++){
    double threshold = 0.005;
    for (int binx=1; binx<nbins-1; binx++){
      if (xx[binx]<90 || xx[binx]>=1000) continue;
      double average = (yy[binx-1] + yy[binx+1])*0.5;

      if (fabs(yy[binx]-average)>threshold*average){
        yy[binx] = average;
//        cout << "Correcting bin " << binx << endl;
      }
    }
  }

  double xx_new[nbins-1];
  double yy_new[nbins-1];

  for (int it=0; it<nIter; it++){
    double threshold = 0.07;
//    if (it==nIter-1) threshold=0;
    if (it>=nIter/2) threshold=0;
    for (int binx=1; binx<nbins-1; binx++){
      TGraph* interpolator;
      TSpline3* spline;
      for (int bb=0; bb<nbins; bb++){
        int newbin = bb;
        if (bb==binx) continue;
        else if (bb>binx) newbin--;
        xx_new[newbin] = xx[bb];
        yy_new[newbin] = yy[bb];
      }
      interpolator = new TGraph(nbins-1, xx_new, yy_new);
      double derivative_first = (yy[1]-yy[0])/(xx[1]-xx[0]);
      double derivative_last = (yy[nbins-1]-yy[nbins-2])/(xx[nbins-1]-xx[nbins-2]);
      spline = new TSpline3("spline", interpolator, "b1e1", derivative_first, derivative_last);

      double val = spline->Eval(xx[binx]);
      if (fabs(yy[binx]-val)>threshold*val) yy[binx] = val;
      delete spline;
      delete interpolator;
    }
  }
  TGraph* tg_fixed = new TGraph(nbins, xx, yy);
  return tg_fixed;
}

void regularizeSignalSystRatio(const int nbins, double xx[], double yy[], double yy_up[], double yy_dn[], int nIter){
  double ratio_up[nbins];
  double ratio_dn[nbins];

  for (int bin=0; bin<nbins; bin++){
    ratio_up[bin] = yy_up[bin] / yy[bin] - 1.;
    ratio_dn[bin] = yy_dn[bin] / yy[bin] - 1.;
  }
  int first_valid = 0;
  if (ratio_up[0]>0 && ratio_dn[0]<0) first_valid = 1;
  if (ratio_up[0]<0 && ratio_dn[0]>0) first_valid = -1;
  int last_valid = 0;
  if (ratio_up[nbins-1]>0 && ratio_dn[nbins-1]<0) last_valid = 1;
  if (ratio_up[nbins-1]<0 && ratio_dn[nbins-1]>0) last_valid = -1;
  if (first_valid==0){
    for (int bin=1; bin<nbins; bin++){
      if (ratio_up[bin]>0 && ratio_dn[bin]<0){
        first_valid = 1;
        ratio_up[0] = ratio_up[bin];
        ratio_dn[0] = ratio_dn[bin];
        break;
      }
      else if (ratio_up[bin]<0 && ratio_dn[bin]>0){
        first_valid = -1;
        ratio_up[0] = ratio_up[bin];
        ratio_dn[0] = ratio_dn[bin];
        break;
      }
      else continue;
    }
  }
  if (last_valid==0){
    for (int bin=0; bin<nbins-1; bin++){
      if (ratio_up[bin]>0 && ratio_dn[bin]<0){
        last_valid = 1;
        ratio_up[nbins-1] = ratio_up[bin];
        ratio_dn[nbins-1] = ratio_dn[bin];
        break;
      }
      else if (ratio_up[bin]<0 && ratio_dn[bin]>0){
        last_valid = -1;
        ratio_up[nbins-1] = ratio_up[bin];
        ratio_dn[nbins-1] = ratio_dn[bin];
        break;
      }
      else continue;
    }
  }
  for (int it=0; it<2; it++){
    double threshold = 0.2;
    for (int binx=1; binx<nbins-1; binx++){
      if ((ratio_up[binx]>0 && ratio_dn[binx]<0) || (ratio_up[binx]<0 && ratio_dn[binx]>0)){
        double average = (ratio_up[binx-1] + ratio_up[binx+1])*0.5;
        if (fabs(ratio_up[binx]-average)>threshold*average){
          ratio_up[binx] = average;
        }
        average = (ratio_dn[binx-1] + ratio_dn[binx+1])*0.5;
        if (fabs(ratio_dn[binx]-average)>threshold*average){
          ratio_dn[binx] = average;
        }
      }
      else{
        ratio_up[binx] = ratio_up[binx-1];
        ratio_dn[binx] = ratio_dn[binx-1];
      }
    }
  }
  
  double xx_new[nbins-1];
  double ratio_up_new[nbins-1];
  double ratio_dn_new[nbins-1];


  for (int it=0; it<nIter; it++){
    double threshold = 0.07;
    if (it>=nIter/2) threshold=0;
    for (int binx=1; binx<nbins-1; binx++){
      TGraph* interpolator;
      TSpline3* spline;
      for (int bb=0; bb<nbins; bb++){
        int newbin = bb;
        if (bb==binx) continue;
        else if (bb>binx) newbin--;
        xx_new[newbin] = xx[bb];
        ratio_up_new[newbin] = ratio_up[bb];
      }
      interpolator = new TGraph(nbins-1, xx_new, ratio_up_new);
      double derivative_first = (ratio_up[1]-ratio_up[0])/(xx[1]-xx[0]);
      double derivative_last = (ratio_up[nbins-1]-ratio_up[nbins-2])/(xx[nbins-1]-xx[nbins-2]);
      spline = new TSpline3("spline", interpolator, "b1e1", derivative_first, derivative_last);

      double val = spline->Eval(xx[binx]);
      if (fabs(ratio_up[binx]-val)>threshold*val){
        if ((val<0 && ratio_up[binx]>0) || (val<0 && ratio_up[binx]>0)) val = -val;
        if (val!=0) ratio_up[binx] = val;
      }
      delete spline;
      delete interpolator;
    }
  }

  for (int it=0; it<nIter; it++){
    double threshold = 0.07;
    if (it>=nIter/2) threshold=0;
    for (int binx=1; binx<nbins-1; binx++){
      TGraph* interpolator;
      TSpline3* spline;
      for (int bb=0; bb<nbins; bb++){
        int newbin = bb;
        if (bb==binx) continue;
        else if (bb>binx) newbin--;
        xx_new[newbin] = xx[bb];
        ratio_dn_new[newbin] = ratio_dn[bb];
      }
      interpolator = new TGraph(nbins-1, xx_new, ratio_dn_new);
      double derivative_first = (ratio_dn[1]-ratio_dn[0])/(xx[1]-xx[0]);
      double derivative_last = (ratio_dn[nbins-1]-ratio_dn[nbins-2])/(xx[nbins-1]-xx[nbins-2]);
      spline = new TSpline3("spline", interpolator, "b1e1", derivative_first, derivative_last);

      double val = spline->Eval(xx[binx]);
      if (fabs(ratio_dn[binx]-val)>threshold*val){
        if ((val<0 && ratio_dn[binx]>0) || (val<0 && ratio_dn[binx]>0)) val = -val;
        if(val!=0) ratio_dn[binx] = val;
      }
      delete spline;
      delete interpolator;
    }
  }

  for (int bin=0; bin<nbins; bin++){
    yy_up[bin] = (ratio_up[bin]+1.)*yy[bin];
    yy_dn[bin] = (ratio_dn[bin]+1.)*yy[bin];
  }
}

