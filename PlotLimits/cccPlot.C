#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH2D.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TPaveText.h"
//#include "RooHZZStyle.C"
//#include "contours.cxx"

#include <sstream>
#include <iostream>

using namespace std;

bool plotInclusive=false;

TLatex *CMSPreliminary(float lumi7TeV=5.1, float lumi8TeV=19.7) {

  stringstream line;
  line << "CMS              #sqrt{s}=7 TeV, L=" << lumi7TeV << " fb^{-1}  #sqrt{s}=8 TeV, L=" << lumi8TeV << " fb^{-1}";
  TLatex* CP = new TLatex(0.15,0.96, line.str().c_str());
  CP->SetNDC(kTRUE);
  CP->SetTextSize(0.032);

  return CP;

}

TPaveText *text(const char *txt, float x1, float y1, float x2, float y2) {
  TPaveText *text = new TPaveText(x1,y1,x2,y2,"brNDC");
  text->AddText(txt);
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.05);
  return text;
}

    
void cccPlot() {

  float fitval[3], fiterrl[3], fiterrh[3];
  for(int cha=0; cha<3; ++cha) {
    fiterrl[cha]=0;
    fiterrh[cha]=0;
  }

  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetLegendFillColor(0);
  
  //TStyle *mystyle = RooHZZStyle("ZZ");
  //mystyle->cd();

  TFile *fituntagged = TFile::Open("125.7_Untagged/higgsCombineTest.MultiDimFit.mH125.7.root");
  TTree *treeuntagged = (TTree*)fituntagged->Get("limit");

  TFile *fitdijet = TFile::Open("125.7_Dijet/higgsCombineTest.MultiDimFit.mH125.7.root");
  TTree *treedijet = (TTree*)fitdijet->Get("limit");

  TFile *fitcomb = TFile::Open("125.7/higgsCombineTest.MultiDimFit.mH125.7.root");
  TTree *treecomb = (TTree*)fitcomb->Get("limit");

  TTree *trees[3]={treeuntagged,treedijet,treecomb};

  for(int cha=0; cha<3; ++cha) {
    
    cout << "Analyzing scan for channel = " << cha << endl;

    float MH;
    float deltaNLL;
    trees[cha]->SetBranchAddress("r", &MH);
    trees[cha]->SetBranchAddress("deltaNLL", &deltaNLL);
    
    float dNLLMinus=1000; float dNLLPlus=1000;
    for(int i=0; i<(int)trees[cha]->GetEntries();++i) {
      trees[cha]->GetEntry(i);
      if(i==0) {fitval[cha]=MH;printf("r0=%f\n",MH);}
      else {
	 if(fabs(2*deltaNLL-1)<dNLLMinus && MH<fitval[cha]) {
	   fiterrl[cha]=MH;
	   dNLLMinus=fabs(2*deltaNLL-1);
	 }
	 if(fabs(2*deltaNLL-1)<dNLLPlus && MH>fitval[cha]) {
	   fiterrh[cha]=MH;
	   dNLLPlus=fabs(2*deltaNLL-1);
	 }
       }
    }
  }

  float min=1000000,max=-1;
  for(int cha=0; cha<3; ++cha) {
    fiterrl[cha]=fitval[cha]-fiterrl[cha];
    fiterrh[cha]=fiterrh[cha]-fitval[cha];
    // patch if the scan arrested too early
    if(fiterrh[cha]==(-fitval[cha])) fiterrh[cha]=fiterrl[cha];
      cout << "Mass for channel " << cha << " = " << fitval[cha] 
	   << " -" << fiterrl[cha] << " +" << fiterrh[cha] << " GeV" << endl;
      if(fitval[cha]-fiterrl[cha]<min)min=fitval[cha]-fiterrl[cha];
      if(fitval[cha]+fiterrh[cha]>max)max=fitval[cha]+fiterrh[cha];
  }
  printf("min %.3f max %.3f\n",min,max);

  //TLatex l; l.SetTextFont(43); l.SetNDC(); l.SetTextSize(25);

    TCanvas *c1 = new TCanvas("c1","c1",550,526);
    //    c1->SetLeftMargin(0.1);
    //    c1->SetBottomMargin(0.15);
    c1->SetGridx(1);

    int nChann = 2;
    if(plotInclusive)nChann++;

    TH2F *frame=new TH2F("frame",";best fit #mu_{S};",1,0.5,2.5,nChann,0,nChann);
    //TH2F *frame=new TH2F("frame",";best fit #mu_{S};",1,0.,2.5,nChann,0,nChann);

    TGraphAsymmErrors points(nChann);
    for (int cha=0; cha<nChann; ++cha) {
      TString channame("");
      if (cha==0) channame+="0/1 jet";
      if (cha==1) channame+="dijet";
      if (cha==2) channame+="Inclusive";
      points.SetPoint(cha,       fitval[cha],  cha+0.5);
      points.SetPointError(cha,  fiterrl[cha], fiterrh[cha], 0, 0);
      frame->GetYaxis()->SetBinLabel(cha+1, channame);
    }
    points.SetLineColor(kRed);
    points.SetLineWidth(3);
    points.SetMarkerStyle(21);
    frame->GetXaxis()->SetNdivisions(8,kFALSE);
    //frame->GetXaxis()->SetNdivisions(5,kFALSE);
    frame->Draw(); //gStyle->SetOptStat(0);
    TBox globalFitBand(fitval[2]-fiterrl[2], 0.0001, fitval[2]+fiterrh[2], (float)nChann-0.01);
    globalFitBand.SetFillStyle(1001);
    globalFitBand.SetFillColor(kGreen);
    globalFitBand.SetLineStyle(0);
    globalFitBand.DrawClone();
    TLine globalFitLine(fitval[2], 0, fitval[2], nChann);
    globalFitLine.SetLineWidth(4);
    globalFitLine.SetLineColor(214);
    globalFitLine.DrawClone();
    points.DrawClone("P SAME");
    
    TPaveText *pt = new TPaveText(0.14,0.956,0.89,0.99,"brNDC");
    pt->SetTextAlign(12);
    pt->SetTextSize(0.03);
    pt->SetTextFont(42);
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    TText *text = pt->AddText(0.01,0.3,"CMS");
    text = pt->AddText(0.3,0.3,"#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}  #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}");
  pt->Draw();


    TLine one(1,0,1,nChann);
    one.Draw("SAME");
    TString outname="plots/bestfit_bycategory";
    if(plotInclusive)outname+="_withIncl";
    TString outnamepdf = outname+".pdf";
    c1->SaveAs(outnamepdf.Data());
    TString  outnamepng=outname+".png";
    c1->SaveAs(outnamepng.Data());
    TString  outnameeps=outname+".eps";
    c1->SaveAs(outnameeps.Data());
    TString  outnameroot=outname+".root";
    c1->SaveAs(outnameroot.Data());
    TString  outnameC=outname+".C";
    c1->SaveAs(outnameC.Data());
}

