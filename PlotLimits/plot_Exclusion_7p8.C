#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <vector>
#include "TCutG.h"
#include "TFile.h"
#include "TH2.h"
#include "TPad.h"

void getLimits(TFile *f, std::vector<double> &v_mh,std::vector<double> &v_mean,std::vector<double> &v_68l,std::vector<double> &v_68h,std::vector<double> &v_95l,std::vector<double> &v_95h,std::vector<double> &v_obs,std::vector<double> &v_obs95);

TGraphAsymmErrors *slidingWindowAverage(TGraphAsymmErrors *input, int slidingWindow);
TGraph *slidingWindowAverage2(TGraph *input, int slidingWindow);
TGraphAsymmErrors * removeGlitches(TGraphAsymmErrors *out);
TGraph * removeGlitches2(TGraph *out);


// --------- Inputs ------- //
TString inputFile =    "results/higgsCombineHZZ4L_ASCLS.E.root";//"results261012_1DoldZshape/higgsCombineHZZ4L_ASCLS.root";
TString inputFileObs = "results/higgsCombineHZZ4L_ASCLS.O.root";//"results261012_1DoldZshape/higgsCombineHZZ4L_ASCLS.root";
const bool addObsLimit = true;
const bool isXSxBR = false;
const bool _DEBUG_ = false;
string method = "FREQ";
Double_t xLow = 99.9;
Double_t xHigh = 1001.0;
Double_t yLow = 0.05;
Double_t yHigh = 20.0;
TString xTitle = "m_{H} (GeV)";
TString yTitle = "#sigma_{95%}/#sigma_{SM}";//"#sigma(H#rightarrow ZZ#rightarrow 4l)_{95% CL}/#sigma(H#rightarrow ZZ#rightarrow 4l)_{SM}";
const bool logy = true;
const bool logx = true;
const bool grid = true;
const bool gridOnTop = true;
const bool points = true;
const bool isTiny = false;
int canvasX = 550;
int canvasY = 526;
//double sqrts = 8.0;
//Double_t lumi = 1.616;
double sqrts = 7.0;
Double_t lumi = 5.1;
std::string plotDir = "plots";
std::string append = "3D_no2l2tau";
// ----------------------- //


using namespace std;

void plot_Exclusion_7p8()
{

  gROOT->ProcessLine(".L ./tdrstyle.C");
  setTDRStyle();
  //tdrStyle->SetPadLeftMargin(0.16);
  //tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetLegendFillColor(0);

  if(isXSxBR)
    {
      gSystem->Load("../CreateDatacards/include/HiggsCSandWidth_cc");
      gSystem->Load("../CreateDatacards/include/HiggsCSandWidthSM4_cc");
      HiggsCSandWidth *myCSW = new HiggsCSandWidth();
      HiggsCSandWidthSM4 *myCSWSM4 = new HiggsCSandWidthSM4();
    }

  TFile *inFile = new TFile(inputFile,"READ");
  TFile *inFileObs = new TFile(inputFileObs,"READ");
  
  // ------------------- Get Values -------------------- //

  vector<double> mH, Val_obs, Val_mean, Val_68h, Val_68l, Val_95h, Val_95l, Val_obs95;
  getLimits(inFile,mH,Val_mean,Val_68l,Val_68h,Val_95l,Val_95h,Val_obs,Val_obs95);
  vector<double> temp5,temp4,temp3,temp2,temp1,mH_Obs;
  getLimits(inFileObs,mH_Obs,temp2,temp3,temp4,temp5,temp1,Val_obs,Val_obs95);
  vector<double> v_masses, v_means, v_lo68, v_hi68, v_lo95, v_hi95, v_obs, v_massesObs;
  vector<double> expExclusion,obsExclusion,expExcl95,obsExcl95;

  for(unsigned int i = 0; i < mH.size(); i++)
    {
      v_masses.push_back( mH[i] );
      v_means.push_back( Val_mean[i] );
      v_lo68.push_back( min( Val_68l[i], Val_68h[i]) );
      v_hi68.push_back( max( Val_68h[i], Val_68l[i]) );
      v_lo95.push_back( min( Val_95l[i], Val_95h[i]) );
      v_hi95.push_back( max( Val_95h[i], Val_95l[i]) );
      if(Val_mean[i] < 1.0) expExclusion.push_back(mH[i]);
      if( max( Val_95h[i], Val_95l[i])< 1.0)expExcl95.push_back(mH[i]);
    }
  for(unsigned int i = 0; i < mH_Obs.size(); i++){
      v_massesObs.push_back( mH_Obs[i] );
      v_obs.push_back(Val_obs[i]);
      if(Val_obs[i] < 1.0 && addObsLimit) obsExclusion.push_back(mH_Obs[i]);
      if(Val_obs95[i]<1.0)obsExcl95.push_back(mH_Obs[i]);
  }

  // ------------------- For XSxBR --------------------- //

  double mHxs;
  double CStot,  CStot_p,  CStot_m;
  double BRHZZ4L, XSBRHZZ4L;
  const int sizeVXSBR = 490;
  
  double a_mhxs[sizeVXSBR], a_xsbr_err_h[sizeVXSBR], a_xsbr_err_l[sizeVXSBR];
  double a_xsbr[sizeVXSBR];
  
  if(isXSxBR)
    {
      int kk = 0;
      for( int k = 110; k < 600; k++)
	{ 
	  mHxs = k;
	  CStot = myCSW->HiggsCS(0,mHxs,sqrts,true);
	  CStot_p = myCSW->HiggsCSErrPlus(1,mHxs,sqrts);
	  CStot_m = myCSW->HiggsCSErrMinus(1,mHxs,sqrts);
	  
	  BRHZZ4L = myCSW->HiggsBR(14,mHxs,true);
	  XSBRHZZ4L = BRHZZ4L*CStot;
	  
	  a_mhxs[kk] = mHxs;
	  a_xsbr_err_h[kk] = XSBRHZZ4L*CStot_p*1000;
	  a_xsbr_err_l[kk] = -1*XSBRHZZ4L*CStot_m*1000;
	  a_xsbr[kk] = XSBRHZZ4L*1000;
	  kk++;
	}
      TGraph *gr_XSBR=new TGraphAsymmErrors(kk,a_mhxs,a_xsbr);
      TGraphAsymmErrors *gr_XSBR68=new TGraphAsymmErrors(kk,a_mhxs,a_xsbr,0,0,a_xsbr_err_l,a_xsbr_err_h);
    }

  // ------------------- Change Values to Arrays -------------------- //

  int nMassEffObs=0,nMassEff=0;
  int nExcluded=0;
  const int sizeV = v_masses.size();
  const int sizeobs = v_massesObs.size();
  double a_masses[sizeV], a_means[sizeV], a_lo68[sizeV], a_hi68[sizeV], a_lo95[sizeV], a_hi95[sizeV], a_obs[sizeobs],a_zero[sizeV], a_massesObs[sizeobs];

  for(unsigned int m = 0; m < v_massesObs.size(); m++){
      a_massesObs[nMassEffObs] = v_massesObs[m];
      a_obs[nMassEffObs] = v_obs[m];
      nMassEffObs++;
  }
  cout<<"Got obs arrays"<<endl;
  for(unsigned int m = 0; m < v_masses.size(); m++)
    {
      if(v_hi68.at(m)>=v_hi95.at(m) || v_lo68.at(m)<=v_lo95.at(m))
	{
	  cout << "Point at M = " << v_masses.at(m) << " excluded" << endl;
	  nExcluded++;
	  continue;
	}
      
      a_masses[nMassEff] = v_masses[m];
      a_means[nMassEff] = v_means[m];
      a_lo68[nMassEff] = v_lo68[m];
      a_hi68[nMassEff] = v_hi68[m];
      a_lo95[nMassEff] = v_lo95[m];
      a_hi95[nMassEff] = v_hi95[m];
      a_zero[nMassEff] = 0;
      if(isXSxBR)
	{
	  double xs = myCSW->HiggsCS(0,a_masses[nMassEff],sqrts,true);
	  double br = myCSW->HiggsBR(14,a_masses[nMassEff],true);
	  double XSxBR_factor = xs*br*1000;
	  cout << XSxBR_factor << endl;

	  a_obs[nMassEff]*=XSxBR_factor;
	  a_means[nMassEff]*=XSxBR_factor;
	  a_hi68[nMassEff]*=XSxBR_factor;
	  a_lo68[nMassEff]*=XSxBR_factor;
	  a_hi95[nMassEff]*=XSxBR_factor;
	  a_lo95[nMassEff]*=XSxBR_factor;
	}
      nMassEff++;
      
    }
  cout << "Excluded " << nExcluded << " sick mass points!" << endl;

  // --------------- Excluded Regions --------------- //
  cout<<"*******EXCLUSION*********"<<endl;
  for(int p = 0; p < expExclusion.size(); p++)
    {
      cout << "Expected Exclusion: " <<  expExclusion[p] << endl;
    }

  for(int s=0;s<expExcl95.size();s++)cout<<"Expected Excluxion 95% "<< expExcl95[s]<<endl;

  if(addObsLimit)
    {
      for(int q = 0; q < obsExclusion.size(); q++)
	{
	  cout << "Observed Exclusion: " <<  obsExclusion[q] << endl;
	}
      for(int t=0;t<obsExcl95.size();t++)cout<<"Observed Excluxion 95% "<< obsExcl95[t]<<endl;
    }

  cout<<endl;
  cout<<"*******LIMITS*********"<<endl;

  for(int r=0;r<nMassEff;r++){
    cout<<"Expected Limit "<<a_masses[r]<<": "<<a_means[r]<<endl;
  }

 if(addObsLimit)
    {
      for(int q2 = 0; q2 < Val_obs.size(); q2++)
	{
	  cout << "Observed limit " <<a_massesObs[q2]<<": "<< Val_obs[q2] << endl;
	}
      for(int t2=0;t2<Val_obs95.size();t2++)cout<<"Observed limit 95% " <<a_massesObs[t2]<<": "<< Val_obs95[t2] << endl;
    }
  
  // ------------------- Draw  -------------------- //



  TCanvas *c = new TCanvas("c","c",canvasX,canvasY);
  TGraphAsymmErrors* gr = new TGraphAsymmErrors(nMassEff, a_masses, a_means);
  TGraphAsymmErrors* grshade_68 = new TGraphAsymmErrors(nMassEff);
  TGraphAsymmErrors* grshade_95 = new TGraphAsymmErrors(nMassEff);
  gr->SetName("Expected");gr->SetTitle("Expected");
  grshade_68->SetName("68");grshade_68->SetTitle("68");
  grshade_95->SetName("95");grshade_95->SetTitle("95");
  TGraph *grObstemp = new TGraph(nMassEffObs, a_massesObs, a_obs);
  grObstemp->Sort();
  double *xtemp=grObstemp->GetX();
  double *ytemp=grObstemp->GetY();
  TGraph *grObs = new TGraph(0);
  for(int igo=0;igo<grObstemp->GetN();igo++){
    if(igo>0)
      if(xtemp[igo]-xtemp[igo-1]<0.5){
	float delta = TMath::Abs(xtemp[igo]-(int)xtemp[igo]);
	bool okd = false;
	if(TMath::Abs(delta-0.5)<0.01)okd=true;
	if(delta<0.01)okd=true;
	if(!okd)continue;
      }
    grObs->Set(grObs->GetN()+1);
    grObs->SetPoint(grObs->GetN()-1,xtemp[igo],ytemp[igo]);
  }
  grObs->SetLineWidth(3);
  grObs->SetLineColor(kBlack);
  grObs->SetMarkerStyle(20);
  grObs->SetMarkerSize(0.7);
 
  cout<<"nMasses: "<<nMassEff<<endl;

   for (int a = 0; a < nMassEff; a++)
     {
       grshade_68->SetPoint(a,a_masses[a],a_means[a]);
       grshade_68->SetPointError(a,0,0,a_means[a]-a_lo68[a],-a_means[a]+a_hi68[a]);
       grshade_95->SetPoint(a,a_masses[a],a_means[a]);
       grshade_95->SetPointError(a,0,0,a_means[a]-a_lo95[a],-a_means[a]+a_hi95[a]);
     }

  grshade_68->Sort();
  grshade_95->Sort();
  gr->SetLineStyle(2);
  gr->SetLineWidth(3);
  gr->SetLineColor(kBlue);
  grshade_68->SetFillColor(kGreen);
  grshade_95->SetFillColor(kYellow);		
  grshade_68->SetLineStyle(2);
  grshade_95->SetLineStyle(2);
  grshade_68->SetLineWidth(3);
  grshade_95->SetLineWidth(3);
  grshade_68->SetLineColor(kGreen);
  grshade_95->SetLineColor(kYellow);
	
	
  char outfileName[192];

  TF1* oneLine = new TF1("oneLine","1",xLow,xHigh);
  oneLine->SetLineColor(kRed);
  oneLine->SetLineWidth(3);

  // --------------- Low Mass Zoom -------------- //
	
  TLegend * box2 = new TLegend(0.5,0.65,0.83,0.9);
  box2->SetFillStyle(0);
  TLegendEntry *cht = box2->AddEntry((TObject*)0,"H#rightarrow ZZ* #rightarrow 4l","");
  cht->SetTextFont(42);
  cht->SetTextSize(0.03);
  if (addObsLimit){ box2->AddEntry(grObs,"Observed","l"); }
  box2->AddEntry(gr,"Expected","l");
  box2->AddEntry(grshade_68,"Expected #pm 1#sigma","lf");
  box2->AddEntry(grshade_95,"Expected #pm 2#sigma","lf");

  double ptLow= 0.28, ptHigh = 0.54;

  TPaveText *pt = new TPaveText(0.14,0.956,0.89,0.99,"brNDC");
  pt->SetTextAlign(12);
  pt->SetTextSize(0.03);
  pt->SetTextFont(42);
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  TText *text = pt->AddText(0.01,0.3,"CMS");
  text = pt->AddText(0.3,0.3,"#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}  #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}");
	
  if(grid) c->SetGrid();
   
  TH1F *hr = c->DrawFrame(105.0,yLow,180.0,yHigh);
	
  gr->Sort();
 
  grshade_95->Draw("e3");
  grshade_68->Draw("e3");
  gr->Draw("L");
  if(isXSxBR)
    {
      gr_XSBR68->Sort();
      gr_XSBR68->SetFillColor(kRed);
      gr_XSBR68->SetFillStyle(3013);
      gr_XSBR68->Draw("3same");
      gr_XSBR->SetLineColor(kRed);
      gr_XSBR->Draw("C");
    }
  if(addObsLimit)
    {
      grObs->Sort();
      if(points)grObs->Draw("CP");
      else grObs->Draw("C");
    }

  if(!isXSxBR)oneLine->Draw("LSAME");


  hr->GetXaxis()->SetTitle(xTitle);
  hr->GetYaxis()->SetTitle(yTitle);
  
  // hr->GetXaxis()->SetLabelFont(42);
  // hr->GetXaxis()->SetLabelOffset(0.007);
  // hr->GetXaxis()->SetLabelSize(0.05);
  // hr->GetXaxis()->SetTitleSize(0.06);
  // hr->GetXaxis()->SetTitleFont(42);
  // hr->GetYaxis()->SetLabelFont(42);
  // hr->GetYaxis()->SetLabelOffset(0.007);
  // hr->GetYaxis()->SetLabelSize(0.05);
  // hr->GetYaxis()->SetTitleSize(0.06);
  // hr->GetYaxis()->SetTitleOffset(1.1);
  // hr->GetYaxis()->SetTitleFont(42);
  // hr->GetZaxis()->SetLabelFont(42);
  // hr->GetZaxis()->SetLabelOffset(0.007);
  // hr->GetZaxis()->SetLabelSize(0.05);
  // hr->GetZaxis()->SetTitleSize(0.06);
  // hr->GetZaxis()->SetTitleFont(42);

  if(logy)gPad->SetLogy();
  
  
  c->Update();
  if(gridOnTop)gPad->RedrawAxis("g");
  pt->Draw("SAME");

  box2->Draw();
  sprintf( outfileName,"%s/UpperLimit_%s_lowMass_%s.eps",plotDir.c_str(),method.c_str(), append.c_str() );
  c->SaveAs(outfileName);
  sprintf( outfileName,"%s/UpperLimit_%s_lowMass_%s.png",plotDir.c_str(),method.c_str(), append.c_str() );	
  c->SaveAs(outfileName);
  sprintf( outfileName,"%s/UpperLimit_%s_lowMass_%s.root",plotDir.c_str(),method.c_str(), append.c_str() );
  c->SaveAs(outfileName);
  sprintf( outfileName,"%s/UpperLimit_%s_lowMass_%s.C",plotDir.c_str(),method.c_str(), append.c_str() );  c->SaveAs(outfileName);

  // --------------- Full Mass Range ---------------- //
	
  TLegend * box3 = new TLegend(0.5,0.7,0.79,0.9);
  box3->SetFillStyle(0);
  TLegendEntry *cht3 = box3->AddEntry((TObject*)0,"H#rightarrow ZZ* #rightarrow 4l","");
  cht3->SetTextFont(42);
  cht3->SetTextSize(0.03);
  if (addObsLimit){ box3->AddEntry(grObs,"Observed","l"); }
  box3->AddEntry(gr,"Expected","l");
  box3->AddEntry(grshade_68,"Expected #pm 1#sigma","lf");
  box3->AddEntry(grshade_95,"Expected #pm 2#sigma","lf");

  TCanvas *cl = new TCanvas("cl","cl",canvasX,canvasY);
  cl->cd();
  if(grid) cl->SetGrid();
  
  TH1F *hrl = cl->DrawFrame(xLow,yLow,xHigh,yHigh-5);
  grshade_95->Draw("e3");
  grshade_68->Draw("e3");
  gr->Draw("C");
  if(isXSxBR)
    {
      gr_XSBR68->SetFillColor(kRed);
      gr_XSBR68->SetFillStyle(3013);
      gr_XSBR68->Draw("3same");
      gr_XSBR->SetLineColor(kRed);
      gr_XSBR->Draw("C");
    }
  if(addObsLimit)
    {
      if(points)grObs->Draw("CP");
      else grObs->Draw("C");
    }

 
  if(!isXSxBR)oneLine->Draw("LSAME");
 

  hrl->GetXaxis()->SetTitle(xTitle);
  hrl->GetYaxis()->SetTitle(yTitle);
   // hrl->GetXaxis()->SetLabelFont(42);
   // hrl->GetXaxis()->SetLabelOffset(0.007);
   // hrl->GetXaxis()->SetLabelSize(0.05);
   // hrl->GetXaxis()->SetTitleSize(0.06);
   // hrl->GetXaxis()->SetTitleFont(42);
   // hrl->GetYaxis()->SetLabelFont(42);
   // hrl->GetYaxis()->SetLabelOffset(0.007);
   // hrl->GetYaxis()->SetLabelSize(0.05);
   // hrl->GetYaxis()->SetTitleSize(0.06);
   // hrl->GetYaxis()->SetTitleOffset(1.1);
   // hrl->GetYaxis()->SetTitleFont(42);
   // hrl->GetZaxis()->SetLabelFont(42);
   // hrl->GetZaxis()->SetLabelOffset(0.007);
   // hrl->GetZaxis()->SetLabelSize(0.05);
   // hrl->GetZaxis()->SetTitleSize(0.06);
   // hrl->GetZaxis()->SetTitleFont(42);

  if(logy)gPad->SetLogy();
  if(logx)
    {
      hrl->GetXaxis()->SetMoreLogLabels();
      hrl->GetXaxis()->SetNoExponent();
      gPad->SetLogx();
      TLine tick; tick.SetLineWidth(1); tick.SetLineColor(1);
      double dyh = (yHigh-5) * 0.08;
      double dyl = yLow * 0.08; //fabs(c1->PixeltoY(c1->VtoPixel(0.95)) - c1->PixeltoY(c1->VtoPixel(0.94)));
      if (gPad->GetLogy() && log((yHigh-5)/yLow) > log(1e6)) { dyh *= 2; dyl *= 2; }
      if (gPad->GetLogy() == 0) { dyh = dyl = 0.01*((yHigh-5)-yLow); }
      if (isTiny) { dyh *= 2; dyl *= 2; }
      for (int j = 100; j < xHigh; j += 10)  {
	  if (j > 400 && j % 20 == 10) continue;
	  tick.DrawLine(j, yLow, j, yLow+(j % 100 == 0 ? 2*dyl : dyl));
	  tick.DrawLine(j, (yHigh-5), j, yHigh-5-(j % 100 == 0 ? 2*dyh : dyh));
      }
    }
  cl->Update();
  hrl->GetXaxis()->SetRangeUser(xLow,xHigh);
  cl->Update();

  if(gridOnTop)gPad->RedrawAxis("g");

  pt->Draw("SAME");
  
  //box3->Draw();
  box2->Draw();
 
  sprintf( outfileName,"%s/UpperLimit_%s_wholeMass_%s.eps",plotDir.c_str(),method.c_str(), append.c_str() );
  cl->SaveAs(outfileName);
  sprintf( outfileName,"%s/UpperLimit_%s_wholeMass_%s.png",plotDir.c_str(),method.c_str(), append.c_str() );	
  cl->SaveAs(outfileName);
  sprintf( outfileName,"%s/UpperLimit_%s_wholeMass_%s.root",plotDir.c_str(),method.c_str(), append.c_str() );
  cl->SaveAs(outfileName);	
  sprintf( outfileName,"%s/UpperLimit_%s_wholeMass_%s.C",plotDir.c_str(),method.c_str(), append.c_str() );  cl->SaveAs(outfileName);

  // ---------------- Root File --------------- //

  TFile outf("testUL_single.root","RECREATE");
  outf.cd();
  gr->Write();
  outf.Write();

	
}



void getLimits(TFile *f, std::vector<double> &v_mh,std::vector<double> &v_mean,std::vector<double> &v_68l,std::vector<double> &v_68h,std::vector<double> &v_95l,std::vector<double> &v_95h,std::vector<double> &v_obs,std::vector<double> &v_obs95)
{

  TTree *tree =(TTree*)f->Get("limit");
  
  double mh,limit,errObs;
  float quant;
  tree->SetBranchAddress("mh",&mh);
  tree->SetBranchAddress("limit",&limit);
  tree->SetBranchAddress("quantileExpected",&quant);
  tree->SetBranchAddress("limitErr",&errObs);

  int nentry=tree->GetEntries();
  cout<<"entries"<<nentry<<endl;
  for(int i=0;i<nentry;i++)
    {
      tree->GetEntry(i);
      if(_DEBUG_)cout << "mH: " << mh << " limit: " << limit << " quantileExpected: " << quant << endl;  
      if(quant>-1.01&&quant<-0.99)
	{
	  v_obs.push_back(limit);
	  v_mh.push_back(mh);
	  v_obs95.push_back(limit+2*errObs/limit);
	}
      else if(quant>0.024 && quant<0.026) {v_95l.push_back(limit); v_mh.push_back(mh);}
      else if(quant>0.15  && quant<0.17 ) v_68l.push_back(limit);
      else if(quant>0.49  && quant<0.51 ) v_mean.push_back(limit);
      else if(quant>0.83  && quant<0.85 ) v_68h.push_back(limit);
      else if(quant>0.974 && quant<0.976) v_95h.push_back(limit);
      else {cout<<"Error! Unknown Quantile =  " << quant << endl;}
      
    }
  
}






TGraphAsymmErrors *slidingWindowAverage(TGraphAsymmErrors *input, int slidingWindow) {

  //cout << 11 << endl;

    TGraphAsymmErrors *out  = (TGraphAsymmErrors *) input->Clone();
    //bool isLogPlot = TString(input->GetName()).Contains("smcls");

    //cout << 22 << endl;

    bool isLogPlot = true;
    for (int i = 0, n = input->GetN(); i < n; ++i) {
        /*if (i >= 1 && i < n-1 && (input->GetX()[i] < input->GetX()[i-1] || input->GetX()[i] > input->GetX()[i+1])) {
            out->GetX()[i] = 0.5*(input->GetX()[i-1] + input->GetX()[i+1]);
            out->GetY()[i] = 0.5*(input->GetY()[i-1] + input->GetY()[i+1]);
            out->GetEYlow()[i] = 0.5*(input->GetEYlow()[i-1] + input->GetEYlow()[i+1]);
            out->GetEYhigh()[i] = 0.5*(input->GetEYhigh()[i-1] + input->GetEYhigh()[i+1]);
            continue;
        }*/

      //cout << 33 << endl;

        double y0 = input->GetY()[i];
        double sum = 0, sumhi = 0, sumlo = 0, sumw = 0;
        for (int j = i-slidingWindow; j <= i+slidingWindow; ++j) 
	  {

            if (j < 0 || j >= n) continue;
            double y = input->GetY()[j], w = 1.0; // /(1.0 + abs(i-j));
	    //cout << j << "  " << y0 << "  " << y << "  " << slidingWindow << endl;
            if (isLogPlot) 
	      {
		if (fabs(log(y0/y)) > log(2)) continue;    
	      }
	    else {
	      if (fabs(y0-y) > 0.1*y0) continue;    
            }
            if (isLogPlot) 
	      {
                sum   += w*log(y);
                sumlo += w*log(input->GetEYlow()[j]);
                sumhi += w*log(input->GetEYhigh()[j]);
	      } 
	    else {
	      sum   += w*y;
	      sumlo += w*input->GetEYlow()[j];
	      sumhi += w*input->GetEYhigh()[j];
            }
            sumw  += w;
	  }
        if (sumw == 0) continue;
        if (isLogPlot) 
	  {
	    out->GetY()[i] = exp(sum/sumw);
	    out->GetEYlow()[i] = exp(sumlo/sumw);
	    out->GetEYhigh()[i] = exp(sumhi/sumw);
	  } 
	else {
	  out->GetY()[i] = sum/sumw;
	  out->GetEYlow()[i] = sumlo/sumw;
	  out->GetEYhigh()[i] = sumhi/sumw;
        }
    }
    return out;
}

TGraph *slidingWindowAverage2(TGraph *input, int slidingWindow) {

  cout << 11 << endl;

    TGraph *out  = (TGraph*) input->Clone();
    //bool isLogPlot = TString(input->GetName()).Contains("smcls");

    cout << 22 << endl;

    bool isLogPlot = true;
    for (int i = 0, n = input->GetN(); i < n; ++i) {
        /*if (i >= 1 && i < n-1 && (input->GetX()[i] < input->GetX()[i-1] || input->GetX()[i] > input->GetX()[i+1])) {
            out->GetX()[i] = 0.5*(input->GetX()[i-1] + input->GetX()[i+1]);
            out->GetY()[i] = 0.5*(input->GetY()[i-1] + input->GetY()[i+1]);
            out->GetEYlow()[i] = 0.5*(input->GetEYlow()[i-1] + input->GetEYlow()[i+1]);
            out->GetEYhigh()[i] = 0.5*(input->GetEYhigh()[i-1] + input->GetEYhigh()[i+1]);
            continue;
        }*/

      cout << 33 << endl;

        double y0 = input->GetY()[i];
        double sum = 0, sumhi = 0, sumlo = 0, sumw = 0;
        for (int j = i-slidingWindow; j <= i+slidingWindow; ++j) 
	  {

            if (j < 0 || j >= n) continue;
            double y = input->GetY()[j], w = 1.0; // /(1.0 + abs(i-j));
	    cout << j << "  " << y0 << "  " << y << "  " << slidingWindow << endl;
            if (isLogPlot) 
	      {
		if (fabs(log(y0/y)) > log(2)) continue;    
	      }
	    else {
	      if (fabs(y0-y) > 0.1*y0) continue;    
            }
            if (isLogPlot) 
	      {
                sum   += w*log(y);
                sumlo += w*log(input->GetEYlow()[j]);
                sumhi += w*log(input->GetEYhigh()[j]);
	      } 
	    else {
	      sum   += w*y;
	      sumlo += w*input->GetEYlow()[j];
	      sumhi += w*input->GetEYhigh()[j];
            }
            sumw  += w;
	  }
        if (sumw == 0) continue;
        if (isLogPlot) 
	  {
	    out->GetY()[i] = exp(sum/sumw);
	    out->GetEYlow()[i] = exp(sumlo/sumw);
	    out->GetEYhigh()[i] = exp(sumhi/sumw);
	  } 
	else {
	  out->GetY()[i] = sum/sumw;
	  out->GetEYlow()[i] = sumlo/sumw;
	  out->GetEYhigh()[i] = sumhi/sumw;
        }
    }
    return out;
}





TGraphAsymmErrors * removeGlitches(TGraphAsymmErrors *out) {
  do {
    int bad = -1; bool hasgood = false;
    for (int i = 0; i < out->GetN(); ++i) {
      if (out->GetEYlow()[i] == 0 && out->GetEYhigh()[i] == 0) {
	bad = i;
      } else {
	hasgood = true;
      }
    }
    if (!hasgood) return out;
    if (bad == -1) return out;
    out->RemovePoint(bad);
  } while (1);
}

TGraph * removeGlitches2(TGraph *out) {
  do {
    int bad = -1; bool hasgood = false;
    for (int i = 0; i < out->GetN(); ++i) {
      if (out->GetEYlow()[i] == 0 && out->GetEYhigh()[i] == 0) {
	bad = i;
      } else {
	hasgood = true;
      }
    }
    if (!hasgood) return out;
    if (bad == -1) return out;
    out->RemovePoint(bad);
  } while (1);
}
