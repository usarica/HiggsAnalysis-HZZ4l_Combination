#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
from ROOT import TH3F
import ROOT
from array import array
from systematicsClass import *
from inputReader import *

## ------------------------------------
##  card and workspace class
## ------------------------------------

class datacardClass_HCTau_1D:

    def __init__(self):
    
        self.ID_4mu = 1
        self.ID_4e  = 2
        self.ID_2e2mu = 3    
        self.isFSR = True

    def loadIncludes(self):
        
        ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
        ROOT.gSystem.AddIncludePath("-Iinclude/")
        ROOT.gROOT.ProcessLine(".L include/tdrstyle.cc")
        ROOT.gSystem.Load("libRooFit")
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
        ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
        ROOT.gSystem.Load("include/HiggsCSandWidthSM4_cc.so")


    # cs x br function 
    def makeXsBrFunction(self,signalProc,rrvMH):
            
        procName = "ggH"
        if(signalProc == 0): procName = "ggH" #dummy, when you sum up all the 5 chans
        if(signalProc == 1): procName = "ggH"
        if(signalProc == 2): procName = "qqH"
        if(signalProc == 3): procName = "WH"
        if(signalProc == 4): procName = "ZH"
        if(signalProc == 5): procName = "ttH"

        
        
        
        channelName = ""
        if (self.channel == self.ID_4mu): channelName = "4mu"
        elif (self.channel == self.ID_4e): channelName = "4e"
        elif (self.channel == self.ID_2e2mu): channelName = "2e2mu"
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)" 

     
        
        myCSWrhf = HiggsCSandWidth()
        
        histXsBr = ROOT.TH1F("hsmxsbr_{0}_{1}".format(procName,channelName),"", 8905, 109.55, 1000.05)
                
        for i in range(1,8906):
            
            mHVal = histXsBr.GetBinCenter(i)
            BR = 0.0 
            if (self.channel == self.ID_2e2mu):
                BR = myCSWrhf.HiggsBR(13,mHVal)
            else:
                BR = myCSWrhf.HiggsBR(12,mHVal)

            if (signalProc == 3 or signalProc == 4 or signalProc == 5):
                #overwrite BR if VH,ttH sample
                #these samples have inclusive Z decay
                BR = myCSWrhf.HiggsBR(11,mHVal)

            if (signalProc==0):
                totXs=0
                for ch in range(1,6):
                    totXs+=myCSWrhf.HiggsCS(ch, mHVal, self.sqrts)
                histXsBr.SetBinContent(i, totXs * BR)
            else:
                histXsBr.SetBinContent(i, myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts) * BR)

            #print '\nmakeXsBrFunction : procName=',procName,'   signalProc=',signalProc,'  mH (input)=',rrvMH.getVal(),
            #print '   CS=',myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts),'   BR=',BR
            
        rdhname = "rdhXsBr_{0}_{1}_{2}".format(procName,self.channel,self.sqrts)
        rdhXsBr = RooDataHist(rdhname,rdhname, ROOT.RooArgList(rrvMH), histXsBr)  
        
        return rdhXsBr
    
    # return trueVar if testStatement else return falseVar
    def getVariable(self,trueVar,falseVar,testStatement):

        if (testStatement): 
            return trueVar
        else:
            return falseVar
    
    # main datacard and workspace function
    def makeCardsWorkspaces(self, theMH, theOutputDir, theInputs,theTemplateDir="templates2D",theDataAppendDir=""):

        ## --------------- SETTINGS AND DECLARATIONS --------------- ##
        DEBUG = False
        self.mH = theMH
        self.SMDsigCut = 1. 
        self.SMDbkgCut = 1. 
        self.lumi = theInputs['lumi']
        self.sqrts = theInputs['sqrts']
        self.channel = theInputs['decayChannel']
        self.bkgMorph = theInputs['useCMS_zz4l_zjet']
        self.outputDir = theOutputDir
        self.templateDir = theTemplateDir
        self.dataAppendDir = theDataAppendDir

        self.ggH_chan = theInputs['ggH']
        self.qqH_chan = theInputs['qqH']
        self.WH_chan = theInputs['WH']
        self.ZH_chan = theInputs['ZH']
        self.ttH_chan = theInputs['ttH']
        self.qqZZ_chan = theInputs['qqZZ']
        self.ggZZ_chan = theInputs['ggZZ']
        self.zjets_chan = theInputs['zjets']
        
        ## ---------------- SET PLOTTING STYLE ---------------- ## 
        ROOT.setTDRStyle(True)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPadLeftMargin(0.16)        

        ## ---------------- VARIABLES FOR LATER --------------- ##
        self.bUseCBnoConvolution = False
        ForXSxBR = False

        myCSW = HiggsCSandWidth()
                
        ## ----------------- WIDTH AND RANGES ----------------- ##
        self.widthHVal =  myCSW.HiggsWidth(0,self.mH)
        if(self.widthHVal < 0.12):
            self.bUseCBnoConvolution = True
        self.isHighMass = False
        if self.mH >= 390:
            if theInputs['useHighMassReweightedShapes']:
                self.isHighMass = True
            else: print "useHighMassReweightedShapes set to FALSE, using non-reweighted shapes!"

            
        print "width: ",self.widthHVal
        
        self.windowVal = max( self.widthHVal, 1.0)
        lowside = 100.0
        highside = 1000.0
        
        if (self.mH >= 275):
            lowside = 180.0
            highside = 650.0
        if (self.mH >= 350):
            lowside = 200.0
            highside = 900.0
        if (self.mH >= 500):
            lowside = 250.0
            highside = 1000.0
        if (self.mH >= 700):
            lowside = 350.0
            highside = 1400.0
        
        self.low_M = max( (self.mH - 20.*self.windowVal), lowside)
        self.high_M = min( (self.mH + 15.*self.windowVal), highside)
       
        if (self.channel == self.ID_4mu): self.appendName = '4mu'
        elif (self.channel == self.ID_4e): self.appendName = '4e'
        elif (self.channel == self.ID_2e2mu): self.appendName = '2e2mu'
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"
            
            
        
        ## ------------------------- SYSTEMATICS CLASSES ----------------------------- ##
    
        systematics = systematicsClass( self.mH, False, self.isFSR, theInputs)

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##
    
        nctau = 41
        bins = 1000
        if(self.bUseCBnoConvolution): bins = 200

        CMS_zz4l_mass_name = "CMS_zz4l_mass"
            
        CMS_zz4l_mass = ROOT.RooRealVar(CMS_zz4l_mass_name,CMS_zz4l_mass_name,self.low_M,self.high_M)    
        CMS_zz4l_mass.setBins(bins)

        x_name = "CMS_zz4l_ctau"
        x_min = 0.
        x_max = 1000.
        x = ROOT.RooRealVar(x_name,x_name,x_min,x_min,x_max)
        x.setBins(nctau-1)

# Variable to morph Txy nominal, up and down systematics
        alphaMorph_Txy_Name =  "CMS_zz4l_TxyMorph_{0:.0f}".format(self.channel)
        alphaMorph_Txy = ROOT.RooRealVar(alphaMorph_Txy_Name,alphaMorph_Txy_Name,0,-10,10)
        alphaMorph_Txy.setConstant(False)
        morphVarList_Txy = ROOT.RooArgList()
        morphVarList_Txy.add(alphaMorph_Txy)


        D1Name = "CMS_zz4l_KD"
        D2Name = "CMS_zz4l_smd"
        D1 = ROOT.RooRealVar(D1Name,D1Name,-1,1) # KD
        D1.setBins(80)
        D2 = ROOT.RooRealVar(D2Name,D2Name,0,1) #SuperMELA
        D2.setBins(40)


        self.LUMI = ROOT.RooRealVar("LUMI_{0:.0f}".format(self.sqrts),"LUMI_{0:.0f}".format(self.sqrts),self.lumi)
        self.LUMI.setConstant(True)
    
        self.MH = ROOT.RooRealVar("MH","MH",self.mH)
        self.MH.setConstant(True)


	# n2, alpha2 are right side parameters of DoubleCB
	# n, alpha are left side parameters of DoubleCB

        n_CB_d = 0.0
        alpha_CB_d = 0.0
        n2_CB_d = 0.0
        alpha2_CB_d = 0.0
        mean_CB_d = 0.0
        sigma_CB_d = 0.0
        mean_BW_d = self.mH
        gamma_BW_d = 0.0
        
        rdhXsBrFuncV_1 = self.makeXsBrFunction(1,self.MH)
        
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ggH",self.channel,self.sqrts)
        rhfXsBrFuncV_1 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_1, 1)
        
        rdhXsBrFuncV_2 = self.makeXsBrFunction(2,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("VBF",self.channel,self.sqrts)
        rhfXsBrFuncV_2 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_2, 1)
        
        rdhXsBrFuncV_3 = self.makeXsBrFunction(3,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("WH",self.channel,self.sqrts)
        rhfXsBrFuncV_3 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_3, 1)
        
        rdhXsBrFuncV_4 = self.makeXsBrFunction(4,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ZH",self.channel,self.sqrts)
        rhfXsBrFuncV_4 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_4, 1)
        
        rdhXsBrFuncV_5 = self.makeXsBrFunction(5,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ttH",self.channel,self.sqrts)
        rhfXsBrFuncV_5 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_5, 1)
        
    
        ## -------- Variable Definitions -------- ##
        name = "CMS_zz4l_mean_e_sig"
        CMS_zz4l_mean_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_mean_e_err_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_mean_e_err = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_err",float(theInputs['CMS_zz4l_mean_e_sig']),-0.99,0.99)
        name = "CMS_zz4l_sigma_e_sig"
        CMS_zz4l_sigma_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_e_sig",3.0,0.0,30.0)
        name = "CMS_zz4l_mean_m_sig"
        CMS_zz4l_mean_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_mean_m_err_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_mean_m_err = ROOT.RooRealVar(name,"CMS_zz4l_mean_m_err",float(theInputs['CMS_zz4l_mean_m_sig']),-0.99,0.99)
        name = "CMS_zz4l_sigma_m_sig"
        CMS_zz4l_sigma_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",3.0,0.0,30.0)
            
        
        name = "CMS_zz4l_alpha2_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_alpha2 = ROOT.RooRealVar(name,"CMS_zz4l_alpha2",1.,-10.,10.)
        name = "CMS_zz4l_n2_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_n2 = ROOT.RooRealVar(name,"CMS_zz4l_n2",2.,-10.,10.)
        name = "CMS_zz4l_alpha_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_alpha = ROOT.RooRealVar(name,"CMS_zz4l_alpha",1.,-10.,10.)
        name = "CMS_zz4l_n_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_n = ROOT.RooRealVar(name,"CMS_zz4l_n",2.,-10.,10.)
        name = "CMS_zz4l_mean_BW_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_mean_BW = ROOT.RooRealVar(name,"CMS_zz4l_mean_BW",self.mH,self.low_M,self.high_M)
        name = "interf_ggH"
        #name = "CMS_zz4l_gamma_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_gamma = ROOT.RooRealVar(name,"CMS_zz4l_gamma",10.,0.001,1000.)
        name = "CMS_zz4l_widthScale_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_widthScale = ROOT.RooRealVar(name,"CMS_zz4l_widthScale",1.0)
            
        one = ROOT.RooRealVar("one","one",1.)
        one.setConstant(True)
    
        CMS_zz4l_mean_BW.setVal( mean_BW_d )
        CMS_zz4l_gamma.setVal(0)
        CMS_zz4l_mean_e_sig.setVal(0)
        CMS_zz4l_mean_e_err.setConstant(kTRUE)
        CMS_zz4l_sigma_e_sig.setVal(0)
        CMS_zz4l_mean_m_sig.setVal(0)
        CMS_zz4l_mean_m_err.setConstant(kTRUE)
        CMS_zz4l_sigma_m_sig.setVal(0)
        CMS_zz4l_alpha.setVal(0)
        CMS_zz4l_n.setVal(0)
        CMS_zz4l_alpha2.setVal(0)
        CMS_zz4l_n2.setVal(0)
    
        CMS_zz4l_widthScale.setConstant(True)
        CMS_zz4l_mean_BW.setConstant(True)

        print "mean_BW ", CMS_zz4l_mean_BW.getVal()
        print "gamma_BW ", CMS_zz4l_gamma.getVal()
        print "mean_e_sig ", CMS_zz4l_mean_e_sig.getVal()
        print "mean_e_err ", CMS_zz4l_mean_e_err.getVal()
        print "sigma_e ", CMS_zz4l_sigma_e_sig.getVal()
        print "mean_m_sig ",CMS_zz4l_mean_m_sig.getVal()
        print "mean_m_err ", CMS_zz4l_mean_m_err.getVal()
        print "sigma_m ", CMS_zz4l_sigma_m_sig.getVal()
        print "alpha ", CMS_zz4l_alpha.getVal()
        print "n ", CMS_zz4l_n.getVal()
        print "alpha2 ", CMS_zz4l_alpha2.getVal()
        print "n2 ", CMS_zz4l_n2.getVal()

                                                                


        ## -------------------- RooFormulaVar's -------------------- ##
        rfv_n_CB = ROOT.RooFormulaVar()
        rfv_alpha_CB = ROOT.RooFormulaVar()
        rfv_n2_CB = ROOT.RooFormulaVar()
        rfv_alpha2_CB = ROOT.RooFormulaVar()
        rfv_mean_CB = ROOT.RooFormulaVar()
        rfv_sigma_CB = ROOT.RooFormulaVar()

        name = "CMS_zz4l_n_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        if self.isHighMass : rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))
        else : rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))
        
        name = "CMS_zz4l_alpha_{0:.0f}_centralValue".format(self.channel)
        if self.isHighMass : rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape_HM'], ROOT.RooArgList(self.MH))
        else : rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape'], ROOT.RooArgList(self.MH))
        
        name = "CMS_zz4l_n2_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        #if self.isHighMass : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
        #else : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
        if self.isHighMass : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape_HM']+")",ROOT.RooArgList(self.MH))
        else : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape']+")",ROOT.RooArgList(self.MH))
        
        name = "CMS_zz4l_alpha2_{0:.0f}_centralValue".format(self.channel)
        if self.isHighMass : rfv_alpha2_CB = ROOT.RooFormulaVar(name,theInputs['alpha2_CB_shape_HM'], ROOT.RooArgList(self.MH))
        else : rfv_alpha2_CB = ROOT.RooFormulaVar(name,theInputs['alpha2_CB_shape'], ROOT.RooArgList(self.MH))
        
        name = "CMS_zz4l_mean_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        
            
        if (self.channel == self.ID_4mu) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_m_err))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_m_err))
        elif (self.channel == self.ID_4e) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_e_sig,CMS_zz4l_mean_e_err))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_e_sig,CMS_zz4l_mean_e_err))
        elif (self.channel == self.ID_2e2mu) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+ (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_e_sig,CMS_zz4l_mean_m_err,CMS_zz4l_mean_e_err))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+ (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_e_sig,CMS_zz4l_mean_m_err,CMS_zz4l_mean_e_err))
        

        
        name = "CMS_zz4l_sigma_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        
        if (self.channel == self.ID_4mu) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig))
        elif (self.channel == self.ID_4e) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_e_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_e_sig))
        elif (self.channel == self.ID_2e2mu) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*TMath::Sqrt((1+@1)*(1+@2))", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*TMath::Sqrt((1+@1)*(1+@2))", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))

        name = "CMS_zz4l_gamma_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        rfv_gamma_BW = ROOT.RooFormulaVar(name,"("+theInputs['gamma_BW_shape_HM']+")"+"*(1+@1*0.05)",ROOT.RooArgList(self.MH,CMS_zz4l_gamma))

        if (DEBUG): print " DEBUG *********  ", theInputs['sigma_CB_shape'] 

        print "n_CB ", rfv_n_CB.getVal()
        print "alpha_CB ", rfv_alpha_CB.getVal()
        print "n2_CB ", rfv_n2_CB.getVal()
        print "alpha2_CB ", rfv_alpha2_CB.getVal()
        print "mean_CB ", rfv_mean_CB.getVal()
        print "sigma_CB ", rfv_sigma_CB.getVal()
        print "gamma_BW ", rfv_gamma_BW.getVal()    

        CMS_zz4l_mean_sig_NoConv = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts),"@0+@1", ROOT.RooArgList(rfv_mean_CB, self.MH))

        print "mean_sig_NoConv ", CMS_zz4l_mean_sig_NoConv.getVal()

        
        
        ## --------------------- SHAPE FUNCTIONS ---------------------- ##
    
        signalCB_ggH = ROOT.RooDoubleCB("signalCB_ggH","signalCB_ggH",CMS_zz4l_mass, self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB, self.bUseCBnoConvolution) , rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ggH = ROOT.RooRelBWUFParam("signalBW_ggH", "signalBW_ggH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ggH =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",CMS_zz4l_mass,signalBW_ggH,signalCB_ggH, 2)
        #High mass pdf
        signalBW_ggH_HM = ROOT.RooRelBWHighMass("signalBW_ggH", "signalBW_ggH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ggH_HM =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",CMS_zz4l_mass,signalBW_ggH_HM,signalCB_ggH, 2)
  
        
        signalCB_VBF = ROOT.RooDoubleCB("signalCB_VBF","signalCB_VBF",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_VBF = ROOT.RooRelBWUFParam("signalBW_VBF", "signalBW_VBF",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_VBF = ROOT.RooFFTConvPdf("sig_VBF","BW (X) CB",CMS_zz4l_mass,signalBW_VBF,signalCB_VBF, 2)
        #High mass pdf
        signalBW_VBF_HM = ROOT.RooRelBWHighMass("signalBW_VBF", "signalBW_VBF",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_VBF_HM = ROOT.RooFFTConvPdf("sig_VBF","BW (X) CB",CMS_zz4l_mass,signalBW_VBF_HM,signalCB_VBF, 2)
                       
        
        signalCB_WH = ROOT.RooDoubleCB("signalCB_WH","signalCB_WH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_WH = ROOT.RooRelBWUFParam("signalBW_WH", "signalBW_WH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_WH = ROOT.RooFFTConvPdf("sig_WH","BW (X) CB",CMS_zz4l_mass,signalBW_WH,signalCB_WH, 2)
        #High mass pdf
        signalBW_WH_HM = ROOT.RooRelBWHighMass("signalBW_WH", "signalBW_WH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_WH_HM = ROOT.RooFFTConvPdf("sig_WH","BW (X) CB",CMS_zz4l_mass,signalBW_WH_HM,signalCB_WH, 2)

        
        signalCB_ZH = ROOT.RooDoubleCB("signalCB_ZH","signalCB_ZH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ZH = ROOT.RooRelBWUFParam("signalBW_ZH", "signalBW_ZH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ZH = ROOT.RooFFTConvPdf("sig_ZH","BW (X) CB",CMS_zz4l_mass,signalBW_ZH,signalCB_ZH, 2)
        #High mass pdf
        signalBW_ZH_HM = ROOT.RooRelBWHighMass("signalBW_ZH", "signalBW_ZH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ZH_HM = ROOT.RooFFTConvPdf("sig_ZH","BW (X) CB",CMS_zz4l_mass,signalBW_ZH_HM,signalCB_ZH, 2)

        
        signalCB_ttH = ROOT.RooDoubleCB("signalCB_ttH","signalCB_ttH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ttH = ROOT.RooRelBWUFParam("signalBW_ttH", "signalBW_ttH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ttH = ROOT.RooFFTConvPdf("sig_ttH","BW (X) CB",CMS_zz4l_mass,signalBW_ttH,signalCB_ttH, 2) 
        #High mass pdf
        signalBW_ttH_HM = ROOT.RooRelBWHighMass("signalBW_ttH", "signalBW_ttH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ttH_HM = ROOT.RooFFTConvPdf("sig_ttH","BW (X) CB",CMS_zz4l_mass,signalBW_ttH_HM,signalCB_ttH, 2)
        
        
        ## Buffer fraction for cyclical behavior
        sig_ggH.setBufferFraction(0.2)
        sig_VBF.setBufferFraction(0.2)
        sig_WH.setBufferFraction(0.2)
        sig_ZH.setBufferFraction(0.2)
        sig_ttH.setBufferFraction(0.2)
        
        sig_ggH_HM.setBufferFraction(0.2)
        sig_VBF_HM.setBufferFraction(0.2)
        sig_WH_HM.setBufferFraction(0.2)
        sig_ZH_HM.setBufferFraction(0.2)
        sig_ttH_HM.setBufferFraction(0.2)

        ## -------------------- 2D SIGNAL SHAPES ------------------------- ##

        print '2D signal shapes'
        mytemplateDir = "{1}/{0:.0f}TeV".format(self.sqrts,self.templateDir)

        signalHistFuncs_Nominal = ROOT.RooArgList()
        signalHistFuncs_TxyUp = ROOT.RooArgList()
        signalHistFuncs_TxyDown = ROOT.RooArgList()

        for tt in range(0,nctau) :
          val_ctau = (x_max-x_min) / (nctau-1) * tt + x_min
          signalTemplates = "{0}_templates_TxyUpDown_CTau{1:.0f}_Modified.root".format(self.appendName,val_ctau)
          templateSigName = "{0}/{1}".format(mytemplateDir,signalTemplates)
          sigTempFile = ROOT.TFile(templateSigName)
          Sig_T = sigTempFile.Get("T_2D_TxyNominal")
          Sig_T.SetName("T_ZZ_{0:.0f}_{1}_KD_{2:.0f}".format(self.sqrts,self.appendName,val_ctau))
          Sig_T_TxyUp = sigTempFile.Get("T_2D_TxyUp")
          Sig_T_TxyUp.SetName("T_ZZ_{0:.0f}_{1}_KD_{2:.0f}_TxyUp".format(self.sqrts,self.appendName,val_ctau))
          Sig_T_TxyDown = sigTempFile.Get("T_2D_TxyDown")
          Sig_T_TxyDown.SetName("T_ZZ_{0:.0f}_{1}_KD_{2:.0f}_TxyDown".format(self.sqrts,self.appendName,val_ctau))

          if t==0:
            dBinsX = Sig_T.GetXaxis().GetNbins()
            dLowX = Sig_T.GetXaxis().GetXmin()
            dHighX = Sig_T.GetXaxis().GetXmax()
            D1.setRange(dLowX,dHighX)
            D1.setBins(dBinsX)
            T_integralName = "normt1_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            T_integral = ROOT.RooConstVar (T_integralName,T_integralName,Sig_T.Integral())
            print "T ",T_integral.getVal()

#            r_fai_pures_norm_Name = "sig_PuresNorm_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#            r_fai_realints_norm_Name = "sig_RealIntsNorm_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#            r_fai_imagints_norm_Name = "sig_ImagIntsNorm_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#            r_fai_pures_norm = ROOT.RooFormulaVar(r_fai_pures_norm_Name,r_fai_pures_norm_Name,"( (1-abs(@0)-abs(@1))*@2+abs(@0)*@3+abs(@1)*@4 )/@2",RooArgList(x,y,T1_integral,T2_integral,T3_integral))
#            r_fai_realints_norm = ROOT.RooFormulaVar(r_fai_realints_norm_Name,r_fai_realints_norm_Name,"( sign(@0)*sqrt(abs(@0)*(1-abs(@0)-abs(@1)))*cos(@2)*@4 + sign(@1)*sqrt (abs(@1)*(1-abs(@0)-abs(@1)))*cos(@3)*@5 + sign(@0*@1)*sqrt(abs(@1)*abs(@0))*cos(@3-@2)*@6 )/@7",RooArgList(x,y,phix,phiy,T4_integral,T5_integral,T6_integral,T1_integral))
#            r_fai_imagints_norm = ROOT.RooFormulaVar(r_fai_imagints_norm_Name,r_fai_imagints_norm_Name,"( sign(@0)*sqrt(abs(@0)*(1-abs(@0)-abs(@1)))*sin(@2)*@4 + sign(@1)*sqrt (abs(@1)*(1-abs(@0)-abs(@1)))*sin(@3)*@5 + sign(@0*@1)*sqrt(abs(@1)*abs(@0))*sin(@3-@2)*@6 )/@7",RooArgList(x,y,phix,phiy,T7_integral,T8_integral,T9_integral,T1_integral))
#            r_fai_norm = ROOT.RooFormulaVar("ggH_norm","ggH_norm","(abs(@3)+abs(@4))>1 ? 0. : TMath::Max((@0+@1+@2)*(1-abs(@5)),0)",RooArgList(r_fai_pures_norm,r_fai_realints_norm,r_fai_imagints_norm,x,y,alpha_zz4l))

          Sig_T_hist = ROOT.RooDataHist ("T_hist_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel,self.sqrts,val_ctau),"", ROOT.RooArgList(D1,D2),Sig_T)
          Sig_T_TxyUp_hist = ROOT.RooDataHist ("T_TxyUp_hist_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel,self.sqrts,val_ctau),"", ROOT.RooArgList(D1,D2),Sig_T_TxyUp)
          Sig_T_TxyDown_hist = ROOT.RooDataHist ("T_TxyDown_hist_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel,self.sqrts,val_ctau),"", ROOT.RooArgList(D1,D2),Sig_T_TxyDown)

          Sig_T_histfunc = ROOT.RooHistFunc ("T_histfunc_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel,self.sqrts,val_ctau),"", ROOT.RooArgSet(D1,D2),Sig_T_hist)
          Sig_T_TxyUp_histfunc = ROOT.RooHistFunc ("T_TxyUp_histfunc_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel,self.sqrts,val_ctau),"", ROOT.RooArgSet(D1,D2),Sig_T_TxyUp_hist)
          Sig_T_TxyDown_histfunc = ROOT.RooHistFunc ("T_TxyDown_histfunc_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel,self.sqrts,val_ctau),"", ROOT.RooArgSet(D1,D2),Sig_T_TxyDown_hist)

          signalHistFuncs_Nominal.add(Sig_T_histfunc)
          signalHistFuncs_TxyUp.add(Sig_T_TxyUp_histfunc)
          signalHistFuncs_TxyDown.add(Sig_T_TxyDown_histfunc)

        ggHpdfName_TxyNominal = "ggH_RooCTauPdf_1D_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggHpdf_TxyNominal = ROOT.HZZ4L_RooCTauPdf_1D(ggHpdfName_TxyNominal,ggHpdfName_TxyNominal,D1,x,signalHistFuncs_Nominal,x_min,x_max)
        ggHpdfName_TxyUp = "ggH_RooCTauPdf_1D_TxyUp_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggHpdf_TxyUp = ROOT.HZZ4L_RooCTauPdf_1D(ggHpdfName_TxyUp,ggHpdfName_TxyUp,D1,x,signalHistFuncs_TxyUp,x_min,x_max)
        ggHpdfName_TxyDown = "ggH_RooCTauPdf_1D_TxyDown_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggHpdf_TxyDown = ROOT.HZZ4L_RooCTauPdf_1D(ggHpdfName_TxyDown,ggHpdfName_TxyDown,D1,x,signalHistFuncs_TxyDown,x_min,x_max)

        funcList_Txy_ggH = ROOT.RooArgList()
        funcList_Txy_ggH.add(ggHpdf_TxyNominal)
        funcList_Txy_ggH.add(ggHpdf_TxyUp)
        funcList_Txy_ggH.add(ggHpdf_TxyDown)
        ggHpdf_Txy_Name = "ggH_Txy_TemplateMorphPdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggH_Txy_pdf = ROOT.FastVerticalInterpHistPdf(ggHpdf_Txy_Name, ggHpdf_Txy_Name, D1,False, funcList_Txy_ggH, morphVarList_Txy,1.0,1)

        signalTemplates_ScaleRes = "{0}_templates_SignalScaleResSyst.root".format(self.appendName)
        templateSigScaleResName = "{0}/{1}".format(mytemplateDir,signalTemplates_ScaleRes)
        sigTempFile_ScaleRes = ROOT.TFile(templateSigScaleResName)
        Sig_T_ScaleResNominal = sigTempFile_ScaleRes.Get("T_2D_ScaleResNominal")
        Sig_T_ScaleResNominal.SetName("T_ZZ_{0:.0f}_{1}_smd".format(self.sqrts,self.appendName))
        Sig_T_ScaleResUp = sigTempFile_ScaleRes.Get("T_2D_ScaleResUp")
        Sig_T_ScaleResUp.SetName("T_ZZ_{0:.0f}_{1}_smd_ScaleResUp".format(self.sqrts,self.appendName))
        Sig_T_ScaleResDown = sigTempFile_ScaleRes.Get("T_2D_ScaleResDown")
        Sig_T_ScaleResDown.SetName("T_ZZ_{0:.0f}_{1}_smd_ScaleResDown".format(self.sqrts,self.appendName))

        dBinsY = Sig_T_ScaleResNominal.GetYaxis().GetNbins()
        dLowY = Sig_T_ScaleResNominal.GetYaxis().GetXmin()
        dHighY = Sig_T_ScaleResNominal.GetYaxis().GetXmax()
        D2.setRange(dLowY,dHighY)
        D2.setBins(dBinsY)

        Sig_T_ScaleResNominal_hist = ROOT.RooDataHist("T_ScaleResNominal_hist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts,"",ROOT.RooArgList(D2),Sig_T_ScaleResNominal)
        Sig_T_ScaleResUp_hist = ROOT.RooDataHist("T_ScaleResUp_hist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts),"",ROOT.RooArgList(D2),Sig_T_ScaleResUp)
        Sig_T_ScaleResDown_hist = ROOT.RooDataHist("T_ScaleResDown_hist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts),"",ROOT.RooArgList(D2),Sig_T_ScaleResDown) 

        Sig_T_ScaleResNominal_pdf = ROOT.RooHistPdf("T_ScaleResNominal_histpdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts),"", ROOT.RooArgSet(D2),Sig_T_ScaleResNominal_hist)
        Sig_T_ScaleResUp_pdf = ROOT.RooHistPdf("T_ScaleResUp_histpdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts),"", ROOT.RooArgSet(D2),Sig_T_ScaleResUp_hist)
        Sig_T_ScaleResDown_pdf = ROOT.RooHistPdf("T_ScaleResDown_histpdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts),"", ROOT.RooArgSet(D2),Sig_T_ScaleResDown_hist)

        ggHpdfName = "ggH_RooProdPdf_Nominal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggHpdfName_systUp = "ggH_RooProdPdf_ScaleResUp_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggHpdfName_systDown = "ggH_RooProdPdf_ScaleResDown_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        ggHpdf = ROOT.RooProdPdf(ggHpdfName,ggHpdfName,ggH_Txy_pdf,Sig_T_ScaleResNominal_pdf)
        ggHpdf_systUp = ROOT.RooProdPdf(ggHpdfName_systUp,ggHpdfName_systUp,ggH_Txy_pdf,Sig_T_ScaleResUp_pdf)
        ggHpdf_systDown = ROOT.RooProdPdf(ggHpdfName_systDown,ggHpdfName_systDown,ggH_Txy_pdf,Sig_T_ScaleResDown_pdf)



        ## ------------------ END 2D SIGNAL SHAPES ------------------------ ##


        ## -------------------------- BACKGROUND SHAPES ---------------------------------- ##
    
        ## qqZZ contribution
        name = "CMS_qqzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a0 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a0",115.3,0.,200.)
        name = "CMS_qqzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a1 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a1",21.96,0.,200.)
        name = "CMS_qqzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a2 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a2",122.8,0.,200.)
        name = "CMS_qqzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a3 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a3",0.03479,0.,1.)
        name = "CMS_qqzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a4 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a4",185.5,0.,200.)
        name = "CMS_qqzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a5 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a5",12.67,0.,200.)
        name = "CMS_qqzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a6 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a6",34.81,0.,100.)
        name = "CMS_qqzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a7 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a7",0.1393,0.,1.)
        name = "CMS_qqzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a8 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a8",66.,0.,200.)
        name = "CMS_qqzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a9 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a9",0.07191,0.,1.)
        name = "CMS_qqzzbkg_a10_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a10 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a10",94.11,0.,200.)
        name = "CMS_qqzzbkg_a11_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a11 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a11",-5.111,-100.,100.)
        name = "CMS_qqzzbkg_a12_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a12 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a12",4834,0.,10000.)
        name = "CMS_qqzzbkg_a13_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a13 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a13",0.2543,0.,1.)
        
        

        if (DEBUG) :
            print "qqZZshape_a0 = ",theInputs['qqZZshape_a0']
            print "qqZZshape_a1 = ",theInputs['qqZZshape_a1']
            print "qqZZshape_a2 = ",theInputs['qqZZshape_a2']
            print "qqZZshape_a3 = ",theInputs['qqZZshape_a3']
            print "qqZZshape_a4 = ",theInputs['qqZZshape_a4']
            print "qqZZshape_a5 = ",theInputs['qqZZshape_a5']
            print "qqZZshape_a6 = ",theInputs['qqZZshape_a6']
            print "qqZZshape_a7 = ",theInputs['qqZZshape_a7']
            print "qqZZshape_a8 = ",theInputs['qqZZshape_a8']
            print "qqZZshape_a9 = ",theInputs['qqZZshape_a9']
            print "qqZZshape_a10 = ",theInputs['qqZZshape_a10']
            print "qqZZshape_a11 = ",theInputs['qqZZshape_a11']
            print "qqZZshape_a12 = ",theInputs['qqZZshape_a12']
            print "qqZZshape_a13 = ",theInputs['qqZZshape_a13']

        
        CMS_qqzzbkg_a0.setVal(theInputs['qqZZshape_a0'])
        CMS_qqzzbkg_a1.setVal(theInputs['qqZZshape_a1'])
        CMS_qqzzbkg_a2.setVal(theInputs['qqZZshape_a2'])
        CMS_qqzzbkg_a3.setVal(theInputs['qqZZshape_a3'])
        CMS_qqzzbkg_a4.setVal(theInputs['qqZZshape_a4'])
        CMS_qqzzbkg_a5.setVal(theInputs['qqZZshape_a5'])
        CMS_qqzzbkg_a6.setVal(theInputs['qqZZshape_a6'])
        CMS_qqzzbkg_a7.setVal(theInputs['qqZZshape_a7'])
        CMS_qqzzbkg_a8.setVal(theInputs['qqZZshape_a8'])
        CMS_qqzzbkg_a9.setVal(theInputs['qqZZshape_a9'])
        CMS_qqzzbkg_a10.setVal(theInputs['qqZZshape_a10'])
        CMS_qqzzbkg_a11.setVal(theInputs['qqZZshape_a11'])
        CMS_qqzzbkg_a12.setVal(theInputs['qqZZshape_a12'])
        CMS_qqzzbkg_a13.setVal(theInputs['qqZZshape_a13'])
        
        CMS_qqzzbkg_a0.setConstant(True)
        CMS_qqzzbkg_a1.setConstant(True)
        CMS_qqzzbkg_a2.setConstant(True)
        CMS_qqzzbkg_a3.setConstant(True)
        CMS_qqzzbkg_a4.setConstant(True)
        CMS_qqzzbkg_a5.setConstant(True)
        CMS_qqzzbkg_a6.setConstant(True)
        CMS_qqzzbkg_a7.setConstant(True)
        CMS_qqzzbkg_a8.setConstant(True)
        CMS_qqzzbkg_a9.setConstant(True)
        CMS_qqzzbkg_a10.setConstant(True)
        CMS_qqzzbkg_a11.setConstant(True)
        CMS_qqzzbkg_a12.setConstant(True)
        CMS_qqzzbkg_a13.setConstant(True)
        
        bkg_qqzz = ROOT.RooqqZZPdf_v2("bkg_qqzzTmp","bkg_qqzzTmp",CMS_zz4l_mass,CMS_qqzzbkg_a0,CMS_qqzzbkg_a1,CMS_qqzzbkg_a2,CMS_qqzzbkg_a3,CMS_qqzzbkg_a4,CMS_qqzzbkg_a5,CMS_qqzzbkg_a6,CMS_qqzzbkg_a7,CMS_qqzzbkg_a8,CMS_qqzzbkg_a9,CMS_qqzzbkg_a10,CMS_qqzzbkg_a11,CMS_qqzzbkg_a12,CMS_qqzzbkg_a13)
        
        ## ggZZ contribution
        name = "CMS_ggzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a0 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a0",115.3,0.,200.)
        name = "CMS_ggzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a1 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a1",21.96,0.,200.)
        name = "CMS_ggzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a2 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a2",122.8,0.,200.)
        name = "CMS_ggzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a3 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a3",0.03479,0.,1.)
        name = "CMS_ggzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
        CMS_ggzzbkg_a4 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a4",185.5,0.,200.)
        name = "CMS_ggzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a5 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a5",12.67,0.,200.)
        name = "CMS_ggzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a6 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a6",34.81,0.,100.)
        name = "CMS_ggzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a7 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a7",0.1393,0.,1.)
        name = "CMS_ggzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a8 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a8",66.,0.,200.)
        name = "CMS_ggzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
        CMS_ggzzbkg_a9 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a9",0.07191,0.,1.)
        
        
        CMS_ggzzbkg_a0.setVal(theInputs['ggZZshape_a0'])
        CMS_ggzzbkg_a1.setVal(theInputs['ggZZshape_a1'])
        CMS_ggzzbkg_a2.setVal(theInputs['ggZZshape_a2'])
        CMS_ggzzbkg_a3.setVal(theInputs['ggZZshape_a3'])
        CMS_ggzzbkg_a4.setVal(theInputs['ggZZshape_a4'])
        CMS_ggzzbkg_a5.setVal(theInputs['ggZZshape_a5'])
        CMS_ggzzbkg_a6.setVal(theInputs['ggZZshape_a6'])
        CMS_ggzzbkg_a7.setVal(theInputs['ggZZshape_a7'])
        CMS_ggzzbkg_a8.setVal(theInputs['ggZZshape_a8'])
        CMS_ggzzbkg_a9.setVal(theInputs['ggZZshape_a9'])
        
        CMS_ggzzbkg_a0.setConstant(True)
        CMS_ggzzbkg_a1.setConstant(True)
        CMS_ggzzbkg_a2.setConstant(True)
        CMS_ggzzbkg_a3.setConstant(True)
        CMS_ggzzbkg_a4.setConstant(True)
        CMS_ggzzbkg_a5.setConstant(True)
        CMS_ggzzbkg_a6.setConstant(True)
        CMS_ggzzbkg_a7.setConstant(True)
        CMS_ggzzbkg_a8.setConstant(True)
        CMS_ggzzbkg_a9.setConstant(True)

        if (DEBUG) :
            print "ggZZshape_a0 = ",theInputs['ggZZshape_a0']
            print "ggZZshape_a1 = ",theInputs['ggZZshape_a1']
            print "ggZZshape_a2 = ",theInputs['ggZZshape_a2']
            print "ggZZshape_a3 = ",theInputs['ggZZshape_a3']
            print "ggZZshape_a4 = ",theInputs['ggZZshape_a4']
            print "ggZZshape_a5 = ",theInputs['ggZZshape_a5']
            print "ggZZshape_a6 = ",theInputs['ggZZshape_a6']
            print "ggZZshape_a7 = ",theInputs['ggZZshape_a7']
            print "ggZZshape_a8 = ",theInputs['ggZZshape_a8']
            print "ggZZshape_a9 = ",theInputs['ggZZshape_a9']
                   
        
        bkg_ggzz = ROOT.RooggZZPdf_v2("bkg_ggzzTmp","bkg_ggzzTmp",CMS_zz4l_mass,CMS_ggzzbkg_a0,CMS_ggzzbkg_a1,CMS_ggzzbkg_a2,CMS_ggzzbkg_a3,CMS_ggzzbkg_a4,CMS_ggzzbkg_a5,CMS_ggzzbkg_a6,CMS_ggzzbkg_a7,CMS_ggzzbkg_a8,CMS_ggzzbkg_a9)
    
        ## Reducible backgrounds
        val_meanL_3P1F = float(theInputs['zjetsShape_mean_3P1F'])
        val_sigmaL_3P1F = float(theInputs['zjetsShape_sigma_3P1F'])
        val_normL_3P1F = float(theInputs['zjetsShape_norm_3P1F'])
        
        val_meanL_2P2F = float(theInputs['zjetsShape_mean_2P2F'])
        val_sigmaL_2P2F = float(theInputs['zjetsShape_sigma_2P2F'])
        val_normL_2P2F = float(theInputs['zjetsShape_norm_2P2F'])
        val_pol0_2P2F = float(theInputs['zjetsShape_pol0_2P2F'])
        val_pol1_2P2F = float(theInputs['zjetsShape_pol1_2P2F'])
        
        val_meanL_2P2F_2 = float(theInputs['zjetsShape_mean_2P2F_2e2mu'])
        val_sigmaL_2P2F_2 = float(theInputs['zjetsShape_sigma_2P2F_2e2mu'])
        val_normL_2P2F_2 = float(theInputs['zjetsShape_norm_2P2F_2e2mu'])


        if (self.channel == self.ID_4mu):
            name = "mlZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet = ROOT.RooRealVar(name,"mean landau Zjet",val_meanL_2P2F)
            name = "slZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet = ROOT.RooRealVar(name,"sigma landau Zjet",val_sigmaL_2P2F)
            print "mean 4mu: ",mlZjet.getVal()
            print "sigma 4mu: ",slZjet.getVal()
            bkg_zjets = ROOT.RooLandau("bkg_zjetsTmp","bkg_zjetsTmp",CMS_zz4l_mass,mlZjet,slZjet)
        elif (self.channel == self.ID_4e):
            name = "mlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet_2p2f = ROOT.RooRealVar(name,"mean landau Zjet 2p2f",val_meanL_2P2F)
            name = "slZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet_2p2f = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f",val_sigmaL_2P2F)
            name = "nlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            nlZjet_2p2f = ROOT.RooRealVar(name,"norm landau Zjet 2p2f",val_normL_2P2F)
            name = "p0Zjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            p0Zjet_2p2f = ROOT.RooRealVar(name,"p0 Zjet 2p2f",val_pol0_2P2F)
            name = "p1Zjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            p1Zjet_2p2f = ROOT.RooRealVar(name,"p1 Zjet 2p2f",val_pol1_2P2F)
            print "mean 2p2f 4e: ",mlZjet_2p2f.getVal()
            print "sigma 2p2f 4e: ",slZjet_2p2f.getVal()
            print "norm 2p2f 4e: ",nlZjet_2p2f.getVal()
            print "pol0 2p2f 4e: ",p0Zjet_2p2f.getVal()
            print "pol1 2p2f 4e: ",p1Zjet_2p2f.getVal()
            bkg_zjets_2p2f = ROOT.RooGenericPdf("bkg_zjetsTmp_2p2f","bkg_zjetsTmp_2p2f","(TMath::Landau(@0,@1,@2))*(1.+ TMath::Exp(@3+@4*@0))",RooArgList(CMS_zz4l_mass,mlZjet_2p2f,slZjet_2p2f,p0Zjet_2p2f,p1Zjet_2p2f))
            
            name = "mlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet_3p1f = ROOT.RooRealVar(name,"mean landau Zjet 3p1f",val_meanL_3P1F)
            name = "slZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet_3p1f = ROOT.RooRealVar(name,"sigma landau Zjet 3p1f",val_sigmaL_3P1F)
            name = "nlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            nlZjet_3p1f = ROOT.RooRealVar(name,"norm landau Zjet 3p1f",val_normL_3P1F)
            print "mean 3p1f 4e: ",mlZjet_3p1f.getVal()
            print "sigma 3p1f 4e: ",slZjet_3p1f.getVal()
            print "norm 3p1f 4e: ",nlZjet_3p1f.getVal()
            bkg_zjets_3p1f = ROOT.RooLandau("bkg_zjetsTmp_3p1f","bkg_zjetsTmp_3p1f",CMS_zz4l_mass,mlZjet_3p1f,slZjet_3p1f)
            
            bkg_zjets = ROOT.RooAddPdf("bkg_zjetsTmp","bkg_zjetsTmp",ROOT.RooArgList(bkg_zjets_2p2f,bkg_zjets_3p1f),ROOT.RooArgList(nlZjet_2p2f,nlZjet_3p1f))
            
        elif (self.channel == self.ID_2e2mu):
            name = "mlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet_2p2f = ROOT.RooRealVar(name,"mean landau Zjet 2p2f",val_meanL_2P2F)
            name = "slZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet_2p2f = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f",val_sigmaL_2P2F)
            name = "nlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            nlZjet_2p2f = ROOT.RooRealVar(name,"norm landau Zjet 2p2f",val_normL_2P2F)
            name = "p0Zjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            p0Zjet_2p2f = ROOT.RooRealVar(name,"p0 Zjet 2p2f",val_pol0_2P2F)
            name = "p1Zjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            p1Zjet_2p2f = ROOT.RooRealVar(name,"p1 Zjet 2p2f",val_pol1_2P2F)
            print "mean 2p2f 2mu2e: ",mlZjet_2p2f.getVal()
            print "sigma 2p2f 2mu2e: ",slZjet_2p2f.getVal()
            print "norm 2p2f 2mu2e: ",nlZjet_2p2f.getVal()
            print "pol0 2p2f 2mu2e: ",p0Zjet_2p2f.getVal()
            print "pol1 2p2f 2mu2e: ",p1Zjet_2p2f.getVal()
            bkg_zjets_2p2f = ROOT.RooGenericPdf("bkg_zjetsTmp_2p2f","bkg_zjetsTmp_2p2f","(TMath::Landau(@0,@1,@2))*(1.+ TMath::Exp(@3+@4*@0))",RooArgList(CMS_zz4l_mass,mlZjet_2p2f,slZjet_2p2f,p0Zjet_2p2f,p1Zjet_2p2f))
            
            name = "mlZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet_2p2f_2 = ROOT.RooRealVar(name,"mean landau Zjet 2p2f 2e2mu",val_meanL_2P2F_2)
            name = "slZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet_2p2f_2 = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f 2e2mu",val_sigmaL_2P2F_2)
            name = "nlZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            nlZjet_2p2f_2 = ROOT.RooRealVar(name,"norm landau Zjet 2p2f 2e2mu",val_normL_2P2F_2)
            print "mean 2p2f 2e2mu: ",mlZjet_2p2f_2.getVal()
            print "sigma 2p2f 2e2mu: ",slZjet_2p2f_2.getVal()
            print "norm 2p2f 2e2mu: ",nlZjet_2p2f_2.getVal()
            bkg_zjets_2p2f_2 = ROOT.RooLandau("bkg_zjetsTmp_2p2f_2","bkg_zjetsTmp_2p2f_2",CMS_zz4l_mass,mlZjet_2p2f_2,slZjet_2p2f_2)
            
            name = "mlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet_3p1f = ROOT.RooRealVar(name,"mean landau Zjet 3p1f",val_meanL_3P1F)
            name = "slZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet_3p1f = ROOT.RooRealVar(name,"sigma landau Zjet 3p1f",val_sigmaL_3P1F)
            name = "nlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            nlZjet_3p1f = ROOT.RooRealVar(name,"norm landau Zjet 3p1f",val_normL_3P1F)
            print "mean 3p1f 2mu2e: ",mlZjet_3p1f.getVal()
            print "sigma 3p1f 2mu2e: ",slZjet_3p1f.getVal()
            print "norm 3p1f 2mu2e: ",nlZjet_3p1f.getVal()
            bkg_zjets_3p1f = ROOT.RooLandau("bkg_zjetsTmp_3p1f","bkg_zjetsTmp_3p1f",CMS_zz4l_mass,mlZjet_3p1f,slZjet_3p1f)
            
            bkg_zjets = ROOT.RooAddPdf("bkg_zjetsTmp","bkg_zjetsTmp",ROOT.RooArgList(bkg_zjets_2p2f,bkg_zjets_3p1f,bkg_zjets_2p2f_2),ROOT.RooArgList(nlZjet_2p2f,nlZjet_3p1f,nlZjet_2p2f_2))

        ## ------------------ 2D BACKGROUND SHAPES ------------------- ##

        print '2D background shapes'
        bkgTemplates = "{0}_templates_Merged_bkg.root".format(self.appendName)
        templateBkgName = "{0}/{1}".format(mytemplateDir,bkgTemplates)
        bkgTempFile = ROOT.TFile(templateBkgName)
        qqZZTemplate = bkgTempFile.Get("template_qqZZ_Dbkg")
        TemplateName = "qqZZTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqZZTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D2),qqZZTemplate)
        PdfName = "qqZZ_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqZZTemplatePdf = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D2),qqZZTempDataHist)

        ggZZTemplate = bkgTempFile.Get("template_ggZZ_Dbkg")
        TemplateName = "ggZZTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D2),ggZZTemplate)
        PdfName = "ggZZ_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZTemplatePdf = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D2),ggZZTempDataHist)

        ZjetsTemplate = bkgTempFile.Get("template_ZX_Dbkg_Nominal")
        TemplateName = "ZjetsTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D2),ZjetsTemplate)
        PdfName = "Zjets_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTemplatePdf = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D2),ZjetsTempDataHist)

        ZjetsTemplateDown = bkgTempFile.Get("template_ZX_Dbkg_ScaleResDown")
        TemplateName = "ZjetsTempDownDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTempDataHistDown = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D2),ZjetsTemplateDown)
        PdfName = "Zjets_TemplateDownPdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTemplatePdfDown = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D2),ZjetsTempDataHistDown)

        ZjetsTemplateUp = bkgTempFile.Get("template_ZX_Dbkg_ScaleResUp")
        TemplateName = "ZjetsTempUpDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTempDataHistUp = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D2),ZjetsTemplateUp)
        PdfName = "Zjets_TemplateUpPdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTemplatePdfUp = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D2),ZjetsTempDataHistUp)

        funcList_zjets = ROOT.RooArgList()
        morphBkgVarName =  "CMS_zz4l_smd_zjets_bkg_{0:.0f}".format(self.channel)
        alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        morphVarListBkg = ROOT.RooArgList()

        if(self.bkgMorph):
            funcList_zjets.add(ZjetsTemplatePdf)
            funcList_zjets.add(ZjetsTemplatePdfUp)
            funcList_zjets.add(ZjetsTemplatePdfDown)
            alphaMorphBkg.setConstant(False)
            morphVarListBkg.add(alphaMorphBkg)
        else:
            funcList_zjets.add(ZjetsTemplatePdf)
            alphaMorphBkg.setConstant(True)
        MorphName = "ZX_TemplateMorphPdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTemplateMorphPdf = ROOT.FastVerticalInterpHistPdf(MorphName,MorphName,D2,False,funcList_zjets,morphVarListBkg,1.0,1)


#        bkgTemplates_TxyNominal = "{0}_templates_bkg_Nominal.root".format(self.appendName)
#        templateBkgName_TxyNominal = "{0}/{1}".format(mytemplateDir,bkgTemplates_TxyNominal)
#        bkgTempFile_TxyNominal = ROOT.TFile(templateBkgName_TxyNominal)
        bkgTempFile_TxyNominal = bkgTempFile
        qqZZTemplate_TxyNominal = bkgTempFile_TxyNominal.Get("template_qqZZ_TxyNominal")
        TemplateName = "qqZZTempDataHist_TxyNominal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqZZTempDataHist_TxyNominal = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1),qqZZTemplate_TxyNominal)
        PdfName = "qqZZ_TemplatePdf_TxyNominal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqZZTemplatePdf_TxyNominal = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1),qqZZTempDataHist_TxyNominal)
        ggZZTemplate_TxyNominal = bkgTempFile_TxyNominal.Get("template_ggZZ_TxyNominal")
        TemplateName = "ggZZTempDataHist_TxyNominal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZTempDataHist_TxyNominal = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1),ggZZTemplate_TxyNominal)
        PdfName = "ggZZ_TemplatePdf_TxyNominal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZTemplatePdf_TxyNominal = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1),ggZZTempDataHist_TxyNominal)
        ZjetsTemplate_TxyNominal = bkgTempFile_TxyNominal.Get("template_ZX_TxyNominal")
        TemplateName = "ZjetsTempDataHist_TxyNominal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTempDataHist_TxyNominal = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1),ZjetsTemplate_TxyNominal)
        PdfName = "Zjets_TemplatePdf_TxyNominal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTemplatePdf_TxyNominal = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1),ZjetsTempDataHist_TxyNominal)

#        bkgTemplates_TxyUp = "{0}_templates_bkg_TxyUp.root".format(self.appendName)
#        templateBkgName_TxyUp = "{0}/{1}".format(mytemplateDir,bkgTemplates_TxyUp)
#        bkgTempFile_TxyUp = ROOT.TFile(templateBkgName_TxyUp)
        bkgTempFile_TxyUp = bkgTempFile_TxyNominal
        qqZZTemplate_TxyUp = bkgTempFile_TxyUp.Get("template_qqZZ_TxyUp")
        TemplateName = "qqZZTempDataHist_TxyUp_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqZZTempDataHist_TxyUp = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1),qqZZTemplate_TxyUp)
        PdfName = "qqZZ_TemplatePdf_TxyUp_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqZZTemplatePdf_TxyUp = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1),qqZZTempDataHist_TxyUp)
        ggZZTemplate_TxyUp = bkgTempFile_TxyUp.Get("template_ggZZ_TxyUp")
        TemplateName = "ggZZTempDataHist_TxyUp_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZTempDataHist_TxyUp = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1),ggZZTemplate_TxyUp)
        PdfName = "ggZZ_TemplatePdf_TxyUp_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZTemplatePdf_TxyUp = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1),ggZZTempDataHist_TxyUp)
        ZjetsTemplate_TxyUp = bkgTempFile_TxyUp.Get("template_ZX_TxyUp")
        TemplateName = "ZjetsTempDataHist_TxyUp_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTempDataHist_TxyUp = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1),ZjetsTemplate_TxyUp)
        PdfName = "Zjets_TemplatePdf_TxyUp_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTemplatePdf_TxyUp = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1),ZjetsTempDataHist_TxyUp)

#        bkgTemplates_TxyDown = "{0}_templates_bkg_TxyDown.root".format(self.appendName)
#        templateBkgName_TxyDown = "{0}/{1}".format(mytemplateDir,bkgTemplates_TxyDown)
#        bkgTempFile_TxyDown = ROOT.TFile(templateBkgName_TxyDown)
        bkgTempFile_TxyDown = bkgTempFile_TxyNominal
        qqZZTemplate_TxyDown = bkgTempFile_TxyDown.Get("template_qqZZ_TxyDown")
        TemplateName = "qqZZTempDataHist_TxyDown_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqZZTempDataHist_TxyDown = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1),qqZZTemplate_TxyDown)
        PdfName = "qqZZ_TemplatePdf_TxyDown_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqZZTemplatePdf_TxyDown = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1),qqZZTempDataHist_TxyDown)
        ggZZTemplate_TxyDown = bkgTempFile_TxyDown.Get("template_ggZZ_TxyDown")
        TemplateName = "ggZZTempDataHist_TxyDown_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZTempDataHist_TxyDown = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1),ggZZTemplate_TxyDown)
        PdfName = "ggZZ_TemplatePdf_TxyDown_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZTemplatePdf_TxyDown = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1),ggZZTempDataHist_TxyDown)
        ZjetsTemplate_TxyDown = bkgTempFile_TxyDown.Get("template_ZX_TxyDown")
        TemplateName = "ZjetsTempDataHist_TxyDown_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTempDataHist_TxyDown = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1),ZjetsTemplate_TxyDown)
        PdfName = "Zjets_TemplatePdf_TxyDown_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ZjetsTemplatePdf_TxyDown = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1),ZjetsTempDataHist_TxyDown)


        funcList_Txy_qqZZ = ROOT.RooArgList()
        funcList_Txy_ggZZ = ROOT.RooArgList()
        funcList_Txy_Zjets = ROOT.RooArgList()
        funcList_Txy_qqZZ.add(qqZZTemplatePdf_TxyNominal)
        funcList_Txy_qqZZ.add(qqZZTemplatePdf_TxyUp)
        funcList_Txy_qqZZ.add(qqZZTemplatePdf_TxyDown)
        funcList_Txy_ggZZ.add(ggZZTemplatePdf_TxyNominal)
        funcList_Txy_ggZZ.add(ggZZTemplatePdf_TxyUp)
        funcList_Txy_ggZZ.add(ggZZTemplatePdf_TxyDown)
        funcList_Txy_Zjets.add(ZjetsTemplatePdf_TxyNominal)
        funcList_Txy_Zjets.add(ZjetsTemplatePdf_TxyUp)
        funcList_Txy_Zjets.add(ZjetsTemplatePdf_TxyDown)

        qqZZpdf_Txy_Name = "qqZZ_Txy_TemplateMorphPdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqZZ_Txy_pdf = ROOT.FastVerticalInterpHistPdf(qqZZpdf_Txy_Name, qqZZpdf_Txy_Name, D1,False, funcList_Txy_qqZZ, morphVarList_Txy,1.0,1)
        qqZZpdfName = "qqZZ_RooProdPdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqZZpdf = ROOT.RooProdPdf(qqZZpdfName,qqZZpdfName,qqZZ_Txy_pdf,qqZZTemplatePdf)

        ggZZpdf_Txy_Name = "ggZZ_Txy_TemplateMorphPdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZ_Txy_pdf = ROOT.FastVerticalInterpHistPdf(ggZZpdf_Txy_Name, ggZZpdf_Txy_Name, D1,False, funcList_Txy_ggZZ, morphVarList_Txy,1.0,1)
        ggZZpdfName = "ggZZ_RooProdPdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZpdf = ROOT.RooProdPdf(ggZZpdfName,ggZZpdfName,ggZZ_Txy_pdf,ggZZTemplatePdf)

        Zjetspdf_Txy_Name = "Zjets_Txy_TemplateMorphPdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        Zjets_Txy_pdf = ROOT.FastVerticalInterpHistPdf(Zjetspdf_Txy_Name, Zjetspdf_Txy_Name, D1,False, funcList_Txy_Zjets, morphVarList_Txy,1.0,1)
        ZjetspdfName = "Zjets_RooProdPdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        Zjetspdf = ROOT.RooProdPdf(ZjetspdfName,ZjetspdfName,Zjets_Txy_pdf,ZjetsTemplateMorphPdf)


        ## ---------------- END 2D BACKGROUND SHAPES ----------------- ##

        ## ----------------------- PLOTS FOR SANITY CHECKS -------------------------- ##

        canv_name = "czz_{0}_{1}".format(self.mH,self.appendName)

        czz = ROOT.TCanvas( canv_name, canv_name, 750, 700 )
        czz.cd()
        zzframe_s = CMS_zz4l_mass.frame(45)
        if self.bUseCBnoConvolution: super(RooDoubleCB,signalCB_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        elif self.isHighMass : super(ROOT.RooFFTConvPdf,sig_ggH_HM).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        else : super(ROOT.RooFFTConvPdf,sig_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        super(ROOT.RooqqZZPdf_v2,bkg_qqzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(4) )
        super(ROOT.RooggZZPdf_v2,bkg_ggzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(6) )
        super(ROOT.RooAbsPdf,bkg_zjets).plotOn(zzframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(6) )
        zzframe_s.Draw()
        figName = "{0}/figs/mzz_{1}_{2}.png".format(self.outputDir, self.mH, self.appendName)
        czz.SaveAs(figName)
        del czz
        
        ## ------------------- LUMI -------------------- ##
        
        rrvLumi = ROOT.RooRealVar("cmshzz4l_lumi","cmshzz4l_lumi",self.lumi)  
        
        ## ----------------------- SIGNAL RATES ----------------------- ##
        
        CMS_zz4l_mass.setRange("shape",self.low_M,self.high_M)
        
        fr_low_M = self.low_M
        fr_high_M = self.high_M        
        if (self.mH >= 450): 
            fr_low_M = 100
            fr_high_M = 1000
        if (self.mH >= 750):
            fr_low_M = 100
            fr_high_M = 1400
            

        CMS_zz4l_mass.setRange("fullrangesignal",fr_low_M,fr_high_M)
        CMS_zz4l_mass.setRange("fullrange",100,1400)
        

        sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        rrva1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a1'])
        sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        rrva2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a2'])
        sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        rrva3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a3'])
        sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        rrva4 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a4'])
        sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        rrvb1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b1'])
        sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        rrvb2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b2'])
        sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        rrvb3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b3'])
        sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
        rrvg1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g1'])
        sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
        rrvg2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g2'])
        sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
        rrvg3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g3'])
        
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        rrva1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa1'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        rrva2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa2'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        rrva3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa3'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        rrva4_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa4'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        rrvb1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb1'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        rrvb2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb2'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        rrvb3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb3'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
        rrvg1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg1'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
        rrvg2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg2'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
        rrvg3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg3'])
        
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        rrva1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa1'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        rrva2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa2'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        rrva3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa3'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        rrva4_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa4'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        rrvb1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb1'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        rrvb2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb2'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        rrvb3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb3'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
        rrvg1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg1'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
        rrvg2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg2'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
        rrvg3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg3'])
        
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        rrva1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa1'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        rrva2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa2'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        rrva3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa3'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        rrva4_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa4'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        rrvb1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb1'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        rrvb2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb2'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        rrvb3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb3'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
        rrvg1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg1'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
        rrvg2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg2'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
        rrvg3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg3'])
        
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        rrva1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa1'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        rrva2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa2'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        rrva3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa3'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        rrva4_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa4'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        rrvb1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb1'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        rrvb2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb2'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        rrvb3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb3'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
        rrvg1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg1'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
        rrvg2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg2'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
        rrvg3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg3'])
        

        if(DEBUG):
            print "sigEff_a1 = ",theInputs['sigEff_a1']
            print "sigEff_a2 = ",theInputs['sigEff_a2']
            print "sigEff_a3 = ",theInputs['sigEff_a3']
            print "sigEff_a4 = ",theInputs['sigEff_a4']
            print "sigEff_b1 = ",theInputs['sigEff_b1']
            print "sigEff_b2 = ",theInputs['sigEff_b2']
            print "sigEff_b3 = ",theInputs['sigEff_b3']
            print "sigEff_g1 = ",theInputs['sigEff_g1']
            print "sigEff_g2 = ",theInputs['sigEff_g2']
            print "sigEff_g3 = ",theInputs['sigEff_g3']

            print "sigEff_qqHa1 = ",theInputs['sigEff_qqHa1']
            print "sigEff_qqHa2 = ",theInputs['sigEff_qqHa2']
            print "sigEff_qqHa3 = ",theInputs['sigEff_qqHa3']
            print "sigEff_qqHa4 = ",theInputs['sigEff_qqHa4']
            print "sigEff_qqHb1 = ",theInputs['sigEff_qqHb1']
            print "sigEff_qqHb2 = ",theInputs['sigEff_qqHb2']
            print "sigEff_qqHb3 = ",theInputs['sigEff_qqHb3']
            print "sigEff_qqHg1 = ",theInputs['sigEff_qqHg1']
            print "sigEff_qqHg2 = ",theInputs['sigEff_qqHg2']
            print "sigEff_qqHg3 = ",theInputs['sigEff_qqHg3']

            print "sigEff_ZHa1 = ",theInputs['sigEff_ZHa1']
            print "sigEff_ZHa2 = ",theInputs['sigEff_ZHa2']
            print "sigEff_ZHa3 = ",theInputs['sigEff_ZHa3']
            print "sigEff_ZHa4 = ",theInputs['sigEff_ZHa4']
            print "sigEff_ZHb1 = ",theInputs['sigEff_ZHb1']
            print "sigEff_ZHb2 = ",theInputs['sigEff_ZHb2']
            print "sigEff_ZHb3 = ",theInputs['sigEff_ZHb3']
            print "sigEff_ZHg1 = ",theInputs['sigEff_ZHg1']
            print "sigEff_ZHg2 = ",theInputs['sigEff_ZHg2']
            print "sigEff_ZHg3 = ",theInputs['sigEff_ZHg3']

            print "sigEff_WHa1 = ",theInputs['sigEff_WHa1']
            print "sigEff_WHa2 = ",theInputs['sigEff_WHa2']
            print "sigEff_WHa3 = ",theInputs['sigEff_WHa3']
            print "sigEff_WHa4 = ",theInputs['sigEff_WHa4']
            print "sigEff_WHb1 = ",theInputs['sigEff_WHb1']
            print "sigEff_WHb2 = ",theInputs['sigEff_WHb2']
            print "sigEff_WHb3 = ",theInputs['sigEff_WHb3']
            print "sigEff_WHg1 = ",theInputs['sigEff_WHg1']
            print "sigEff_WHg2 = ",theInputs['sigEff_WHg2']
            print "sigEff_WHg3 = ",theInputs['sigEff_WHg3']

            print "sigEff_ttHa1 = ",theInputs['sigEff_ttHa1']
            print "sigEff_ttHa2 = ",theInputs['sigEff_ttHa2']
            print "sigEff_ttHa3 = ",theInputs['sigEff_ttHa3']
            print "sigEff_ttHa4 = ",theInputs['sigEff_ttHa4']
            print "sigEff_ttHb1 = ",theInputs['sigEff_ttHb1']
            print "sigEff_ttHb2 = ",theInputs['sigEff_ttHb2']
            print "sigEff_ttHb3 = ",theInputs['sigEff_ttHb3']
            print "sigEff_ttHg1 = ",theInputs['sigEff_ttHg1']
            print "sigEff_ttHg2 = ",theInputs['sigEff_ttHg2']
            print "sigEff_ttHg3 = ",theInputs['sigEff_ttHg3']
           

        sigEffName_ggH = "hzz4lggHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigEffName_qqH = "hzz4lqqHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigEffName_WH = "hzz4lWHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigEffName_ZH = "hzz4lZHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigEffName_ttH = "hzz4lttHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        listEff = ROOT.RooArgList(rrva1,rrva2,rrva3,rrva4,rrvb1,rrvb2,rrvb3,self.MH)
        listEff.add(rrvg1)
        listEff.add(rrvg2)
        listEff.add(rrvg3)

        listEff_qqh = ROOT.RooArgList(rrva1_qqh,rrva2_qqh,rrva3_qqh,rrva4_qqh,rrvb1_qqh,rrvb2_qqh,rrvb3_qqh,self.MH)
        listEff_qqh.add(rrvg1_qqh)
        listEff_qqh.add(rrvg2_qqh)
        listEff_qqh.add(rrvg3_qqh)

        listEff_wh = ROOT.RooArgList(rrva1_wh,rrva2_wh,rrva3_wh,rrva4_wh,rrvb1_wh,rrvb2_wh,rrvb3_wh,self.MH)
        listEff_wh.add(rrvg1_wh)
        listEff_wh.add(rrvg2_wh)
        listEff_wh.add(rrvg3_wh)

        listEff_zh = ROOT.RooArgList(rrva1_zh,rrva2_zh,rrva3_zh,rrva4_zh,rrvb1_zh,rrvb2_zh,rrvb3_zh,self.MH)
        listEff_zh.add(rrvg1_zh)
        listEff_zh.add(rrvg2_zh)
        listEff_zh.add(rrvg3_zh)

        listEff_tth = ROOT.RooArgList(rrva1_tth,rrva2_tth,rrva3_tth,rrva4_tth,rrvb1_tth,rrvb2_tth,rrvb3_tth,self.MH)
        listEff_tth.add(rrvg1_tth)
        listEff_tth.add(rrvg2_tth)
        listEff_tth.add(rrvg3_tth)
        
        rfvSigEff_ggH = ROOT.RooFormulaVar(sigEffName_ggH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",listEff) #ROOT.RooArgList(rrva1,rrva2,rrva3,rrva4,rrvb1,rrvb2,rrvb3,self.MH,rrvg1,rrvg2,rrvg3))

        rfvSigEff_qqH = ROOT.RooFormulaVar(sigEffName_qqH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",listEff_qqh)
        rfvSigEff_ZH = ROOT.RooFormulaVar(sigEffName_ZH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",listEff_zh)
        rfvSigEff_WH = ROOT.RooFormulaVar(sigEffName_WH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",listEff_wh)
        rfvSigEff_ttH = ROOT.RooFormulaVar(sigEffName_ttH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",listEff_tth)
        #from TF1 *polyFunc= new TF1("polyFunc","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x)+[7]*TMath::Gaus(x,[8],[9])", 110., xMax);
        
        ## following printout is needed ,  dont remove it
        print " @@@@@@@@ ggHeff ",rfvSigEff_ggH.getVal()
        print " @@@@@@@@ qqHeff ",rfvSigEff_qqH.getVal()
        print " @@@@@@@@ ZHeff ",rfvSigEff_ZH.getVal()
        print " @@@@@@@@ WHeff ",rfvSigEff_WH.getVal()
        print " @@@@@@@@ ttHeff ",rfvSigEff_ttH.getVal()
    
        CS_ggH = myCSW.HiggsCS(1,self.mH,self.sqrts)
        CS_VBF = myCSW.HiggsCS(2,self.mH,self.sqrts)
        CS_WH = myCSW.HiggsCS(3,self.mH,self.sqrts)
        CS_ZH = myCSW.HiggsCS(4,self.mH,self.sqrts)
        CS_ttH = myCSW.HiggsCS(5,self.mH,self.sqrts)
    
        BRH2e2mu = myCSW.HiggsBR(13,self.mH)
        BRH4mu = myCSW.HiggsBR(12,self.mH)
        BRH4e = myCSW.HiggsBR(12,self.mH)
        BR = 0.0
        if( self.channel == self.ID_4mu ): BR = BRH4mu
        if( self.channel == self.ID_4e ): BR = BRH4e
        if( self.channel == self.ID_2e2mu ): BR = BRH2e2mu

        #HZZ Branching ratio for ZH,WH,ttH samples
        BRZZ = myCSW.HiggsBR(11,self.mH)
    
        sigEfficiency_ggH = rfvSigEff_ggH.getVal()
        sigEfficiency_qqH = rfvSigEff_qqH.getVal()
        sigEfficiency_ZH = rfvSigEff_ZH.getVal()
        sigEfficiency_WH = rfvSigEff_WH.getVal()
        sigEfficiency_ttH = rfvSigEff_ttH.getVal()

        if(DEBUG):
            print "CS_ggH: ",CS_ggH,", CS_VBF: ",CS_VBF,", CS_WH: ",CS_WH,", CS_ZH: ",CS_ZH
            print ", CS_ttH: ",CS_ttH,", BRH2e2mu: ",BRH2e2mu,", BRH4mu: ",BRH4mu,", BRH4e: ",BRH4e,", BRZZ: ",BRZZ

        print "CS_ggH: ",CS_ggH,", CS_VBF: ",CS_VBF,", CS_WH: ",CS_WH,", CS_ZH: ",CS_ZH
        print ", CS_ttH: ",CS_ttH,", BRH2e2mu: ",BRH2e2mu,", BRH4mu: ",BRH4mu,", BRH4e: ",BRH4e,", BRZZ: ",BRZZ


        
        ## SIG YIELDS
        sigRate_ggH = CS_ggH*BR*sigEfficiency_ggH*1000.*self.lumi
        sigRate_VBF = CS_VBF*BR*sigEfficiency_qqH*1000.*self.lumi
        sigRate_WH = CS_WH*BRZZ*sigEfficiency_WH*1000.*self.lumi
        sigRate_ZH = CS_ZH*BRZZ*sigEfficiency_ZH*1000.*self.lumi
        sigRate_ttH = CS_ttH*BRZZ*sigEfficiency_ttH*1000.*self.lumi

        rfvSMD_Ratio_ggH = ROOT.RooFormulaVar()
        rfvSMD_Ratio_qqH = ROOT.RooFormulaVar()
        rfvSMD_Ratio_WH = ROOT.RooFormulaVar()
        rfvSMD_Ratio_ZH = ROOT.RooFormulaVar()
        rfvSMD_Ratio_ttH = ROOT.RooFormulaVar()

        tag_Ratio_Name = "hzz4l_SMD_ratio_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rfvSMD_Ratio_ggH = ROOT.RooRealVar(tag_Ratio_Name,tag_Ratio_Name,self.SMDsigCut)
        tag_Ratio_Name = "hzz4l_SMD_ratio_qqH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rfvSMD_Ratio_qqH = ROOT.RooRealVar(tag_Ratio_Name,tag_Ratio_Name,self.SMDsigCut)
        tag_Ratio_Name = "hzz4l_SMD_ratio_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rfvSMD_Ratio_WH = ROOT.RooRealVar(tag_Ratio_Name,tag_Ratio_Name,self.SMDsigCut)
        tag_Ratio_Name = "hzz4l_SMD_ratio_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rfvSMD_Ratio_ZH = ROOT.RooRealVar(tag_Ratio_Name,tag_Ratio_Name,self.SMDsigCut)
        tag_Ratio_Name = "hzz4l_SMD_ratio_ttH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rfvSMD_Ratio_ttH = ROOT.RooRealVar(tag_Ratio_Name,tag_Ratio_Name,self.SMDsigCut)

        print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvSMD_Ratio_ggH.getVal()
        print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvSMD_Ratio_qqH.getVal()
        print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvSMD_Ratio_WH.getVal()
        print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvSMD_Ratio_ZH.getVal()
        print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvSMD_Ratio_ttH.getVal()
        sigRate_ggH *= rfvSMD_Ratio_ggH.getVal()
        sigRate_VBF *= rfvSMD_Ratio_qqH.getVal()
        sigRate_WH *= rfvSMD_Ratio_WH.getVal()
        sigRate_ZH *= rfvSMD_Ratio_ZH.getVal()
        sigRate_ttH *= rfvSMD_Ratio_ttH.getVal()
       
        tmpNormSigNoConv = signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        tmpNormSigConv = sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        tmpNormSigHM   = sig_ggH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
      
        normalizationSignal = 0.0
        if self.isHighMass : normalizationSignal = tmpNormSigHM
        else : normalizationSignal = self.getVariable(tmpNormSigNoConv,tmpNormSigConv,self.bUseCBnoConvolution)
            
        print "#################### ",signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        print "#################### ",signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        print "#################### ",sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        print "#################### norm Signal",normalizationSignal
        
        sclFactorSig_ggH = sigRate_ggH/normalizationSignal
        sclFactorSig_VBF = sigRate_VBF/normalizationSignal
        sclFactorSig_WH = sigRate_WH/normalizationSignal
        sclFactorSig_ZH = sigRate_ZH/normalizationSignal
        sclFactorSig_ttH = sigRate_ttH/normalizationSignal

        integral_ggH = 0.0
        integral_VBF = 0.0
        integral_WH  = 0.0
        integral_ZH  = 0.0
        integral_ttH = 0.0

        if self.isHighMass : integral_ggH = sig_ggH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_ggH = self.getVariable(signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass : integral_VBF = sig_VBF_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_VBF = self.getVariable(signalCB_VBF.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_VBF.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass : integral_WH = sig_WH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_WH = self.getVariable(signalCB_WH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_WH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)
        

        if self.isHighMass : integral_ZH = sig_ZH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_ZH = self.getVariable(signalCB_ZH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ZH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass : integral_ttH = sig_ttH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_ttH = self.getVariable(signalCB_ttH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ttH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)
        
        sigRate_ggH_Shape = sclFactorSig_ggH*integral_ggH
        sigRate_VBF_Shape = sclFactorSig_VBF*integral_VBF
        sigRate_WH_Shape = sclFactorSig_WH*integral_WH
        sigRate_ZH_Shape = sclFactorSig_ZH*integral_ZH
        sigRate_ttH_Shape = sclFactorSig_ttH*integral_ttH
        

        normSigName = "cmshzz4l_normalizationSignal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rrvNormSig = ROOT.RooRealVar()

        

        if self.isHighMass :
            rrvNormSig = ROOT.RooRealVar(normSigName,normSigName, sig_ggH_HM.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal())
        else :
            rrvNormSig = ROOT.RooRealVar(normSigName,normSigName, self.getVariable(signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),self.bUseCBnoConvolution))
        rrvNormSig.setConstant(True)
        print "!!!%%%*** ",rrvNormSig.getVal()
        print "!!!%%%*** ",integral_ggH
        

        #rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),self.getVariable(signalCB_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)),ROOT.RooArgList(rfvSigEff_ggH, rhfXsBrFuncV_1))

        rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_ggH),ROOT.RooArgList(rfvSigEff_ggH, rhfXsBrFuncV_1))
        
        print "Compare integrals: integral_ggH=",integral_ggH,"  ; calculated=",self.getVariable(signalCB_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)
        
        rfvSigRate_VBF = ROOT.RooFormulaVar("qqH_norm","@0*@1*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_VBF),ROOT.RooArgList(rfvSigEff_qqH, rhfXsBrFuncV_2))
        
        
        rfvSigRate_WH = ROOT.RooFormulaVar("WH_norm","@0*@1*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_WH),ROOT.RooArgList(rfvSigEff_WH, rhfXsBrFuncV_3))
        
        
        rfvSigRate_ZH = ROOT.RooFormulaVar("ZH_norm","@0*@1*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_ZH),ROOT.RooArgList(rfvSigEff_ZH, rhfXsBrFuncV_4))
        
        
        rfvSigRate_ttH = ROOT.RooFormulaVar("ttH_norm","@0*@1*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_ttH),ROOT.RooArgList(rfvSigEff_ttH, rhfXsBrFuncV_5))
        

        print signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),"   ",sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal()
        print signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),"   ",sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal()
        
        print " @@@@@@@ norm sig = ",rrvNormSig.getVal()
        print " @@@@@@@ rfvSigRate_ggH = ",rfvSigRate_ggH.getVal()
        print " sigRate_ggH_Shape=",sigRate_ggH_Shape
        print " @@@@@@@ rfvSigRate_VBF = ",rfvSigRate_VBF.getVal()
        print " sigRate_VBF_Shape=",sigRate_VBF_Shape
        print " @@@@@@@ rfvSigRate_WH = ",rfvSigRate_WH.getVal()
        print " sigRate_WH_Shape=",sigRate_WH_Shape
        print " @@@@@@@ rfvSigRate_ZH = ",rfvSigRate_ZH.getVal()
        print " sigRate_ZH_Shape=",sigRate_ZH_Shape
        print " @@@@@@@ rfvSigRate_ttH = ",rfvSigRate_ttH.getVal()
        print " sigRate_ttH_Shape=",sigRate_ttH_Shape
        sigRate_Total_Shape_analytical = sigRate_ggH_Shape+sigRate_VBF_Shape+sigRate_WH_Shape+sigRate_ZH_Shape+sigRate_ttH_Shape
        print "Sum of analytical sigRate_XYZ_Shape=",sigRate_Total_Shape_analytical
        ## SET RATES TO 1 
        ## DC RATES WILL BE MULTIPLIED
        ## BY RATES IMPORTED TO WS
        #sigRate_ggH_Shape = 1
        #sigRate_VBF_Shape = 1
        #sigRate_WH_Shape = 1
        #sigRate_ZH_Shape = 1
        #sigRate_ttH_Shape = 1

        sigRate_ggH_input = theInputs['ggH_rate']
        if sigRate_ggH_input < 0:
            sigRate_ggH_input=sigRate_ggH_Shape
	else:
	    print "ggH Rate: ",sigRate_ggH_input
	    sigRate_ggH_Shape=sigRate_ggH_input

        eff_qqH_input = theInputs['qqH_eff']
        if eff_qqH_input >= 0:
	    print "qqH Custom Efficiency: ",eff_qqH_input
	    sigRate_VBF_Shape=eff_qqH_input*CS_VBF*BR*1000.*self.lumi
	    print "sigRate_VBF_Shape after custom eficiency: ",sigRate_VBF_Shape

        sigRate_Total_Shape = sigRate_ggH_Shape+sigRate_VBF_Shape+sigRate_WH_Shape+sigRate_ZH_Shape+sigRate_ttH_Shape
        sigRate_ggH_Shape=sigRate_Total_Shape
	print "Total yield: ",sigRate_ggH_Shape

             
        ## ----------------------- BACKGROUND RATES ----------------------- ##

        ## rates per lumi for scaling
        bkgRate_qqzz = theInputs['qqZZ_rate']/theInputs['qqZZ_lumi']
        bkgRate_ggzz = theInputs['ggZZ_rate']/theInputs['qqZZ_lumi']
        bkgRate_zjets = theInputs['zjets_rate']/theInputs['zjets_lumi']
        
        ## Get Normalizations
        normalizationBackground_qqzz = bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_ggzz = bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_zjets = bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()

        print "channel: "+self.appendName
        print "fullrange zjets: ",normalizationBackground_zjets
        
        sclFactorBkg_qqzz = self.lumi*bkgRate_qqzz/normalizationBackground_qqzz
        sclFactorBkg_ggzz = self.lumi*bkgRate_ggzz/normalizationBackground_ggzz
        sclFactorBkg_zjets = self.lumi*bkgRate_zjets/normalizationBackground_zjets
               
        bkgRate_qqzz_Shape = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_ggzz_Shape = sclFactorBkg_ggzz * bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_zjets_Shape = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()

        rfvSMD_Ratio_qqZZ = ROOT.RooFormulaVar()
        rfvSMD_Ratio_ggZZ = ROOT.RooFormulaVar()
        rfvSMD_Ratio_Zjets = ROOT.RooFormulaVar()

        tag_Ratio_Name = "hzz4l_SMD_ratio_qqZZ_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rfvSMD_Ratio_qqZZ = ROOT.RooRealVar(tag_Ratio_Name,tag_Ratio_Name,self.SMDbkgCut)
        tag_Ratio_Name = "hzz4l_SMD_ratio_ggZZ_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rfvSMD_Ratio_ggZZ = ROOT.RooRealVar(tag_Ratio_Name,tag_Ratio_Name,self.SMDbkgCut)
        tag_Ratio_Name = "hzz4l_SMD_ratio_Zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rfvSMD_Ratio_Zjets = ROOT.RooRealVar(tag_Ratio_Name,tag_Ratio_Name,self.SMDbkgCut)

        print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvSMD_Ratio_qqZZ.getVal()
        print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvSMD_Ratio_ggZZ.getVal()
        print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvSMD_Ratio_Zjets.getVal()
        
        bkgRate_qqzz_Shape *= rfvSMD_Ratio_qqZZ.getVal()
        bkgRate_ggzz_Shape *= rfvSMD_Ratio_ggZZ.getVal()
        bkgRate_zjets_Shape *= rfvSMD_Ratio_Zjets.getVal()
        
        if(DEBUG):
            print "Shape signal rate: ",sigRate_ggH_Shape,", background rate: ",bkgRate_qqzz_Shape,", ",bkgRate_zjets_Shape," in ",low_M," - ",high_M
            CMS_zz4l_mass.setRange("lowmassregion",100.,160.)
            bkgRate_qqzz_lowmassregion = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_ggzz_lowmassregion = sclFactorBkg_ggzz * bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_zjets_lowmassregion = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            lowmassyield = bkgRate_qqzz_lowmassregion + bkgRate_ggzz_lowmassregion + bkgRate_zjets_lowmassregion
            print "low mass yield: ",lowmassyield
        
        ## --------------------------- DATASET --------------------------- ##

        dataFileDir = "CMSdata"
        dataTreeName = "data_obs"
        if (self.dataAppendDir == ''):
            dataFileName = "{0}/hzz{1}_{2}.root".format(dataFileDir,self.appendName,self.lumi)
        else:
            dataFileName = "{0}_{1}/hzz{2}_{3}.root".format(dataFileDir,self.dataAppendDir,self.appendName,self.lumi)
        if (DEBUG): print dataFileName," ",dataTreeName 
        data_obs_file = ROOT.TFile(dataFileName)

        print data_obs_file.Get(dataTreeName)
        
        if not (data_obs_file.Get(dataTreeName)):
            print "File, \"",dataFileName,"\", or tree, \"",dataTreeName,"\", not found" 
            print "Exiting..."
            sys.exit()
        
        data_obs_tree = data_obs_file.Get(dataTreeName)
        data_obs = ROOT.RooDataSet()
        datasetName = "data_obs_{0}".format(self.appendName)
        

        data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,D1,D2,D3))

            
        ## --------------------------- WORKSPACE -------------------------- ##

        endsInP5 = False
        tmpMH = self.mH
        if ( math.fabs(math.floor(tmpMH)-self.mH) > 0.001): endsInP5 = True
        if (DEBUG): print "ENDS IN P5  ",endsInP5

        name_Shape = ""
        name_ShapeWS = ""
        name_ShapeWS2 = ""

        if (endsInP5): name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        
        if (endsInP5): name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else: name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)

        if(DEBUG): print name_Shape,"  ",name_ShapeWS2
        
        w = ROOT.RooWorkspace("w","w")
        
        w.importClassCode(RooqqZZPdf_v2.Class(),True)
        w.importClassCode(RooggZZPdf_v2.Class(),True)
        w.importClassCode(HZZ4L_RooCTauPdf.Class(),True)
        w.importClassCode(RooFormulaVar.Class(),True)
        if self.isHighMass :
            w.importClassCode(RooRelBWHighMass.Class(),True)
            
                
                
        getattr(w,'import')(data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?
#        getattr(w,'import')(r_fai_norm) ### Should this be renamed?


        ggHpdf.SetNameTitle("ggH","ggH")
        getattr(w,'import')(ggHpdf, ROOT.RooFit.RecycleConflictNodes())
        ggHpdf_systUp.SetNameTitle("ggH_Res{0}Up".format(self.appendName),"ggH_Res{0}Up".format(self.appendName))
        getattr(w,'import')(ggHpdf_systUp, ROOT.RooFit.RecycleConflictNodes())
        ggHpdf_systDown.SetNameTitle("ggH_Res{0}Down".format(self.appendName),"ggH_Res{0}Down".format(self.appendName))
        getattr(w,'import')(ggHpdf_systDown, ROOT.RooFit.RecycleConflictNodes())

        qqZZpdf.SetNameTitle("bkg_qqzz","bkg_qqzz")
        ggZZpdf.SetNameTitle("bkg_ggzz","bkg_ggzz")
        Zjetspdf.SetNameTitle("bkg_zjets","bkg_zjets")
        getattr(w,'import')(qqZZpdf, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(ggZZpdf, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(Zjetspdf, ROOT.RooFit.RecycleConflictNodes())
  
        w.writeToFile(name_ShapeWS)
        
        ## --------------------------- DATACARDS -------------------------- ##

        systematics.setSystematics(bkgRate_qqzz_Shape, bkgRate_ggzz_Shape, bkgRate_zjets_Shape)

        ## If the channel is not declared in inputs, set rate = 0
        if not self.ggH_chan:  sigRate_ggH_Shape = 0
        if not self.qqH_chan:  sigRate_VBF_Shape = 0
        if not self.WH_chan:   sigRate_WH_Shape = 0
        if not self.ZH_chan:   sigRate_ZH_Shape = 0
        if not self.ttH_chan:  sigRate_ttH_Shape = 0

        if not self.qqZZ_chan:  bkgRate_qqzz_Shape = 0
        if not self.ggZZ_chan:  bkgRate_ggzz_Shape = 0
        if not self.zjets_chan: bkgRate_zjets_Shape = 0

        rates = {}
        rates['ggH'] = sigRate_ggH_Shape
        rates['qqH'] = sigRate_VBF_Shape
        rates['WH']  = sigRate_WH_Shape
        rates['ZH']  = sigRate_ZH_Shape
        rates['ttH'] = sigRate_ttH_Shape

        rates['qqZZ']  = bkgRate_qqzz_Shape
        rates['ggZZ']  = bkgRate_ggzz_Shape
        rates['zjets'] = bkgRate_zjets_Shape
        rates['ttbar'] = 0
        rates['zbb']   = 0
        

        ## Write Datacards
        fo = open( name_Shape, "wb")
        self.WriteDatacard(fo,theInputs, name_ShapeWS2, rates, data_obs.numEntries())
        
        systematics.WriteSystematics(fo, theInputs)
        systematics.WriteShapeSystematics(fo,theInputs)

        
        fo.close()


    def WriteDatacard(self,file,theInputs,nameWS,theRates,obsEvents):

        numberSig = self.numberOfSigChan(theInputs)
        numberBg  = self.numberOfBgChan(theInputs)
        
        file.write("imax 1\n")
        file.write("jmax {0}\n".format(numberSig+numberBg-1))
        file.write("kmax *\n")
        
        file.write("------------\n")
        file.write("shapes * * {0} w:$PROCESS w:$PROCESS_$SYSTEMATIC\n".format(nameWS))
        file.write("------------\n")
        

        file.write("bin a{0} \n".format(self.channel))
        file.write("observation {0} \n".format(obsEvents))
        
        file.write("------------\n")
        file.write("## mass window [{0},{1}] \n".format(self.low_M,self.high_M))
        file.write("bin ")        

        channelList=['ggH','qqH','WH','ZH','ttH','qqZZ','ggZZ','zjets','ttbar','zbb']

        channelName=['ggH','qqH','WH','ZH','ttH','bkg_qqzz','bkg_ggzz','bkg_zjets','bkg_ttbar','bkg_zbb']
         
        for chan in channelList:
            if theInputs[chan]:
                file.write("a{0} ".format(self.channel))
        file.write("\n")
                                        
        file.write("process ")

        i=0

        for chan in channelList:
            #print 'checking if ',chan,' is in the list of to-do'
            #print "{0} ".format(channelName[i])
            if theInputs[chan]:
                file.write("{0} ".format(channelName[i]))
                #print 'writing in card index=',i,'  chan=',chan
                #print "{0} ".format(channelName[i])
            i+=1

        
        file.write("\n")
            
        processLine = "process "

        for x in range(-numberSig+1,1):
            processLine += "{0} ".format(x)

        for y in range(1,numberBg+1):
            processLine += "{0} ".format(y)

        file.write(processLine)
        file.write("\n")
            
        file.write("rate ")
        for chan in channelList:
            if theInputs[chan]:
                file.write("{0:.4f} ".format(theRates[chan]))
        file.write("\n")
        file.write("------------\n")


        
    def numberOfSigChan(self,inputs):

        counter=0

        if inputs['ggH']: counter+=1
        if inputs['qqH']: counter+=1
        if inputs['WH']:  counter+=1
        if inputs['ZH']:  counter+=1
        if inputs['ttH']: counter+=1
        
        return counter

    def numberOfBgChan(self,inputs):

        counter=0

        if inputs['qqZZ']:  counter+=1
        if inputs['ggZZ']:  counter+=1
        if inputs['zjets']: counter+=1
        if inputs['ttbar']: counter+=1
        if inputs['zbb']:   counter+=1
        
        return counter

