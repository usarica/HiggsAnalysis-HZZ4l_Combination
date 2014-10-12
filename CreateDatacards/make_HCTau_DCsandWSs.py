#!/usr/bin/python
#-----------------------------------------------
# Latest update: 2012.08.30
# by Matt Snowball
#-----------------------------------------------
import sys, os, pwd, commands
import optparse, shlex, re
import math
from ROOT import *
import ROOT
from array import array
from datacardClass_HCTau_1D import *
from inputReader import *

#define function for parsing options
def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    
    parser.add_option('-i', '--input', dest='inputDir', type='string', default="",    help='inputs directory')
    parser.add_option('-a', '--append', dest='appendName', type='string', default="",    help='append name for cards dir')
    parser.add_option('-b', action='store_true', dest='noX', default=True ,help='no X11 windows')
    parser.add_option('-t', '--templateDir', type='string', dest='templateDir', default="templates2D" ,help='directory with 2D template histos')
#    parser.add_option('-m', '--model', type='string', dest='model', default="1D" ,help='model: 1D, phase, 2D ')
    parser.add_option('-r', '--datadir', type='string', dest='dataDirAppend', default="" ,help='dataDirAppend: Reference CMSdata folder per measurement')

    
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()


    if (opt.appendName == ''):
        print 'Please pass an append name for the cards directory! Exiting...'
        sys.exit()

    if (opt.inputDir == ''):
        print 'Please pass an input directory! Exiting...'
        sys.exit()



    
# define make directory function
def makeDirectory(subDirName):
    if (not os.path.exists(subDirName)):
        cmd = 'mkdir -p '+subDirName
        status, output = commands.getstatusoutput(cmd)
        if status !=0:
            print 'Error in creating submission dir '+subDirName+'. Exiting...'
            sys.exit()
    else:
        print 'Directory '+subDirName+' already exists. Exiting...'
    #    sys.exit()


#define function for processing of os command
def processCmd(cmd):
#    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        sys.exit()



def creationLoop(directory):
    global opt, args

    startMass=[ 125.6]
    stepSizes=[ 0.1]
    endVal=[     1]

    myClass = datacardClass_HCTau_1D()

    myClass.loadIncludes()
    
    myReader4e = inputReader(opt.inputDir+"/inputs_4e.txt")
    myReader4e.readInputs()
    theInputs4e = myReader4e.getInputs()
    
    myReader4mu = inputReader(opt.inputDir+"/inputs_4mu.txt")
    myReader4mu.readInputs()
    theInputs4mu = myReader4mu.getInputs()
    
    myReader2e2mu = inputReader(opt.inputDir+"/inputs_2e2mu.txt")
    myReader2e2mu.readInputs()
    theInputs2e2mu = myReader2e2mu.getInputs()

    
    a=0
    while (a < len(startMass) ):
	
	c = 0
        while (c < endVal[a] ): 
            
            mStart = startMass[a]
            step = stepSizes[a]
            mh = mStart + ( step * c ) 
            mhs = str(mh).replace('.0','')

            print mh

            makeDirectory(directory+'/HCG/'+mhs)
            myClass.makeCardsWorkspaces(mh,directory,theInputs4e,opt.templateDir,opt.dataDirAppend)
            myClass.makeCardsWorkspaces(mh,directory,theInputs4mu,opt.templateDir,opt.dataDirAppend)
            myClass.makeCardsWorkspaces(mh,directory,theInputs2e2mu,opt.templateDir,opt.dataDirAppend)
                
            c += 1
            

	a += 1






# the main procedure
def make_HCTau_DCsandWSs():
    
    # parse the arguments and options
    global opt, args
    parseOptions()

    if (opt.appendName != ''):
        dirName = 'cards_'+opt.appendName
    

    subdir = ['HCG','HCG_XSxBR','figs']

    for d in subdir:
        makeDirectory(dirName+'/'+d)
        

    creationLoop(dirName)
    

    sys.exit()



# run the create_RM_cfg() as main()
if __name__ == "__main__":
    make_HCTau_DCsandWSs()


