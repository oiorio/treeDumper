#!/usr/bin/env python
"""
This is a small script that submits a config over many datasets
"""
import os, commands
from optparse import OptionParser

def getOptions() :
    """
    Parse and return the arguments provided by the user.
    """
    usage = ('usage: python submit_all.py -c CONFIG -d DIR ')

    parser = OptionParser(usage=usage)    
    parser.add_option("-c", "--config", dest="config",
        help=("The crab script you want to submit "),
        metavar="CONFIG")
    parser.add_option("-d", "--dir", dest="dir",
        help=("The crab directory you want to use "),
        metavar="DIR")
    parser.add_option("-f", "--datasets", dest="datasets",
        help=("File listing datasets to run over"),
        metavar="FILE")
    parser.add_option( "--jecVersion", dest="jecVersion",
        default="",
        help=("wildcard for the jecs to use"),
        metavar="JECV")
    parser.add_option( "--isData", dest="isData",
                       default=False, action="store_true",
                       help=("Is it data?") )
                      
    parser.add_option("-n", "--dryrun", dest="dryrun",
        help=("Dry run"),
        metavar="DRYRUN")

    (options, args) = parser.parse_args()


    if options.config == None or options.dir == None:
        parser.error(usage)
    
    return options
    

def main():

    options = getOptions()
    print "IsData? ", options.isData
    datasetsFile = open( options.datasets )
    jobsLines = datasetsFile.readlines()
    jobs = []

    from CRABClient.UserUtilities import config
    config = config()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    config.section_("General")
    config.General.workArea = options.dir
   
    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = options.config
    config.JobType.allowUndistributedCMSSW = True
#    config.JobType.pyCfgParams = ['isData=' + str(options.isData), 'changeJECs=True', 'channel=wzjets']
    config.JobType.pyCfgParams = ['isData=' + str(options.isData), 'changeJECs=True','channel=wprime']

    #config.JobType.inputFiles = ['Fall15_25nsV2_DATA.db', 'Fall15_25nsV2_MC.db']
    #config.JobType.inputFiles = ["Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt", "Fall15_25nsV2_MC_L1RC_AK4PFchs.txt","Fall15_25nsV2_MC_L2Relative_AK4PFchs.txt", "Fall15_25nsV2_MC_L3Absolute_AK4PFchs.txt","Fall15_25nsV2_MC_L2L3Residual_AK4PFchs.txt", "Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt","Fall15_25nsV2_DATA_L1RC_AK4PFchs.txt","Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt","Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt",  "Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt"]
    
    #Input part:

    
    wildcard= options.jecVersion+"*txt"
    lscmd = "ls JECs/"+wildcard 
    outs=commands.getstatusoutput(lscmd)
    print "status: ", outs[0]," result: ",outs[1]
    inputs=[]
    for l in outs[1].split("\n"):
        inputs.append(l)

    inputs.append('Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt')
    inputs.append('DeepCSV_94XSF_V3_B_F.csv')
    inputs.append('CSVv2_Moriond17_B_H.csv')
    inputs.append('CSVv2_ichep.csv')
    inputs.append('cMVAv2_Moriond17_B_H.csv')
    inputs.append('cMVAv2_ichep.csv')
    inputs.append('btagging_cmva.root')

    
    config.JobType.inputFiles = inputs
    print "inputs are", config.JobType.inputFiles
    oldinputs = [
        'Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK8PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_L1RC_AK4PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_L1RC_AK8PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_L2Residual_AK4PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK8PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_UncertaintySources_AK4PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_UncertaintySources_AK8PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFchs.txt',
        'Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK8PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_L1FastJet_AK8PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_L1RC_AK4PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_L1RC_AK8PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK8PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_L2Relative_AK8PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_L2Residual_AK4PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_L3Absolute_AK8PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_UncertaintySources_AK4PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_UncertaintySources_AK8PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_Uncertainty_AK4PFchs.txt',
        'Summer16_23Sep2016EFV4_DATA_Uncertainty_AK8PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_L1FastJet_AK8PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_L1RC_AK4PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_L1RC_AK8PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_L2L3Residual_AK8PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_L2Relative_AK8PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_L2Residual_AK4PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_L3Absolute_AK8PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_UncertaintySources_AK4PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_UncertaintySources_AK8PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_Uncertainty_AK4PFchs.txt',
        'Summer16_23Sep2016GV4_DATA_Uncertainty_AK8PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_L1FastJet_AK8PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_L1RC_AK4PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_L1RC_AK8PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_L2L3Residual_AK8PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_L2Relative_AK8PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_L2Residual_AK4PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_L3Absolute_AK8PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_UncertaintySources_AK4PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_UncertaintySources_AK8PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PFchs.txt',
        'Summer16_23Sep2016HV4_DATA_Uncertainty_AK8PFchs.txt',
        'Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt',
        'Summer16_23Sep2016V4_MC_L1FastJet_AK8PFchs.txt',
        'Summer16_23Sep2016V4_MC_L1RC_AK4PFchs.txt',
        'Summer16_23Sep2016V4_MC_L1RC_AK8PFchs.txt',
        'Summer16_23Sep2016V4_MC_L2L3Residual_AK4PFchs.txt',
        'Summer16_23Sep2016V4_MC_L2L3Residual_AK8PFchs.txt',
        'Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt',
        'Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt',
        'Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt',
        'Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt',
        'Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt',
        'Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt',
        ]
    config.section_("Data")
    config.Data.ignoreLocality = True
    config.Data.inputDataset = None
    config.Data.inputDBS = 'phys03'
    config.Data.splitting = 'FileBased' 
#    config.Data.splitting = 'Automatic' 
    config.Data.unitsPerJob = 5
    config.Data.publication = False    
#    config.Data.outLFNDirBase = '/store/user/cgiuglia/trees/May12/'
#    config.Data.outLFNDirBase = '/store/user/oiorio/ttDM/trees/2018/May28/'
#    config.Data.outLFNDirBase = '/store/user/oiorio/Wprime/2018/June/June13/'
    config.Data.outLFNDirBase = '/store/user/oiorio/Tprime/trees/2018/Oct3/'

    config.section_("Site")
#    config.Site.storageSite = 'T2_IT_Pisa'
    config.Site.storageSite = 'T2_CH_CSCS'
    config.Site.whitelist = ['T2_CH_CERN','T2_IT_*','T2_DE_*','T2_CH_*']
    #config.Site.blacklist = ['T2_DE_RWTH']

    print 'Using config ' + options.config
    print 'Writing to directory ' + options.dir

    
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print 'Cannot execute command'
            print hte.headers

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    datasetsFile = open( options.datasets )
    jobsLines = datasetsFile.readlines()
    jobs = []


    for ijob in jobsLines :
        s = ijob.rstrip()
        jobs.append( s )
        print '  --> added ' + s

        if (ijob.rstrip()).find('2016B') :
            run = 'B'
        if (ijob.rstrip()).find('2016C') :                                                                                                            
            run = 'C'                                                                                                                               
        if (ijob.rstrip()).find('2016D') :                                                                                                       
            run = 'D'                                                                                                                                 
        if (ijob.rstrip()).find('2016E') :                                                                                                                
            run= 'E'                                                                                                                                      
        if (ijob.rstrip()).find('2016F') :                                                                                                                  
            run= 'F'                                                                                                                                         
        if (ijob.rstrip()).find('2016G') :                                                                                                                
            run= 'G'

    for ijob, job in enumerate(jobs) :

        eraLabel=''
        print "**** Check job rstrip: ", job.rstrip()        
        if ('2016B' in job.rstrip()) :
            eraLabel = 'BCD'
        elif ('2016C' in job.rstrip()) :   
            eraLabel = 'BCD'
        elif ('2016D' in job.rstrip()) :                                                                                                       
            eraLabel = 'BCD'
        elif ('2016E' in job.rstrip()) :
            eraLabel = 'EF'
        elif ('2016F' in job.rstrip()) :
            eraLabel = 'EF'
        elif ('2016G' in job.rstrip()) :                                                                                                    
            eraLabel = 'G'
        elif ('2016H' in job.rstrip()) :                                                                                                    
            eraLabel = 'H'
        else:
            eraLabel = 'BCD'
        #        print "-------> ERA: ", eraLabel

        if(options.isData):config.JobType.pyCfgParams = ['isData=' + str(options.isData), 'changeJECs=True', "EraLabel="+eraLabel]
        print "====> Config: ", config.JobType.pyCfgParams

        ptbin = job.split('/')[1]
        cond = job.split('/')[2]
        runs = cond.split('-')[1:-1]
        run = '-'.join(runs)
        print "===> PtBin: ", ptbin
        if(options.isData):
            ptbin = ptbin.split('-')[:1]
            ptbin = '-'.join(ptbin)

        print 'ptbin: ',ptbin
        config.General.requestName = ptbin +'_'+options.dir 
        print "1st Request name: ", config.General.requestName        
        if(options.isData): 
            print "I am in the wrong hole -> isData is: ", options.isData
            config.General.requestName = ptbin + run + "_" + opt.dir  
        print "2nd Request name: ", config.General.requestName        
        config.Data.inputDataset = job
        config.Data.outputDatasetTag = ptbin+"_"+options.dir
        if(options.isData): config.Data.outputDatasetTag = ptbin + run + "_"+options.dir
        print 'Submitting ' + config.General.requestName + ', dataset = ' + job
        print 'Configuration :'
        print "Request name: ", config.General.requestName        
#print config
        try :
            from multiprocessing import Process
            print options.dryrun
            if options.dryrun == None:
                p = Process(target=submit, args=(config,))
                p.start()
                p.join()
            #submit(config)
        except :
            print 'Not submitted.'
        

if __name__ == '__main__':
    main()            
