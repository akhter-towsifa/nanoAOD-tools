#!/usr/bin/env python
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *
from PhysicsTools.NanoAODTools.postprocessing.trigger.triggerFilter import triggerFilter

from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import countHistogramsProducer
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import runsAndLumis

import argparse

parser = argparse.ArgumentParser()
#parser.add_argument( '--isMC', default=False, type=bool)
parser.add_argument( '--isMC', default="False")
parser.add_argument('--crab', default=False, action='store_true')
parser.add_argument('--outfile', default='.')
parser.add_argument('--maxEntries', default=None)
parser.add_argument("infile")
parser.add_argument("dataYear")
parser.add_argument("runPeriod")
parser.add_argument("triggers")
parser.add_argument("btagWP", type=float)
parser.add_argument("btag_type")
parser.add_argument("selector")
parser.add_argument("keep_and_drop")
args = parser.parse_args()

isMC = args.isMC=="True"
print("isMC", args.isMC, isMC)
infile = args.infile
outfile = args.outfile
if args.maxEntries is not None: maxEntries = int(args.maxEntries)
else: maxEntries = args.maxEntries
dataYear = args.dataYear
runPeriod = args.runPeriod
print(args.triggers)
triggers = eval(args.triggers)
btagWP = float(args.btagWP)
btag_type = args.btag_type
selector = args.selector
keep_and_drop = args.keep_and_drop
crab = args.crab

print('Arguments:')
print('\t {}'.format(args))

#isMC=False
#infile=''
#outfile='.'
#maxEntries = None
#dataYear =  '2017'
#runPeriod = ''
#triggers = ['HLT_DoubleEle33_CaloIdL_MW']
#btagWP = 0.4941
#btag_type = 'deepcsv'
#selector = 'minseok'
#crab = True

#select right module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.heepV72018PromptProducer import heepV72018PromptProducer
if isMC:
    from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSF
    from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import createJMECorrector
    if dataYear=='2016':
        from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import lepSF2016 as lepSF
        from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2016 as puWeight
    if dataYear=='2017':
        from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import lepSF2017 as lepSF
        from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2017 as puWeight
    if dataYear=='2018':
        from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import lepSF2018 as lepSF
        from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2018 as puWeight
if dataYear=='2016':
    from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import muonScaleRes2016 as muonScaleRes
if dataYear=='2017':
    from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import muonScaleRes2017 as muonScaleRes
if dataYear=='2018':
    from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import muonScaleRes2018 as muonScaleRes

#different preselection producerts
if selector=='inclusive':
    from PhysicsTools.NanoAODTools.postprocessing.bff.bffInclusive_preselectionModule import bffInclusivePreselProducer as preselectorProducer
if selector=='bff':
    from PhysicsTools.NanoAODTools.postprocessing.bff.bffPreselModule import bffPreselProducer as preselectorProducer
if selector=='bff_eff':
    from PhysicsTools.NanoAODTools.postprocessing.bff.bffBtagEff import bffBtagEffProducer as preselectorProducer    
if selector=='minseok':
    from PhysicsTools.NanoAODTools.postprocessing.bff.bffPreselModule_minseok import bffPreselProducer as preselectorProducer

#prepare file if glob
if not crab:
    print("using input file", infile)
    import glob
    infile = glob.glob(infile)
else:
    from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles
    print("get input file (crab)")
    infile = inputFiles()
print("infile", infile)
#2017 EE fix
conditions_dict = {
    '2016':{'metBranchName':'MET', 'heepBranchName': 'cutBased_HEEP'},
    '2017':{'metBranchName':'METFixEE2017', 'heepBranchName': 'cutBased_HEEP'},
    '2018':{'metBranchName':'MET', 'heepBranchName': 'cutBased_HEEPV7p0_2018Prompt'}
}

#set up process for mc and data
if isMC:
    jmeCorrections = createJMECorrector(
            isMC=isMC,
            dataYear=dataYear,
            runPeriod=runPeriod,
            jesUncert="Total",
            metBranchName=conditions_dict[dataYear]['metBranchName'],
            applySmearing=True,
            jetType="AK4PFchs",
            noGroom=False
        )
    modules=[
            countHistogramsProducer(),
            #triggerFilter(triggers),
            btagSF(dataYear, algo=btag_type),
            jmeCorrections(),
            puWeight(),
            muonScaleRes(),
            lepSF(),
            heepV72018PromptProducer(),
            preselectorProducer(btagWP, triggers, isMC=isMC, btag_type=btag_type, **conditions_dict[dataYear])
        ]
    p = PostProcessor(outfile,
            infile,
            modules=modules,
            provenance=True,
            fwkJobReport=True,
            outputbranchsel=keep_and_drop,
            maxEntries=maxEntries
        )
else:
    modules=[
            countHistogramsProducer(),
            #triggerFilter(triggers),
            muonScaleRes(),
            heepV72018PromptProducer(),
            preselectorProducer(btagWP, triggers, isMC=isMC, btag_type=btag_type, **conditions_dict[dataYear])
        ]

    p = PostProcessor(outfile,
            infile,
            modules=modules,
            provenance=True,
            fwkJobReport=True,
            outputbranchsel=keep_and_drop,
            jsonInput=runsAndLumis(),
            maxEntries=maxEntries
        )

p.run()

print("DONE")
