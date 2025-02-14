sample_list = [
{"era": 2016, "ismc":1, "type": "TT", "name": "TT", "das": "/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "ST", "name": "ST_anti", "das": "/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "ST", "name": "ST", "das": "/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DB", "name": "WW", "das": "/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DB", "name": "WZ", "das": "/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DB", "name": "ZZ", "das": "/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DY", "name": "ZToEE_120_200", "das": "/ZToEE_NNPDF30_13TeV-powheg_M_120_200/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DY", "name": "ZToEE_200_400", "das": "/ZToEE_NNPDF30_13TeV-powheg_M_200_400/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DY", "name": "ZToEE_400_800", "das": "/ZToEE_NNPDF30_13TeV-powheg_M_400_800/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DY", "name": "ZToEE_50_120", "das": "/ZToEE_NNPDF30_13TeV-powheg_M_50_120/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DY", "name": "ZToEE_800_1400", "das": "/ZToEE_NNPDF30_13TeV-powheg_M_800_1400/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DY", "name": "ZToMuMu_120_200", "das": "/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DY", "name": "ZToMuMu_200_400", "das": "/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DY", "name": "ZToMuMu_400_800", "das": "/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DY", "name": "ZToMuMu_50_120", "das": "/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "DY", "name": "ZToMuMu_800_1400", "das": "/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_200_dbs0p04", "das": "/BFFZprimeToMuMu_M_200_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_200_dbs0p5", "das": "/BFFZprimeToMuMu_M_200_dbs0p5_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_200_dbs1p0", "das": "/BFFZprimeToMuMu_M_200_dbs1p0_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_250_dbs0p04", "das": "/BFFZprimeToMuMu_M_250_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_300_dbs0p04", "das": "/BFFZprimeToMuMu_M_300_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_350_dbs0p04", "das": "/BFFZprimeToMuMu_M_350_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_350_dbs0p5", "das": "/BFFZprimeToMuMu_M_350_dbs0p5_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_350_dbs1p0", "das": "/BFFZprimeToMuMu_M_350_dbs1p0_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_400_dbs0p04", "das": "/BFFZprimeToMuMu_M_400_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_500_dbs0p04", "das": "/BFFZprimeToMuMu_M_500_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_500_dbs0p5", "das": "/BFFZprimeToMuMu_M_500_dbs0p5_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
{"era": 2016, "ismc":1, "type": "BFF", "name": "BFF_500_dbs1p0", "das": "/BFFZprimeToMuMu_M_500_dbs1p0_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM"},
    
    
{"das": "/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "ST_tW_antitop_5f_inclusiveDecays_13TeV","type": "ST"},
{"das": "/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "ST_tW_top_5f_inclusiveDecays_13TeV","type": "ST"},
{"das": "/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "ST_t_anti","type": "ST"},
{"das": "/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v2/NANOAODSIM", "era": 2016, "ismc": 1, "name": "ST_t","type": "ST"},
{"das": "/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "ST_s","type": "ST"},
{"das": "/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "WWZ_TuneCUETP8M1_13TeV","type": "TB"},
{"das": "/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "ZZZ_TuneCUETP8M1_13TeV","type": "TB"},
{"das": "/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "WZZ_TuneCUETP8M1_13TeV","type": "TB"},
{"das": "/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "WWW_4F_TuneCUETP8M1_13TeV","type": "TB"},
{"das": "/WJetsToQQ_HT400to600_qc19_3j_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "WJetsToQQ_HT400to600_qc19_3j_TuneCUETP8M1_13TeV","type": "WJets"},
{"das": "/WJetsToQQ_HT-800toInf_qc19_3j_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "WJetsToQQ_HT_800_inf","type": "WJets"},
{"das": "/WJetsToQQ_HT600to800_qc19_3j_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "WJetsToQQ_HT600to800_qc19_3j_TuneCUETP8M1_13TeV","type": "WJets"},
{"das": "/WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "WJetsToQQ_HT","type": "WJets"},
{"das": "/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "WJetsToLNu_TuneCUETP8M1_13TeV","type": "WJets"},
{"das": "/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "TTWJetsToQQ_TuneCUETP8M1_13TeV","type": "TTV"},
{"das": "/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "TTWJetsToLNu_TuneCUETP8M1_13TeV","type": "TTV"},
{"das": "/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "TTZToQQ_TuneCUETP8M1_13TeV","type": "TTV"},
{"das": "/TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "TTZToLL_M","type": "TTV"},
{"das": "/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM", "era": 2016, "ismc": 1, "name": "TTZToLLNuNu_M","type": "TTV"},
    
{"era": 2016, "ismc":1, "type": "Higgs", "name": "ggH", "das": "/GluGluHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM" },
{"era": 2016, "ismc":1, "type": "Higgs", "name": "VBF", "das": "/VBFHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnlo_pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM" },
{"era": 2016, "ismc":1, "type": "Higgs", "name": "Wp", "das": "/WplusH_HToMuMu_WToAll_M125_TuneCP5_PSweights_13TeV_powheg_pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM" },
{"era": 2016, "ismc":1, "type": "Higgs", "name": "Wm", "das": "/WminusH_HToMuMu_WToAll_M125_TuneCP5_PSweights_13TeV_powheg_pyt-hia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM" },
{"era": 2016, "ismc":1, "type": "Higgs", "name": "ZH", "das": "/ZH_HToMuMu_ZToAll_M125_TuneCP5_PSweights_13TeV_powheg_pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM" },
{"era": 2016, "ismc":1, "type": "Higgs", "name": "ttH", "das": "/ttHToMuMu_M125_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM" },
    
    
{"era": 2016, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_50_nlo", "das": "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/NANOAODSIM"},
{"era": 2016, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_50_mad_v2", "das": "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM"},
{"era": 2016, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_50_mad", "das": "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/NANOAODSIM"},
{"era": 2016, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_100_200", "das": "/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM"},
{"era": 2016, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_200_400", "das": "/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/NANOAODSIM"},
{"era": 2016, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_400_500", "das": "/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM"},
{"era": 2016, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_500_700", "das": "/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM"},
{"era": 2016, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_700_800", "das": "/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM"},
{"era": 2016, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_800_1000", "das": "/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM"},

{"era": 2017, "ismc":1, "type": "TT", "name": "TT", "das": "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_new_pmx_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "ST", "name": "ST_anti", "das": "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "ST", "name": "ST", "das": "/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DB", "name": "WW", "das": "/WW_TuneCP5_13TeV-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DB", "name": "WZ", "das": "/WZ_TuneCP5_13TeV-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DB", "name": "ZZ", "das": "/ZZ_TuneCP5_13TeV-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_200_dbs0p04", "das": "/BFFZprimeToMuMu_M_200/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_200_dbs0p5", "das": "/BFFZprimeToMuMu_M_200_dbs0p5/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_200_dbs1p0", "das": "/BFFZprimeToMuMu_M_200_dbs1p0/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_250_dbs0p04", "das": "/BFFZprimeToMuMu_M_250/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_300_dbs0p04", "das": "/BFFZprimeToMuMu_M_300/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_350_dbs0p04", "das": "/BFFZprimeToMuMu_M_350/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_350_dbs0p5", "das": "/BFFZprimeToMuMu_M_350_dbs0p5/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_350_dbs1p0", "das": "/BFFZprimeToMuMu_M_350_dbs1p0/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_400_dbs0p04", "das": "/BFFZprimeToMuMu_M_400/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_450_dbs0p04", "das": "/BFFZprimeToMuMu_M_450/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_500_dbs0p04", "das": "/BFFZprimeToMuMu_M_500/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_500_dbs0p5", "das": "/BFFZprimeToMuMu_M_500_dbs0p5/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "BFF", "name": "BFF_500_dbs1p0", "das": "/BFFZprimeToMuMu_M_500_dbs1p0/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DY", "name": "ZToEE_120_200", "das": "/ZToEE_NNPDF31_13TeV-powheg_M_120_200/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DY", "name": "ZToEE_200_400", "das": "/ZToEE_NNPDF31_13TeV-powheg_M_200_400/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DY", "name": "ZToEE_400_800", "das": "/ZToEE_NNPDF31_13TeV-powheg_M_400_800/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DY", "name": "ZToEE_50_120", "das": "/ZToEE_NNPDF31_13TeV-powheg_M_50_120/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DY", "name": "ZToEE_800_1400", "das": "/ZToEE_NNPDF31_13TeV-powheg_M_800_1400/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DY", "name": "ZToMuMu_120_200", "das": "/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DY", "name": "ZToMuMu_200_400", "das": "/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DY", "name": "ZToMuMu_400_800", "das": "/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DY", "name": "ZToMuMu_50_120", "das": "/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc":1, "type": "DY", "name": "ZToMuMu_800_1400", "das": "/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},


{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_200_dbs1p0", "das": "/BFFZprimeToMuMu_M-200-dbs1p0_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_200_dbs0p04", "das": "/BFFZprimeToMuMu_M-200_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_200_dbs0p5", "das": "/BFFZprimeToMuMu_M-200_dbs0p5_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_250_dbs0p04", "das": "/BFFZprimeToMuMu_M-250_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_300_dbs0p04", "das": "/BFFZprimeToMuMu_M-300_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_350_dbs0p04", "das": "/BFFZprimeToMuMu_M-350_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_350_dbs0p5", "das": "/BFFZprimeToMuMu_M-350_dbs0p5_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_350_dbs1p0", "das": "/BFFZprimeToMuMu_M-350_dbs1p0_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_400_dbs0p04", "das": "/BFFZprimeToMuMu_M-400_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_450_dbs0p04", "das": "/BFFZprimeToMuMu_M-450_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_500_dbs0p04", "das": "/BFFZprimeToMuMu_M-500_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_500_dbs0p5", "das": "/BFFZprimeToMuMu_M-500_dbs0p5_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "BFF", "name": "BFF_500_dbs1p0", "das": "/BFFZprimeToMuMu_M-500_dbs1p0_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DY", "name": "ZToEE_120_200", "das": "/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_120_200/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DY", "name": "ZToEE_200_400", "das": "/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_200_400/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DY", "name": "ZToEE_400_800", "das": "/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_400_800/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DY", "name": "ZToEE_50_120", "das": "/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_50_120/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DY", "name": "ZToEE_800_1400", "das": "/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_800_1400/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DY", "name": "ZToMuMu_120_200", "das": "/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DY", "name": "ZToMuMu_200_400", "das": "/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DY", "name": "ZToMuMu_400_800", "das": "/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DY", "name": "ZToMuMu_50_120", "das": "/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DY", "name": "ZToMuMu_800_1400", "das": "/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "TT", "name": "TT", "das": "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "ST", "name": "ST_anti", "das": "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "ST", "name": "ST", "das": "/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DB", "name": "WW", "das": "/WW_TuneCP5_13TeV-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DB", "name": "WZ", "das": "/WZ_TuneCP5_13TeV-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc":1, "type": "DB", "name": "ZZ", "das": "/ZZ_TuneCP5_13TeV-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},


    
{"era": 2016, "ismc":0, "type": "data", "name": "SingleMuonB", "das": "/SingleMuon/Run2016B-02Apr2020_ver1-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "SingleMuonB", "das": "/SingleMuon/Run2016B-02Apr2020_ver2-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "SingleMuonC", "das": "/SingleMuon/Run2016C-02Apr2020-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "SingleMuonD", "das": "/SingleMuon/Run2016D-02Apr2020-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "SingleMuonE", "das": "/SingleMuon/Run2016E-02Apr2020-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "SingleMuonF", "das": "/SingleMuon/Run2016F-02Apr2020-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "SingleMuonG", "das": "/SingleMuon/Run2016G-02Apr2020-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "SingleMuonH", "das": "/SingleMuon/Run2016H-02Apr2020-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "DoubleEGB", "das": "/DoubleEG/Run2016B-02Apr2020_ver1-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "DoubleEGB", "das": "/DoubleEG/Run2016B-02Apr2020_ver2-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "DoubleEGC", "das": "/DoubleEG/Run2016C-02Apr2020-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "DoubleEGD", "das": "/DoubleEG/Run2016D-02Apr2020-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "DoubleEGE", "das": "/DoubleEG/Run2016E-02Apr2020-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "DoubleEGF", "das": "/DoubleEG/Run2016F-02Apr2020-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "DoubleEGG", "das": "/DoubleEG/Run2016G-02Apr2020-v1/NANOAOD"},
{"era": 2016, "ismc":0, "type": "data", "name": "DoubleEGH", "das": "/DoubleEG/Run2016H-02Apr2020-v1/NANOAOD"},

{"era": 2017, "ismc":0, "type": "data", "name": "SingleMuonB", "das": "/SingleMuon/Run2017B-02Apr2020-v1/NANOAOD"},
{"era": 2017, "ismc":0, "type": "data", "name": "SingleMuonC", "das": "/SingleMuon/Run2017C-02Apr2020-v1/NANOAOD"},
{"era": 2017, "ismc":0, "type": "data", "name": "SingleMuonD", "das": "/SingleMuon/Run2017D-02Apr2020-v1/NANOAOD"},
{"era": 2017, "ismc":0, "type": "data", "name": "SingleMuonE", "das": "/SingleMuon/Run2017E-02Apr2020-v1/NANOAOD"},
{"era": 2017, "ismc":0, "type": "data", "name": "SingleMuonF", "das": "/SingleMuon/Run2017F-02Apr2020-v1/NANOAOD"},
{"era": 2017, "ismc":0, "type": "data", "name": "DoubleEGB", "das": "/DoubleEG/Run2017B-02Apr2020-v1/NANOAOD"},
{"era": 2017, "ismc":0, "type": "data", "name": "DoubleEGC", "das": "/DoubleEG/Run2017C-02Apr2020-v1/NANOAOD"},
{"era": 2017, "ismc":0, "type": "data", "name": "DoubleEGD", "das": "/DoubleEG/Run2017D-02Apr2020-v1/NANOAOD"},
{"era": 2017, "ismc":0, "type": "data", "name": "DoubleEGE", "das": "/DoubleEG/Run2017E-02Apr2020-v1/NANOAOD"},
{"era": 2017, "ismc":0, "type": "data", "name": "DoubleEGF", "das": "/DoubleEG/Run2017F-02Apr2020-v1/NANOAOD"},

{"era": 2018, "ismc":0, "type": "data", "name": "SingleMuonA", "das": "/SingleMuon/Run2018A-02Apr2020-v1/NANOAOD"},
{"era": 2018, "ismc":0, "type": "data", "name": "SingleMuonB", "das": "/SingleMuon/Run2018B-02Apr2020-v1/NANOAOD"},
{"era": 2018, "ismc":0, "type": "data", "name": "SingleMuonC", "das": "/SingleMuon/Run2018C-02Apr2020-v1/NANOAOD"},
{"era": 2018, "ismc":0, "type": "data", "name": "SingleMuonD", "das": "/SingleMuon/Run2018D-02Apr2020-v1/NANOAOD"},
{"era": 2018, "ismc":0, "type": "data", "name": "EGammaA", "das": "/EGamma/Run2018A-02Apr2020-v1/NANOAOD"},
{"era": 2018, "ismc":0, "type": "data", "name": "EGammaB", "das": "/EGamma/Run2018B-02Apr2020-v1/NANOAOD"},
{"era": 2018, "ismc":0, "type": "data", "name": "EGammaC", "das": "/EGamma/Run2018C-02Apr2020-v1/NANOAOD"},
{"era": 2018, "ismc":0, "type": "data", "name": "EGammaD", "das": "/EGamma/Run2018D-02Apr2020-v1/NANOAOD"},     


#dy samples for stitching
{"era": 2017, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_50_nlo", "das": "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_50_nlo_v2", "das": "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8_ext3-v1/NANOAODSIM"},
{"era": 2017, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_50_mad_v2", "das": "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv7-PU2017RECOSIMstep_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8_ext1-v1/NANOAODSIM"},
{"era": 2017, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_50_mad", "das": "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv7-PU2017RECOSIMstep_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_100_200", "das": "/DYJetsToLL_M-100to200_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_200_400", "das": "/DYJetsToLL_M-200to400_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_400_500", "das": "/DYJetsToLL_M-400to500_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_500_700", "das": "/DYJetsToLL_M-500to700_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_700_800", "das": "/DYJetsToLL_M-700to800_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},
{"era": 2017, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_800_1000", "das": "/DYJetsToLL_M-800to1000_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM"},


{"era": 2018, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_50_nlo", "das": "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_50_nlo_v2", "das": "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21_ext2-v1/NANOAODSIM"},
{"era": 2018, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_50_mad", "das": "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_100_200", "das": "/DYJetsToLL_M-100to200_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_200_400", "das": "/DYJetsToLL_M-200to400_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_400_500", "das": "/DYJetsToLL_M-400to500_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_500_700", "das": "/DYJetsToLL_M-500to700_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_700_800", "das": "/DYJetsToLL_M-700to800_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},
{"era": 2018, "ismc": 1.0, "type": "DYJets", "name": "DYJLL_M_800_1000", "das": "/DYJetsToLL_M-800to1000_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"},


]

import pandas as pd

sample_df = pd.DataFrame(sample_list)
