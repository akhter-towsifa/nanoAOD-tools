import ROOT
import os
import numpy as np
import json
import copy
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class lepSFProducer(Module):
    def __init__(self, muonSelectionTag, electronSelectionTag):
        self.mu_f_name = ["trigger","ID","ISO"]
        self.flip_eta = [0,0,0]
        self.not_error = [0,0,0]
        self.error = [1,1,1]
        #identifies which f is for the trigger
        self.trigger_index = 0
        self.mu_f_sys_name = ["triggerUp", "triggerDown", "ID","ISO"]
        self.sys_flip_eta = [1,1,0,0]
        self.sys_error = [0,0,1,1]
        #2016 muon legacy
        if muonSelectionTag=="muonSF_2016_GH_legacy":
            self.mu_weights = [1]
            self.mu_f=[[
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/2016_RunGH_SF_ID.root",
                "muon_legacy/2016_RunGH_SF_ISO.root"
            ]]
            self.mu_h = [
                "SF_2016_var",
                "NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt_stat",
                 "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt_stat"
            ]
            self.mu_sys_f=[[
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/2016_RunGH_SF_ID.root",
                "muon_legacy/2016_RunGH_SF_ISO.root"
            ]]
            self.mu_sys_h = [
                "SF_2016_errorUpper",
                "SF_2016_errorUpper"
                "NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt_syst",
                "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt_syst"
            ]

        if muonSelectionTag=="muonSF_2016_BCDEF_legacy":
            self.mu_weights = [1]
            self.mu_f=[[
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/2016_RunBCDEF_SF_ID.root",
                "muon_legacy/2016_RunBCDEF_SF_ISO.root"
            ]]
            self.mu_h = [
                "SF_2016_errorUpper",
                "SF_2016_errorUpper"
                "NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt_stat",
                "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt_stat"
            ]
            self.mu_sys_f=[[  
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/2016_RunBCDEF_SF_ID.root",
                "muon_legacy/2016_RunBCDEF_SF_ISO.root"
            ]]
            self.mu_sys_h = [
                "SF_2016_errorUpper",
                "SF_2016_errorUpper", 
                "NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt_syst",
                "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt_syst"
            ]


        if muonSelectionTag=="muonSF_2016_weighted_legacy":
            luminosity_proportion = 0.55
            self.mu_weights = [luminosity_proportion,1-luminosity_proportion] #luminositiy weighted for different trigger menus
            self.mu_f=[
                          [
                              "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                              "muon_legacy/2016_RunBCDEF_SF_ID.root",
                              "muon_legacy/2016_RunBCDEF_SF_ISO.root"
                          ],
                          [
                              "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                              "muon_legacy/2016_RunGH_SF_ID.root",
                              "muon_legacy/2016_RunGH_SF_ISO.root"
                          ]
                  ]
            self.mu_h = [
                "SF_2016_var",
                "NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt_stat",
                "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt_stat"
            ]
            self.mu_sys_f=[[
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/2016_RunBCDEF_SF_ID.root",
                "muon_legacy/2016_RunBCDEF_SF_ISO.root"
            ],
            [
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/2016_RunGH_SF_ID.root",
                "muon_legacy/2016_RunGH_SF_ISO.root"
            ],
            ]
            self.mu_sys_h = [
                "SF_2016_errorUpper",
                "SF_2016_errorUpper",  
                "NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt_syst",
                "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt_syst"
            ]
        #2017 muon
        if muonSelectionTag=="muonSF_2017":
            self.mu_weights = [1]
            self.mu_f=[[
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/2017_RunBCDEF_SF_ID_syst.root",
                "muon_legacy/2017_RunBCDEF_SF_ISO_syst.root"
            ]]
            self.mu_h = [
                "SF_2017_var",
                "NUM_HighPtID_DEN_genTracks_pair_newTuneP_probe_pt_abseta_stat",
                "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_pair_newTuneP_probe_pt_abseta_stat"
            ]
            self.mu_sys_f=[[
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/2017_RunBCDEF_SF_ID_syst.root",
                "muon_legacy/2017_RunBCDEF_SF_ISO_syst.root"
            ]]
            self.mu_sys_h = [
                "SF_2017_errorUpper",
                "SF_2017_errorUpper", 
                "NUM_HighPtID_DEN_genTracks_pair_newTuneP_probe_pt_abseta_syst",
                "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_pair_newTuneP_probe_pt_abseta_syst"
            ]
        #2018 muon
        if muonSelectionTag=="muonSF_2018":
            self.mu_weights = [1]
            self.mu_f=[[
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/2018_RunABCD_SF_ID.root",
                "muon_legacy/2018_RunABCD_SF_ISO.root"
            ]]
            self.mu_h = [
                "SF_2018_var",
                "NUM_HighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta_stat",
                "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_pair_newTuneP_probe_pt_abseta_stat"
            ]
            self.mu_sys_f=[[
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root",
                "muon_legacy/2018_RunABCD_SF_ID.root",
                "muon_legacy/2018_RunABCD_SF_ISO.root"
            ]]
            self.mu_sys_h = [
                "SF_2018_errorUpper",
                "SF_2018_errorUpper", 
                "NUM_HighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta_syst",
                "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_pair_newTuneP_probe_pt_abseta_syst"
            ]

        #egamma
        if electronSelectionTag=="egamma_2016_legacy":
            self.el_SF_EB = 0.983
            self.el_SF_EE = 0.991
            self.el_stat_EB = 0.001
            self.el_stat_EE = 0.001
            self.el_sys_EB = [.01,0.0022,.03]
            self.el_sys_EE = [.01,0.0143,.04]
        if electronSelectionTag=="egamma_2017":
            self.el_SF_EB = 0.967
            self.el_SF_EE = 0.973
            self.el_stat_EE = 0.001
            self.el_stat_EB = 0.002
            self.el_sys_EB = [.01,0.0022,.03]
            self.el_sys_EE = [.02,0.0143,.05]
        if electronSelectionTag=="egamma_2018":
            self.el_SF_EB = 0.969
            self.el_SF_EE = 0.984
            self.el_stat_EE = 0.000
            self.el_stat_EB = 0.001
            self.el_sys_EB = [.01,0.0022,.03]
            self.el_sys_EE = [.02,0.0143,.05]

        self.mu_sys_f = [["%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/leptonSF/"
            % os.environ['CMSSW_BASE'] + f for f in f_list] for f_list in self.mu_sys_f]
        self.mu_f = [["%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/leptonSF/"
            % os.environ['CMSSW_BASE'] + f for f in f_list] for f_list in self.mu_f]

       # if "/LeptonEfficiencyCorrector_cc.so" not in ROOT.gSystem.GetLibraries(
       # ):
       #     print("Load C++ Worker")
       #     ROOT.gROOT.ProcessLine(
       #         ".L %s/src/PhysicsTools/NanoAODTools/python/postprocessing/helpers/LeptonEfficiencyCorrector.cc+"

       #         % os.environ['CMSSW_BASE'])
        for library in [ "libCondFormatsJetMETObjects", "libPhysicsToolsNanoAODTools" ]:
            if library not in ROOT.gSystem.GetLibraries():
                print("Load Library '%s'" % library.replace("lib", ""))
                ROOT.gSystem.Load(library) 
    def beginJob(self):
        self._worker_mu = [[] for i in range(len(self.mu_f))]
        for i in range(len(self.mu_f)):
            for j, f_list in enumerate(self.mu_f[i]):
                print(f_list)
                roo_vec_file = ROOT.std.vector(str)(1)
                roo_vec_hist = ROOT.std.vector(str)(1)
                roo_vec_file[0] = f_list
                roo_vec_hist[0] = self.mu_h[j]
                self._worker_mu[i].append( ROOT.LeptonEfficiencyCorrectorCppWorker(roo_vec_file,roo_vec_hist) )
        self._worker_sys_mu = [[] for i in range(len(self.mu_f))]
        for i in range(len(self.mu_sys_f)):
            for j, f_list in enumerate(self.mu_sys_f[i]):
                roo_vec_file = ROOT.std.vector(str)(1)
                roo_vec_hist = ROOT.std.vector(str)(1)
                roo_vec_file[0] = f_list
                roo_vec_hist[0] = self.mu_sys_h[j]
                self._worker_sys_mu[i].append( ROOT.LeptonEfficiencyCorrectorCppWorker(roo_vec_file,roo_vec_hist) )
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for type_SF in self.mu_f_name:
            self.out.branch("Muon_effSF_{}".format(type_SF), "F", lenVar="nMuon")
            self.out.branch("Muon_effSF_stat_{}".format(type_SF), "F", lenVar="nMuon")
        for type_SF in self.mu_f_sys_name:
            self.out.branch("Muon_effSF_sys_{}".format(type_SF), "F", lenVar="nMuon")
        self.out.branch("Electron_effSF", "F", lenVar="nElectron")
        self.out.branch("Electron_effSF_stat", "F", lenVar="nElectron")
        self.out.branch("Electron_effSF_sys", "F", lenVar="nElectron")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def computeSFlist(self, worker, hist_names, muons, errors, flips):
            sf_mu = [[] for i in range(len(worker))]
            for i in range(len(worker)):
                for error, flip, j in zip(errors, flips, range(len(worker[i]))):
                    LEC_object = worker[i][j]
                    if "abseta" in hist_names[j]: pdgId = 13
                    else: pdgId = 99999
                    if not error: sf_mu[i].append([LEC_object.getSF(pdgId, *flip_list([mu.pt, mu.eta], flip)) for mu in muons ])
                    else: sf_mu[i].append([LEC_object.getSFErr(pdgId, *flip_list([mu.pt, mu.eta], flip)) for mu in muons ])
            averaged_sf_mu = np.multiply(self.mu_weights[0], sf_mu[0])
            for i in range(len(sf_mu)-1):
                averaged_sf_mu = np.add(averaged_sf_mu,np.multiply(self.mu_weights[i+1],sf_mu[i+1]))
            return averaged_sf_mu
    def computeElSF(self, el):
        eT = el.p4().Et()
        eta = el.eta
        if abs(eta) > 1.566 and abs(eta) < 2.5: SF = self.el_SF_EE
        elif abs(eta) < 1.4442: SF = self.el_SF_EB
        else: return 0
        return SF
    def computeElSys(self, el):
        eT = el.p4().Et()
        eta = el.eta
        if abs(eta) > 1.566 and abs(eta) < 2.5: sys_variables = self.el_sys_EE
        elif abs(eta) < 1.4442: sys_variables = self.el_sys_EB
        else: return -1
        if eT < 90: return sys_variables[0]
        else: return min((1+(eT-90)*sys_variables[1])/100,sys_variables[2])
    def computeElStat(self, el):
        eta = el.eta
        if abs(eta) > 1.566 and abs(eta) < 2.5: stat_variables = self.el_stat_EE
        elif abs(eta) < 1.4442: stat_variables = self.el_stat_EE
        else: return -0
        return stat_variables
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        trigObjs = Collection(event, "TrigObj")
        muonTrigObjs = list(filter(muon_selection, trigObjs))
        averaged_sf_mu = self.computeSFlist(self._worker_mu,self.mu_h,muons, self.not_error, self.flip_eta)
        averaged_sf_mu_stat = self.computeSFlist(self._worker_mu,self.mu_h,muons, self.error, self.flip_eta)
        # trigger histogram error and center value is switched, so flip them
        trigger_center_value = copy.deepcopy(averaged_sf_mu[0])
        trigger_error = copy.deepcopy(averaged_sf_mu_stat[0])
        averaged_sf_mu[0] = trigger_error
        averaged_sf_mu_stat[0] = trigger_center_value
        averaged_sf_mu_sys = self.computeSFlist(self._worker_sys_mu,self.mu_sys_h,muons, self.sys_error, self.sys_flip_eta)
        #set trigger sf to 1 if not trigger object
        for i, muon in enumerate(muons):
            muonIsTriggerObj = False
            for muonTrigObj in muonTrigObjs:
                if deltaR(muonTrigObj, muon) < 0.4:
                    muonIsTriggerObj = True
                    break
            if not muonIsTriggerObj:
                averaged_sf_mu[self.trigger_index][i] = 1
                averaged_sf_mu_stat[self.trigger_index][i] = 0
                averaged_sf_mu_sys[0][i] = 0
                averaged_sf_mu_sys[1][i] = 0       
        for i, type_SF in enumerate(self.mu_f_name):
            self.out.fillBranch("Muon_effSF_{}".format(type_SF), averaged_sf_mu[i])
            self.out.fillBranch("Muon_effSF_stat_{}".format(type_SF), averaged_sf_mu_stat[i])
        for i, type_SF in enumerate(self.mu_f_sys_name):
            self.out.fillBranch("Muon_effSF_sys_{}".format(type_SF), averaged_sf_mu_sys[i])

        sf_el =  map(self.computeElSF, electrons)
        el_sys =  map(self.computeElSys, electrons)
        el_sys = np.multiply(el_sys,sf_el)
        el_stat = map(self.computeElStat, electrons)

        self.out.fillBranch("Electron_effSF", sf_el)
        self.out.fillBranch("Electron_effSF_stat", el_stat)
        self.out.fillBranch("Electron_effSF_sys", el_sys)
        return True
                      
def deltaR(obj1, obj2):
    return ((obj1.phi-obj2.phi)**2+(obj1.eta-obj2.eta)**2)**.5
def isKthBitSet(n, k): return bool(n & (1 << (k-1)))
def muon_selection(trigObj):
    return (trigObj.id==13) & isKthBitSet(trigObj.filterBits, 11)
def flip_list(arr, flip):
    if flip: arr.reverse()
    return arr
        
# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

lepSF2016 = lambda : lepSFProducer( "muonSF_2016_weighted_legacy", "egamma_2016_legacy")
lepSF2016_GH_legacy = lambda : lepSFProducer( "muonSF_2016_GH_legacy", "egamma_2016_legacy")
lepSF2016_BtoF_legacy = lambda : lepSFProducer( "muonSF_2016_BCDEF_legacy", "egamma_2016_legacy")

lepSF2017 = lambda : lepSFProducer( "muonSF_2017", "egamma_2017")

lepSF2018 = lambda : lepSFProducer( "muonSF_2018", "egamma_2018")
