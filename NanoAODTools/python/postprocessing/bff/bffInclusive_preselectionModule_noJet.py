'''
This module should be able to create a new branch storing diLepton
mass given the following pass-- a) 2 muons with opposite signs,
b) no additional jets
'''

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.bff import bffBtagEff #contains jet information
import numpy as np
from root_numpy import tree2array


class bffInclusivePreselProducer_noJet(Module):
  def __init__(self):
    self.muSel = lambda x,pt: ((x.pt_corrected > pt) & (abs(x.eta) < 2.4) & (x.tightId > 0) 
                              & (x.pfRelIso04_all < 0.25))
    self.eleSel = lambda x,pt: ((x.pt > pt) & (abs(x.eta) < 2.4) & x.cutBased_HEEP > 0)
    self.diLepMass = -1
    self.diLepPt = -1
    self.lep_1 = ROOT.TLorentzVector()
    self.lep_2 = ROOT.TLorentzVector()
    self.alljetSel = lambda j: ((j.pt > 20) & (abs(j.eta) < 2.4) & ((j.jetId >> 1) & 1) & ((j.puId & 1) | (j.pt>50)))
    self.isMC = True
    self.btagWP = btagWP
    def deepcsv(jet):
        return jet.btagDeepB > self.btagWP
    def deepflavour(jet):
        return jet.btagDeepFlavB > self.btagWP
        #set right filtering function
    if btag_type=="deepcsv":
        self.select_btag = deepcsv
    elif btag_type=="deepflavour":
        self.select_btag = deepflavour

    pass
  
  def beginJob(self):
    pass
  def endJob(self):
    pass
  
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    if wrappedOutputTree.tree().GetListOfBranches().FindObject("Jet_btagSF"):
        self.isMC = True
    else: 
        self.isMC = False

    self.out = wrappedOutputTree
    self.out.branch("DiLepMass", "F")
    self.out.branch("DiLepPt", "F")
    self.out.branch("DiLepMass_noJet", "F") #new mass branch created
    self.out.branch("DiLepPt_noJet", "F")   #new pt branch created
    self.out.branch("IncMumu", "I")
    self.out.branch("IncEe", "I")
    self.out.branch("IncMumu_noJet", "I")   #new branch
    self.out.branch("IncEmu", "I")
    self.out.branch("IncMumuLowPt", "I")
    self.out.branch("IncEeLowPt", "I")
    self.out.branch("IncMumuLowPt_noJet", "I")  #new branch
    self.out.branch("IncEmuLowPt", "I")
    self.Pass = ROOT.std.vector(ROOT.std.vector('int'))()
    self.Total= ROOT.std.vector(ROOT.std.vector('int'))()
    self.out._tree.Branch("PassBtag", self.Pass)
    self.out._tree.Branch("TotalBtag",self.Total)

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass
  
  def selectDiMu(self, electrons, muons):
    if len(muons) != 2:
      return False
    if len(electrons) != 0:
      return False
    if (muons[0].charge+muons[1].charge) != 0:
      return False
    self.lep_1 = muons[0].p4()*(muons[0].pt_corrected/muons[0].pt)
    self.lep_2 = muons[1].p4()*(muons[1].pt_corrected/muons[1].pt)
    diLep = self.lep_1 + self.lep_2
    self.diLepMass = diLep.M()
    self.diLepPt = diLep.Pt()
    self.out.fillBranch("DiLepMass", self.diLepMass)
    self.out.fillBranch("DiLepPt", self.diLepPt)
    return True

  def selectDiEle(self, electrons, muons):
    if len(electrons) != 2:
      return False
    if len(muons) != 0:
      return False
    if (electrons[0].charge+electrons[1].charge) != 0:
      return False
    self.lep_1 = electrons[0].p4()
    self.lep_2 = electrons[1].p4()
    diLep = self.lep_1 + self.lep_2
    self.diLepMass = diLep.M()
    self.diLepPt = diLep.Pt()
    self.out.fillBranch("DiLepMass", self.diLepMass)
    self.out.fillBranch("DiLepPt", self.diLepPt)
    return True
  
  def selectDiMu_noJet(self, electrons, muons):
    if len(muons) != 2:
      return False
    if len(electrons) != 0:
      return False
    if (muons[0].charge+muons[1].charge) != 0:
      return False
    if len(Jet) != 0:
      return False    #requires no additional jets
    self.lep_1 = muons[0].p4()*(muons[0].pt_corrected/muons[0].pt)
    self.lep_2 = muons[1].p4()*(muons[1].pt_corrected/muons[1].pt)
    diLep = self.lep_1 + self.lep_2
    self.diLepMass = diLep.M()
    self.diLepPt = diLep.Pt()
    self.out.fillBranch("DiLepMass_noJet", self.diLepMass)
    self.out.fillBranch("DiLepPt_noJet", self.diLepPt)
    return True

  def selectEleMu(self, electrons, muons):
    if len(electrons) != 1:
      return False
    if len(muons) != 1:
      return False
    if (electrons[0].charge+muons[0].charge) != 0:
      return False
    self.lep_1 = electrons[0].p4()
    self.lep_2 = muons[0].p4()*(muons[0].pt_corrected/muons[0].pt)
    diLep = self.lep_1 + self.lep_2
    self.diLepMass = diLep.M()
    self.diLepPt = diLep.Pt()
    self.out.fillBranch("DiLepMass", self.diLepMass)
    self.out.fillBranch("DiLepPt", self.diLepPt)
    return True
  
  def analyze(self, event):
    '''process event, return True (go to next module) or False (fail, go to next event)'''
    electrons = sorted(filter(lambda x: self.eleSel(x,53), Collection(event, "Electron")), key=lambda x: x.pt)
    muons = sorted(filter(lambda x: self.muSel(x,53), Collection(event, "Muon")), key=lambda x: x.pt_corrected)

    electronsLowPt = sorted(filter(lambda x: self.eleSel(x,24), Collection(event, "Electron")), key=lambda x: x.pt)
    muonsLowPt = sorted(filter(lambda x: self.muSel(x,24), Collection(event, "Muon")), key=lambda x: x.pt_corrected)

    isDiMu = self.selectDiMu(electrons, muons)
    isDiEle = self.selectDiEle(electrons, muons)
    isDiMu_noJet = self.selectDiMu_noJet(electrons, muons)
    isEleMu = self.selectEleMu(electrons, muons)

    nLowPtLep = len(electronsLowPt)+len(muonsLowPt)

    self.out.fillBranch("IncMumu", isDiMu)
    self.out.fillBranch("IncEe", isDiEle)
    self.out.fillBranch("IncMumu_noJet", isDiMu_noJet)
    self.out.fillBranch("IncEmu", isEleMu)

    self.out.fillBranch("IncMumuLowPt", isDiMu and nLowPtLep<3)
    self.out.fillBranch("IncEeLowPt", isDiEle and nLowPtLep<3)
    self.out.fillBranch("IncMumuLowPt_noJet", isDiMu_noJet and nLowPtLep<3)
    self.out.fillBranch("IncEmuLowPt", isEleMu and nLowPtLep<3)

    return (isDiMu or isDiEle or isDiMu_noJet or isEleMu)

# define modules using the syntax 'name = lambda: constructor' to avoid having them loaded when not needed
bffInclusivePreselModuleConstr = lambda: bffInclusivePreselProducer_noJet()