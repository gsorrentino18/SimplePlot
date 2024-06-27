import math
import numpy as np # TODO is importing this everywhere slowing things down? does python have IFDEF commands?

from calculate_functions import calculate_mt
from branch_functions import add_trigger_branches, add_DeepTau_branches, add_gen_branches

def isInsideRcone(eleEta, elePhi, tauEta, tauPhi, Rcone):
    deltaR=math.sqrt( (eleEta-tauEta)**2 + (elePhi-tauPhi)**2 )
    if (deltaR < Rcone): 
       return True
    else:
       return False

def genmatching(event_dictionary, event_number, recoEta, recoPhi, DEBUG):
  unpack_gen = ["nGenPart", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_pdgId", "GenPart_status", "GenPart_statusFlags"]
  unpack_gen = (event_dictionary.get(key) for key in unpack_gen)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_gen]
  isGM=False
  for i, nGen, genPt, genEta, genPhi, genPdgId, genStatus, genFlag in zip(*to_check):
     #if (i!=event_number): continue
     for n in range(nGen):
        if (genPdgId[n] == 11 and genStatus[n]==1 and genFlag[n] == 13884):
           isGM = isGenMatch(recoEta, recoPhi, genEta[n], genPhi[n], 0.4)
           if (DEBUG and isGM): print("n: ", n, " genStatus: ", genStatus[n], " genPdgId:", genPdgId[n], " genStatusFlag: ", genFlag[n], " isGM:", isGM)
           #if (isGM): break   
  return isGM

def make_etau_cut(event_dictionary, DeepTau_version, skip_DeepTau=False):
  '''
  Works similarly to 'make_ditau_cut'. 
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_etau = ["Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_iso",
                 "Electron_dxy", "Electron_dz", "Electron_charge", "Tau_dxy", "Tau_dz", "Tau_charge", "Tau_decayMode",
                 "PuppiMET_pt", "PuppiMET_phi", "Lepton_tauIdx", "Lepton_elIdx", "l1_indices", "l2_indices",
                 "Electron_mvaNoIso_WP90", "Electron_pfRelIso03_all"]
  unpack_etau = add_DeepTau_branches(unpack_etau, DeepTau_version)
  unpack_etau = add_trigger_branches(unpack_etau, final_state_mode="etau")
  unpack_etau = (event_dictionary.get(key) for key in unpack_etau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_etau]

  pass_cuts, FS_mt = [], []
  FS_el_pt, FS_el_eta, FS_el_phi, FS_el_iso, FS_el_dxy, FS_el_dz, FS_el_chg = [], [], [], [], [], [], []
  FS_tau_pt, FS_tau_eta, FS_tau_phi, FS_tau_dxy, FS_tau_dz, FS_tau_chg = [], [], [], [], [], []
  xtrigger_flag = []
 
  for i, lep_pt, lep_eta, lep_phi, lep_iso,\
      el_dxy, el_dz, el_chg, tau_dxy, tau_dz, tau_chg, tau_decayMode,\
      MET_pt, MET_phi, tau_idx, el_idx,l1_idx, l2_idx, \
      noIso, el_iso, \
      vJet, vMu, vEle, trg30el, trg32el, trg35el, crosstrg in zip(*to_check):

    # some handling to figure out which FS index applies to what lepton
    # note for the DeepTauID we use the tau branch index directly instead of the lepton branch
    # (for tau branches we need the tau_idx, for lepton branches we can simply use the l1_idx, l2_idx)
    tauFSLoc, tauBranchLoc, elFSLoc, elBranchLoc = 999, 999, 999, 999
    if (tau_idx[l1_idx] != -1 and el_idx[l2_idx] != -1):
      tauFSLoc = l1_idx
      tauBranchLoc = tau_idx[l1_idx]
      elFSLoc  = l2_idx
      elBranchLoc = el_idx[l2_idx]
    elif (tau_idx[l2_idx] != -1 and el_idx[l1_idx] != -1):
      tauFSLoc = l2_idx
      tauBranchLoc = tau_idx[l2_idx]
      elFSLoc  = l1_idx
      elBranchLoc = el_idx[l1_idx]
    else:
      print("Should not print :)")

    elPtVal    = lep_pt[elFSLoc] 
    elEtaVal   = lep_eta[elFSLoc]
    elPhiVal   = lep_phi[elFSLoc]
    elIsoVal   = lep_iso[elFSLoc]
    elDxyVal   = abs(el_dxy[elBranchLoc])
    elDzVal    = el_dz[elBranchLoc]
    elChgVal   = el_chg[elBranchLoc]
    tauPtVal   = lep_pt[tauFSLoc] 
    tauEtaVal  = lep_eta[tauFSLoc]
    tauPhiVal  = lep_phi[tauFSLoc]
    tauDxyVal  = abs(tau_dxy[tauBranchLoc])
    tauDzVal   = tau_dz[tauBranchLoc]
    tauChgVal  = tau_chg[tauBranchLoc]
    mtVal      = calculate_mt(elPtVal, elPhiVal, MET_pt, MET_phi)
    #passMT     = (mtVal < 60.0) #TAU CQM
    passMT     = (mtVal < 50.0)

    pfRelIso03_all = el_iso[elBranchLoc]

    #ROOTmtVal  = calculate_mt_pyROOT(muPtVal, muEtaVal, muPhiVal, mu_M[muLoc], MET_pt, MET_phi)
    #passROOTMT = (ROOTmtVal < 50.0)

    passTauPtAndEta  = ((tauPtVal > 30.0) and (abs(tauEtaVal) < 2.3))
    #passTauPtAndEta  = ((tauPtVal > 20.0) and (abs(tauEtaVal) < 2.5)) #TAU CQM
    
    pass31ElPt   = ((trg30el) and (elPtVal > 31.0) and (abs(elEtaVal) < 2.1)) 
    #pass33ElPt   = ((trg32el) and (elPtVal > 33.0) and (abs(elEtaVal) < 2.1) and (pfRelIso03_all<0.10)) #TAU CQM 
    pass33ElPt   = ((trg32el) and (elPtVal > 33.0) and (abs(elEtaVal) < 2.1)) 
    pass36ElPt   = ((trg35el) and (elPtVal > 36.0) and (abs(elEtaVal) < 2.1)) 
    # upper bound on cross trigger will change if lower single electron trigger included
    # HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1
    passElPtCrossTrigger = ((crosstrg) and ((25.0 < elPtVal < 31.0) and (abs(elEtaVal) < 2.1))
                                       and ((tauPtVal > 35.0)       and (abs(tauEtaVal) < 2.1)) ) 

    # Medium v Jet, VLoose v Muon, Tight v Ele
    passTauDT  = ((vJet[tauBranchLoc] >= 5) and (vMu[tauBranchLoc] >= 4) and (vEle[tauBranchLoc] >= 6)) #ele Tight >=6

    crosstrg_only = False
    if ((crosstrg==1) and (trg30el==0) and (trg32el==0) and (trg35el==0)): crosstrg_only = True

    #if ((passMT and passTauPtAndEta and (pass31ElPt or pass33ElPt or pass36ElPt) and (crosstrg_only==False) and passTauDT)
    #   or (passMT and passTauPtAndEta and (pass31ElPt or pass33ElPt or pass36ElPt) and (crosstrg_only==False) and skip_DeepTau)):

    if ((passMT and passTauPtAndEta and (pass31ElPt) and passTauDT)
       or (passMT and passTauPtAndEta and (pass31ElPt) and skip_DeepTau)):

    #if ((passMT and passTauPtAndEta and (pass31ElPt or pass33ElPt or pass36ElPt or passElPtCrossTrigger) and passTauDT)
    #   or (passMT and passTauPtAndEta and (pass31ElPt or pass33ElPt or pass36ElPt or passElPtCrossTrigger) and skip_DeepTau)):
      pass_cuts.append(i)
      FS_el_pt.append(elPtVal)
      FS_el_eta.append(elEtaVal)
      FS_el_phi.append(elPhiVal)
      FS_el_iso.append(elIsoVal)
      FS_el_dxy.append(elDxyVal)
      FS_el_dz.append(elDzVal)
      FS_el_chg.append(elChgVal)
      FS_tau_pt.append(tauPtVal)
      FS_tau_eta.append(tauEtaVal)
      FS_tau_phi.append(tauPhiVal)
      FS_tau_dxy.append(tauDxyVal)
      FS_tau_dz.append(tauDzVal)
      FS_tau_chg.append(tauChgVal)
      FS_mt.append(mtVal)
      xtrigger_flag.append(crosstrg)

  event_dictionary["pass_cuts"]  = np.array(pass_cuts)
  event_dictionary["FS_el_pt"]   = np.array(FS_el_pt)
  event_dictionary["FS_el_eta"]  = np.array(FS_el_eta)
  event_dictionary["FS_el_phi"]  = np.array(FS_el_phi)
  event_dictionary["FS_el_iso"]  = np.array(FS_el_iso)
  event_dictionary["FS_el_dxy"]  = np.array(FS_el_dxy)
  event_dictionary["FS_el_dz"]   = np.array(FS_el_dz)
  event_dictionary["FS_el_chg"]  = np.array(FS_el_chg)
  event_dictionary["FS_tau_pt"]  = np.array(FS_tau_pt)
  event_dictionary["FS_tau_eta"] = np.array(FS_tau_eta)
  event_dictionary["FS_tau_phi"] = np.array(FS_tau_phi)
  event_dictionary["FS_tau_dxy"] = np.array(FS_tau_dxy)
  event_dictionary["FS_tau_dz"]  = np.array(FS_tau_dz)
  event_dictionary["FS_tau_chg"] = np.array(FS_tau_chg)
  event_dictionary["FS_mt"]      = np.array(FS_mt)
  event_dictionary["xtrigger_flag"] = np.array(xtrigger_flag)

  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after etau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary

def make_etau_AR_cut(event_dictionary, DeepTau_version): 
  unpack_etau_AR_vars = ["event", "Lepton_tauIdx", "Lepton_elIdx", "Lepton_iso", "l1_indices", "l2_indices"]
  unpack_etau_AR_vars = add_DeepTau_branches(unpack_etau_AR_vars, DeepTau_version)
  unpack_etau_AR_vars = (event_dictionary.get(key) for key in unpack_etau_AR_vars)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_etau_AR_vars]
  pass_AR_cuts = []
  for i, event, tau_idx, ele_idx, lep_iso, l1_idx, l2_idx, vJet, _, _ in zip(*to_check):
    # keep indices where tau fails and muon passes iso 
    ele_lep_idx = l1_idx if ele_idx[l1_idx] != -1 else l2_idx
    ele_iso = lep_iso[ele_lep_idx]
    tau_branchIdx  = tau_idx[l1_idx] + tau_idx[l2_idx] + 1
    if ((vJet[tau_branchIdx] < 5) and (ele_iso<0.15)):
      pass_AR_cuts.append(i)
    
  event_dictionary["pass_AR_cuts"] = np.array(pass_AR_cuts)
  return event_dictionary

