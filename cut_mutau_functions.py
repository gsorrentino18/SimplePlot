import numpy as np # TODO is importing this everywhere slowing things down? does python have IFDEF commands?

from calculate_functions import calculate_mt
from branch_functions import add_trigger_branches, add_DeepTau_branches

def make_mutau_cut(event_dictionary, DeepTau_version):
  '''
  Works similarly to 'make_ditau_cut'. 
  Notably, the mutau cuts are more complicated, but it is simple to 
  extend the existing methods as long as one can stomach the line breaks.
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_mutau = ["Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_iso",
                  "Muon_dxy", "Muon_dz", "Muon_charge", "Tau_dxy", "Tau_dz", "Tau_decayMode",
                  #"MET_pt", "MET_phi", "Muon_phi",
                  "PuppiMET_pt", "PuppiMET_phi",
                  "Lepton_tauIdx", "Lepton_muIdx", "l1_indices", "l2_indices"]
  #TODO add this "CleanJet_btagWP" (no effect in August skims since it was always 1)
  unpack_mutau = add_DeepTau_branches(unpack_mutau, DeepTau_version)
  unpack_mutau = add_trigger_branches(unpack_mutau, final_state_mode="mutau")
  unpack_mutau = (event_dictionary.get(key) for key in unpack_mutau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_mutau] # "*" unpacks a tuple
  pass_cuts, FS_mt = [], []
  FS_mu_pt, FS_mu_eta, FS_mu_phi, FS_mu_iso, FS_mu_dxy, FS_mu_dz, FS_mu_chg = [], [], [], [], [], [], []
  FS_tau_pt, FS_tau_eta, FS_tau_phi, FS_tau_dxy, FS_tau_dz = [], [], [], [], []
  # note these are in the same order as the variables in the first line of this function :)
      #MET_pt, MET_phi, muon_phi, tau_idx, mu_idx,\
  for i, lep_pt, lep_eta, lep_phi, lep_iso,\
      mu_dxy, mu_dz, mu_chg, tau_dxy, tau_dz, tau_decayMode,\
      MET_pt, MET_phi, tau_idx, mu_idx,\
      l1_idx, l2_idx, vJet, vMu, vEle, trg24mu, trg27mu, crosstrg in zip(*to_check):

    # some handling to figure out which FS index applies to what lepton
    # note for the DeepTauID we use the tau branch index directly instead of the lepton branch
    # (for tau branches we need the tau_idx, for lepton branches we can simply use the l1_idx, l2_idx)
    tauFSLoc, tauBranchLoc, muFSLoc, muBranchLoc = 999, 999, 999, 999
    if (tau_idx[l1_idx] != -1 and mu_idx[l2_idx] != -1):
      tauFSLoc = l1_idx
      tauBranchLoc = tau_idx[l1_idx]
      muLoc  = l2_idx
      muBranchLoc = mu_idx[l2_idx]
    elif (tau_idx[l2_idx] != -1 and mu_idx[l1_idx] != -1):
      tauFSLoc = l2_idx
      tauBranchLoc = tau_idx[l2_idx]
      muLoc  = l1_idx
      muBranchLoc = mu_idx[l1_idx]
    else:
      print("Should not print :)")

    muPtVal    = lep_pt[muLoc] 
    muEtaVal   = lep_eta[muLoc]
    muPhiVal   = lep_phi[muLoc]
    #muPhiVal   = muon_phi[muBranchLoc]
    muIsoVal   = lep_iso[muLoc]
    muDxyVal   = abs(mu_dxy[muBranchLoc])
    muDzVal    = mu_dz[muBranchLoc]
    muChgVal   = mu_chg[muBranchLoc]
    tauPtVal   = lep_pt[tauFSLoc] 
    tauEtaVal  = lep_eta[tauFSLoc]
    tauPhiVal  = lep_phi[tauFSLoc]
    tauDxyVal  = abs(tau_dxy[tauBranchLoc])
    tauDzVal   = tau_dz[tauBranchLoc]
    mtVal      = calculate_mt(muPtVal, muPhiVal, MET_pt, MET_phi)
    passMT     = (mtVal < 65.0) #(mtVal < 50.0) #mine
    #ROOTmtVal  = calculate_mt_pyROOT(muPtVal, muEtaVal, muPhiVal, mu_M[muLoc], MET_pt, MET_phi)
    #passROOTMT = (ROOTmtVal < 50.0)

    #passTauPtAndEta  = ((tauPtVal > 30.0) and (abs(tauEtaVal) < 2.3)) # mine
    passTauPtAndEta  = ((tauPtVal > 20.0) and (abs(tauEtaVal) < 2.5))
    #pass25MuPt   = ((trg24mu) and (muPtVal > 25.0) and (abs(muEtaVal) < 2.4)) #mine
    #pass28MuPt   = ((trg27mu) and (muPtVal > 28.0) and (abs(muEtaVal) < 2.4)) #mine
    pass25MuPt    = (trg24mu or trg27mu) and (muPtVal > 25.0) and (abs(muEtaVal) < 2.4)
    # HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1
    passMuPtCrossTrigger = ((crosstrg) and ((21.0 < muPtVal < 25.0) and (abs(muEtaVal) < 2.1))
                                       and ((tauPtVal > 32.0)       and (abs(tauEtaVal) < 2.1)) ) 

    # Medium v Jet, Tight v Muon, VVVLoose v Ele
    #passTauDT  = ((vJet[tauBranchLoc] >= 5) and (vMu[tauBranchLoc] >= 4) and (vEle[tauBranchLoc] >= 1))
    passTauDT  = ((vJet[tauBranchLoc] >= 5) and (vMu[tauBranchLoc] >= 4) and (vEle[tauBranchLoc] >= 2))

    skip_DM2 = (tau_decayMode[tauBranchLoc] != 2)
    #restrict_tau_decayMode = (tau_decayMode[tauBranchLoc] == 0)

    #if (passMT and (passTauPtAndEta and (pass25MuPt or pass28MuPt or passMuPtCrossTrigger)) and passTauDT): #mine
    #if (passMT and (passTauPtAndEta and pass25MuPt and passTauDT) and skip_DM2 and restrict_tau_decayMode):
    if (passMT and (passTauPtAndEta and pass25MuPt and passTauDT) and skip_DM2):
      pass_cuts.append(i)
      FS_tau_pt.append(tauPtVal)
      FS_tau_eta.append(tauEtaVal)
      FS_tau_phi.append(tauPhiVal)
      FS_tau_dxy.append(tauDxyVal)
      FS_tau_dz.append(tauDzVal)
      FS_mu_pt.append(muPtVal)
      FS_mu_eta.append(muEtaVal)
      FS_mu_phi.append(muPhiVal)
      FS_mu_iso.append(muIsoVal)
      FS_mu_dxy.append(muDxyVal)
      FS_mu_dz.append(muDzVal)
      FS_mu_chg.append(muChgVal)
      FS_mt.append(mtVal)

  event_dictionary["pass_cuts"] = np.array(pass_cuts)
  event_dictionary["FS_mu_pt"]  = np.array(FS_mu_pt)
  event_dictionary["FS_mu_eta"] = np.array(FS_mu_eta)
  event_dictionary["FS_mu_phi"] = np.array(FS_mu_phi)
  event_dictionary["FS_mu_iso"] = np.array(FS_mu_iso)
  event_dictionary["FS_mu_dxy"] = np.array(FS_mu_dxy)
  event_dictionary["FS_mu_dz"]  = np.array(FS_mu_dz)
  event_dictionary["FS_mu_chg"] = np.array(FS_mu_chg)
  event_dictionary["FS_tau_pt"]  = np.array(FS_tau_pt)
  event_dictionary["FS_tau_eta"] = np.array(FS_tau_eta)
  event_dictionary["FS_tau_phi"] = np.array(FS_tau_phi)
  event_dictionary["FS_tau_dxy"] = np.array(FS_tau_dxy)
  event_dictionary["FS_tau_dz"]  = np.array(FS_tau_dz)
  event_dictionary["FS_mt"]    = np.array(FS_mt)
  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after mutau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary
