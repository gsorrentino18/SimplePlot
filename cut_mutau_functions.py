import numpy as np # TODO is importing this everywhere slowing things down? does python have IFDEF commands?

from calculate_functions import calculate_mt
from branch_functions import add_trigger_branches, add_DeepTau_branches

def make_mutau_cut(event_dictionary, DeepTau_version, skip_DeepTau=False):
  '''
  Works similarly to 'make_ditau_cut'. 
  Notably, the mutau cuts are more complicated, but it is simple to 
  extend the existing methods as long as one can stomach the line breaks.
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_mutau = ["Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_iso",
                  "Muon_dxy", "Muon_dz", "Muon_charge", "Tau_dxy", "Tau_dz", "Tau_charge", "Tau_decayMode",
                  "PuppiMET_pt", "PuppiMET_phi",
                  "Lepton_tauIdx", "Lepton_muIdx", "l1_indices", "l2_indices", 
                  #"Tau_rawPNetVSjet", "Tau_rawPNetVSmu", "Tau_rawPNetVSe"
                  ]
  #TODO add this "CleanJet_btagWP" (no effect in August skims since it was always 1)
  unpack_mutau = add_DeepTau_branches(unpack_mutau, DeepTau_version)
  unpack_mutau = add_trigger_branches(unpack_mutau, final_state_mode="mutau")
  unpack_mutau = (event_dictionary.get(key) for key in unpack_mutau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_mutau] # "*" unpacks a tuple
  pass_cuts, FS_mt = [], []
  FS_mu_pt, FS_mu_eta, FS_mu_phi, FS_mu_iso, FS_mu_dxy, FS_mu_dz, FS_mu_chg = [], [], [], [], [], [], []
  FS_tau_pt, FS_tau_eta, FS_tau_phi, FS_tau_dxy, FS_tau_dz, FS_tau_chg, FS_tau_DM = [], [], [], [], [], [], []
  #FS_tau_PNet_v_jet, FS_tau_PNet_v_mu, FS_tau_PNet_v_e = [], [], []
  # note these are in the same order as the variables in the first line of this function :)
  # goes after l1_idx, l2_idx,
      #PNetvJet, PNetvMu, PNetvE,\
  for i, lep_pt, lep_eta, lep_phi, lep_iso,\
      mu_dxy, mu_dz, mu_chg, tau_dxy, tau_dz, tau_chg, tau_decayMode,\
      MET_pt, MET_phi, tau_idx, mu_idx,\
      l1_idx, l2_idx,\
      vJet, vMu, vEle, trg24mu, trg27mu, crosstrg in zip(*to_check):

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
    muIsoVal   = lep_iso[muLoc]
    muDxyVal   = abs(mu_dxy[muBranchLoc])
    muDzVal    = abs(mu_dz[muBranchLoc])
    muChgVal   = mu_chg[muBranchLoc]
    tauPtVal   = lep_pt[tauFSLoc] 
    tauEtaVal  = lep_eta[tauFSLoc]
    tauPhiVal  = lep_phi[tauFSLoc]
    tauDxyVal  = abs(tau_dxy[tauBranchLoc])
    tauDzVal   = abs(tau_dz[tauBranchLoc])
    tauChgVal  = tau_chg[tauBranchLoc]
    mtVal      = calculate_mt(muPtVal, muPhiVal, MET_pt, MET_phi)
    #passMT     = (mtVal < 65.0) # SF groups
    passMT     = (mtVal < 50.0) # mine
    #passMT     = True # passthrough for mt check

    #tauPNetvJetVal = PNetvJet[tauBranchLoc]
    #tauPNetvMuVal  = PNetvMu[tauBranchLoc]
    #tauPNetvEVal   = PNetvE[tauBranchLoc]

    # my selection
    passTauPtAndEta  = ((tauPtVal > 30.0) and (abs(tauEtaVal) < 2.3)) # mine
    pass25MuPt   = ((trg24mu) and (muPtVal > 25.0) and (abs(muEtaVal) < 2.4)) #mine
    pass28MuPt   = ((trg27mu) and (muPtVal > 28.0) and (abs(muEtaVal) < 2.4)) #mine
    # HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1
    passMuPtCrossTrigger = ((crosstrg) and ((21.0 < muPtVal < 25.0) and (abs(muEtaVal) < 2.1))
                                       and ((tauPtVal > 32.0)       and (abs(tauEtaVal) < 2.1)) ) 
    passMuPtCrossTrigger = False # dummy to turn off crosstrg

    # Medium v Jet, Tight v Muon, VVVLoose v Ele
    passTauDT  = ((vJet[tauBranchLoc] >= 5) and (vMu[tauBranchLoc] >= 4) and (vEle[tauBranchLoc] >= 1))
    # end my selection

    # SF group's selection
    #passTauPtAndEta  = ((tauPtVal > 20.0) and (abs(tauEtaVal) < 2.5))
    #pass25MuPt    = (trg24mu or trg27mu) and (muPtVal > 25.0) and (abs(muEtaVal) < 2.4)
 
    # Medium v Jet, Tight v Muon, VVLoose v Ele
    #passTauDT  = ((vJet[tauBranchLoc] >= 5) and (vMu[tauBranchLoc] >= 4) and (vEle[tauBranchLoc] >= 2))
    #skip_DM2 = (tau_decayMode[tauBranchLoc] != 2)
    # end SF group's selection

    if skip_DeepTau: passTauDT = False; # makes OR below mutually exclusive # for ARregion

    #restrict_tau_decayMode = (tau_decayMode[tauBranchLoc] == 0)

    if ( (passMT and (passTauPtAndEta and (pass25MuPt or pass28MuPt or passMuPtCrossTrigger)) and passTauDT) # mine
      or (passMT and (passTauPtAndEta and (pass25MuPt or pass28MuPt or passMuPtCrossTrigger)) and skip_DeepTau) ):
    #if ( (passMT and (passTauPtAndEta and pass25MuPt and passTauDT) and skip_DM2) # SF groups
    #  or (passMT and (passTauPtAndEta and pass25MuPt and skip_DeepTau) and skip_DM2) ):
      pass_cuts.append(i)
      FS_mu_pt.append(muPtVal)
      FS_mu_eta.append(muEtaVal)
      FS_mu_phi.append(muPhiVal)
      FS_mu_iso.append(muIsoVal)
      FS_mu_dxy.append(muDxyVal)
      FS_mu_dz.append(muDzVal)
      FS_mu_chg.append(muChgVal)

      FS_tau_pt.append(tauPtVal)
      FS_tau_eta.append(tauEtaVal)
      FS_tau_phi.append(tauPhiVal)
      FS_tau_dxy.append(tauDxyVal)
      FS_tau_dz.append(tauDzVal)
      FS_tau_chg.append(tauChgVal)
      FS_tau_DM.append(tau_decayMode[tauBranchLoc])

      #FS_tau_PNet_v_jet.append(tauPNetvJetVal)
      #FS_tau_PNet_v_mu.append(tauPNetvMuVal)
      #FS_tau_PNet_v_e.append(tauPNetvEVal)

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
  event_dictionary["FS_tau_chg"] = np.array(FS_tau_chg)
  event_dictionary["FS_tau_DM"]  = np.array(FS_tau_DM)
  event_dictionary["FS_mt"]    = np.array(FS_mt)
  #event_dictionary["FS_tau_rawPNetVSjet"] = np.array(FS_tau_PNet_v_jet)
  #event_dictionary["FS_tau_rawPNetVSmu"]  = np.array(FS_tau_PNet_v_mu)
  #event_dictionary["FS_tau_rawPNetVSe"]   = np.array(FS_tau_PNet_v_e)
  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after mutau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary


def make_mutau_AR_cut(event_dictionary, DeepTau_version):
  unpack_mutau_AR_vars = ["event", "Lepton_tauIdx", "Lepton_muIdx", "Lepton_iso", "l1_indices", "l2_indices"]
  unpack_mutau_AR_vars = add_DeepTau_branches(unpack_mutau_AR_vars, DeepTau_version)
  unpack_mutau_AR_vars = (event_dictionary.get(key) for key in unpack_mutau_AR_vars)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_mutau_AR_vars]
  pass_AR_cuts = []
  for i, event, tau_idx, mu_idx, lep_iso, l1_idx, l2_idx, vJet, _, _ in zip(*to_check):
    # keep indices where tau fails and muon passes iso 
    mu_lep_idx = l1_idx if mu_idx[l1_idx] != -1 else l2_idx
    muon_iso = lep_iso[mu_lep_idx]
    tau_branchIdx  = tau_idx[l1_idx] + tau_idx[l2_idx] + 1
    if ((vJet[tau_branchIdx] < 5) and (muon_iso<0.15)):
      pass_AR_cuts.append(i)
  
  event_dictionary["pass_AR_cuts"] = np.array(pass_AR_cuts)
  return event_dictionary

def make_mutau_TnP_cut(event_dictionary, DeepTau_version, skip_DeepTau=False, requireProbe=False):
  '''
  Works similarly to 'make_ditau_cut'. 
  Notably, the mutau cuts are more complicated, but it is simple to 
  extend the existing methods as long as one can stomach the line breaks.
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_mutau = ["Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_iso",
                  "Muon_dxy", "Muon_dz", "Muon_charge", "Tau_dxy", "Tau_dz", "Tau_charge", "Tau_decayMode",
                  "PuppiMET_pt", "PuppiMET_phi", "HTT_m_vis",
                  "Lepton_tauIdx", "Lepton_muIdx", "l1_indices", "l2_indices"]
  #TODO add this "CleanJet_btagWP" (no effect in August skims since it was always 1)
  unpack_mutau = add_DeepTau_branches(unpack_mutau, DeepTau_version)
  unpack_mutau = add_trigger_branches(unpack_mutau, final_state_mode="mutau")
  unpack_mutau = (event_dictionary.get(key) for key in unpack_mutau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_mutau] # "*" unpacks a tuple
  pass_cuts, FS_mt = [], []
  FS_mu_pt, FS_mu_eta, FS_mu_phi, FS_mu_iso, FS_mu_dxy, FS_mu_dz, FS_mu_chg = [], [], [], [], [], [], []
  FS_tau_pt, FS_tau_eta, FS_tau_phi, FS_tau_dxy, FS_tau_dz, FS_tau_chg = [], [], [], [], [], []
  pass_tag, pass_probe = [], []
  # note these are in the same order as the variables in the first line of this function :)
  for i, lep_pt, lep_eta, lep_phi, lep_iso,\
      mu_dxy, mu_dz, mu_chg, tau_dxy, tau_dz, tau_chg, tau_decayMode,\
      MET_pt, MET_phi, mvis, tau_idx, mu_idx,\
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
    muIsoVal   = lep_iso[muLoc]
    muDxyVal   = abs(mu_dxy[muBranchLoc])
    muDzVal    = abs(mu_dz[muBranchLoc])
    muChgVal   = mu_chg[muBranchLoc]
    tauPtVal   = lep_pt[tauFSLoc] 
    tauEtaVal  = lep_eta[tauFSLoc]
    tauPhiVal  = lep_phi[tauFSLoc]
    tauDxyVal  = abs(tau_dxy[tauBranchLoc])
    tauDzVal   = abs(tau_dz[tauBranchLoc])
    tauChgVal  = tau_chg[tauBranchLoc]

    mtVal      = calculate_mt(muPtVal, muPhiVal, MET_pt, MET_phi)
    passMT     = (mtVal < 30.0) #(mtVal < 50.0) #mine

    zmasswindow = (40.0 < mvis < 80.0)

    passTag   = trg24mu  # HLT_IsoMu24
    passProbe = crosstrg # HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1

    passOfflineTau = ((tauPtVal > 20.0) and (abs(tauEtaVal) < 2.3))
    passOfflineMu  = (trg24mu and (muPtVal > 27.0) and (abs(muEtaVal < 2.3)))

    skip_DMs = ((tau_decayMode[tauBranchLoc] != 5) and (tau_decayMode[tauBranchLoc] != 6))

    passTauDT  = (vJet[tauBranchLoc] >= 5)

    # TODO : actually figure out what you want
    # numerator shouldn't have tau pt restrictions, but denom should
    # could duplicate data, run it through two different selections, then compare those
    # sounds fine
    if ((passMT and zmasswindow and passTag and (passOfflineTau and passOfflineMu) and skip_DMs and passTauDT) or
        (requireProbe and 
        (passMT and zmasswindow and passTag and passProbe and passOfflineMu and skip_DMs and passTauDT)) ):
      pass_cuts.append(i)
      FS_mu_pt.append(muPtVal)
      FS_mu_eta.append(muEtaVal)
      FS_mu_phi.append(muPhiVal)
      FS_mu_iso.append(muIsoVal)
      FS_mu_dxy.append(muDxyVal)
      FS_mu_dz.append(muDzVal)
      FS_mu_chg.append(muChgVal)

      FS_tau_pt.append(tauPtVal)
      FS_tau_eta.append(tauEtaVal)
      FS_tau_phi.append(tauPhiVal)
      FS_tau_dxy.append(tauDxyVal)
      FS_tau_dz.append(tauDzVal)
      FS_tau_chg.append(tauChgVal)

      FS_mt.append(mtVal)

      pass_tag.append(passTag)
      pass_probe.append(passProbe)

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
  event_dictionary["FS_tau_chg"] = np.array(FS_tau_chg)
  event_dictionary["FS_mt"]    = np.array(FS_mt)
  event_dictionary["pass_tag"]   = np.array(pass_tag)
  event_dictionary["pass_probe"] = np.array(pass_probe)
  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after mutau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary



