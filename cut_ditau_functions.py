import numpy as np

from branch_functions import add_trigger_branches, add_DeepTau_branches

def make_ditau_cut(event_dictionary, DeepTau_version, skip_DeepTau=False):
  '''
  Use a minimal set of branches to define selection criteria and identify events which pass.
  A separate function uses the generated branch "pass_cuts" to remove the info from the
  loaded samples.
  Note: the zip method in python is a row-scanner, so the for loop below looks like this
  Events | pt | eta | tau_idx
  ###########################
       1 | 27 | 0.5 | 1
       2 | 35 | 1.5 | 0
       3 | 40 | 2.1 | 0
  i.e. we see variables of events sequentially.
  With this info, we make a simple check and store relevant variables.
  Note: stored variable branches are appended to other functions so that cutting
  events works properly
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_ditau = ["Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_tauIdx", 
                  "Tau_dxy", "Tau_dz", "Tau_decayMode", "Tau_charge", "l1_indices", "l2_indices",
                  #"Tau_rawPNetVSjet", "Tau_rawPNetVSmu", "Tau_rawPNetVSe"
                  ]
  unpack_ditau = add_DeepTau_branches(unpack_ditau, DeepTau_version)
  unpack_ditau = add_trigger_branches(unpack_ditau, final_state_mode="ditau")
  unpack_ditau = (event_dictionary.get(key) for key in unpack_ditau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_ditau] # "*" unpacks a tuple
  pass_cuts = []
  FS_t1_pt, FS_t1_eta, FS_t1_phi, FS_t1_dxy, FS_t1_dz, FS_t1_chg, FS_t1_DM = [], [], [], [], [], [], []
  FS_t2_pt, FS_t2_eta, FS_t2_phi, FS_t2_dxy, FS_t2_dz, FS_t2_chg, FS_t2_DM = [], [], [], [], [], [], []
  #FS_t1_PNet_v_jet, FS_t1_PNet_v_mu, FS_t1_PNet_v_e = [], [], []
  #FS_t2_PNet_v_jet, FS_t2_PNet_v_mu, FS_t2_PNet_v_e = [], [], []
  FS_t1_DeepTau_v_jet, FS_t1_DeepTau_v_mu, FS_t1_DeepTau_v_e = [], [], []
  FS_t2_DeepTau_v_jet, FS_t2_DeepTau_v_mu, FS_t2_DeepTau_v_e = [], [], []
  # note these are in the same order as the variables in the first line of this function :)
  # TODO: double-check counts with/without trigger :)
      #PNetvJet, PNetvMu, PNetvE, vJet, vMu, vEle,\
  for i, lep_pt, lep_eta, lep_phi, tau_idx,\
      tau_dxy, tau_dz, tau_decayMode, tau_chg, l1_idx, l2_idx,\
      vJet, vMu, vEle,\
      ditau_trig, ditau_jet_low_trig, ditau_jet_high_trig,\
      ditau_VBFRun2_trig, ditau_VBFRun3_trig in zip(*to_check):
    # assign object pts and etas
    t1_pt  = lep_pt[l1_idx]
    t2_pt  = lep_pt[l2_idx]
    t1_eta = lep_eta[l1_idx]
    t2_eta = lep_eta[l2_idx]
    # need nCleanJets variable to know what to fill
    # also, can happen that trigger with jets fires but nCleanJets says those are bad jets
    # would be good info to log/study/dump to terminal one time
    #j1_pt  = ???
    #j2_pt  = ???
    #j1_eta = ???
    #j2_eta = ???
    triggers = [ditau_trig, ditau_jet_low_trig, ditau_jet_high_trig,\
                ditau_VBFRun2_trig, ditau_VBFRun3_trig]
    #triggers = [ditau_trig, 0, 0, ditau_VBFRun2_trig, 0] # testing
    passKinems = pass_kinems_by_trigger(triggers, t1_pt, t2_pt, t1_eta, t2_eta)
         # jet cuts need to be made before FS cuts
    #passKinems = (lep_pt[l1_idx] >= 40.0 and lep_pt[l2_idx] >= 40.0)
    # Medium v Jet, VLoose v Muon, VVVLoose v Ele
    t1passDT   = (vJet[tau_idx[l1_idx]] >= 5 and vMu[tau_idx[l1_idx]] >= 1 and vEle[tau_idx[l1_idx]] >= 2)
    t2passDT   = (vJet[tau_idx[l2_idx]] >= 5 and vMu[tau_idx[l2_idx]] >= 1 and vEle[tau_idx[l2_idx]] >= 2)
    t1_decayMode = tau_decayMode[tau_idx[l1_idx]]
    t2_decayMode = tau_decayMode[tau_idx[l2_idx]]
    #good_tau_decayMode = ((t1_decayMode == 11) and (t2_decayMode == 11))
    good_tau_decayMode = True

    t1_chg = tau_chg[tau_idx[l1_idx]]
    t2_chg = tau_chg[tau_idx[l2_idx]]
    pass_charge = False
    if (t1_chg*t2_chg == -1): pass_charge = True
  
    if ( (pass_charge and passKinems and t1passDT and t2passDT and good_tau_decayMode) or 
         (passKinems and skip_DeepTau and t2passDT and good_tau_decayMode) ):
      pass_cuts.append(i)
      # TODO update me with the variable names used earlier
      FS_t1_pt.append(lep_pt[l1_idx])
      FS_t1_eta.append(lep_eta[l1_idx])
      FS_t1_phi.append(lep_phi[l1_idx])
      FS_t1_dxy.append(abs(tau_dxy[tau_idx[l1_idx]]))
      FS_t1_dz.append(abs(tau_dz[tau_idx[l1_idx]]))
      FS_t1_chg.append(t1_chg)
      FS_t1_DM.append(t1_decayMode)
      #FS_t1_PNet_v_jet.append(PNetvJet[tau_idx[l1_idx]])
      #FS_t1_PNet_v_mu.append(PNetvMu[tau_idx[l1_idx]])
      #FS_t1_PNet_v_e.append(PNetvE[tau_idx[l1_idx]])
      FS_t1_DeepTau_v_jet.append(vJet[tau_idx[l1_idx]])
      FS_t1_DeepTau_v_mu.append(vMu[tau_idx[l1_idx]])
      FS_t1_DeepTau_v_e.append(vEle[tau_idx[l1_idx]])
      FS_t2_pt.append(lep_pt[l2_idx])
      FS_t2_eta.append(lep_eta[l2_idx])
      FS_t2_phi.append(lep_phi[l2_idx])
      FS_t2_dxy.append(abs(tau_dxy[tau_idx[l2_idx]]))
      FS_t2_dz.append(abs(tau_dz[tau_idx[l2_idx]]))
      FS_t2_chg.append(t2_chg)
      FS_t2_DM.append(t2_decayMode)
      #FS_t2_PNet_v_jet.append(PNetvJet[tau_idx[l2_idx]])
      #FS_t2_PNet_v_mu.append(PNetvMu[tau_idx[l2_idx]])
      #FS_t2_PNet_v_e.append(PNetvE[tau_idx[l2_idx]])
      FS_t2_DeepTau_v_jet.append(vJet[tau_idx[l2_idx]])
      FS_t2_DeepTau_v_mu.append(vMu[tau_idx[l2_idx]])
      FS_t2_DeepTau_v_e.append(vEle[tau_idx[l2_idx]])

  event_dictionary["pass_cuts"] = np.array(pass_cuts)
  event_dictionary["FS_t1_pt"]  = np.array(FS_t1_pt)
  event_dictionary["FS_t1_eta"] = np.array(FS_t1_eta)
  event_dictionary["FS_t1_phi"] = np.array(FS_t1_phi)
  event_dictionary["FS_t1_dxy"] = np.array(FS_t1_dxy)
  event_dictionary["FS_t1_dz"]  = np.array(FS_t1_dz)
  event_dictionary["FS_t1_chg"] = np.array(FS_t1_chg)
  event_dictionary["FS_t1_DM"] = np.array(FS_t1_DM)
  #event_dictionary["FS_t1_rawPNetVSjet"] = np.array(FS_t1_PNet_v_jet)
  #event_dictionary["FS_t1_rawPNetVSmu"]  = np.array(FS_t1_PNet_v_mu)
  #event_dictionary["FS_t1_rawPNetVSe"]   = np.array(FS_t1_PNet_v_e)
  event_dictionary["FS_t1_DeepTauVSjet"] = np.array(FS_t1_DeepTau_v_jet)
  event_dictionary["FS_t1_DeepTauVSmu"]  = np.array(FS_t1_DeepTau_v_mu)
  event_dictionary["FS_t1_DeepTauVSe"]   = np.array(FS_t1_DeepTau_v_e)
  event_dictionary["FS_t2_pt"]  = np.array(FS_t2_pt)
  event_dictionary["FS_t2_eta"] = np.array(FS_t2_eta)
  event_dictionary["FS_t2_phi"] = np.array(FS_t2_phi)
  event_dictionary["FS_t2_dxy"] = np.array(FS_t2_dxy)
  event_dictionary["FS_t2_dz"]  = np.array(FS_t2_dz)
  event_dictionary["FS_t2_chg"] = np.array(FS_t2_chg)
  event_dictionary["FS_t2_DM"] = np.array(FS_t2_DM)
  #event_dictionary["FS_t2_rawPNetVSjet"] = np.array(FS_t2_PNet_v_jet)
  #event_dictionary["FS_t2_rawPNetVSmu"]  = np.array(FS_t2_PNet_v_mu)
  #event_dictionary["FS_t2_rawPNetVSe"]   = np.array(FS_t2_PNet_v_e)
  event_dictionary["FS_t2_DeepTauVSjet"] = np.array(FS_t2_DeepTau_v_jet)
  event_dictionary["FS_t2_DeepTauVSmu"]  = np.array(FS_t2_DeepTau_v_mu)
  event_dictionary["FS_t2_DeepTauVSe"]   = np.array(FS_t2_DeepTau_v_e)

  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after ditau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary


def make_ditau_cut_FF(event_dictionary, DeepTau_version):
  '''
  Use a minimal set of branches to define selection criteria and identify events which pass.
  A separate function uses the generated branch "pass_cuts" to remove the info from the
  loaded samples.
  Note: the zip method in python is a row-scanner, so the for loop below looks like this
  Events | pt | eta | tau_idx
  ###########################
       1 | 27 | 0.5 | 1
       2 | 35 | 1.5 | 0
       3 | 40 | 2.1 | 0
  i.e. we see variables of events sequentially.
  With this info, we make a simple check and store relevant variables.
  Note: stored variable branches are appended to other functions so that cutting
  events works properly
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_ditau = ["Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_tauIdx", 
                  "Tau_dxy", "Tau_dz", "Tau_decayMode", "Tau_charge", "l1_indices", "l2_indices",
                  #"Tau_rawPNetVSjet", "Tau_rawPNetVSmu", "Tau_rawPNetVSe"
                  ]
  unpack_ditau = add_DeepTau_branches(unpack_ditau, DeepTau_version)
  unpack_ditau = add_trigger_branches(unpack_ditau, final_state_mode="ditau")
  unpack_ditau = (event_dictionary.get(key) for key in unpack_ditau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_ditau] # "*" unpacks a tuple
  pass_cuts = []
  FS_t1_pt, FS_t1_eta, FS_t1_phi, FS_t1_dxy, FS_t1_dz, FS_t1_chg, FS_t1_DM = [], [], [], [], [], [], []
  FS_t2_pt, FS_t2_eta, FS_t2_phi, FS_t2_dxy, FS_t2_dz, FS_t2_chg, FS_t2_DM = [], [], [], [], [], [], []
  #FS_t1_PNet_v_jet, FS_t1_PNet_v_mu, FS_t1_PNet_v_e = [], [], []
  #FS_t2_PNet_v_jet, FS_t2_PNet_v_mu, FS_t2_PNet_v_e = [], [], []
  FS_t1_DeepTau_v_jet, FS_t1_DeepTau_v_mu, FS_t1_DeepTau_v_e = [], [], []
  FS_t2_DeepTau_v_jet, FS_t2_DeepTau_v_mu, FS_t2_DeepTau_v_e = [], [], []
  # note these are in the same order as the variables in the first line of this function :)
  # TODO: double-check counts with/without trigger :)
      #PNetvJet, PNetvMu, PNetvE, vJet, vMu, vEle,\
  for i, lep_pt, lep_eta, lep_phi, tau_idx,\
      tau_dxy, tau_dz, tau_decayMode, tau_chg, l1_idx, l2_idx,\
      vJet, vMu, vEle,\
      ditau_trig, ditau_jet_low_trig, ditau_jet_high_trig,\
      ditau_VBFRun2_trig, ditau_VBFRun3_trig in zip(*to_check):
    # assign object pts and etas
    t1_pt  = lep_pt[l1_idx]
    t2_pt  = lep_pt[l2_idx]
    t1_eta = lep_eta[l1_idx]
    t2_eta = lep_eta[l2_idx]
    # need nCleanJets variable to know what to fill
    # also, can happen that trigger with jets fires but nCleanJets says those are bad jets
    # would be good info to log/study/dump to terminal one time
    #j1_pt  = ???
    #j2_pt  = ???
    #j1_eta = ???
    #j2_eta = ???
    triggers = [ditau_trig, ditau_jet_low_trig, ditau_jet_high_trig,\
                ditau_VBFRun2_trig, ditau_VBFRun3_trig]
    #triggers = [ditau_trig, 0, 0, ditau_VBFRun2_trig, 0] # testing
    passKinems = pass_kinems_by_trigger(triggers, t1_pt, t2_pt, t1_eta, t2_eta)
         # jet cuts need to be made before FS cuts
    #passKinems = (lep_pt[l1_idx] >= 40.0 and lep_pt[l2_idx] >= 40.0)
    # Medium v Jet, VLoose v Muon, VVVLoose v Ele
    t1_decayMode = tau_decayMode[tau_idx[l1_idx]]
    t2_decayMode = tau_decayMode[tau_idx[l2_idx]]
    #good_tau_decayMode = ((t1_decayMode == 11) and (t2_decayMode == 11))
    good_tau_decayMode = True

    t1passDT   = (vMu[tau_idx[l1_idx]] >= 1 and vEle[tau_idx[l1_idx]] >= 2)
    t2passDT   = (vMu[tau_idx[l2_idx]] >= 1 and vEle[tau_idx[l2_idx]] >= 2)

    t1_chg = tau_chg[tau_idx[l1_idx]]
    t2_chg = tau_chg[tau_idx[l2_idx]]
  
    if ( (passKinems and good_tau_decayMode and t1passDT and t2passDT) ):
      pass_cuts.append(i)
      # TODO update me with the variable names used earlier
      FS_t1_pt.append(lep_pt[l1_idx])
      FS_t1_eta.append(lep_eta[l1_idx])
      FS_t1_phi.append(lep_phi[l1_idx])
      FS_t1_dxy.append(abs(tau_dxy[tau_idx[l1_idx]]))
      FS_t1_dz.append(abs(tau_dz[tau_idx[l1_idx]]))
      FS_t1_chg.append(t1_chg)
      FS_t1_DM.append(t1_decayMode)
      #FS_t1_PNet_v_jet.append(PNetvJet[tau_idx[l1_idx]])
      #FS_t1_PNet_v_mu.append(PNetvMu[tau_idx[l1_idx]])
      #FS_t1_PNet_v_e.append(PNetvE[tau_idx[l1_idx]])
      FS_t1_DeepTau_v_jet.append(vJet[tau_idx[l1_idx]])
      FS_t1_DeepTau_v_mu.append(vMu[tau_idx[l1_idx]])
      FS_t1_DeepTau_v_e.append(vEle[tau_idx[l1_idx]])
      FS_t2_pt.append(lep_pt[l2_idx])
      FS_t2_eta.append(lep_eta[l2_idx])
      FS_t2_phi.append(lep_phi[l2_idx])
      FS_t2_dxy.append(abs(tau_dxy[tau_idx[l2_idx]]))
      FS_t2_dz.append(abs(tau_dz[tau_idx[l2_idx]]))
      FS_t2_chg.append(t2_chg)
      FS_t2_DM.append(t2_decayMode)
      #FS_t2_PNet_v_jet.append(PNetvJet[tau_idx[l2_idx]])
      #FS_t2_PNet_v_mu.append(PNetvMu[tau_idx[l2_idx]])
      #FS_t2_PNet_v_e.append(PNetvE[tau_idx[l2_idx]])
      FS_t2_DeepTau_v_jet.append(vJet[tau_idx[l2_idx]])
      FS_t2_DeepTau_v_mu.append(vMu[tau_idx[l2_idx]])
      FS_t2_DeepTau_v_e.append(vEle[tau_idx[l2_idx]])

  event_dictionary["pass_cuts"] = np.array(pass_cuts)
  event_dictionary["FS_t1_pt"]  = np.array(FS_t1_pt)
  event_dictionary["FS_t1_eta"] = np.array(FS_t1_eta)
  event_dictionary["FS_t1_phi"] = np.array(FS_t1_phi)
  event_dictionary["FS_t1_dxy"] = np.array(FS_t1_dxy)
  event_dictionary["FS_t1_dz"]  = np.array(FS_t1_dz)
  event_dictionary["FS_t1_chg"] = np.array(FS_t1_chg)
  event_dictionary["FS_t1_DM"] = np.array(FS_t1_DM)
  #event_dictionary["FS_t1_rawPNetVSjet"] = np.array(FS_t1_PNet_v_jet)
  #event_dictionary["FS_t1_rawPNetVSmu"]  = np.array(FS_t1_PNet_v_mu)
  #event_dictionary["FS_t1_rawPNetVSe"]   = np.array(FS_t1_PNet_v_e)
  event_dictionary["FS_t1_DeepTauVSjet"] = np.array(FS_t1_DeepTau_v_jet)
  event_dictionary["FS_t1_DeepTauVSmu"]  = np.array(FS_t1_DeepTau_v_mu)
  event_dictionary["FS_t1_DeepTauVSe"]   = np.array(FS_t1_DeepTau_v_e)
  event_dictionary["FS_t2_pt"]  = np.array(FS_t2_pt)
  event_dictionary["FS_t2_eta"] = np.array(FS_t2_eta)
  event_dictionary["FS_t2_phi"] = np.array(FS_t2_phi)
  event_dictionary["FS_t2_dxy"] = np.array(FS_t2_dxy)
  event_dictionary["FS_t2_dz"]  = np.array(FS_t2_dz)
  event_dictionary["FS_t2_chg"] = np.array(FS_t2_chg)
  event_dictionary["FS_t2_DM"] = np.array(FS_t2_DM)
  #event_dictionary["FS_t2_rawPNetVSjet"] = np.array(FS_t2_PNet_v_jet)
  #event_dictionary["FS_t2_rawPNetVSmu"]  = np.array(FS_t2_PNet_v_mu)
  #event_dictionary["FS_t2_rawPNetVSe"]   = np.array(FS_t2_PNet_v_e)
  event_dictionary["FS_t2_DeepTauVSjet"] = np.array(FS_t2_DeepTau_v_jet)
  event_dictionary["FS_t2_DeepTauVSmu"]  = np.array(FS_t2_DeepTau_v_mu)
  event_dictionary["FS_t2_DeepTauVSe"]   = np.array(FS_t2_DeepTau_v_e)

  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after ditau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary



def pass_kinems_by_trigger(triggers, t1_pt, t2_pt, t1_eta, t2_eta, 
                           j1_pt=-999, j2_pt=-999, j1_eta=-999, j2_eta=-999, mjj=-999):
  '''
  Helper function to apply different object kinematic criteria depending on trigger used
  '''
  # first determine loosest trigger passed
  # ditau > ditau+onejet > ditau+Run3 > ditau+Run2 (last two unsure of ordering)
  # then apply cuts to offline objects based on loosest trigger
  # -- will show ditau+onejet contributing to >1 jet topologies
  # -- okay. ditau shows this behavior anyways and it is expected
  # -- this is to Konstantin's point about gain in VBF specifically,
  # -- i.e. is it justified with those other triggers in place?
  passTrig = False
  passTauKinems = False
  passJetKinems = False
  ditau_trig, ditau_jet_low_trig, ditau_jet_high_trig, ditau_VBFRun2_trig, ditau_VBFRun3_trig = triggers
  if ditau_trig:
    passTrig = True # could put trigger name here, then sort later by what passes
    passTauKinems = (t1_pt > 40 and t2_pt > 40 and abs(t1_eta) < 2.1 and abs(t2_eta) < 2.1)
    passJetKinems = True # 30 30 600 for 2jet, but will be global cut earlier and just checked here, unsure global 1jet
  if (ditau_jet_low_trig or ditau_jet_high_trig) and not (ditau_trig):
    passTrig = True
    passTauKinems = (t1_pt > 35 and t2_pt > 35 and abs(t1_eta) < 2.1 and abs(t2_eta) < 2.1)
    #passJetKinems = (j1_pt > 65 and abs(j1_eta) < 4.7) # 2jet 30 600 # TODO: assuming higher trigger wasn't used, check
    passJetKinems = True
  if (ditau_VBFRun2_trig) and not (ditau_trig or ditau_jet_low_trig or ditau_jet_high_trig or ditau_VBFRun3_trig):
    passTrig = True
    passTauKinems = (t1_pt > 25 and t2_pt > 25 and abs(t1_eta) < 2.1 and abs(t2_eta) < 2.1)
    #passJetKinems = (j1_pt > 120 and j2_pt > 40 and mjj > 700 and abs(j1_pt) < 4.7 and abs(j2_pt) < 4.7)
    passJetKinems = True
  if (ditau_VBFRun3_trig) and not (ditau_trig or ditau_jet_low_trig or ditau_jet_high_trig):
    passTrig = True
    passTauKinems = (t1_pt > 50 and t2_pt > 25 and abs(t1_eta) < 2.1 and abs(t2_eta) < 2.1)
    #passJetKinems = (j1_pt > 45 and j2_pt > 45 and mjj > 600 and abs(j1_pt) < 4.7 and abs(j2_pt) < 4.7)
    passJetKinems = True

  passKinems = (passTrig and passTauKinems and passJetKinems)
  return passKinems

#TODO make into one function
def make_ditau_AR_cut(event_dictionary, DeepTau_version):
  unpack_ditau_AR_vars = ["Lepton_tauIdx", "l1_indices", "l2_indices"]
  unpack_ditau_AR_vars = add_DeepTau_branches(unpack_ditau_AR_vars, DeepTau_version)
  unpack_ditau_AR_vars = (event_dictionary.get(key) for key in unpack_ditau_AR_vars)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_ditau_AR_vars]
  pass_AR_cuts = []
  for i, tau_idx, l1_idx, l2_idx, vJet, _, _ in zip(*to_check):
    # keep indices where first tau fails and 2nd passes
    if (vJet[tau_idx[l1_idx]] < 5) and (vJet[tau_idx[l2_idx]] >= 5):
      pass_AR_cuts.append(i)
  
  event_dictionary["pass_AR_cuts"] = np.array(pass_AR_cuts)
  return event_dictionary

def make_ditau_SR_cut(event_dictionary, DeepTau_version):
  unpack_ditau_SR_vars = ["Lepton_tauIdx", "l1_indices", "l2_indices"]
  unpack_ditau_SR_vars = add_DeepTau_branches(unpack_ditau_SR_vars, DeepTau_version)
  unpack_ditau_SR_vars = (event_dictionary.get(key) for key in unpack_ditau_SR_vars)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_ditau_SR_vars]
  pass_SR_cuts = []
  for i, tau_idx, l1_idx, l2_idx, vJet, _, _ in zip(*to_check):
    # keep indices where first tau fails and 2nd passes
    if (vJet[tau_idx[l1_idx]] >= 5) and (vJet[tau_idx[l2_idx]] >= 5):
      pass_SR_cuts.append(i)
  
  event_dictionary["pass_SR_cuts"] = np.array(pass_SR_cuts)
  return event_dictionary

def make_ditau_AR_aiso_cut(event_dictionary, DeepTau_version):
  unpack_ditau_AR_aiso_vars = ["Lepton_tauIdx", "l1_indices", "l2_indices"]
  unpack_ditau_AR_aiso_vars = add_DeepTau_branches(unpack_ditau_AR_aiso_vars, DeepTau_version)
  unpack_ditau_AR_aiso_vars = (event_dictionary.get(key) for key in unpack_ditau_AR_aiso_vars)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_ditau_AR_aiso_vars]
  pass_AR_aiso_cuts = []
  for i, tau_idx, l1_idx, l2_idx, vJet, _, _ in zip(*to_check):
    # keep indices where first tau fails and 2nd passes
    if (vJet[tau_idx[l1_idx]] < 5) and (vJet[tau_idx[l2_idx]] < 5):
      pass_AR_aiso_cuts.append(i)
  
  event_dictionary["pass_AR_aiso_cuts"] = np.array(pass_AR_aiso_cuts)
  return event_dictionary

def make_ditau_SR_aiso_cut(event_dictionary, DeepTau_version):
  unpack_ditau_SR_aiso_vars = ["Lepton_tauIdx", "l1_indices", "l2_indices"]
  unpack_ditau_SR_aiso_vars = add_DeepTau_branches(unpack_ditau_SR_aiso_vars, DeepTau_version)
  unpack_ditau_SR_aiso_vars = (event_dictionary.get(key) for key in unpack_ditau_SR_aiso_vars)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_ditau_SR_aiso_vars]
  pass_SR_aiso_cuts = []
  for i, tau_idx, l1_idx, l2_idx, vJet, _, _ in zip(*to_check):
    # keep indices where first tau fails and 2nd passes
    if (vJet[tau_idx[l1_idx]] >= 5) and (vJet[tau_idx[l2_idx]] < 5):
      pass_SR_aiso_cuts.append(i)
  
  event_dictionary["pass_SR_aiso_cuts"] = np.array(pass_SR_aiso_cuts)
  return event_dictionary




