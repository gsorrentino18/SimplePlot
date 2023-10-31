import numpy as np

### README
# this file contains functions to perform cuts and self-contained studies

from calculate_functions   import calculate_mt
from triggers_dictionary import triggers_dictionary
from utility_functions import time_print


def append_lepton_indices(event_dictionary):
  '''
  Read the entries of "FSLeptons" and extract the values to place in separate branches.
  It was easier to do this once when the data is first loaded than to do it every time
  that it is needed. 
  '''
  FSLeptons = event_dictionary["FSLeptons"]
  l1_indices, l2_indices = [], []
  for event in FSLeptons:
    if len(event)>2: print(f"More than one FS pair: {event}")
    l1_indices.append(event[0])
    l2_indices.append(event[1])
  event_dictionary["l1_indices"] = np.array(l1_indices)
  event_dictionary["l2_indices"] = np.array(l2_indices)
  return event_dictionary


def make_ditau_cut(event_dictionary, DeepTau_version):
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
  unpack_ditau = ["Lepton_pt", "Lepton_eta", "Lepton_tauIdx", "l1_indices", "l2_indices"]
  unpack_ditau = add_DeepTau_branches(unpack_ditau, DeepTau_version)
  unpack_ditau = (event_dictionary.get(key) for key in unpack_ditau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_ditau] # "*" unpacks a tuple
  pass_cuts, FS_t1_pt, FS_t2_pt, FS_t1_eta, FS_t2_eta = [], [], [], [], []
  # note these are in the same order as the variables in the first line of this function :)
  # TODO: double-check counts with/without trigger :)
  for i, lep_pt, lep_eta, tau_idx, l1_idx, l2_idx, vJet, vMu, vEle in zip(*to_check):
    passKinems = (lep_pt[l1_idx] >= 40.0 and lep_pt[l2_idx] >= 40.0)
    # Medium v Jet, VLoose v Muon, VVVLoose v Ele
    t1passDT   = (vJet[tau_idx[l1_idx]] >= 5 and vMu[tau_idx[l1_idx]] >= 1 and vEle[tau_idx[l1_idx]] >= 1)
    t2passDT   = (vJet[tau_idx[l2_idx]] >= 5 and vMu[tau_idx[l2_idx]] >= 1 and vEle[tau_idx[l2_idx]] >= 1)
    if (passKinems and t1passDT and t2passDT):
      pass_cuts.append(i)
      FS_t1_pt.append(lep_pt[l1_idx])
      FS_t2_pt.append(lep_pt[l2_idx])
      FS_t1_eta.append(lep_eta[l1_idx])
      FS_t2_eta.append(lep_eta[l2_idx])

  event_dictionary["pass_cuts"] = np.array(pass_cuts)
  event_dictionary["FS_t1_pt"]  = np.array(FS_t1_pt)
  event_dictionary["FS_t2_pt"]  = np.array(FS_t2_pt)
  event_dictionary["FS_t1_eta"] = np.array(FS_t1_eta)
  event_dictionary["FS_t2_eta"] = np.array(FS_t2_eta)
  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after ditau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary


def make_mutau_cut(event_dictionary, DeepTau_version):
  '''
  Works similarly to 'make_ditau_cut'. 
  Notably, the mutau cuts are more complicated, but it is simple to 
  extend the existing methods as long as one can stomach the line breaks.
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_mutau = ["Tau_pt", "Tau_eta", 
                  "Muon_pt", "Muon_eta", "Muon_phi", 
                  "PuppiMET_pt", "PuppiMET_phi",
                  "Lepton_tauIdx", "Lepton_muIdx", "l1_indices", "l2_indices"]
  #TODO add this "CleanJet_btagWP" (no effect in August skims since it was always 1)
  unpack_mutau = add_DeepTau_branches(unpack_mutau, DeepTau_version)
  unpack_mutau = add_trigger_branches(unpack_mutau, final_state_mode="mutau")
  unpack_mutau = (event_dictionary.get(key) for key in unpack_mutau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_mutau] # "*" unpacks a tuple
  pass_cuts, FS_mu_pt, FS_tau_pt, FS_mu_eta, FS_tau_eta, FS_mt = [], [], [], [], [], []
  # note these are in the same order as the variables in the first line of this function :)
  for i, tau_pt, tau_eta, mu_pt, mu_eta, mu_phi, MET_pt, MET_phi, tau_idx, mu_idx,\
      l1_idx, l2_idx, vJet, vMu, vEle, trg24mu, trg27mu, crosstrg in zip(*to_check):

    tauLoc     = tau_idx[l1_idx] + tau_idx[l2_idx] + 1
    muLoc      = mu_idx[l1_idx]  + mu_idx[l2_idx]  + 1
    tauEtaVal  = tau_eta[tauLoc]
    tauPtVal   = tau_pt[tauLoc] 
    muPtVal    = mu_pt[muLoc] 
    muEtaVal   = mu_eta[muLoc]
    muPhiVal   = mu_phi[muLoc]
    mtVal      = calculate_mt(muPtVal, muPhiVal, MET_pt, MET_phi)
    passMT     = (mtVal < 50.0)
    #ROOTmtVal  = calculate_mt_pyROOT(muPtVal, muEtaVal, muPhiVal, mu_M[muLoc], MET_pt, MET_phi)
    #passROOTMT = (ROOTmtVal < 50.0)

    passTauPtAndEta  = ((tauPtVal > 30.0) and (abs(tauEtaVal) < 2.3))
    pass25MuPt   = ((trg24mu) and (muPtVal > 25.0) and (abs(muEtaVal) < 2.4))
    pass28MuPt   = ((trg27mu) and (muPtVal > 28.0) and (abs(muEtaVal) < 2.4))
    # HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1
    passMuPtCrossTrigger = ((crosstrg) and ((21.0 < muPtVal < 25.0) and (abs(muEtaVal) < 2.1))
                                       and ((tauPtVal > 32.0)       and (abs(tauEtaVal) < 2.1)) ) 

    # Medium v Jet, Tight v Muon, VVVLoose v Ele
    passTauDT  = ((vJet[tauLoc] >= 5) and (vMu[tauLoc] >= 4) and (vEle[tauLoc] >= 1))

    if (passMT and (passTauPtAndEta and (pass25MuPt or pass28MuPt or passMuPtCrossTrigger)) and passTauDT):
      pass_cuts.append(i)
      FS_mu_pt.append(muPtVal)
      FS_tau_pt.append(tauPtVal)
      FS_mu_eta.append(muEtaVal)
      FS_tau_eta.append(tauEtaVal)
      FS_mt.append(mtVal)

  event_dictionary["pass_cuts"] = np.array(pass_cuts)
  event_dictionary["FS_mu_pt"]  = np.array(FS_mu_pt)
  event_dictionary["FS_tau_pt"] = np.array(FS_tau_pt)
  event_dictionary["FS_mu_eta"] = np.array(FS_mu_eta)
  event_dictionary["FS_tau_eta"] = np.array(FS_tau_eta)
  event_dictionary["FS_mt"]    = np.array(FS_mt)
  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after mutau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary


def make_etau_cut(event_dictionary, DeepTau_version):
  '''
  Works similarly to 'make_ditau_cut'. 
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  # TODO shouldn't i be using Lepton_pt/eta?
  unpack_etau = ["Tau_pt", "Tau_eta", 
                 "Electron_pt", "Electron_eta", "Electron_phi", 
                 "PuppiMET_pt", "PuppiMET_phi",
                 "Lepton_tauIdx", "Lepton_elIdx", "l1_indices", "l2_indices"]
  unpack_etau = add_DeepTau_branches(unpack_etau, DeepTau_version)
  unpack_etau = add_trigger_branches(unpack_etau, final_state_mode="etau")
  unpack_etau = (event_dictionary.get(key) for key in unpack_etau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_etau]
  pass_cuts, FS_el_pt, FS_tau_pt, FS_el_eta, FS_tau_eta, FS_mt = [], [], [], [], [], []
  for i, tau_pt, tau_eta, el_pt, el_eta, el_phi, MET_pt, MET_phi,\
      tau_idx, el_idx, l1_idx, l2_idx,\
      vJet, vMu, vEle,\
      trg32el, trg35el, crosstrg in zip(*to_check):

    tauLoc     = tau_idx[l1_idx] + tau_idx[l2_idx] + 1
    elLoc      = el_idx[l1_idx]  + el_idx[l2_idx]  + 1
    tauEtaVal  = tau_eta[tauLoc]
    tauPtVal   = tau_pt[tauLoc] 
    elPtVal    = el_pt[elLoc] 
    elEtaVal   = el_eta[elLoc]
    elPhiVal   = el_phi[elLoc]
    mtVal      = calculate_mt(elPtVal, elPhiVal, MET_pt, MET_phi)
    passMT     = (mtVal < 50.0)
    #ROOTmtVal  = calculate_mt_pyROOT(muPtVal, muEtaVal, muPhiVal, mu_M[muLoc], MET_pt, MET_phi)
    #passROOTMT = (ROOTmtVal < 50.0)

    passTauPtAndEta  = ((tauPtVal > 30.0) and (abs(tauEtaVal) < 2.3))
    pass33ElPt   = ((trg32el) and (elPtVal > 33.0) and (abs(elEtaVal) < 2.1))
    pass36ElPt   = ((trg35el) and (elPtVal > 36.0) and (abs(elEtaVal) < 2.1))
    # upper bound on cross trigger will change if lower single electron trigger included
    # HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1
    passElPtCrossTrigger = ((crosstrg) and ((25.0 < elPtVal < 33.0) and (abs(elEtaVal) < 2.1))
                                       and ((tauPtVal > 35.0)       and (abs(tauEtaVal) < 2.1)) ) 

    # Medium v Jet, VLoose v Muon, Tight v Ele
    passTauDT  = ((vJet[tauLoc] >= 5) and (vMu[tauLoc] >= 1) and (vEle[tauLoc] >= 6))

    if (passMT and passTauPtAndEta and (pass33ElPt or pass36ElPt or passElPtCrossTrigger) and passTauDT):
      pass_cuts.append(i)
      FS_el_pt.append(elPtVal)
      FS_tau_pt.append(tauPtVal)
      FS_el_eta.append(elEtaVal)
      FS_tau_eta.append(tauEtaVal)
      FS_mt.append(mtVal)

  event_dictionary["pass_cuts"]  = np.array(pass_cuts)
  event_dictionary["FS_el_pt"]   = np.array(FS_el_pt)
  event_dictionary["FS_tau_pt"]  = np.array(FS_tau_pt)
  event_dictionary["FS_el_eta"]  = np.array(FS_el_eta)
  event_dictionary["FS_tau_eta"] = np.array(FS_tau_eta)
  event_dictionary["FS_mt"]      = np.array(FS_mt)
  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after ditau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary


def make_jet_cut(event_dictionary):
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_jetVars = ["nCleanJet", "CleanJet_pt", "CleanJet_eta"]
  unpack_jetVars = (event_dictionary.get(key) for key in unpack_jetVars)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_jetVars] # "*" unpacks a tuple
  pass_zero_jet_cuts, pass_one_jet_cuts, pass_two_jet_cuts, pass_two_or_more_jet_cuts = [], [], [], []
  nCleanJetGT30 = []
  for i, nJet, jet_pt, jet_eta in zip(*to_check):
    passingJets = 0
    for ijet in range(0, nJet):
      if (jet_pt[ijet] > 30.0) and (jet_eta[ijet] < 4.7):
        passingJets += 1
    nCleanJetGT30.append(passingJets)

    if   passingJets == 0: pass_zero_jet_cuts.append(i)
    elif passingJets == 1: pass_one_jet_cuts.append(i)
    elif passingJets == 2: pass_two_jet_cuts.append(i)
    elif passingJets >= 2: pass_two_or_more_jet_cuts.append(i)

  event_dictionary["pass_zero_jet_cuts"]        = np.array(pass_zero_jet_cuts)
  event_dictionary["pass_one_jet_cuts"]         = np.array(pass_one_jet_cuts)
  event_dictionary["pass_two_jet_cuts"]         = np.array(pass_two_jet_cuts)
  event_dictionary["pass_two_or_more_jet_cuts"] = np.array(pass_two_or_more_jet_cuts)
  event_dictionary["nCleanJetGT30"] = np.array(nCleanJetGT30)
  #nEvents_postzerojetcut      = len(np.array(pass_zero_jet_cuts))
  #nEvents_postonejetcut       = len(np.array(pass_one_jet_cuts))
  #nEvents_posttwojetcut       = len(np.array(pass_two_jet_cuts))
  #nEvents_posttwoormorejetcut = len(np.array(pass_two_or_more_jet_cuts))
  print(f"nEvents before and after 0  jet cuts = {nEvents_precut}, {len(np.array(pass_zero_jet_cuts))}")
  print(f"nEvents before and after 1  jet cuts = {nEvents_precut}, {len(np.array(pass_one_jet_cuts))}")
  print(f"nEvents before and after 2  jet cuts = {nEvents_precut}, {len(np.array(pass_two_jet_cuts))}")
  print(f"nEvents before and after ≥2 jet cuts = {nEvents_precut}, {len(np.array(pass_two_or_more_jet_cuts))}")
  return event_dictionary


def manual_dimuon_lepton_veto(event_dictionary):
  '''
  Works similarly to 'make_ditau_cut' except the branch "pass_manual_lepton_veto"
  is made specifically for the dimuon final state. Some special handling is required
  due to the way events are selected in step2 of the NanoTauFramework
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_veto = ["Lepton_pdgId", "Lepton_iso"]
  unpack_veto = (event_dictionary.get(key) for key in unpack_veto)
  to_check    = [range(len(event_dictionary["Lepton_pt"])), *unpack_veto]
  pass_manual_lepton_veto = []
  for i, lep_pdgId_array, lep_iso_array in zip(*to_check):
    for pdgId, iso in zip(lep_pdgId_array, lep_iso_array):
      nIsoEle, nIsoMu = 0, 0
      if abs(pdgId) == 11:
        nIsoEle += 1 if (iso < 0.3) else 0
      elif abs(pdgId) == 13:
        nIsoMu  += 1 if (iso < 0.3) else 0
      if nIsoEle < 1 and nIsoMu < 3:
        pass_manual_lepton_veto.append(i)

  event_dictionary["pass_manual_lepton_veto"] = np.array(pass_manual_lepton_veto)
  print(f"events before and after manual dimuon lepton veto = {nEvents_precut}, {len(np.array(pass_manual_lepton_veto))}")
  return event_dictionary


def make_dimuon_cut(event_dictionary):
  '''
  Works similarly to 'make_ditau_cut'. 
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_dimuon = ["Lepton_pt", "Lepton_eta", "Lepton_iso", "HTT_m_vis", "HTT_dR", "l1_indices", "l2_indices"]
  unpack_dimuon = (event_dictionary.get(key) for key in unpack_dimuon)
  to_check      = [range(len(event_dictionary["Lepton_pt"])), *unpack_dimuon]
  pass_cuts, FS_m1_pt, FS_m2_pt, FS_m1_eta, FS_m2_eta = [], [], [], [], []
  for i, pt, eta, iso, mvis, dR, l1_idx, l2_idx in zip(*to_check):
    passKinematics = (pt[l1_idx] > 26 and pt[l2_idx] > 20 and mvis > 20 and dR > 0.5)
    passIso        = (iso[l1_idx] < 0.15 and iso[l2_idx] < 0.15)
    if (passKinematics and passIso):
      pass_cuts.append(i)
      FS_m1_pt.append(pt[l1_idx])
      FS_m2_pt.append(pt[l2_idx])
      FS_m1_eta.append(eta[l1_idx])
      FS_m2_eta.append(eta[l2_idx])

  event_dictionary["pass_cuts"] = np.array(pass_cuts)
  event_dictionary["FS_m1_pt"]  = np.array(FS_m1_pt)
  event_dictionary["FS_m2_pt"]  = np.array(FS_m2_pt)
  event_dictionary["FS_m1_eta"] = np.array(FS_m1_eta)
  event_dictionary["FS_m2_eta"] = np.array(FS_m2_eta)
  print(f"events before and after dimuon cuts = {nEvents_precut}, {len(np.array(pass_cuts))}")
  return event_dictionary


def make_run_cut(event_dictionary, good_runs):
  '''
  Given a set of runs, create a branch of events belonging to that set.
  The branch is later used to reject all other events.
  '''
  good_runs = np.sort(good_runs)
  first_run, last_run = good_runs[0], good_runs[-1]
  print(f"first run {first_run}, last run {last_run}")
  # check if it's within the range, then check if it's in the list
  pass_run_cut = []
  for i, run in enumerate(event_dictionary["run"]):
    if first_run <= run <= last_run:
      if run in good_runs:
        pass_run_cut.append(i) 

  event_dictionary["pass_run_cut"] = np.array(pass_run_cut)
  return event_dictionary


def apply_cut(event_dictionary, cut_branch, protected_branches=[]):
  '''
  Remove all entries in 'event_dictionary' not in 'cut_branch' using the numpy 'take' method.
  Branches that are added during previous cut steps are added here because their entries
  already pass cuts by construction.
  The returned event_dictionary now only contains events passing all cuts.

  If all events are removed by cut, print a message to alert the user.
  The deletion is actually handled in the main body when the size of the dictionary is checked.
  '''
  delete_sample = False
  if len(event_dictionary[cut_branch]) == 0:
    print("All events removed, sample deleted")
    delete_sample = True

  for branch in event_dictionary:
    if delete_sample:
      pass
    elif ((branch != cut_branch) and (branch not in protected_branches)):
      event_dictionary[branch] = np.take(event_dictionary[branch], event_dictionary[cut_branch])

  return event_dictionary


def Era_F_trigger_study(data_events, final_state_mode):
  '''
  Compact function for 2022 era F trigger study, where ChargedIsoTau
  triggers were briefly enabled for Run2-Run3 Tau trigger studies. 
  '''
  FS_triggers = triggers_dictionary[final_state_mode]
  for trigger in FS_triggers:
    print(f" {trigger} has {np.sum(data_events[trigger])} events")

  good_runs = [361971, 361989, 361990, 361994, 362058, 362059, 362060, 
               362061, 362062, 362063, 362064, 362087, 362091, 362104, 
               362105, 362106, 362107, 362148, 362153, 362154, 362159, 
               362161, 362163, 362166, 362167]
  data_events = make_run_cut(data_events, good_runs)
  data_events = apply_cut(data_events, "pass_run_cut") # will break if used

  print("after reducing run range")
  for trigger in FS_triggers:
    print(f" {trigger} has {np.sum(data_events[trigger])} events")
  
  return data_events


def study_triggers():
  '''
  Template function for returning ORs/ANDs of HLT triggers in an organized way.
  Will be extended at an opportune moment.
  '''
  Run2OR, Run2AND, Run3OR, Run3AND = 0, 0, 0, 0

  mutau_triggers = [data_events[trigger] for trigger in add_trigger_branches([], "mutau")]
  for HLT_single1, HLT_single2, HLT_crossRun2, HLT_crossRun3 in zip(*mutau_triggers):
    if HLT_single1 or HLT_single2 or HLT_crossRun2:
      Run2OR  += 1
    if HLT_single1 or HLT_single2 or HLT_crossRun3:
      Run3OR  += 1
    if (HLT_single1 or HLT_single2) and HLT_crossRun2:
      Run2AND += 1
    if (HLT_single1 or HLT_single2) and HLT_crossRun3:
      Run3AND += 1
 
  print(f"Run2 OR/AND: {Run2OR}\t{Run2AND}")
  print(f"Run3 OR/AND: {Run3OR}\t{Run3AND}")


def apply_final_state_cut(event_dictionary, final_state_mode, DeepTau_version):
  '''
  Organizational function that generalizes call to a (set of) cuts based on the
  final cut. Importantly, the function that rejects events, 'apply_cut',
  is called elsewhere
  '''
  protected_branches = set_protected_branches(final_state_mode=final_state_mode)
  if final_state_mode == "ditau":
    event_dictionary = make_ditau_cut(event_dictionary, DeepTau_version)
    event_dictionary = apply_cut(event_dictionary, "pass_cuts", protected_branches)
  elif final_state_mode == "mutau":
    event_dictionary = make_mutau_cut(event_dictionary, DeepTau_version)
    event_dictionary = apply_cut(event_dictionary, "pass_cuts", protected_branches)
  elif final_state_mode == "etau":
    event_dictionary = make_etau_cut(event_dictionary, DeepTau_version)
    event_dictionary = apply_cut(event_dictionary, "pass_cuts", protected_branches)
  elif final_state_mode == "dimuon":
    event_dictionary = manual_dimuon_lepton_veto(event_dictionary)
    event_dictionary = apply_cut(process_events, "pass_manual_lepton_veto")
    event_dictionary = make_dimuon_cut(event_dictionary)
    event_dictionary = apply_cut(event_dictionary, "pass_cuts")
  else:
    print(f"No cuts to apply for {final_state_mode} final state.")
  return event_dictionary


def apply_jet_cut(event_dictionary, jet_mode):
  '''
  Organizational function to reduce event_dictionary to contain only
  events with jets passing certain criteria. Enables plotting of jet objects
  jet_mode can be "pass_zero_jet_cut", "pass_one_jet_cut", "pass_two_jet_cut"
                 "pass_one_or_more_jet_cut", "pass_two_or_more_jet_cut",
                 "pass_inclusive_ptGT30_etaLT4p7_cut"
  '''
  event_dictionary   = make_jet_cut(event_dictionary)
  protected_branches = set_protected_branches(jet_mode=jet_mode)
  event_dictionary   = apply_cut(event_dictionary, jet_mode, protected_branches)
  return event_dictionary


def apply_cuts_to_process(process, process_dictionary, final_state_mode, DeepTau_version="2p5"):
  '''
  Organizational function to hold two function calls and empty list handling that
  is performed for all loaded datasets in our framework.
  Can be extended to hold additional standard cuts (i.e. jets) or the returned
  value can be cut on as needed.
  '''
  time_print(f"Processing {process}")
  process_events = process_dictionary[process]["info"]
  if len(process_events["run"])==0: return None

  process_events = append_lepton_indices(process_events)
  cut_events = apply_final_state_cut(process_events, final_state_mode, DeepTau_version)
  if len(cut_events["pass_cuts"])==0: return None 

  # add jet mode cut, default to no cut (inclusive jets)
  # also add proper branches to jet plotting

  return cut_events


def set_branches(final_state_mode, DeepTau_version="2p5"):
  common_branches = [
    "run", "luminosityBlock", "event", "Generator_weight",
    "FSLeptons", "Lepton_pt", "Lepton_eta",
    "nCleanJet", "CleanJet_pt", "CleanJet_eta",
    "HTT_m_vis", "HTT_dR",
  ]
  branches = common_branches
  branches = add_final_state_branches(branches, final_state_mode)
  branches = add_DeepTau_branches(branches, DeepTau_version)
  branches = add_trigger_branches(branches, final_state_mode)
  return branches


def add_final_state_branches(branches_, final_state_mode):
  '''
  Helper function to add only relevant branches to loaded branches based on final state.
  '''
  final_state_branches = {
    "ditau"  : ["Lepton_tauIdx"],

    "mutau"  : ["Tau_pt", "Tau_eta", 
                "Muon_pt", "Muon_eta", "Muon_phi",
                "Lepton_tauIdx", "Lepton_muIdx", 
                "PuppiMET_pt", "PuppiMET_phi"],

    "etau"   : ["Tau_pt", "Tau_eta",
                "Electron_pt", "Electron_eta", "Electron_phi",
                "Lepton_tauIdx", "Lepton_elIdx",
                "PuppiMET_pt", "PuppiMET_phi"],

    "dimuon" : ["Muon_pt", "Muon_eta",
                "HTT_m_vis", "HTT_dR",
                "Lepton_iso", "Lepton_pdgId", "Lepton_muIdx"],
  }

  branch_to_add = final_state_branches[final_state_mode]
  for new_branch in branch_to_add:
    branches_.append(new_branch)
  
  return branches_


def add_DeepTau_branches(branches_, DeepTauVersion):
  '''
  Helper function to add DeepTauID branches
  '''
  if DeepTauVersion == "2p1":
    for DeepTau_v2p1_branch in ["Tau_idDeepTau2017v2p1VSjet", "Tau_idDeepTau2017v2p1VSmu", "Tau_idDeepTau2017v2p1VSe"]:
      branches_.append(DeepTau_v2p1_branch)

  elif DeepTauVersion == "2p5":
    for DeepTau_v2p5_branch in ["Tau_idDeepTau2018v2p5VSjet", "Tau_idDeepTau2018v2p5VSmu", "Tau_idDeepTau2018v2p5VSe"]:
      branches_.append(DeepTau_v2p5_branch)

  else:
    print(f"no branches added with argument {DeepTauVersion}. Try 2p1 or 2p5.")

  return branches_


def add_trigger_branches(branches_, final_state_mode):
  '''
  Helper function to add HLT branches used by a given final state
  '''
  for trigger in triggers_dictionary[final_state_mode]:
    branches_.append(trigger)
  return branches_


def set_vars_to_plot(final_state_mode, jet_mode="none"):
  '''
  Helper function to keep plotting variables organized
  '''
  jet_plotting_vars = {
    "none": [],
    "pass_zero_jet_cuts"  : [],
    "pass_one_jet_cuts"   : ["CleanJetGT30_pt_1", "CleanJetGT30_eta_1"],

    "pass_two_jet_cuts"   : ["CleanJetGT30_pt_1", "CleanJetGT30_eta_1",
                             "CleanJetGT30_pt_2", "CleanJetGT30_eta_2"],

    "pass_two_or_more_jet_cuts" : [
      "CleanJetGT30_pt_1", "CleanJetGT30_eta_1",
      "CleanJetGT30_pt_2", "CleanJetGT30_eta_2",
      "CleanJetGT30_pt_3", "CleanJetGT30_eta_3"],
  }

  final_state_plotting_vars = {
    "none"   : [],
    "ditau"  : ["FS_t1_pt", "FS_t2_pt", "FS_t1_eta", "FS_t2_eta"],
    "mutau"  : ["FS_mu_pt", "FS_tau_pt", "FS_mu_eta", "FS_tau_eta",
                "FS_mt", "PuppiMET_pt"],
    "etau"   : ["FS_el_pt", "FS_tau_pt", "FS_el_eta", "FS_tau_eta",
                "FS_mt", "PuppiMET_pt"],
    "dimuon" : ["FS_m1_pt", "FS_m2_pt", "FS_m1_eta", "FS_m2_eta"],
  }

  vars_to_plot = ["HTT_m_vis", "HTT_dR"] # common to all final states
  FS_vars_to_add = final_state_plotting_vars[final_state_mode]
  for var in FS_vars_to_add:
    vars_to_plot.append(var)
  jet_vars_to_add = jet_plotting_vars[jet_mode]
  for jet_var in jet_vars_to_add:
    vars_to_plot.append(jet_vars_to_add)

  return vars_to_plot

def set_protected_branches(final_state_mode="none", jet_mode="none", DeepTau_version="none"):
  '''
  Set branches to be protected (i.e. not cut on) when using "apply_cut."
  Generally, you should protect any branches introduced by a cut.
  '''
  if final_state_mode == "none": # no mode given, assume jet cut
    initial_branches   = ["HTT_m_vis", "HTT_dR"]
    protected_branches = ["pass_zero_jet_cuts", "pass_one_jet_cuts", "pass_two_jet_cuts",
                          "pass_two_or_more_jet_cuts"]

  else:
    initial_branches = set_branches(final_state_mode, DeepTau_version="2p5")
    vars_to_trim     = set_vars_to_plot(final_state_mode, jet_mode)
    protected_branches = [var for var in vars_to_trim if var not in initial_branches]

  return protected_branches




