import numpy as np

from get_and_set_functions import add_DeepTau_branches, add_trigger_branches
from calculate_functions   import calculate_mt

def append_lepton_indices(event_dictionary):
  '''
  Read the entries of "FSLeptons" and extract the values to place in separate branches.
  It was easier to do this once when the data is first loaded than to do it every time
  that it is needed. 
  '''
  FSLeptons = event_dictionary["FSLeptons"]
  l1_indices, l2_indices = [], []
  for event in FSLeptons:
    l1_indices.append(event[0])
    l2_indices.append(event[1])
  event_dictionary["l1_indices"] = np.array(l1_indices)
  event_dictionary["l2_indices"] = np.array(l2_indices)
  return event_dictionary


def make_ditau_cut(event_dictionary, DeepTauVersion):
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
  unpack_ditau = add_DeepTau_branches(unpack_ditau, DeepTauVersion)
  unpack_ditau = (event_dictionary.get(key) for key in unpack_ditau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_ditau] # "*" unpacks a tuple
  pass_cuts, FS_t1_pt, FS_t2_pt, FS_t1_eta, FS_t2_eta = [], [], [], [], []
  # note these are in the same order as the variables in the first line of this function :)
  for i, lep_pt, lep_eta, tau_idx, l1_idx, l2_idx, vJet, vMu, vEle in zip(*to_check):
    passKinems = (lep_pt[l1_idx] >= 40 and lep_pt[l2_idx] >= 40)
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

def make_mutau_cut(event_dictionary, DeepTauVersion):
  '''
  Works similarly to 'make_ditau_cut'. 
  Notably, the mutau cuts are more complicated, but it is simple to 
  extend the existing methods as long as one can stomach the line breaks.
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_mutau = ["Tau_pt", "Tau_eta", "Muon_pt", "Muon_eta", "Muon_phi", "PuppiMET_pt", "PuppiMET_phi",
                  "Lepton_tauIdx", "Lepton_muIdx", "l1_indices", "l2_indices"]
  #TODO add this "CleanJet_btagWP" (no effect in August skims since it was always 1)
  #unpack_mutau.append("HTT_Lep_pt") # TODO delete after testing
  #unpack_mutau.append("HTT_Tau_pt") # TODO delete after testing
  unpack_mutau = add_DeepTau_branches(unpack_mutau, DeepTauVersion)
  unpack_mutau = add_trigger_branches(unpack_mutau, final_state_mode="mutau")
  unpack_mutau = (event_dictionary.get(key) for key in unpack_mutau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_mutau] # "*" unpacks a tuple
  pass_cuts, FS_mu_pt, FS_tau_pt, FS_mu_eta, FS_tau_eta, HTT_mt = [], [], [], [], [], []
  # note these are in the same order as the variables in the first line of this function :)
  for i, tau_pt, tau_eta, mu_pt, mu_eta, mu_phi, MET_pt, MET_phi, tau_idx, mu_idx,\
      l1_idx, l2_idx, vJet, vMu, vEle, trg24mu, trg27mu, crosstrg, _ in zip(*to_check):
      #l1_idx, l2_idx, HTT_Lep_pt, HTT_Tau_pt, vJet, vMu, vEle, trg24mu, trg27mu, crosstrg, _ in zip(*to_check):

    tauLoc     = tau_idx[l1_idx] + tau_idx[l2_idx] + 1
    muLoc      = mu_idx[l1_idx]  + mu_idx[l2_idx]  + 1
    tauEtaVal  = tau_eta[tauLoc]
    tauPtVal   = tau_pt[tauLoc] 
    muPtVal    = mu_pt[muLoc] 
    muEtaVal   = mu_eta[muLoc]
    muPhiVal   = mu_phi[muLoc]
    mtVal      = calculate_mt(muPtVal, muPhiVal, MET_pt, MET_phi)

    #goodMuonsAndTausCut  = "HTT_Tau_pt > 30 && ( (HLT_IsoMu24 && HTT_Lep_pt > 25.) || \
    #                      (!HLT_IsoMu24 && HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 && \
    #                       HTT_Lep_pt > 21. && Tau_eta[Lepton_tauIdx[FSLeptons[0]] + Lepton_tauIdx[FSLeptons[1]] + 1]) )"

    #dumbTestTau = (HTT_Tau_pt > 30.)
    #dumbTestMuTrig = (trg24mu and HTT_Lep_pt > 25.)
    #dumbTestCrossTrig = (not (trg24mu or trg27mu) and crosstrg and HTT_Lep_pt > 21 and tauEtaVal < 2.1)
    dumbTestMuTrig = (trg24mu and muPtVal > 25.)
    dumbTestCrossTrig = (not (trg24mu or trg27mu) and crosstrg and muPtVal > 21 and tauEtaVal < 2.1)

    passMT     = (mtVal < 50.)
    passTauPt  = (tauPtVal > 30.)
    pass25MuPt   = (trg24mu and muPtVal > 25.)
    pass28MuPt   = (trg27mu and muPtVal > 28.)
    passMuPtCrossTrigger = (crosstrg and (21. < muPtVal < 25.) and abs(tauEtaVal) < 2.1)
    passTauDT  = (vJet[tauLoc] >= 5 and vMu[tauLoc] >= 4 and vEle[tauLoc] >= 1)
    ## for tt test
    ##passMuPtCrossTrigger = ( not (pass25MuPt or pass28MuPt) and crosstrg and (muPtVal > 21.) and tauEtaVal < 2.1)
    ##if ( passMT and (passTauPt and (pass25MuPt or passMuPtCrossTrigger)) and passTauDT): 

    #if ( passMT and (passTauPt and (pass25MuPt or pass28MuPt or passMuPtCrossTrigger)) and passTauDT):
    noCut = True
    if noCut:
    #if (dumbTestTau and (dumbTestMuTrig or dumbTestCrossTrig) and passTauDT and passMT):
    #if (passTauPt and (dumbTestMuTrig or dumbTestCrossTrig) and passTauDT and passMT):
      pass_cuts.append(i)
      FS_mu_pt.append(muPtVal)
      FS_tau_pt.append(tauPtVal)
      FS_mu_eta.append(muEtaVal)
      FS_tau_eta.append(tauEtaVal)
      HTT_mt.append(mtVal)

  event_dictionary["pass_cuts"] = np.array(pass_cuts)
  event_dictionary["FS_mu_pt"]  = np.array(FS_mu_pt)
  event_dictionary["FS_tau_pt"] = np.array(FS_tau_pt)
  event_dictionary["FS_mu_eta"] = np.array(FS_mu_eta)
  event_dictionary["FS_tau_eta"] = np.array(FS_tau_eta)
  event_dictionary["HTT_mt"]    = np.array(HTT_mt)
  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after mutau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary

def manual_dimuon_lepton_veto(event_dictionary):
  '''
  Works similarly to 'make_ditau_cut' except the branch "pass_manual_lepton_veto"
  is made specifically for the dimuon final state. Some special handling is required
  due to the way events are selected in step2 of the NanoTauFramework
  '''
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
  print(f"events passing manual dimuon lepton veto = {len(np.array(pass_manual_lepton_veto))}")
  return event_dictionary


def make_dimuon_cut(event_dictionary):
  '''
  Works similarly to 'make_ditau_cut'. 
  '''
  unpack_dimuon = ["Lepton_pt", "Lepton_iso", "HTT_m_vis", "HTT_dR", "l1_indices", "l2_indices"]
  unpack_dimuon = (event_dictionary.get(key) for key in unpack_dimuon)
  to_check      = [range(len(event_dictionary["Lepton_pt"])), *unpack_dimuon]
  pass_cuts = []
  #TODO : extend function to add dimuon variables for plotting
  for i, pt, iso, mvis, dR, l1_idx, l2_idx in zip(*to_check):
    passKinematics = (pt[l1_idx] > 26 and pt[l2_idx] > 20 and mvis > 20 and dR > 0.5)
    passIso        = (iso[l1_idx] < 0.15 and iso[l2_idx] < 0.15)
    if (passKinematics and passIso):
      pass_cuts.append(i)

  event_dictionary["pass_cuts"] = np.array(pass_cuts)
  print(f"events passing dimuon cuts = {len(np.array(pass_cuts))}")
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


def apply_cut(event_dictionary, cut_branch):
  '''
  Remove all entries in 'event_dictionary' not in 'cut_branch' using the numpy 'take' method.
  Branches that are added during previous cut steps are added here because their entries
  already pass cuts by construction.
  The returned event_dictionary now only contains events passing all cuts.
  '''
  branches_added_during_cut = ["FS_t1_pt", "FS_t2_pt", "FS_t1_eta", "FS_t2_eta",
                               "FS_mu_pt", "FS_tau_pt", "FS_mu_eta", "FS_tau_eta",
                               "HTT_mt"]
  for branch in event_dictionary:
    if branch != cut_branch:
      if branch in branches_added_during_cut:
        pass
      else:
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
  data_events = apply_cut(data_events, "pass_run_cut")

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


def make_final_state_cut(event_dictionary, useDeepTauVersion, final_state_mode):
  '''
  Organizational function that generalizes call to a (set of) cuts based on the
  final cut. Importantly, the function that rejects events, 'apply_cut',
  is called elsewhere
  '''
  if final_state_mode == "ditau":
    event_dictionary = make_ditau_cut(event_dictionary, useDeepTauVersion)
  elif final_state_mode == "mutau":
    event_dictionary = make_mutau_cut(event_dictionary, useDeepTauVersion)
  elif final_state_mode == "etau":
    event_dictionary = make_etau_cut(event_dictionary, useDeepTauVersion)
  elif final_state_mode == "dimuon":
    event_dictionary = manual_dimuon_lepton_veto(event_dictionary)
    event_dictionary = apply_cut(process_events, "pass_manual_lepton_veto")
    event_dictionary = make_dimuon_cut(event_dictionary)
  else:
    print(f"No cuts to apply for {final_state_mode} final state.")
  return event_dictionary


