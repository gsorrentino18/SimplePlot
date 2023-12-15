import numpy as np

### README
# this file contains functions to perform cuts and self-contained studies

from calculate_functions import calculate_mt
from utility_functions   import time_print

from cut_ditau_functions import make_ditau_cut # is this something where using import * is justified?
from cut_mutau_functions import make_mutau_cut
from branch_functions    import add_trigger_branches, add_DeepTau_branches


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


def append_flavor_indices(event_dictionary, final_state_mode, keep_fakes=False):
  unpack_flav = ["l1_indices", "l2_indices", "Lepton_tauIdx", "Tau_genPartFlav"]
  unpack_flav = (event_dictionary.get(key) for key in unpack_flav)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_flav]
  FS_t1_flav, FS_t2_flav = [], []
  pass_gen_cuts, genuine_events, jet_fake_events, lep_fake_events = [], [], [], []
  for i, l1_idx, l2_idx, tau_idx, tau_flav in zip(*to_check):
    genuine, lep_fake, jet_fake = False, False, False
    t1_flav = -1
    t2_flav = -1
    if final_state_mode == "ditau":
      t1_flav = tau_flav[tau_idx[l1_idx]]
      t2_flav = tau_flav[tau_idx[l2_idx]]
      if (t1_flav == 5) and (t2_flav == 5):
        # genuine tau --> both taus are taus at gen level
        genuine = True
        genuine_events.append(i)
      elif (t1_flav == 0) or (t2_flav == 0):
        # jet fake --> one tau is faked by jet
        jet_fake = True
        jet_fake_events.append(i)
      elif (t1_flav < 5 and t1_flav > 0) or (t2_flav < 5 and t1_flav > 0):
        # lep fake --> both taus are faked by lepton
        # event with one tau faking jet enters category above first due to ordering
        # implies also the case where both are faked but one is faked by lepton 
        # is added to jet fakes, which i think is fine
        lep_fake = True
        lep_fake_events.append(i)
    elif ((final_state_mode == "mutau") or (final_state_mode == "etau")):
      t1_flav = tau_flav[tau_idx[l1_idx] + tau_idx[l2_idx] + 1] # update with NanoAODv12 samples
      if (t1_flav == 5):
        genuine = True
        genuine_events.append(i)
      elif (t1_flav == 0):
        jet_fake = True
        jet_fake_events.append(i)
      elif (t1_flav < 5 and t1_flav > 0):
        lep_fake = True
        lep_fake_events.append(i)

    else:
      print(f"No gen matching for that final state ({final_state_mode}), crashing...")
      return None
  
    if (keep_fakes==False) and ((genuine) or (lep_fake)):
      # save genuine background events and lep_fakes, remove jet fakes with gen matching
      FS_t1_flav.append(t1_flav)
      FS_t2_flav.append(t2_flav)
      pass_gen_cuts.append(i)

    if (keep_fakes==True) and ((genuine) or (lep_fake) or (jet_fake)):
      # save all events and their flavors, even if they are jet fakes
      FS_t1_flav.append(t1_flav)
      FS_t2_flav.append(t2_flav)
      pass_gen_cuts.append(i)

  event_dictionary["FS_t1_flav"] = np.array(FS_t1_flav)
  event_dictionary["FS_t2_flav"] = np.array(FS_t2_flav)
  event_dictionary["pass_gen_cuts"] = np.array(pass_gen_cuts)
  return event_dictionary


def make_TnP_cut(event_dictionary, DeepTau_version, numerator=False):
  '''
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_TnP = ["Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_tauIdx",
                "l1_indices", "l2_indices"]
  return event_dictionary

'''
def make_mutau_cut(event_dictionary, DeepTau_version):
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_mutau = ["Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_iso",
                  "Muon_dxy", "Muon_dz", "Tau_dxy", "Tau_dz", "Tau_decayMode",
                  #"MET_pt", "MET_phi", "Muon_phi",
                  "PuppiMET_pt", "PuppiMET_phi",
                  "Lepton_tauIdx", "Lepton_muIdx", "l1_indices", "l2_indices"]
  #TODO add this "CleanJet_btagWP" (no effect in August skims since it was always 1)
  unpack_mutau = add_DeepTau_branches(unpack_mutau, DeepTau_version)
  unpack_mutau = add_trigger_branches(unpack_mutau, final_state_mode="mutau")
  unpack_mutau = (event_dictionary.get(key) for key in unpack_mutau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_mutau] # "*" unpacks a tuple
  pass_cuts, FS_mt = [], []
  FS_mu_pt, FS_mu_eta, FS_mu_phi, FS_mu_iso, FS_mu_dxy, FS_mu_dz = [], [], [], [], [], []
  FS_tau_pt, FS_tau_eta, FS_tau_phi, FS_tau_dxy, FS_tau_dz = [], [], [], [], []
  # note these are in the same order as the variables in the first line of this function :)
      #MET_pt, MET_phi, muon_phi, tau_idx, mu_idx,\
  for i, lep_pt, lep_eta, lep_phi, lep_iso,\
      mu_dxy, mu_dz, tau_dxy, tau_dz, tau_decayMode,\
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
    tauPtVal   = lep_pt[tauFSLoc] 
    tauEtaVal  = lep_eta[tauFSLoc]
    tauPhiVal  = lep_phi[tauFSLoc]
    tauDxyVal  = abs(tau_dxy[tauBranchLoc])
    tauDzVal   = tau_dz[tauBranchLoc]
    mtVal      = calculate_mt(muPtVal, muPhiVal, MET_pt, MET_phi)
    passMT     = True #(mtVal < 65.0) #(mtVal < 50.0) #mine
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
      FS_mt.append(mtVal)

  event_dictionary["pass_cuts"] = np.array(pass_cuts)
  event_dictionary["FS_mu_pt"]  = np.array(FS_mu_pt)
  event_dictionary["FS_mu_eta"] = np.array(FS_mu_eta)
  event_dictionary["FS_mu_phi"] = np.array(FS_mu_phi)
  event_dictionary["FS_mu_iso"] = np.array(FS_mu_iso)
  event_dictionary["FS_mu_dxy"] = np.array(FS_mu_dxy)
  event_dictionary["FS_mu_dz"]  = np.array(FS_mu_dz)
  event_dictionary["FS_tau_pt"]  = np.array(FS_tau_pt)
  event_dictionary["FS_tau_eta"] = np.array(FS_tau_eta)
  event_dictionary["FS_tau_phi"] = np.array(FS_tau_phi)
  event_dictionary["FS_tau_dxy"] = np.array(FS_tau_dxy)
  event_dictionary["FS_tau_dz"]  = np.array(FS_tau_dz)
  event_dictionary["FS_mt"]    = np.array(FS_mt)
  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after mutau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary
'''

def make_etau_cut(event_dictionary, DeepTau_version):
  '''
  Works similarly to 'make_ditau_cut'. 
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_etau = ["Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_iso",
                 "Electron_dxy", "Electron_dz", "Tau_dxy", "Tau_dz",
                 "PuppiMET_pt", "PuppiMET_phi",
                 "Lepton_tauIdx", "Lepton_elIdx", "l1_indices", "l2_indices"]
  unpack_etau = add_DeepTau_branches(unpack_etau, DeepTau_version)
  unpack_etau = add_trigger_branches(unpack_etau, final_state_mode="etau")
  unpack_etau = (event_dictionary.get(key) for key in unpack_etau)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_etau]
  pass_cuts, FS_mt = [], []
  FS_el_pt, FS_el_eta, FS_el_phi, FS_el_iso, FS_el_dxy, FS_el_dz = [], [], [], [], [], []
  FS_tau_pt, FS_tau_eta, FS_tau_phi, FS_tau_dxy, FS_tau_dz = [], [], [], [], []
  for i, lep_pt, lep_eta, lep_phi, lep_iso,\
      el_dxy, el_dz, tau_dxy, tau_dz,\
      MET_pt, MET_phi, tau_idx, el_idx,\
      l1_idx, l2_idx, vJet, vMu, vEle, trg32el, trg35el, crosstrg in zip(*to_check):

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
    tauPtVal   = lep_pt[tauFSLoc] 
    tauEtaVal  = lep_eta[tauFSLoc]
    tauPhiVal  = lep_phi[tauFSLoc]
    tauDxyVal  = abs(tau_dxy[tauBranchLoc])
    tauDzVal   = tau_dz[tauBranchLoc]
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
    passTauDT  = ((vJet[tauBranchLoc] >= 5) and (vMu[tauBranchLoc] >= 1) and (vEle[tauBranchLoc] >= 6))

    if (passMT and passTauPtAndEta and (pass33ElPt or pass36ElPt or passElPtCrossTrigger) and passTauDT):
      pass_cuts.append(i)
      FS_el_pt.append(elPtVal)
      FS_el_eta.append(elEtaVal)
      FS_el_phi.append(elPhiVal)
      FS_el_iso.append(elIsoVal)
      FS_el_dxy.append(elDxyVal)
      FS_el_dz.append(elDzVal)
      FS_tau_pt.append(tauPtVal)
      FS_tau_eta.append(tauEtaVal)
      FS_tau_phi.append(tauPhiVal)
      FS_tau_dxy.append(tauDxyVal)
      FS_tau_dz.append(tauDzVal)
      FS_mt.append(mtVal)

  event_dictionary["pass_cuts"]  = np.array(pass_cuts)
  event_dictionary["FS_el_pt"]   = np.array(FS_el_pt)
  event_dictionary["FS_el_eta"]  = np.array(FS_el_eta)
  event_dictionary["FS_el_phi"]  = np.array(FS_el_phi)
  event_dictionary["FS_el_iso"]  = np.array(FS_el_iso)
  event_dictionary["FS_el_dxy"]  = np.array(FS_el_dxy)
  event_dictionary["FS_el_dz"]   = np.array(FS_el_dz)
  event_dictionary["FS_tau_pt"]  = np.array(FS_tau_pt)
  event_dictionary["FS_tau_eta"] = np.array(FS_tau_eta)
  event_dictionary["FS_tau_phi"] = np.array(FS_tau_phi)
  event_dictionary["FS_tau_dxy"] = np.array(FS_tau_dxy)
  event_dictionary["FS_tau_dz"]  = np.array(FS_tau_dz)
  event_dictionary["FS_mt"]      = np.array(FS_mt)
  nEvents_postcut = len(np.array(pass_cuts))
  print(f"nEvents before and after ditau cuts = {nEvents_precut}, {nEvents_postcut}")
  return event_dictionary


def set_FF_values(final_state_mode, jet_mode_and_DeepTau_version):
  '''
  '''
  # should have aiso/iso as well
  FF_values = {
    # FS : { "jet_mode" : [intercept, slope] }  
    "ditau" : { 
      "custom_0j_2p5_FF" :    [0.273506,-0.000930995],#[0.273506, -0.000930995], # quick check with cutoff
      "custom_1j_2p5_FF" :    [0.243313,-0.000929758],#[0.222382, -0.000994264],
      "custom_GTE2j_2p5_FF" : [0.222382,-0.000994264],#[0.230973, -0.000759698],

      "custom_0j_2p5_CH" :    [1.14356 ,-0.000254023],
      "custom_1j_2p5_CH" :    [1.101 ,4.88618e-05],
      "custom_GTE2j_2p5_CH" : [1.11098,7.22662e-05],

      "custom_0j_2p5_check_FF"   : [0.315051, -0.00227768, 1.12693e-05],
      "custom_0j_2p5_check_Clos" : [1.90974, -0.0257612, 0.000156502],
      "custom_0j_2p5_check_CH"   : [1.15254, -0.00342661, 2.08671e-05],

      "custom_1j_2p5_check_FF"   : [0.330999, -0.00374668, 2.47092e-05],
      "custom_1j_2p5_check_Clos" : [1.56709, -0.0145289, 3.08194e-05],
      "custom_1j_2p5_check_CH"   : [1.29089, -0.00589395, 3.75957e-05],

      "custom_GTE2j_2p5_check_FF"   : [0.304739, -0.00343086, 1.88761e-05],
      "custom_GTE2j_2p5_check_Clos" : [1.99696, -0.0322061, 0.000239253],
      "custom_GTE2j_2p5_check_CH"   : [1.15858, -0.00220756, 1.75213e-05],

      "0j_2p1"      : [0.409537, -0.00166789],
      "1j_2p1"      : [0.338192, -0.00114901],
      "GTE2j_2p1"   : [0.274382, -0.000810031],
      "Closure_2p1" : [1, -0.001], # dummy values
      "OSSS_Bias_2p1" : [1, -0.001], # dummy values
      "0j_2p5"      : [0.277831, -0.000975272],
      "1j_2p5"      : [0.264218, -0.00121849],
      "GTE2j_2p5"   : [0.2398, -0.00124643],
      "Closure_2p5" : [1.0459, -0.00102224],
      "OSSS_Bias_2p5" : [1, 0], # dummy values
    },
    "mutau" : {  # wrong 2p5
      "0j_2p5"     : [0.037884, 0.000648851],
      "1j_2p5"     : [0.0348384, 0.000630731],
      "GTE2j_2p5"  : [0.0342287, 0.000358899],
    },
    "etau"  : {#Dummy values
      "0j_2p5"     : [1, 1], 
      "1j_2p5"     : [1, 1],
      "GTE2j_2p5"  : [1, 1],
    },
  } 
  intercept = FF_values[final_state_mode][jet_mode_and_DeepTau_version][0]
  slope     = FF_values[final_state_mode][jet_mode_and_DeepTau_version][1]
  #slope2    = FF_values[final_state_mode][jet_mode_and_DeepTau_version][2]

  #return intercept, slope, slope2
  return intercept, slope


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


def add_FF_weights(event_dictionary, jet_mode, DeepTau_version):
  unpack_FFVars = ["Lepton_pt", "HTT_m_vis", "l1_indices", "l2_indices"]
  unpack_FFVars = (event_dictionary.get(key) for key in unpack_FFVars)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_FFVars]
  FF_weights = []
  bins = [50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300]
  ditau_weight_map = {
    # these are all very close to 0.999, which means in the most pessimistic case
    # you lose 10 events in every 10000 events.
    # since we have at most 90000 events in the 0j case, the biggest error we can have is
    # 90 events. This is not big enough to effect plots dramatically, so as long as the
    # probs are close to 0.999, they don't need to be updated every time
    # SHOULD update them the final time
    "0j_2p5" : [bins,
          [0.998869, # < 50
           0.999233, 0.999735, 0.999866, 0.999873, 0.999901, 0.999944, 0.999944, 0.999939, 
           0.999930, 0.999924, 0.999921, 0.999913, 0.999900, 0.999899, 0.999885, 
           0.999871, 0.999826, 0.999691]], # > 200
    "1j_2p5" : [bins,
          [0.999676, 
           0.999512, 0.999607, 0.999771, 0.999810, 0.999827, 0.999864, 0.999851, 0.999826, 
           0.999807, 0.999795, 0.999776, 0.999775, 0.999743, 0.999740, 0.999726, 
           0.999700, 0.999630, 0.999492]],
    "GTE2j_2p5" : [bins,
          [0.999411, 
           0.999046, 0.999038, 0.999246, 0.999336, 0.999382, 0.999425, 0.999362, 0.999296, 
           0.999252, 0.999181, 0.999143, 0.999101, 0.999101, 0.999038, 0.999003, 
           0.998991, 0.998904, 0.998680]],
    "0j_2p1" : [bins,
          [0.996032, 
           0.997514, 0.999023, 0.999554, 0.999771, 0.999772, 0.999850, 0.999858, 
           0.999850, 0.999843, 0.999824, 0.999811, 0.999792, 0.999783, 0.999781, 0.999735, 
           0.999683, 0.999594, 0.999428]],
    "1j_2p1" : [bins,
          [0.998659, 
           0.998212, 0.998626, 0.999254, 0.999635, 0.999633, 0.999658, 0.999630, 
           0.999600, 0.999549, 0.999526, 0.999469, 0.999428, 0.999394, 0.999348, 0.999345, 
           0.999223, 0.999137, 0.998901]],
    "GTE2j_2p1" : [bins,
          [0.997985, 
           0.997208, 0.997390, 0.998144, 0.998912, 0.998832, 0.998663, 0.998491, 
           0.998358, 0.998206, 0.998065, 0.998010, 0.997843, 0.997769, 0.997668, 0.997625, 
           0.997477, 0.997151, 0.997048]],
  }
  
  FF_key = jet_mode + "_" + DeepTau_version

  intercept, slope = set_FF_values("ditau", "custom_"+jet_mode+"_2p5_FF")
  OSSS_bias_intercept, OSSS_bias_slope = set_FF_values("ditau", "custom_"+jet_mode+"_2p5_CH")

  #intercept, slope = set_FF_values("ditau", FF_key)
  #closure_intercept, closure_slope = set_FF_values("ditau", "Closure_" + DeepTau_version)
  #OSSS_bias_intercept, OSSS_bias_slope = set_FF_values("ditau", "OSSS_Bias_" + DeepTau_version)

  #intercept, slope, slope2 = set_FF_values("ditau", "custom_"+jet_mode+"_2p5_check_FF")
  #closure_intercept, closure_slope, closure_slope2 = set_FF_values("ditau", "custom_"+jet_mode+"_2p5_check_Clos")
  #OSSS_bias_intercept, OSSS_bias_slope, OSSS_bias_slope2 = set_FF_values("ditau", "custom_"+jet_mode+"_2p5_check_CH")

  for i, lep_pt, m_vis, l1_idx, l2_idx in zip(*to_check):
    if m_vis < bins[0]: # 50
      one_minus_MC_over_data_weight = ditau_weight_map[FF_key][1][0] # first weight
    elif m_vis > bins[-3]: # > 200
      if m_vis > bins[-1]: # > 300
        one_minus_MC_over_data_weight = ditau_weight_map[FF_key][1][-1] # last weight
      elif bins[-2] < m_vis < bins[-1]: # between 250 and 300
        one_minus_MC_over_data_weight = ditau_weight_map[FF_key][1][-2]
      elif bins[-3] < m_vis < bins[-2]: # between 200 and 250
        one_minus_MC_over_data_weight = ditau_weight_map[FF_key][1][-3]
    else: # mvis between 50 and 200
      m_vis_idx = int(m_vis // 10) - 5 # makes 50 bin zero idx
      m_vis_weight_idx = m_vis_idx + 1 # 0 in weights is < 50 weight
      one_minus_MC_over_data_weight = ditau_weight_map[FF_key][1][m_vis_weight_idx]

    l1_pt = lep_pt[l1_idx] if lep_pt[l1_idx] < 120.0 else 120.0
    l2_pt = lep_pt[l2_idx] if lep_pt[l2_idx] < 200.0 else 200.0
    m_vis = m_vis if m_vis < 350.0 else 350.0
    FF_weight = one_minus_MC_over_data_weight*(intercept + l1_pt * slope)
    #FF_weight *= (closure_intercept + lep_pt[l2_idx] * closure_slope)
    FF_weight *= (OSSS_bias_intercept + m_vis * OSSS_bias_slope)

    #FF_weight = one_minus_MC_over_data_weight*(intercept + l1_pt * slope + l1_pt*l1_pt*slope2)
    #FF_weight *= (closure_intercept + l2_pt * closure_slope + l2_pt*l2_pt*closure_slope2)
    #FF_weight *= (OSSS_bias_intercept + m_vis * OSSS_bias_slope + m_vis*m_vis*OSSS_bias_slope2)

    FF_weights.append(FF_weight)
  event_dictionary["FF_weight"] = np.array(FF_weights)
  return event_dictionary


def make_jet_cut(event_dictionary, jet_mode):
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_jetVars = ["nCleanJet", "CleanJet_pt", "CleanJet_eta"]
  unpack_jetVars = (event_dictionary.get(key) for key in unpack_jetVars)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_jetVars] # "*" unpacks a tuple
  nCleanJetGT30, pass_0j_cuts, pass_1j_cuts, pass_2j_cuts, pass_3j_cuts, pass_GTE2j_cuts = [], [], [], [], [], []
  CleanJetGT30_pt_1, CleanJetGT30_pt_2, CleanJetGT30_pt_3    = [], [], []
  CleanJetGT30_eta_1, CleanJetGT30_eta_2, CleanJetGT30_eta_3 = [], [], []
  for i, nJet, jet_pt, jet_eta in zip(*to_check):
    passingJets = 0
    passingJetsPt, passingJetsEta = [], []
    for ijet in range(0, nJet):
      if (jet_pt[ijet] > 30.0) and (jet_eta[ijet] < 4.7):
        passingJets += 1
        passingJetsPt.append(jet_pt[ijet])
        passingJetsEta.append(jet_eta[ijet])
    nCleanJetGT30.append(passingJets)

    if passingJets == 0: 
      pass_0j_cuts.append(i)

    if (passingJets == 1) and (jet_mode == "Inclusive" or jet_mode == "1j"): 
      # for GTE2j, fill this in block below
      pass_1j_cuts.append(i)
      CleanJetGT30_pt_1.append(passingJetsPt[0])
      CleanJetGT30_eta_1.append(passingJetsEta[0])

    if (passingJets == 2) and (jet_mode == "Inclusive" or jet_mode == "2j"):
      pass_2j_cuts.append(i)
      CleanJetGT30_pt_1.append(passingJetsPt[0])
      CleanJetGT30_pt_2.append(passingJetsPt[1])
      CleanJetGT30_eta_1.append(passingJetsEta[0])
      CleanJetGT30_eta_2.append(passingJetsEta[1])

    if (passingJets >= 2) and (jet_mode == "GTE2j"): 
      pass_GTE2j_cuts.append(i)
      CleanJetGT30_pt_1.append(passingJetsPt[0])
      CleanJetGT30_pt_2.append(passingJetsPt[1])
      CleanJetGT30_eta_1.append(passingJetsEta[0])
      CleanJetGT30_eta_2.append(passingJetsEta[1])
  #    elif passingJets >= 2 and (jet_mode == "GTE2j"):
  #      pass_2j_cuts.append(i)
  #      CleanJetGT30_pt_1.append(passingJetsPt[0])
  #      CleanJetGT30_pt_2.append(passingJetsPt[1])
  #      CleanJetGT30_eta_1.append(passingJetsEta[0])
  #      CleanJetGT30_eta_2.append(passingJetsEta[1])
  #      if passingJets == 3:
  #        pass_3j_cuts.append(i)
  #        CleanJetGT30_pt_3.append(passingJetsPt[2])
  #        CleanJetGT30_eta_3.append(passingJetsEta[2])
  #      elif passingJets > 3:
  #        pass
  #        #print("More than 3 FS jets")

  event_dictionary["nCleanJetGT30"]   = np.array(nCleanJetGT30)

  if jet_mode == "Inclusive":
    # fill branches like above
    event_dictionary["pass_0j_cuts"]    = np.array(pass_0j_cuts)
    event_dictionary["pass_1j_cuts"]    = np.array(pass_1j_cuts)
    event_dictionary["pass_2j_cuts"]    = np.array(pass_2j_cuts)
    event_dictionary["pass_3j_cuts"]    = np.array(pass_3j_cuts)
    event_dictionary["pass_GTE2j_cuts"] = np.array(pass_GTE2j_cuts)

    event_dictionary["CleanJetGT30_pt_1"]  = np.array(CleanJetGT30_pt_1)
    event_dictionary["CleanJetGT30_pt_2"]  = np.array(CleanJetGT30_pt_2)
    event_dictionary["CleanJetGT30_pt_3"]  = np.array(CleanJetGT30_pt_3)
    event_dictionary["CleanJetGT30_eta_1"] = np.array(CleanJetGT30_eta_1)
    event_dictionary["CleanJetGT30_eta_2"] = np.array(CleanJetGT30_eta_2)
    event_dictionary["CleanJetGT30_eta_3"] = np.array(CleanJetGT30_eta_3)
  
  elif jet_mode == "0j":
    # literally don't do any of the above
    event_dictionary["pass_0j_cuts"]    = np.array(pass_0j_cuts)

  elif jet_mode == "1j":
    # only do the 1j things
    event_dictionary["pass_1j_cuts"]    = np.array(pass_1j_cuts)
    event_dictionary["CleanJetGT30_pt_1"]  = np.array(CleanJetGT30_pt_1)
    event_dictionary["CleanJetGT30_eta_1"] = np.array(CleanJetGT30_eta_1)

  elif jet_mode == "2j":
    # only do the 2j things
    event_dictionary["pass_2j_cuts"]    = np.array(pass_2j_cuts)
    event_dictionary["CleanJetGT30_pt_1"]  = np.array(CleanJetGT30_pt_1)
    event_dictionary["CleanJetGT30_pt_2"]  = np.array(CleanJetGT30_pt_2)
    event_dictionary["CleanJetGT30_eta_1"] = np.array(CleanJetGT30_eta_1)
    event_dictionary["CleanJetGT30_eta_2"] = np.array(CleanJetGT30_eta_2)

  elif jet_mode == "3j" or jet_mode == "GTE2j":
    # importantly different from inclusive
    #event_dictionary["pass_2j_cuts"]    = np.array(pass_2j_cuts)
    #event_dictionary["pass_3j_cuts"]    = np.array(pass_3j_cuts)
    event_dictionary["pass_GTE2j_cuts"]    = np.array(pass_GTE2j_cuts)
    event_dictionary["CleanJetGT30_pt_1"]  = np.array(CleanJetGT30_pt_1)
    event_dictionary["CleanJetGT30_pt_2"]  = np.array(CleanJetGT30_pt_2)
    #event_dictionary["CleanJetGT30_pt_3"]  = np.array(CleanJetGT30_pt_3)
    event_dictionary["CleanJetGT30_eta_1"] = np.array(CleanJetGT30_eta_1)
    event_dictionary["CleanJetGT30_eta_2"] = np.array(CleanJetGT30_eta_2)
    #event_dictionary["CleanJetGT30_eta_3"] = np.array(CleanJetGT30_eta_3)

  # can only do this if inclusive
  if jet_mode == "Inclusive":
    print("nEvents with exactly 0,1,2,3 jets and ≥2 jets")
    print(f"{len(np.array(pass_0j_cuts))}, {len(np.array(pass_1j_cuts))}, {len(np.array(pass_2j_cuts))}, {len(np.array(pass_3j_cuts))}, {len(np.array(pass_GTE2j_cuts))}")

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
    event_passes_manual_lepton_veto = False
    nIsoEle, nIsoMu = 0, 0 # there are many pdgId=15 particles, but we assume those are fake taus
    for pdgId, iso in zip(lep_pdgId_array, lep_iso_array):
      if (abs(pdgId) == 11) and (iso < 0.3):
        nIsoEle += 1
      elif (abs(pdgId) == 13) and (iso < 0.3):
        nIsoMu  += 1
      else:
        pass

      if nIsoEle > 0:
        event_passes_manual_lepton_veto = False
      elif nIsoMu > 2:
        event_passes_manual_lepton_veto = False
      else:
        event_passes_manual_lepton_veto = True

    if event_passes_manual_lepton_veto:
      pass_manual_lepton_veto.append(i)

  event_dictionary["pass_manual_lepton_veto"] = np.array(pass_manual_lepton_veto)
  print(f"events before and after manual dimuon lepton veto = {nEvents_precut}, {len(np.array(pass_manual_lepton_veto))}")
  return event_dictionary


def make_dimuon_cut(event_dictionary, useMiniIso=False):
  '''
  Works similarly to 'make_ditau_cut'. 
  '''
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_dimuon = ["Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_iso", 
                   "Lepton_muIdx", "Muon_dxy", "Muon_dz",
                   "HTT_m_vis", "HTT_dR", "l1_indices", "l2_indices"]
  unpack_dimuon = (event_dictionary.get(key) for key in unpack_dimuon)
  to_check      = [range(len(event_dictionary["Lepton_pt"])), *unpack_dimuon]
  pass_cuts = []
  FS_m1_pt, FS_m1_eta, FS_m1_phi, FS_m1_iso, FS_m1_dxy, FS_m1_dz = [], [], [], [], [], []
  FS_m2_pt, FS_m2_eta, FS_m2_phi, FS_m2_iso, FS_m2_dxy, FS_m2_dz = [], [], [], [], [], []
  for i, pt, eta, phi, iso, muIdx, mu_dxy, mu_dz, mvis, dR, l1_idx, l2_idx in zip(*to_check):
    # removed (dR > 0.5) and changed (mvis > 20) cut. Our minimum dR is 0.3 from skim level
    #passKinematics = (pt[l1_idx] > 26 and pt[l2_idx] > 20 and (mvis > 20) and (dR > 0.5)
    passKinematics = (pt[l1_idx] > 26 and pt[l2_idx] > 20 and (70 < mvis < 130))
    if (useMiniIso == False):
      passIso      = (iso[l1_idx] < 0.25 and iso[l2_idx] < 0.25) # for PFRelIso, Loose 25, Medium 20, Tight 15
    if (useMiniIso == True):
      passIso      = (iso[l1_idx] < 0.40 and iso[l2_idx] < 0.40) # for MiniIso, Loose 40, Medium 20, Tight 10
    if (passKinematics and passIso):
      pass_cuts.append(i)
      FS_m1_pt.append(pt[l1_idx])
      FS_m1_eta.append(eta[l1_idx])
      FS_m1_phi.append(phi[l1_idx])
      FS_m1_iso.append(iso[l1_idx])
      FS_m1_dxy.append(abs(mu_dxy[muIdx[l1_idx]]))
      FS_m1_dz.append(mu_dz[muIdx[l1_idx]])
      FS_m2_pt.append(pt[l2_idx])
      FS_m2_eta.append(eta[l2_idx])
      FS_m2_phi.append(phi[l2_idx])
      FS_m2_iso.append(iso[l2_idx])
      FS_m2_dxy.append(abs(mu_dxy[muIdx[l2_idx]]))
      FS_m2_dz.append(mu_dz[muIdx[l2_idx]])

  event_dictionary["pass_cuts"] = np.array(pass_cuts)
  event_dictionary["FS_m1_pt"]  = np.array(FS_m1_pt)
  event_dictionary["FS_m1_eta"] = np.array(FS_m1_eta)
  event_dictionary["FS_m1_phi"] = np.array(FS_m1_phi)
  event_dictionary["FS_m1_iso"] = np.array(FS_m1_iso)
  event_dictionary["FS_m1_dxy"] = np.array(FS_m1_dxy)
  event_dictionary["FS_m1_dz"] = np.array(FS_m1_dz)
  event_dictionary["FS_m2_pt"]  = np.array(FS_m2_pt)
  event_dictionary["FS_m2_eta"] = np.array(FS_m2_eta)
  event_dictionary["FS_m2_phi"] = np.array(FS_m2_phi)
  event_dictionary["FS_m2_iso"] = np.array(FS_m2_iso)
  event_dictionary["FS_m2_dxy"] = np.array(FS_m2_dxy)
  event_dictionary["FS_m2_dz"] = np.array(FS_m2_dz)
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
    print("All events removed, sample will be deleted")
    delete_sample = True
    return None
 
  #print(f"cut branch: {cut_branch}") # DEBUG
  #print(f"protected branches: {protected_branches}") # DEBUG
  for branch in event_dictionary:
    if delete_sample:
      pass

    # special handling, will need to be adjusted by hand for excatly 2j or 3j studies # DEBUG
    # this only works for GTE2j, not Inclusive because the "apply_cut" method for jets is never called there # DEBUG
    if (("pass_GTE2j_cuts" in event_dictionary) and
        (branch == "HTT_DiJet_dEta_fromHighestMjj" or branch == "HTT_DiJet_MassInv_fromHighestMjj")):
      #print("very special GTE2j handling underway") # DEBUG
      event_dictionary[branch] = np.take(event_dictionary[branch], event_dictionary["pass_GTE2j_cuts"])
      #if (branch == "CleanJetGT30_pt_3" or branch == "CleanJetGT30_eta_3"):
      #  event_dictionary[branch] = np.take(event_dictionary[branch], event_dictionary["pass_3j_cuts"])

    elif ((branch != cut_branch) and (branch not in protected_branches)):
      #print(f"going to cut {branch}, {len(event_dictionary[branch])}") # DEBUG
      event_dictionary[branch] = np.take(event_dictionary[branch], event_dictionary[cut_branch])

  return event_dictionary


def Era_F_trigger_study(data_events, final_state_mode):
  '''
  Compact function for 2022 era F trigger study, where ChargedIsoTau
  triggers were briefly enabled for Run2-Run3 Tau trigger studies. 
  '''
  from triggers_dictionary import triggers_dictionary
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


def apply_final_state_cut(event_dictionary, final_state_mode, DeepTau_version, useMiniIso=False):
  '''
  Organizational function that generalizes call to a (set of) cuts based on the
  final cut. Importantly, the function that rejects events, 'apply_cut',
  is called elsewhere
  '''
  # setting inclusive in the jet_mode includes all jet branches in protected branches
  # this is okay because in the current ordering (FS cut then jet cut), no jet branches
  # are event created yet.
  protected_branches = set_protected_branches(final_state_mode=final_state_mode, jet_mode="Inclusive")
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
    # old samples need manual lepton veto
    if (useMiniIso == False):
      event_dictionary = manual_dimuon_lepton_veto(event_dictionary)
      event_dictionary = apply_cut(event_dictionary, "pass_manual_lepton_veto")
      event_dictionary = make_dimuon_cut(event_dictionary)
      event_dictionary = apply_cut(event_dictionary, "pass_cuts", protected_branches)
    # new samples don't and they use a different iso
    if (useMiniIso == True):
      event_dictionary = make_dimuon_cut(event_dictionary, useMiniIso==True)
      event_dictionary = apply_cut(event_dictionary, "pass_cuts", protected_branches)
  else:
    print(f"No cuts to apply for {final_state_mode} final state.")
  return event_dictionary


def apply_AR_cut(event_dictionary, final_state_mode, jet_mode, DeepTau_version):
  '''
  Organizational function
  added 'free_pass_AR' so that FS branches are populated without any selection applied
  added 'skip_DeepTau' to apply a partial selection (all but leading tau deeptau reqs)
  '''
  protected_branches = ["None"]
  event_dictionary = append_lepton_indices(event_dictionary)
  if (final_state_mode == "ditau") and (jet_mode != "Inclusive"):
    event_dictionary = make_ditau_AR_cut(event_dictionary, DeepTau_version)
    event_dictionary = apply_cut(event_dictionary, "pass_AR_cuts", protected_branches)
    event_dictionary = apply_jet_cut(event_dictionary, jet_mode)
    # non-standard FS cut
    event_dictionary   = make_ditau_cut(event_dictionary, DeepTau_version, free_pass_AR=False, skip_DeepTau=True)
    protected_branches = set_protected_branches(final_state_mode="ditau", jet_mode="none")
    event_dictionary   = apply_cut(event_dictionary, "pass_cuts", protected_branches)
    #
    event_dictionary = add_FF_weights(event_dictionary, jet_mode, DeepTau_version)
  else:
    print(f"{final_state_mode} : {jet_mode} not possible. Continuing without AR or FF method applied.")
  return event_dictionary


def apply_jet_cut(event_dictionary, jet_mode):
  '''
  Organizational function to reduce event_dictionary to contain only
  events with jets passing certain criteria. Enables plotting of jet objects
  jet_mode can be "Inclusive", "0j", "1j", "2j", "3j", "GTE2j",
  '''
  jet_cut_branch = {
    "Inclusive" : "Inclusive",
    "0j" : "pass_0j_cuts",
    "1j" : "pass_1j_cuts",
    "2j" : "pass_2j_cuts",
    "3j" : "pass_3j_cuts",
    "GTE2j" : "pass_GTE2j_cuts",
  }
  event_dictionary   = make_jet_cut(event_dictionary, jet_mode)
  protected_branches = set_protected_branches(final_state_mode="none", jet_mode=jet_mode)
  if jet_mode == "Inclusive":
    print("jet mode is Inclusive, no jet cut performed")
  else:
    event_dictionary = apply_cut(event_dictionary, jet_cut_branch[jet_mode], protected_branches)
  return event_dictionary


def apply_HTT_FS_cuts_to_process(process, process_dictionary, 
                                 final_state_mode, jet_mode="Inclusive", 
                                 DeepTau_version="2p5", useMiniIso=False):
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
  if ("Data" not in process) and (final_state_mode=="ditau"): # unify with other DMs later
    process_events = append_flavor_indices(process_events, final_state_mode)
    process_events = apply_cut(process_events, "pass_gen_cuts", protected_branches=["FS_t1_flav", "FS_t2_flav", "pass_gen_cuts"])
    if (process_events==None or len(process_events["run"])==0): return None

  if ("Data" not in process) and (final_state_mode=="mutau"): # unify with other DMs later
    if ("TT" in process) or ("WJ" in process):
      process_events = append_flavor_indices(process_events, final_state_mode, keep_fakes=True)
      process_events = apply_cut(process_events, "pass_gen_cuts", protected_branches=["FS_t1_flav", "FS_t2_flav", "pass_gen_cuts"])
      if (process_events==None or len(process_events["run"])==0): return None
    else:
      process_events = append_flavor_indices(process_events, final_state_mode)
      process_events = apply_cut(process_events, "pass_gen_cuts", protected_branches=["FS_t1_flav", "FS_t2_flav", "pass_gen_cuts"])
      if (process_events==None or len(process_events["run"])==0): return None

  FS_cut_events = apply_final_state_cut(process_events, final_state_mode, DeepTau_version, useMiniIso=useMiniIso)
  if (FS_cut_events==None or len(FS_cut_events["run"])==0): return None 
  cut_events = apply_jet_cut(FS_cut_events, jet_mode)
  if (cut_events==None or len(cut_events["run"])==0): return None

  # TODO : want to move to this
  #jet_cut_events = apply_jet_cut(process_events, jet_mode)
  #if len(jet_cut_events["run"])==0: return None
  #FS_cut_events = apply_final_state_cut(jet_cut_events, final_state_mode, DeepTau_version, useMiniIso=useMiniIso)
  #if len(FS_cut_events["run"])==0: return None  

  return FS_cut_events


def set_good_events(final_state_mode, disable_triggers=False, useMiniIso=False):
  '''
  Return a string defining a 'good_events' flag used by uproot to preskim input events
  to only those passing these simple requirements. 'good_events' changes based on
  final_state_mode, and the trigger condition is removed if a trigger study is 
  being conducted (since requiring the trigger biases the study).
  '''
  good_events = ""
  if disable_triggers: print("*"*20 + " removed trigger requirement " + "*"*20)

  # relevant definitions from NanoTauAnalysis /// modules/TauPairSelector.py
  # HTT_SRevent and HTT_ARevent require opposite sign objects
  # HTT_SRevent = ((pdgIdPair < 0) 
  #            and ( ((LeptonIso < 0.2) and (abs(pdgIdPair)==11*13)) or (LeptonIso < 0.15)) 
  #            and TauPassVsJet and (self.leptons[finalpair[1]].pt > 15))
  # HTT_ARevent = ((pdgIdPair < 0) 
  #            and ( ((LeptonIso < 0.2) and (abs(pdgIdPair)==11*13)) or (LeptonIso < 0.15)) 
  #            and (not TauPassVsJet) and (self.leptons[finalpair[1]].pt > 15))
  #     # All SR requirements besides TauPassVsJet
  # HTT_SSevent = ((pdgIdPair > 0) 
  #            and ( ((LeptonIso < 0.2) and (abs(pdgIdPair)==11*13)) or (LeptonIso < 0.15)) 
  #            and TauPassVsJet and (self.leptons[finalpair[1]].pt > 15)) 
  #     # All SR requirements besides opposite sign
  
  # apply FS cut separately so it can be used with reject_duplicate_events
  if final_state_mode == "ditau":
    #good_events = "(HTT_SRevent) & (METfilters) & (LeptonVeto==0) & (abs(HTT_pdgId)==15*15) & (Trigger_ditau)"
    #triggers = "(HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1\
    #           | HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60\
    #           | HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75\
    #           | HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1\
    #           | HLT_DoublePFJets40_Mass500_MediumDeepTauPFTauHPS45_L2NN_MediumDeepTauPFTauHPS20_eta2p1)"
    triggers = "(HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1)"

    good_events = "(HTT_SRevent) & (METfilters) & (LeptonVeto==0) & (abs(HTT_pdgId)==15*15) & " + triggers
    if disable_triggers: good_events = good_events.replace(" & (Trigger_ditau)", "")

  elif final_state_mode == "mutau":
    good_events = "(HTT_SRevent) & (METfilters) & (LeptonVeto==0) & (abs(HTT_pdgId)==13*15) & (Trigger_mutau)"
    #good_events = "(HTT_SRevent) & (METfilters) & (abs(HTT_pdgId)==13*15) & (Trigger_mutau)"
    if disable_triggers: good_events = good_events.replace(" & (Trigger_mutau)", "")

  elif final_state_mode == "etau":
    good_events = "(HTT_SRevent) & (METfilters) & (LeptonVeto==0) & (abs(HTT_pdgId)==11*15) & (Trigger_etau)"
    if disable_triggers: good_events = good_events.replace(" & (Trigger_etau)", "")

  elif final_state_mode == "dimuon":
    # lepton veto must be applied manually for this final state
    if (useMiniIso == False):
      good_events = "(METfilters) & (HTT_pdgId==-13*13) & (HLT_IsoMu24)"
    if (useMiniIso == True):
      good_events = "(METfilters) & (LeptonVeto==0) & (HTT_pdgId==-13*13) & (HLT_IsoMu24)"
    if disable_triggers: good_events = good_events.replace(" & (HLT_IsoMu24)", "")

  return good_events


def set_branches(final_state_mode, DeepTau_version):
  common_branches = [
    "run", "luminosityBlock", "event", "Generator_weight",
    "FSLeptons", "Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_iso",
    "Tau_genPartFlav", "Tau_decayMode",
    "nCleanJet", "CleanJet_pt", "CleanJet_eta",
    "HTT_m_vis", "HTT_dR",
    #"HTT_DiJet_dEta_fromHighestMjj", "HTT_DiJet_MassInv_fromHighestMjj",
  ]
  branches = common_branches
  branches = add_final_state_branches(branches, final_state_mode)
  if final_state_mode != "dimuon": branches = add_DeepTau_branches(branches, DeepTau_version)
  branches = add_trigger_branches(branches, final_state_mode)
  return branches


def add_final_state_branches(branches_, final_state_mode):
  '''
  Helper function to add only relevant branches to loaded branches based on final state.
  '''
  final_state_branches = {
    "ditau"  : ["Lepton_tauIdx", "Tau_dxy", "Tau_dz"],

    "mutau"  : ["Muon_dxy", "Muon_dz",
                "Tau_dxy", "Tau_dz",
                "Lepton_tauIdx", "Lepton_muIdx", 
                #"MET_pt", "MET_phi", "Muon_phi"],
                "PuppiMET_pt", "PuppiMET_phi"],

    "etau"   : ["Electron_dxy", "Electron_dz",
                "Tau_dxy", "Tau_dz",
                "Lepton_tauIdx", "Lepton_elIdx",
                "PuppiMET_pt", "PuppiMET_phi"],

    "dimuon" : ["Lepton_pdgId", "Lepton_muIdx",
                "Muon_dxy", "Muon_dz"],
  }

  branch_to_add = final_state_branches[final_state_mode]
  for new_branch in branch_to_add:
    branches_.append(new_branch)
  
  return branches_



# this is ugly and bad and i am only doing this out of desperation
clean_jet_vars = {
    "Inclusive" : ["nCleanJetGT30",
      #"CleanJetGT30_pt_1", "CleanJetGT30_eta_1",
      #"CleanJetGT30_pt_2", "CleanJetGT30_eta_2",
      #"CleanJetGT30_pt_3", "CleanJetGT30_eta_3",
    ],

    "0j" : ["nCleanJetGT30"],
    "1j" : ["nCleanJetGT30", "CleanJetGT30_pt_1", "CleanJetGT30_eta_1"],
    "2j" : ["nCleanJetGT30",
      "CleanJetGT30_pt_1", "CleanJetGT30_eta_1",
      "CleanJetGT30_pt_2", "CleanJetGT30_eta_2",
    ],
    "3j" : ["nCleanJetGT30",
      #"CleanJetGT30_pt_1", "CleanJetGT30_eta_1",
      #"CleanJetGT30_pt_2", "CleanJetGT30_eta_2",
      #"CleanJetGT30_pt_3", "CleanJetGT30_eta_3",
    ],
    "GTE2j" : ["nCleanJetGT30",
      "CleanJetGT30_pt_1", "CleanJetGT30_eta_1",
      "CleanJetGT30_pt_2", "CleanJetGT30_eta_2",
      #"CleanJetGT30_pt_3", "CleanJetGT30_eta_3",
    ],
}

final_state_vars = {
    # can't put nanoaod branches here because this dictionary is used to protect branches created internally
    "none"   : [],
    "ditau"  : ["FS_t1_pt", "FS_t1_eta", "FS_t1_phi", "FS_t1_dxy", "FS_t1_dz",
                "FS_t2_pt", "FS_t2_eta", "FS_t2_phi", "FS_t2_dxy", "FS_t2_dz",
                "FS_t1_flav", "FS_t2_flav"],

    "mutau"  : ["FS_mu_pt", "FS_mu_eta", "FS_mu_phi", "FS_mu_iso", "FS_mu_dxy", "FS_mu_dz",
                "FS_tau_pt", "FS_tau_eta", "FS_tau_phi", "FS_tau_dxy", "FS_tau_dz",
                "FS_mt", "FS_t1_flav", "FS_t2_flav"],#, "PuppiMET_pt"],

    "etau"   : ["FS_el_pt", "FS_el_eta", "FS_el_phi", "FS_el_iso", "FS_el_dxy", "FS_el_dz",
                "FS_tau_pt", "FS_tau_eta", "FS_tau_phi", "FS_tau_dxy", "FS_tau_dz",
                "FS_mt", "FS_t1_flav", "FS_t2_flav"], #"PuppiMET_pt"],

    "dimuon" : ["FS_m1_pt", "FS_m1_eta", "FS_m1_phi", "FS_m1_iso", "FS_m1_dxy", "FS_m1_dz",
                "FS_m2_pt", "FS_m2_eta", "FS_m2_phi", "FS_m2_iso", "FS_m2_dxy", "FS_m2_dz"],
}

def set_vars_to_plot(final_state_mode, jet_mode="none"):
  '''
  Helper function to keep plotting variables organized
  Shouldn't this be in  plotting functions?
  '''
  vars_to_plot = ["HTT_m_vis", "HTT_dR", "PuppiMET_pt"] # common to all final states # TODO add MET here, add Tau_decayMode
  FS_vars_to_add = final_state_vars[final_state_mode]
  for var in FS_vars_to_add:
    vars_to_plot.append(var)

  jet_vars_to_add = clean_jet_vars[jet_mode]
  #if (jet_mode=="Inclusive") or (jet_mode=="GTE2j"):
  #  jet_vars_to_add += ["HTT_DiJet_dEta_fromHighestMjj", "HTT_DiJet_MassInv_fromHighestMjj"]
  for jet_var in jet_vars_to_add:
    vars_to_plot.append(jet_var)

  return vars_to_plot

# TODO fix this function and make it more straightforward
# way too easy to get confused with it currently
def set_protected_branches(final_state_mode, jet_mode, DeepTau_version="none"):
  '''
  Set branches to be protected (i.e. not cut on) when using "apply_cut."
  Generally, you should protect any branches introduced by a cut.

  protect all "FS" branches for FS cuts
  protect all "pass_xj_cuts" and "JetGT30_" branches for jet cuts
  '''

  if final_state_mode != "none": # not cutting FS branches
    protected_branches = final_state_vars[final_state_mode]
    # all "HTT_" branches automatically handled, just protecting "FS_" branches which were introduced by a cut
  
  elif final_state_mode == "none":
    if jet_mode == "Inclusive": # cutting FS branches, but not the jet branches
      protected_branches = ["pass_0j_cuts", "pass_1j_cuts", "pass_2j_cuts", "pass_3j_cuts", "pass_GTE2j_cuts"]
      # should fromHighestMjj branches be protected? it seems not
      protected_branches += clean_jet_vars[jet_mode]
      protected_branches = [var for var in protected_branches if var != "nCleanJetGT30"] # unprotect one branch

    elif jet_mode == "0j": # cutting FS branches, protecting just one jet branch
      protected_branches = ["pass_0j_cuts"]
      protected_branches += clean_jet_vars[jet_mode]
      protected_branches = [var for var in protected_branches if var != "nCleanJetGT30"] # unprotect one branch

    elif jet_mode == "1j":
      protected_branches = ["pass_0j_cuts", "pass_1j_cuts"]
      protected_branches += clean_jet_vars[jet_mode]
      protected_branches = [var for var in protected_branches if var != "nCleanJetGT30"] # unprotect one branch

    elif jet_mode == "2j":
      protected_branches = ["pass_0j_cuts", "pass_1j_cuts", "pass_2j_cuts"]
      protected_branches += clean_jet_vars[jet_mode]
      protected_branches = [var for var in protected_branches if var != "nCleanJetGT30"] # unprotect one branch

    elif jet_mode == "3j":
      protected_branches = ["pass_0j_cuts", "pass_1j_cuts", "pass_2j_cuts", "pass_3j_cuts"]
      protected_branches += clean_jet_vars[jet_mode]
      protected_branches = [var for var in protected_branches if var != "nCleanJetGT30"] # unprotect one branch

    elif jet_mode == "GTE2j":
      protected_branches = ["pass_0j_cuts", "pass_1j_cuts", "pass_2j_cuts", "pass_3j_cuts", "pass_GTE2j_cuts"]
      protected_branches += clean_jet_vars[jet_mode]
      protected_branches = [var for var in protected_branches if var != "nCleanJetGT30"] # unprotect one branch

  else:
    print("final state mode must be specified as 'none' or a valid final state to properly protect your branches")

  return protected_branches



