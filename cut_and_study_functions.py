import numpy as np

### README
# this file contains functions to perform cuts and self-contained studies

from calculate_functions import calculate_mt, hasbit, getBin, highest_mjj_pair
from utility_functions   import time_print, text_options, log_print

from cut_ditau_functions import make_ditau_cut, make_ditau_AR_cut
from cut_mutau_functions import make_mutau_cut, make_mutau_AR_cut, make_mutau_TnP_cut
from cut_etau_functions  import make_etau_cut, make_etau_AR_cut
from branch_functions    import add_trigger_branches, add_DeepTau_branches, add_Zpt_branches

# TODO : consider putting this function in a different file and importing it here
from MC_dictionary import MC_dictionary
from XSec import XSecRun3 as XSec
def load_and_store_NWEvents(process, event_dictionary):
  '''
  Read the NWEvents value for a sample and store it in the MC_dictionary,
  overriding the hardcoded values from V11 samples. Delete the NWEvents branch after.
  '''
  MC_dictionary[process]["NWEvents"] = event_dictionary["NWEvents"][0]
  MC_dictionary[process]["XSecMCweight"] = event_dictionary["XSecMCweight"][0]
  #print(process, MC_dictionary[process]["NWEvents"]) # DEBUG
  event_dictionary.pop("NWEvents")
  event_dictionary.pop("XSecMCweight")

def customize_DY(process, final_state_mode):
  for DYtype in ["DYGen", "DYLep", "DYJet"]:
    MC_dictionary[DYtype]["XSecMCweight"] = MC_dictionary[process]["XSecMCweight"]
    MC_dictionary[DYtype]["NWEvents"] = MC_dictionary[process]["NWEvents"]
  if (process == "DYIncNLO"): # double-check 
    # overwrite DYGen, DYLep, DYJet values with NLO values
    for subprocess in ["DYGen", "DYLep", "DYJet"]:
      MC_dictionary[subprocess]["XSec"]         = XSec["DYJetsToLL_M-50"]
      MC_dictionary[subprocess]["NWEvents"]     = MC_dictionary["DYIncNLO"]["NWEvents"]
      MC_dictionary[subprocess]["plot_scaling"] = 1  # override kfactor
  label_text = { "ditau" : r"$Z{\rightarrow}{\tau_h}{\tau_h}$",
                 "mutau" : r"$Z{\rightarrow}{\mu}{\tau_h}$",
                 "etau"  : r"$Z{\rightarrow}{e}{\tau_h}$",
                 "emu"   : r"$Z{\rightarrow}{e}{\mu}$",}
  MC_dictionary["DYGen"]["label"] = label_text[final_state_mode]


def append_Zpt_weight(event_dictionary):
  unpack_Zpt = [
    "nGenPart", "GenPart_pdgId", "GenPart_status", "GenPart_statusFlags",
    "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass",
  ]
  unpack_Zpt = (event_dictionary.get(key) for key in unpack_Zpt)
  Gen_Zpt, Gen_Z_mass, Gen_Zpt_weight = [], [], []

  # could make our own weights like this with a little effort
  # load 2D ROOT hist from local file
  from ROOT import TLorentzVector, TFile, TH2
  zptroot = TFile("SFs/zpt_reweighting_LO_2022.root", "open")
  zpthist = zptroot.Get("zptmass_histo")
  for nGen, pdgId, status, statusFlags, pt, eta, phi, mass in zip(*unpack_Zpt):
    good_lep_vecs = []
    for iparticle in range(nGen):
      pdgId_part  = abs(pdgId[iparticle])
      status_part = status[iparticle]
      flags_part  = statusFlags[iparticle]
      if ( ((pdgId_part==11 or pdgId_part==13) and status_part==1 and hasbit(flags_part, 8))
        or (pdgId_part==15 and status_part==2 and hasbit(flags_part, 8)) ): # 8 : fromHardProcess
        lep_vec = TLorentzVector() # surprisingly, you can't combine this with the following line
        lep_vec.SetPtEtaPhiM(pt[iparticle], eta[iparticle], phi[iparticle], mass[iparticle])
        good_lep_vecs.append(lep_vec)
    # end loop over particles in event
    #print(f"Z boson lep cands in event: {len(good_lep_vecs)}") # always 2
    zmass, zpt = 0.0, 0.0
    if (len(good_lep_vecs) == 2):
      zboson = good_lep_vecs[0] + good_lep_vecs[-1] # adding only cands in the list
      zmass = zboson.M()
      zpt   = zboson.Pt()

    zptweight = 1.0
    if not (zmass==0.0 and zpt==0.0):
      xbin = getBin(zmass, zpthist.GetXaxis())
      ybin = getBin(zpt, zpthist.GetYaxis())
      zptweight = zpthist.GetBinContent(xbin, ybin)
      if zptweight<=0.0: zptweight=1.0
    Gen_Zpt_weight.append(zptweight)

  event_dictionary["Weight_DY_Zpt_by_hand"] = np.array(Gen_Zpt_weight)
  return event_dictionary


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


def append_electron_flavor_indices(event_dictionary, final_state_mode):
  unpack_flav = ["l1_indices", "l2_indices", "Lepton_elIdx", "Lepton_tauIdx", "Electron_genPartFlav"]
  unpack_flav = (event_dictionary.get(key) for key in unpack_flav)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_flav]
  FS_ele_flav= []
  true_ele = []
  for i, l1_idx, l2_idx, ele_idx, tau_idx, ele_flav in zip(*to_check):
    ele = -1
    if (tau_idx[l1_idx] != -1 and ele_idx[l2_idx] != -1):
      ele = ele_flav[ele_idx[l2_idx]]
    elif (tau_idx[l2_idx] != -1 and ele_idx[l1_idx] != -1):
      ele = ele_flav[ele_idx[l1_idx]]
    else:
      print("Should not print :)")

    if (ele ==1 or ele==15): #prompt ele, or ele coming from tau
      FS_ele_flav.append(ele)
      true_ele.append(i)
    elif (ele!=1 or ele!=15):
      FS_ele_flav.append(ele)
    else:
      print(f"No gen matching for that final state ({final_state_mode}), crashing...")
      return None

  event_dictionary["FS_ele_flav"] = np.array(FS_ele_flav)
  event_dictionary["true_ele"]  = np.array(true_ele)

  return event_dictionary


def append_flavor_indices(event_dictionary, final_state_mode, keep_fakes=False):
  unpack_flav = ["l1_indices", "l2_indices", "Lepton_tauIdx", "Tau_genPartFlav"]
  unpack_flav = (event_dictionary.get(key) for key in unpack_flav)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_flav]
  FS_t1_flav, FS_t2_flav = [], []
  pass_gen_cuts, event_flavor = [], []
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
        event_flavor.append("G")
      elif (t1_flav == 0) or (t2_flav == 0):
        # jet fake --> one tau is faked by jet
        jet_fake = True
        event_flavor.append("J")
      elif (t1_flav < 5 and t1_flav > 0) or (t2_flav < 5 and t1_flav > 0):
        # lep fake --> both taus are faked by lepton
        # event with one tau faking jet enters category above first due to ordering
        # implies also the case where both are faked but one is faked by lepton 
        # is added to jet fakes, which i think is fine
        lep_fake = True
        event_flavor.append("L")
    elif ((final_state_mode == "mutau") or (final_state_mode == "etau")):
      t1_flav = tau_flav[tau_idx[l1_idx] + tau_idx[l2_idx] + 1] # update with NanoAODv12 samples -->is this ok?
      if (t1_flav == 5):
        genuine = True
        event_flavor.append("G")
      elif (t1_flav == 0):
        jet_fake = True
        event_flavor.append("J")
      elif (t1_flav < 5 and t1_flav > 0):
      #elif (t1_flav == 1):
        lep_fake = True
        #print("t1_flav:", t1_flav, "t2_flav:", t2_flav)
        event_flavor.append("L")

    else:
      print(f"No gen matching for that final state ({final_state_mode}), crashing...")
      return None
  
    if (keep_fakes==False) and ((genuine) or (lep_fake)):
      # save genuine background events and lep_fakes, remove jet fakes with gen matching
      # used in all categories because fakes are estimated with FF method
      FS_t1_flav.append(t1_flav)
      FS_t2_flav.append(t2_flav)
      pass_gen_cuts.append(i)

    if (keep_fakes==True) and ((genuine) or (lep_fake) or (jet_fake)):
      # save all events and their flavors, even if they are jet fakes
      # used to split DY to genuine, lep fakes, and jet fakes in all categories
      FS_t1_flav.append(t1_flav)
      FS_t2_flav.append(t2_flav)
      pass_gen_cuts.append(i)

  event_dictionary["FS_t1_flav"] = np.array(FS_t1_flav)
  event_dictionary["FS_t2_flav"] = np.array(FS_t2_flav)
  event_dictionary["pass_gen_cuts"] = np.array(pass_gen_cuts)
  event_dictionary["event_flavor"]  = np.array(event_flavor)
  return event_dictionary


def set_FF_values(final_state_mode, jet_mode_and_DeepTau_version, determining_FF):
  '''
  '''
  # should have aiso/iso as well
  FF_values = {
    # FS : { "jet_mode" : [intercept, slope] }  
    "etau"  : {

      "0j_2p5"     : [1, 1], 
      "1j_2p5"     : [1, 1],
      "GTE2j_2p5"  : [1, 1],

      #PostEE
      "custom_0j_2p5_FF" : [0.009, 0.0018], #fake factors
      "custom_1j_2p5_FF" : [0.033, 0.00097],
      "custom_GTE2j_2p5_FF" : [0.038, 0.00058],
      
      "custom_0j_2p5_mvis"    : [1.15, -0.0015], # bias correction as function of mvis
      "custom_1j_2p5_mvis"    : [1.15, -0.0015],
      "custom_GTE2j_2p5_mvis" : [1.15, -0.0015],
      
      "custom_0j_2p5_iso"    : [1.09, -0.4746], # bias correction as function of isolation
      "custom_1j_2p5_iso"    : [1.09, -0.4746],
      "custom_GTE2j_2p5_iso" : [1.09, -0.4746],
    
      "custom_0j_2p5_closure"    : [1.47, -0.0068], # closure correction as a function of pT
      "custom_1j_2p5_closure"    : [1.47, -0.0068],
      "custom_GTE2j_2p5_closure" : [1.47, -0.0068],
      
      #PreEE
      #"custom_0j_2p5_FF" : [-0.015, 0.0023], #fake factors
      #"custom_1j_2p5_FF" : [0.00034, 0.0016],
      #"custom_GTE2j_2p5_FF" : [0.019,  0.00099],

      #"custom_0j_2p5_mvis"    : [1.20, -0.0019], # bias correction as function of mvis
      #"custom_1j_2p5_mvis"    : [1.20, -0.0019],
      #"custom_GTE2j_2p5_mvis" : [1.20, -0.0019],

      #"custom_0j_2p5_iso"    : [1.04, -0.467], # bias correction as function of isolation
      #"custom_1j_2p5_iso"    : [1.04, -0.467],
      #"custom_GTE2j_2p5_iso" : [1.04, -0.467],

      #"custom_0j_2p5_closure"    : [1.67, -0.008], # closure correction as a function of pT
      #"custom_1j_2p5_closure"    : [1.67, -0.008],
      #"custom_GTE2j_2p5_closure" : [1.67, -0.008],

    },
  } 
  if (determining_FF == True):
    print("determining FF dictionary, setting dummy values and not plotting QCD")
    intercept, slope = 0.2, 1.1
  else:
    intercept = FF_values[final_state_mode][jet_mode_and_DeepTau_version][0]
    slope     = FF_values[final_state_mode][jet_mode_and_DeepTau_version][1]
    #slope2    = FF_values[final_state_mode][jet_mode_and_DeepTau_version][2]

  #return intercept, slope, slope2
  return intercept, slope


def add_FF_weights(event_dictionary, final_state_mode, jet_mode, DeepTau_version, determining_FF=False):
  unpack_FFVars = ["Lepton_pt", "HTT_m_vis", "HTT_Lep_iso", "l1_indices", "l2_indices", "Lepton_tauIdx", "Lepton_elIdx"]
  unpack_FFVars = (event_dictionary.get(key) for key in unpack_FFVars)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_FFVars]
  FF_weights = []
  bins = [40, 60, 80, 100, 120, 140, 160, 180, 999999]

  FF_mvis_weights = {
  # FS_mode : jet_mode : semilep_mode : mvis fake fraction [0, 10, 20, ..., 300]
  "etau" : {
    "0j" : {
      "WJ" : [ 0.0, 0.44948888, 0.48709053, 0.67524474, 0.50798577,
       0.40759256, 0.41038405, 0.4246403 , 0.42892607, 0.43369386],
      "QCD" : [0.183483, 0.044200, 0.292518, 0.221791, 0.216506, 0.252396, 0.221537, 0.179232, 0.331143]
    },
    "1j" : {
       "WJ" : [ 0.0, 0.35919842, 0.38588244, 0.43133431, 0.46609478,
       0.47824467, 0.56399762, 0.51112266, 0.40956644, 0.40627854,
       0.45590943, 0.62126431, 0.57938354, 0.54590267, 0.42468222],
       "QCD" : [0.269993, 0.110099, 0.110525, 0.181433, 0.158913, 0.244922, 0.117961, 0.333551, 0.307554]
    },
    "GTE2j" : {
       "WJ" : [ 0.0, 0.35919842, 0.38588244, 0.43133431, 0.46609478,
       0.47824467, 0.56399762, 0.51112266, 0.40956644, 0.40627854,
       0.45590943, 0.62126431, 0.57938354, 0.54590267, 0.42468222],
       "QCD" : [0.382559, 0.288983, 0.308179, 0.364518, 0.335701, 0.317249, 0.297538, 0.347137, 0.422986]
    },
  }
  }
  
  intercept, slope = set_FF_values(final_state_mode, "custom_"+jet_mode+"_2p5_FF", determining_FF)
  OSSS_bias_intercept, OSSS_bias_slope = set_FF_values(final_state_mode, "custom_"+jet_mode+"_2p5_mvis", determining_FF)
  closure_intercept, closure_slope = set_FF_values(final_state_mode, "custom_"+jet_mode+"_2p5_closure", determining_FF)
  iso_intercept, iso_slope = set_FF_values(final_state_mode, "custom_"+jet_mode+"_2p5_iso", determining_FF)

  one_minus_MC_over_data_weight = -999999
  for i, lep_pt, m_vis, lep_iso, l1_idx, l2_idx, tau_idx, el_idx in zip(*to_check):

    if m_vis < bins[0]: # 40
       one_minus_MC_over_data_weight = FF_mvis_weights[final_state_mode][jet_mode]['QCD'][0]
    else:
       for bin_index in range(1,9):
          if (m_vis > bins[bin_index-1] and m_vis < bins[bin_index]):
             one_minus_MC_over_data_weight = FF_mvis_weights[final_state_mode][jet_mode]['QCD'][bin_index]
             break
          else: continue

    if (tau_idx[l1_idx] != -1 and el_idx[l2_idx] != -1):
      tauFSLoc = l1_idx
      elFSLoc  = l2_idx
    elif (tau_idx[l2_idx] != -1 and el_idx[l1_idx] != -1):
      tauFSLoc = l2_idx
      elFSLoc  = l1_idx

    l1_pt = lep_pt[tauFSLoc] if lep_pt[tauFSLoc] < 100.0 else 100.0
    l2_pt = lep_pt[elFSLoc] if lep_pt[elFSLoc] < 150.0 else 150.0
    #l1_pt = lep_pt[l1_idx] if lep_pt[l1_idx] < 120.0 else 120.0
    #l2_pt = lep_pt[l2_idx] if lep_pt[l2_idx] < 200.0 else 200.0
    m_vis = m_vis if m_vis < 350.0 else 350.0

    if (l1_pt) < 35:
       FF_weight = 0.01  + l1_pt * 0.00175591 
    else:
       FF_weight = (intercept + l1_pt * slope)

    #FF_weight *= one_minus_MC_over_data_weight*(intercept + l1_pt * slope)
    FF_weight *= (OSSS_bias_intercept + m_vis * OSSS_bias_slope)*(closure_intercept + l2_pt * closure_slope)*(iso_intercept + lep_iso * iso_slope)

    FF_weights.append(FF_weight)
  event_dictionary["FF_weight"] = np.array(FF_weights)
  return event_dictionary


def make_jet_cut(event_dictionary, jet_mode):
  nEvents_precut = len(event_dictionary["Lepton_pt"])
  unpack_jetVars = ["nCleanJet", "CleanJet_pt", "CleanJet_eta", "CleanJet_phi", "CleanJet_mass"]
  unpack_jetVars = (event_dictionary.get(key) for key in unpack_jetVars)
  to_check = [range(len(event_dictionary["Lepton_pt"])), *unpack_jetVars] # "*" unpacks a tuple
  nCleanJetGT30, pass_0j_cuts, pass_1j_cuts, pass_2j_cuts, pass_3j_cuts, pass_GTE2j_cuts = [], [], [], [], [], []
  CleanJetGT30_pt_1, CleanJetGT30_pt_2, CleanJetGT30_pt_3    = [], [], []
  CleanJetGT30_eta_1, CleanJetGT30_eta_2, CleanJetGT30_eta_3 = [], [], []
  CleanJetGT30_phi_1, CleanJetGT30_phi_2, CleanJetGT30_phi_3 = [], [], []
  mjj_array, detajj_array = [], []
  from ROOT import TLorentzVector 
  for i, nJet, jet_pt, jet_eta, jet_phi, jet_mass in zip(*to_check):
    passingJets = 0
    passingJetsPt, passingJetsEta, passingJetsPhi, passingJetsMass = [], [], [], []
    for ijet in range(0, nJet):
      if (jet_pt[ijet] > 30.0) and (jet_eta[ijet] < 4.7):
        passingJets += 1
        passingJetsPt.append(jet_pt[ijet])
        passingJetsEta.append(jet_eta[ijet])
        passingJetsPhi.append(jet_phi[ijet])
        passingJetsMass.append(jet_mass[ijet])
    nCleanJetGT30.append(passingJets)

    if passingJets == 0: 
      pass_0j_cuts.append(i)

    if (passingJets == 1) and (jet_mode == "Inclusive" or jet_mode == "1j"): 
      pass_1j_cuts.append(i)
      CleanJetGT30_pt_1.append(passingJetsPt[0])
      CleanJetGT30_eta_1.append(passingJetsEta[0])
      CleanJetGT30_phi_1.append(passingJetsPhi[0])

    if (passingJets == 2) and (jet_mode == "Inclusive" or jet_mode == "2j"):
      pass_2j_cuts.append(i)
      CleanJetGT30_pt_1.append(passingJetsPt[0])
      CleanJetGT30_pt_2.append(passingJetsPt[1])
      CleanJetGT30_eta_1.append(passingJetsEta[0])
      CleanJetGT30_eta_2.append(passingJetsEta[1])
      CleanJetGT30_phi_1.append(passingJetsPhi[0])
      CleanJetGT30_phi_2.append(passingJetsPhi[1])
      j1_vec, j2_vec = TLorentzVector(), TLorentzVector() # surprisingly, you can't combine this with the following line
      j1_vec.SetPtEtaPhiM(passingJetsPt[0], passingJetsEta[0], passingJetsPhi[0], passingJetsMass[0])
      j2_vec.SetPtEtaPhiM(passingJetsPt[1], passingJetsEta[1], passingJetsPhi[1], passingJetsMass[1])
      mjj_array.append((j1_vec + j2_vec).M())
      detajj_array.append(abs(j1_vec.Eta() - j2_vec.Eta()))

    if (passingJets >= 2) and (jet_mode == "GTE2j"): 
      pass_GTE2j_cuts.append(i)
      TLorentzVector_Jets = []
      for i in range(passingJets):
        temp_jet_vec = TLorentzVector()
        temp_jet_vec.SetPtEtaPhiM(passingJetsPt[i], passingJetsEta[i], passingJetsPhi[i], passingJetsMass[i])
        TLorentzVector_Jets.append(temp_jet_vec)
      j1_TVec, j2_TVec = highest_mjj_pair(TLorentzVector_Jets)
      CleanJetGT30_pt_1.append(passingJetsPt[0])
      CleanJetGT30_pt_2.append(passingJetsPt[1])
      CleanJetGT30_eta_1.append(passingJetsEta[0])
      CleanJetGT30_eta_2.append(passingJetsEta[1])
      CleanJetGT30_phi_1.append(passingJetsPhi[0])
      CleanJetGT30_phi_2.append(passingJetsPhi[1])
      mjj_array.append((j1_TVec+j2_TVec).M())
      detajj_array.append(abs(j1_TVec.Eta()-j2_TVec.Eta()))

  event_dictionary["nCleanJetGT30"]   = np.array(nCleanJetGT30)

  if jet_mode == "pass":
    print("debug jet mode, only filling nCleanJetGT30")

  elif jet_mode == "Inclusive":
    pass
    # fill branches like above
    #event_dictionary["pass_0j_cuts"]    = np.array(pass_0j_cuts)
    #event_dictionary["pass_1j_cuts"]    = np.array(pass_1j_cuts)
    #event_dictionary["pass_2j_cuts"]    = np.array(pass_2j_cuts)
    #event_dictionary["pass_3j_cuts"]    = np.array(pass_3j_cuts)
    #event_dictionary["pass_GTE2j_cuts"] = np.array(pass_GTE2j_cuts)

    #event_dictionary["CleanJetGT30_pt_1"]  = np.array(CleanJetGT30_pt_1)
    #event_dictionary["CleanJetGT30_pt_2"]  = np.array(CleanJetGT30_pt_2)
    #event_dictionary["CleanJetGT30_pt_3"]  = np.array(CleanJetGT30_pt_3)
    #event_dictionary["CleanJetGT30_eta_1"] = np.array(CleanJetGT30_eta_1)
    #event_dictionary["CleanJetGT30_eta_2"] = np.array(CleanJetGT30_eta_2)
    #event_dictionary["CleanJetGT30_eta_3"] = np.array(CleanJetGT30_eta_3)
  
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
    event_dictionary["FS_mjj"] = np.array(mjj_array)
    event_dictionary["FS_detajj"] = np.array(detajj_array)

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
    event_dictionary["CleanJetGT30_phi_1"] = np.array(CleanJetGT30_phi_1)
    event_dictionary["CleanJetGT30_phi_2"] = np.array(CleanJetGT30_phi_2)
    #event_dictionary["CleanJetGT30_phi_3"] = np.array(CleanJetGT30_phi_3)
    event_dictionary["FS_mjj"] = np.array(mjj_array)
    event_dictionary["FS_detajj"] = np.array(detajj_array)

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
  DEBUG = False # set this to true to show print output from this function
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
    print(text_options["red"] + "ALL EVENTS REMOVED! SAMPLE WILL BE DELETED! " + text_options["reset"])
    delete_sample = True
    return None
 
  if DEBUG: print(f"cut branch: {cut_branch}")
  if DEBUG: print(f"protected branches: {protected_branches}")
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
      if DEBUG: print(f"going to cut {branch}, {len(event_dictionary[branch])}")
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


def apply_final_state_cut(event_dictionary, final_state_mode, DeepTau_version, useMiniIso=False, GENMATCHING=False):
  '''
  Organizational function that generalizes call to a (set of) cuts based on the
  final cut. Importantly, the function that rejects events, 'apply_cut',
  is called elsewhere
  '''
  # setting inclusive in the jet_mode includes all jet branches in protected branches
  # this is okay because in the current ordering (FS cut then jet cut), no jet branches
  # are event created yet.
  #if (final_state_mode == "mutau_TnP"):
  #  protected_branches = set_protected_branches(final_state_mode="mutau_TnP", jet_mode="Inclusive")
  #else:
  protected_branches = set_protected_branches(final_state_mode=final_state_mode, jet_mode="Inclusive")
  if final_state_mode == "ditau":
    event_dictionary = make_ditau_cut(event_dictionary, DeepTau_version)
    event_dictionary = apply_cut(event_dictionary, "pass_cuts", protected_branches)
  elif final_state_mode == "mutau":
    event_dictionary = make_mutau_cut(event_dictionary, DeepTau_version)
    event_dictionary = apply_cut(event_dictionary, "pass_cuts", protected_branches)
  elif final_state_mode == "mutau_TnP": # special mode for Tau TRG studies
    event_dictionary = make_mutau_TnP_cut(event_dictionary, DeepTau_version)
    event_dictionary = apply_cut(event_dictionary, "pass_cuts", protected_branches)
  elif final_state_mode == "etau":
    if (GENMATCHING==False):
      event_dictionary = make_etau_cut(event_dictionary, DeepTau_version)
      event_dictionary = apply_cut(event_dictionary, "pass_cuts", protected_branches)
    else: print("Customized event selection")
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


def apply_flavor_cut(event_dictionary):
  # get list of event indices with events matching flavor key
  event_flavor_array = event_dictionary["Cuts"]["event_flavor"]
  # cut out other events
  event_dictionary = apply_cut(event_dictionary, "pass_flav_cut") # no protected branches
  return event_dictionary


def apply_AR_cut(process, event_dictionary, final_state_mode, jet_mode, DeepTau_version, determining_FF):
  '''
  Organizational function
  added 'skip_DeepTau' to apply a partial selection (all but leading tau deeptau reqs)
  '''
  protected_branches = ["None"]
  event_dictionary = append_lepton_indices(event_dictionary)
  if ("Data" not in process):
    load_and_store_NWEvents(process, event_dictionary)
    if ("DY" in process): customize_DY(process, final_state_mode)
  if (final_state_mode != "dimuon"):
    # non-standard FS cut
    if (final_state_mode == "ditau"):
      event_dictionary = make_ditau_AR_cut(event_dictionary, DeepTau_version)
      event_dictionary = apply_cut(event_dictionary, "pass_AR_cuts", protected_branches)
      event_dictionary = apply_jet_cut(event_dictionary, jet_mode)
      event_dictionary = make_ditau_cut(event_dictionary, DeepTau_version, skip_DeepTau=True)
    if (final_state_mode == "mutau"):
      event_dictionary = make_mutau_AR_cut(event_dictionary, DeepTau_version)
      event_dictionary = apply_cut(event_dictionary, "pass_AR_cuts", protected_branches)
      event_dictionary = apply_jet_cut(event_dictionary, jet_mode)
      event_dictionary = make_mutau_cut(event_dictionary, DeepTau_version, skip_DeepTau=True)
    if (final_state_mode == "etau"):
      event_dictionary = make_etau_AR_cut(event_dictionary, DeepTau_version)
      event_dictionary = apply_cut(event_dictionary, "pass_AR_cuts", protected_branches)
      event_dictionary = apply_jet_cut(event_dictionary, jet_mode)
      event_dictionary = make_etau_cut(event_dictionary, DeepTau_version, skip_DeepTau=True)
    protected_branches = set_protected_branches(final_state_mode=final_state_mode, jet_mode="none")
    event_dictionary   = apply_cut(event_dictionary, "pass_cuts", protected_branches)
    event_dictionary = add_FF_weights(event_dictionary, final_state_mode, jet_mode, DeepTau_version,
                                      determining_FF=determining_FF)
    determining_FF = False # paranoid thing to prevent leaking the variable
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
    "pass"      : "Inclusive", #DEBUG
    "0j" : "pass_0j_cuts",
    "1j" : "pass_1j_cuts",
    "2j" : "pass_2j_cuts",
    "3j" : "pass_3j_cuts",
    "GTE2j" : "pass_GTE2j_cuts",
  }
  event_dictionary   = make_jet_cut(event_dictionary, jet_mode)
  protected_branches = set_protected_branches(final_state_mode="none", jet_mode=jet_mode)
  if jet_mode == "Inclusive" or jet_mode == "pass":
    print("jet mode is Inclusive, no jet cut performed")
  else:
    event_dictionary = apply_cut(event_dictionary, jet_cut_branch[jet_mode], protected_branches)
  return event_dictionary


def apply_HTT_FS_cuts_to_process(process, process_dictionary, log_file,
                                 final_state_mode, jet_mode="Inclusive", 
                                 DeepTau_version="2p5", useMiniIso=False):
  '''
  Organizational function to hold two function calls and empty list handling that
  is performed for all loaded datasets in our framework.
  Can be extended to hold additional standard cuts (i.e. jets) or the returned
  value can be cut on as needed.
  '''
  log_print(f"Processing {process}", log_file)
  process_events = process_dictionary[process]["info"]
  if len(process_events["run"])==0: return None

  process_events = append_lepton_indices(process_events)
  protected_branches = ["FS_t1_flav", "FS_t2_flav","pass_gen_cuts", "event_flavor"]

  if ("Data" not in process):
    load_and_store_NWEvents(process, process_events)
    if ("DY" in process): customize_DY(process, final_state_mode) #append_Zpt_weight(process_events)
    #keep_fakes = True  #remove after SF check
    #process_events = append_flavor_indices(process_events, final_state_mode, keep_fakes=keep_fakes)
    #process_events = apply_cut(process_events, "pass_gen_cuts", protected_branches=protected_branches)
    #if (process_events==None or len(process_events["run"])==0): return None

    keep_fakes = False
    if ((("TT" in process) or ("WJ" in process) or ("DY" in process)) and (final_state_mode=="mutau" or final_state_mode=="etau")):
      # when FF method is finished/improved no longer need to keep TT and WJ fakes
      keep_fakes = True
    if (("DY" in process) and (final_state_mode=="ditau")):
      keep_fakes = True
    process_events = append_flavor_indices(process_events, final_state_mode, keep_fakes=keep_fakes)
    process_events = apply_cut(process_events, "pass_gen_cuts", protected_branches=protected_branches)
    if ((final_state_mode=="etau") and ("DY" in process)):
      process_events = append_electron_flavor_indices(process_events, final_state_mode)
      process_events = apply_cut(process_events, "true_ele")
    if (process_events==None or len(process_events["run"])==0): return None

  FS_cut_events = apply_final_state_cut(process_events, final_state_mode, DeepTau_version, useMiniIso=useMiniIso)
  if (FS_cut_events==None or len(FS_cut_events["run"])==0): return None 
  cut_events = apply_jet_cut(FS_cut_events, jet_mode)
  if (cut_events==None or len(cut_events["run"])==0): return None

  # TODO : want to move to this
  # re TODO actually want to move to a splitting/copying paradigm instead of modification in place
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

  good_events = "(HTT_SRevent) & (METfilters) & (LeptonVeto==0) & (JetMapVeto_EE_30GeV) & (JetMapVeto_HotCold_30GeV)"

  if final_state_mode == "ditau":
    triggers = "(HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1\
               | HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60\
               | HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75\
               | HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1\
               | HLT_DoublePFJets40_Mass500_MediumDeepTauPFTauHPS45_L2NN_MediumDeepTauPFTauHPS20_eta2p1)"
    #triggers = "(HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1)"

    good_events += " & (abs(HTT_pdgId)==15*15) & " + triggers
    if disable_triggers: good_events = good_events.replace(" & (Trigger_ditau)", "")

  elif final_state_mode == "mutau":
    good_events += " & (abs(HTT_pdgId)==13*15) & (Trigger_mutau)"
    if disable_triggers: good_events = good_events.replace(" & (Trigger_mutau)", "")

  elif final_state_mode == "etau":
    #good_events += " & (abs(HTT_pdgId)==11*15) & (HLT_Ele30_WPTight_Gsf | HLT_Ele32_WPTight_Gsf | HLT_Ele35_WPTight_Gsf) & (HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 == 0)"
    good_events += " & (abs(HTT_pdgId)==11*15) & (HLT_Ele30_WPTight_Gsf | HLT_Ele32_WPTight_Gsf | HLT_Ele35_WPTight_Gsf | HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1)"
    if disable_triggers: good_events = good_events.replace(" & (Trigger_etau)", "")

  # non-HTT FS modes
  elif final_state_mode == "mutau_TnP":
    good_events = "(METfilters) & (abs(HTT_pdgId)==13*15)"

  elif final_state_mode == "dimuon":
    # lepton veto must be applied manually for this final state
    if (useMiniIso == False):
      good_events = "(METfilters) & (HTT_pdgId==-13*13) & (HLT_IsoMu24)"
    if (useMiniIso == True):
      good_events = "(METfilters) & (LeptonVeto==0) & (HTT_pdgId==-13*13) & (HLT_IsoMu24)"
    if disable_triggers: good_events = good_events.replace(" & (HLT_IsoMu24)", "")

  return good_events


def set_branches(final_state_mode, DeepTau_version, process="None"):
  common_branches = [
    "run", "luminosityBlock", "event", "Generator_weight", "NWEvents", "XSecMCweight",
    "TauSFweight", "MuSFweight", "ElSFweight","ElSFweight_Trig", "ElSFweight_Reco", "ElSFweight_ID", "Weight_DY_Zpt_LO", "Weight_DY_Zpt_NLO", "PUweight", "Weight_TTbar_NNLO",
    "FSLeptons", "Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_iso",
    "Electron_genPartFlav",
    "Tau_genPartFlav", "Tau_decayMode",
    "nCleanJet", "CleanJet_pt", "CleanJet_eta", "CleanJet_phi", "CleanJet_mass",
    "HTT_m_vis", "HTT_dR", "HTT_pT_l1l2", "FastMTT_PUPPIMET_mT", "FastMTT_PUPPIMET_mass",
    #"Tau_rawPNetVSjet", "Tau_rawPNetVSmu", "Tau_rawPNetVSe",
    "PV_npvs", "Pileup_nPU",
    "Electron_mvaNoIso_WP90", "Electron_pfRelIso03_all",
    "Electron_hoe", "Electron_lostHits", "Electron_mass", "Electron_r9", "Electron_sieie", "Electron_pfRelIso03_chg",
    "Electron_dr03TkSumPt", "Electron_dr03EcalRecHitSumEt", "Electron_dr03HcalDepth1TowerSumEt"
  ]
  branches = common_branches
  branches = add_final_state_branches(branches, final_state_mode)
  if final_state_mode != "dimuon": branches = add_DeepTau_branches(branches, DeepTau_version)
  branches = add_trigger_branches(branches, final_state_mode)

  if (process == "DY"): branches = add_Zpt_branches(branches)
  return branches


def add_final_state_branches(branches_, final_state_mode):
  '''
  Helper function to add only relevant branches to loaded branches based on final state.
  '''
  final_state_branches = {
    "ditau"  : ["Lepton_tauIdx", "Tau_dxy", "Tau_dz", "Tau_charge", "PuppiMET_pt"],

    "mutau"  : ["Muon_dxy", "Muon_dz", "Muon_charge",
                "Tau_dxy", "Tau_dz", "Tau_charge",
                "Lepton_tauIdx", "Lepton_muIdx",
                "PuppiMET_pt", "PuppiMET_phi", "HTT_Lep_iso"],

    "mutau_TnP"  : ["Muon_dxy", "Muon_dz", "Muon_charge",
                "Tau_dxy", "Tau_dz", "Tau_charge",
                "Lepton_tauIdx", "Lepton_muIdx",
                "PuppiMET_pt", "PuppiMET_phi"],

    "etau"   : ["Electron_dxy", "Electron_dz", "Electron_charge", 
                "Tau_dxy", "Tau_dz", "Tau_charge", 
                "Lepton_tauIdx", "Lepton_elIdx",
                "PuppiMET_pt", "PuppiMET_phi","HTT_Lep_iso"],

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
    "GTE2j" : ["nCleanJetGT30", 
               "CleanJetGT30_pt_1", "CleanJetGT30_eta_1", "CleanJetGT30_phi_1",
               "CleanJetGT30_pt_2", "CleanJetGT30_eta_2", "CleanJetGT30_phi_2",
               "FS_mjj", "FS_detajj",
              ],
}

final_state_vars = {
    # can't put nanoaod branches here because this dictionary is used to protect branches created internally
    "none"   : [],
    "ditau"  : ["FS_t1_pt", "FS_t1_eta", "FS_t1_phi", "FS_t1_dxy", "FS_t1_dz", "FS_t1_chg", "FS_t1_DM",
                "FS_t2_pt", "FS_t2_eta", "FS_t2_phi", "FS_t2_dxy", "FS_t2_dz", "FS_t2_chg", "FS_t2_DM",
                "FS_t1_flav", "FS_t2_flav", 
                #"FS_t1_rawPNetVSjet", "FS_t1_rawPNetVSmu", "FS_t1_rawPNetVSe",
                #"FS_t2_rawPNetVSjet", "FS_t2_rawPNetVSmu", "FS_t2_rawPNetVSe",
                "FS_t1_DeepTauVSjet", "FS_t1_DeepTauVSmu", "FS_t1_DeepTauVSe", 
                "FS_t2_DeepTauVSjet", "FS_t2_DeepTauVSmu", "FS_t2_DeepTauVSe", 
                ],

    "mutau"  : ["FS_mu_pt", "FS_mu_eta", "FS_mu_phi", "FS_mu_iso", "FS_mu_dxy", "FS_mu_dz", "FS_mu_chg",
                "FS_tau_pt", "FS_tau_eta", "FS_tau_phi", "FS_tau_dxy", "FS_tau_dz", "FS_tau_chg", "FS_tau_DM",
                "FS_mt", "FS_t1_flav", "FS_t2_flav", 
                #"FS_tau_rawPNetVSjet", "FS_tau_rawPNetVSmu", "FS_tau_rawPNetVSe"
                ],

    "mutau_TnP"  : ["FS_mu_pt", "FS_mu_eta", "FS_mu_phi", "FS_mu_iso", "FS_mu_dxy", "FS_mu_dz", "FS_mu_chg",
                "FS_tau_pt", "FS_tau_eta", "FS_tau_phi", "FS_tau_dxy", "FS_tau_dz", "FS_tau_chg",
                "FS_mt", "FS_t1_flav", "FS_t2_flav", "pass_tag", "pass_probe"],

    "etau"   : ["FS_el_pt", "FS_el_eta", "FS_el_phi", "FS_el_iso", "FS_el_dxy", "FS_el_dz", "FS_el_chg",
                "FS_tau_pt", "FS_tau_eta", "FS_tau_phi", "FS_tau_dxy", "FS_tau_dz", "FS_tau_chg",
                "FS_mt", "FS_t1_flav", "FS_t2_flav", "xtrigger_flag",
                "FS_el_hoe", "FS_el_lostHits", "FS_el_mass", "FS_el_r9", "FS_el_sieie", "FS_el_isochg", 
                "FS_el_dr03TkSumPt", "FS_el_dr03EcalRecHitSumEt", "FS_el_dr03HcalDepth1TowerSumEt"],

    "dimuon" : ["FS_m1_pt", "FS_m1_eta", "FS_m1_phi", "FS_m1_iso", "FS_m1_dxy", "FS_m1_dz",
                "FS_m2_pt", "FS_m2_eta", "FS_m2_phi", "FS_m2_iso", "FS_m2_dxy", "FS_m2_dz"],
}

def set_vars_to_plot(final_state_mode, jet_mode="none"):
  '''
  Helper function to keep plotting variables organized
  Shouldn't this be in  plotting functions?
  '''
  vars_to_plot = ["HTT_m_vis", "HTT_dR", "HTT_Lep_iso"]
                  #"HTT_DiJet_MassInv_fromHighestMjj", "HTT_DiJet_dEta_fromHighestMjj"] 
                  # common to all final states # TODO add MET here, add Tau_decayMode
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
    if jet_mode == "Inclusive" or jet_mode=="pass": # cutting FS branches, but not the jet branches
      jet_mode = "Inclusive"
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



