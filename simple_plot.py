# I'll make my own plotter, with blackjack and hookers!
# Authored by Braden Allmond, Sep 11, 2023
 
# libraries
import uproot
import numpy as np
import sys
from os import getlogin, path
import matplotlib.pyplot as plt
from datetime import datetime, timezone

# user files
from binning_dictionary   import binning_dictionary
from MC_dictionary        import MC_dictionary
from triggers_dictionary  import triggers_dictionary

from file_map      import testing_file_map, full_file_map
from text_options  import text_options

# fb ^ -1
luminosities = {
  "2022 F"   : 17.61,
  "2022 G"   :  3.06,
  "2022 F&G" : 20.67, # calculate by hand when new quantities are relevant
}

def get_midpoints(input_bins):
  midpoints = []
  for i, ibin in enumerate(input_bins):
    if (i+1 != len(input_bins)):
      midpoints.append( ibin + (input_bins[i+1] - ibin)/2 )
  midpoints = np.array(midpoints)
  return midpoints


def setup_ratio_plot():
  gs = gridspec_kw = {'height_ratios': [4, 1], 'hspace': 0}
  fig, (upper_ax, lower_ax) = plt.subplots(nrows=2, sharex=True, gridspec_kw=gridspec_kw)
  return (upper_ax, lower_ax)


def make_bins(variable_name):
  nbins, xmin, xmax = binning_dictionary[variable_name]
  check_uniformity = (xmax-xmin)/nbins
  if (check_uniformity % 1 != 0 and check_uniformity % 0.1 != 0):
    print(f"nbins, xmin, xmax : {nbins}, {xmin}, {xmax}")
    print(f"(xmax-xmin)/nbins = {check_uniformity}, results in bad bin edges and centers")
  xbins = np.linspace(xmin, xmax, nbins)
  return xbins


def plot_data(histogram_axis, xbins, data_info, luminosity):
  # TODO: understand errors for poisson statistics
  # TODO: understand all of statistics
  #stat_error = np.array([1000/np.sqrt(entry) if entry > 0 else 0 for entry in data_info]) #wrong but scaled up...
  #stat_error  = sum_of_data**2 * np.ones(np.shape(data_info)) #wrong
  sum_of_data = np.sum(data_info)
  stat_error  = np.array([np.sqrt(entry * (1 - entry/sum_of_data)) if entry > 0 else 0 for entry in data_info]) # 
  #stat_error = np.array([np.sqrt(entry) if entry > 0 else 0 for entry in data_info]) # unsure if correct?
  midpoints   = get_midpoints(xbins) # len(midpoints) = len(xbins) - 1
  label = f"Data [{np.sum(data_info):>.0f}]"
  histogram_axis.errorbar(midpoints, data_info, yerr=stat_error, 
                          color="black", marker="o", linestyle='none', markersize=2, label=label)
  #histogram_axis.plot(midpoints, data_info, color="black", marker="o", linestyle='none', markersize=2, label="Data")
  # above plots without error bars


def add_CMS_preliminary(axis):
  CMS_text = "CMS"
  axis.text(0.01, 1.02, CMS_text, transform=hist_ax.transAxes, fontsize=16, weight='bold')
  preliminary_text = "Preliminary"
  axis.text(0.12, 1.02, preliminary_text, transform=hist_ax.transAxes, fontsize=16, style='italic')


def spruce_up_plot(histogram_axis, ratio_plot_axis, variable_name, luminosity):
  add_CMS_preliminary(histogram_axis)
  # reverse dictionary search to get correct era title from luminosity
  title = [key for key in luminosities.items() if key[1] == luminosity][0][0]
  title_string = f"{title}, {luminosity}" + r"$fb^{-1}$"
  histogram_axis.set_title(title_string, loc='right')
  histogram_axis.set_ylabel("Events / Bin")
  histogram_axis.minorticks_on()

  yticks = histogram_axis.yaxis.get_major_ticks()
  yticks[0].label1.set_visible(False) # hides a zero that overlaps with the upper plot
  plt.minorticks_on()

  ratio_plot_axis.set_ylim([0.6, 1.4]) # 0.0, 2.0 also make sense
  ratio_plot_axis.set_xlabel(variable_name) # shared axis label
  ratio_plot_axis.set_ylabel("Obs. / Exp.")
  ratio_plot_axis.axhline(y=1, color='grey', linestyle='--')
  

# TODO : delete after commiting
#def update_legend_labels(histogram_axis, yields, display_yields=True):
  # TODO set yields in plotters individually
  #handles, old_labels = hist_ax.get_legend_handles_labels()
  #new_labels = []
  # reorganize artists and labels by hand
  #handles.append(handles[0])
  #handles.pop(0)
  #old_labels.append(old_labels[0])
  #old_labels.pop(0)
  #
  #for i, label in enumerate(old_labels):
  #  label += f" [{yields[i]:>.0f}]"
  #  new_labels.append(label)
  #if display_yields:
  #  hist_ax.legend(handles, new_labels)
  #else:
  #  hist_ax.legend(handles, old_labels)
  #  print(f"labels and yields : {new_labels}") # put labels in terminal if not on plot
#  pass


def set_MC_process_info(process, luminosity, scaling=False, signal=False):
  color = MC_dictionary[process]["color"]
  label = MC_dictionary[process]["label"]
  if scaling:
  # factor of 1000 comes from lumi and XSec units of fb^-1 = 10E15 b^-1 and pb = 10E-12 b respectively
    plot_scaling = MC_dictionary[process]["plot_scaling"] # 1 for all non-signal processes by default
    scaling = 1000. * plot_scaling * luminosity * MC_dictionary[process]["XSec"] / MC_dictionary[process]["NWevents"]
  if signal:
    label += " x" + str(plot_scaling)
  return (color, label, scaling)


def calculate_underoverflow(events, xbins, weights):
  count_bin_values = [-999999., xbins[0], xbins[-1], 999999.]
  values, bins = np.histogram(events, count_bin_values, weights=weights)
  underflow_value, overflow_value = values[0], values[-1]
  return underflow_value, overflow_value


def get_binned_info(process_name, process_variable, xbins, process_weights, luminosity):
  scaling = 1 if "Data" in process_name else set_MC_process_info(process_name, luminosity, scaling=True)[2]
  weights = scaling*process_weights
  underflow, overflow = calculate_underoverflow(process_variable, xbins, weights)
  binned_values, _    = np.histogram(process_variable, xbins, weights=weights)
  binned_values[0]   += underflow
  binned_values[-1]  += overflow
  return binned_values
 

def accumulate_MC_subprocesses(parent_process, process_dictionary):
  accumulated_values = 0
  for MC_process in process_dictionary:
    if get_parent_process(MC_process) == parent_process:
      accumulated_values += process_dictionary[MC_process]["BinnedEvents"]
  return accumulated_values

  
def plot_signal(histogram_axis, xbins, signal_dictionary, luminosity):
  for signal in signal_dictionary:
    color, label, _ = set_MC_process_info(signal, luminosity, scaling=True, signal=True)
    current_hist = signal_dictionary[signal]["BinnedEvents"]
    label += f" [{np.sum(current_hist):>.0f}]"
    stairs = histogram_axis.stairs(current_hist, xbins, color=color, label=label, fill=False)


def plot_MC(histogram_axis, xbins, stack_dictionary, luminosity):
  bottom = 0
  for MC_process in stack_dictionary:
    color, label, _ = set_MC_process_info(MC_process, luminosity)
    current_hist = stack_dictionary[MC_process]["BinnedEvents"]
    # TODO handle variable binning with list of differences
    label += f" [{np.sum(current_hist):>.0f}]"
    bars = histogram_axis.bar(xbins[0:-1], current_hist, width=xbins[1]-xbins[0],
                            color=color, label=label, 
                            bottom=bottom, fill=True, align='edge')
    bottom += current_hist # Stack, set top of current graph as bottom of next graph

 
def calculate_yields(data, backgrounds, signals):
  yields = [np.sum(data)]
  total_background, total_signal = 0, 0
  for process in backgrounds: 
    process_yield = np.sum(backgrounds[process]["BinnedEvents"])
    yields.append(process_yield)
    total_background += process_yield
  for signal in signals:
    signal_yield = np.sum(backgrounds[process]["BinnedEvents"])
    yields.append(signal_yield)
    total_signal += signal_yield

  print("signal-to-background information")
  print(f"S/B      : {total_signal/total_background:.3f}")
  print(f"S/(S+B)  : {total_signal/(total_signal+total_background):.3f}")
  print(f"S/√(B)   : {total_signal/np.sqrt(total_background):.3f}")
  print(f"S/√(S+B) : {total_signal/np.sqrt(total_signal+total_background):.3f}")
  
  return yields


def make_ratio_plot(ratio_axis, xbins, numerator_data, denominator_data):
  # TODO fix error calculation, should be np.sqrt(entry / sum? )
  numerator_statistical_error   = [1/np.sqrt(entry) if entry > 0 else 0 for entry in numerator_data]
  denominator_statistical_error = [1/np.sqrt(entry) if entry > 0 else 0 for entry in denominator_data]
  combined_statistical_error    = np.add(numerator_statistical_error, denominator_statistical_error)

  midpoints = get_midpoints(xbins)
  ratio_axis.errorbar(midpoints, numerator_data/denominator_data, yerr=combined_statistical_error,
                    color="black", marker="o", linestyle='none', markersize=2)


def append_lepton_indices(event_dictionary):
  FSLeptons = event_dictionary["FSLeptons"]
  l1_indices, l2_indices = [], []
  for event in FSLeptons:
    l1_indices.append(event[0])
    l2_indices.append(event[1])
  event_dictionary["l1_indices"] = np.array(l1_indices)
  event_dictionary["l2_indices"] = np.array(l2_indices)
  return event_dictionary


def make_ditau_cut(event_dictionary, DeepTauVersion):
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
  unpack_dimuon = ["Lepton_pt", "Lepton_iso", "HTT_m_vis", "HTT_dR", "l1_indices", "l2_indices"]
  unpack_dimuon = (event_dictionary.get(key) for key in unpack_dimuon)
  to_check      = [range(len(event_dictionary["Lepton_pt"])), *unpack_dimuon]
  pass_cuts = []
  for i, pt, iso, mvis, dR, l1_idx, l2_idx in zip(*to_check):
    passKinematics = (pt[l1_idx] > 26 and pt[l2_idx] > 20 and mvis > 20 and dR > 0.5)
    passIso        = (iso[l1_idx] < 0.15 and iso[l2_idx] < 0.15)
    if (passKinematics and passIso):
      pass_cuts.append(i)

  event_dictionary["pass_cuts"] = np.array(pass_cuts)
  print(f"events passing dimuon cuts = {len(np.array(pass_cuts))}")
  return event_dictionary


def calculate_mt(lep_pt, lep_phi, MET_pt, MET_phi):
  # useful also for etau, emu
  delta_phi = phi_mpi_pi(lep_phi - MET_phi)
  mt = np.sqrt(2 * lep_pt * MET_pt * (1 - np.cos(delta_phi) ) ) 
  #sum_pt_2  = (lep_pt + MET_pt)**2
  #sum_ptx_2 = (lep_pt*np.cos(lep_phi) + MET_pt*np.cos(MET_phi))**2
  #sum_pty_2 = (lep_pt*np.sin(lep_phi) + MET_pt*np.sin(MET_phi))**2
  #mt = sum_pt_2 - sum_ptx_2 - sum_pty_2 # alternate calculation, same quantity
  return mt
  

def make_run_cut(event_dictionary, good_runs):
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


def add_DeepTau_branches(branches_, DeepTauVersion):
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
  for trigger in triggers_dictionary[final_state_mode]:
    branches_.append(trigger)
  return branches_


def Era_F_trigger_study(data_events, final_state_mode):
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
  # template, adjust when used
  # TODO general trigger study algorithm for ORs/ANDs of given trigger set
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


def center(input_string):
  spacer = "-"
  screen = 76
  center = (screen - (len(input_string)))//2
  return spacer*center + input_string + spacer*center
  

def attention(input_string):
  print(center("THE FINAL STATE MODE IS"))
  if getlogin() == "ballmond":
    center_val = (76 - 3*len(input_string))//2
    s_1 = text_options["bold_italic_blink"] + text_options["green"] + input_string + text_options["reset"]
    s_2 = text_options["bold_italic_blink"] + text_options["yellow"] + input_string + text_options["reset"]
    s_3 = text_options["bold_italic_blink"] + text_options["purple"] + input_string + text_options["reset"]
    s_full = center_val*" " + s_1+s_2+s_3 + center_val*" " 
    print(s_full)
  print(center(input_string))


def set_good_events(final_state_mode, trigger_study=False):
  good_events = ""
  if trigger_study: print("*"*20 + " removed trigger cut for yield study " + "*"*20)
  if final_state_mode == "ditau":
    good_events = "(HTT_SRevent) & (abs(HTT_pdgId)==15*15) & (METfilters) & (LeptonVeto==0)"
    if not trigger_study: good_events += " & (Trigger_ditau)"
  elif final_state_mode == "mutau":
    good_events = "(HTT_SRevent) & (abs(HTT_pdgId)==15*13) & (METfilters) & (LeptonVeto==0)"
    if not trigger_study: good_events += " & (Trigger_mutau) & (Trigger_ditau==0)"
  elif final_state_mode == "dimuon":
    good_events = "(HTT_pdgId==-13*13) & (METfilters) & (HLT_IsoMu24)" # lepton veto must be applied manually for this final state

  print(f"good events pass : {good_events}")
  return good_events


def add_final_state_branches(branches_, final_state_mode):
  branch_to_add = []
  if final_state_mode == "ditau":
    pass

  elif final_state_mode == "mutau" or final_state_mode == "etau":
    branch_to_add = ["Tau_pt", "Tau_eta", "MET_pt", "MET_phi", "PuppiMET_pt", "PuppiMET_phi"]
    if final_state_mode == "mutau":
      branch_to_add += ["Muon_pt", "Muon_eta", "Muon_phi"]
    elif final_state_mode == "etau":
      branch_to_add += ["Electron_pt", "Electron_eta", "Electron_phi"]

  elif final_state_mode == "dimuon":
      branch_to_add += ["Muon_pt", "Lepton_pdgId", "Lepton_iso", "HTT_m_vis", "HTT_dR"]

  else:
    print("Hey, that's not a valid final state mode. Try ditau, mutau, etau, or mumu.")

  for new_branch in branch_to_add:
    branches_.append(new_branch)
  
  return branches_


def make_final_state_cut(event_dictionary, useDeepTauVersion, final_state_mode):
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


def fill_process_list(process_list, file_directory, branches, good_events, final_state_mode, testing=False):
  file_map = testing_file_map if testing else full_file_map
  if final_state_mode == "ditau": del(file_map["DataMuon"])
  if final_state_mode == "mutau": del(file_map["DataTau"])
  for process in file_map:
    time_print(f"Loading {file_map[process]}")
    file_string = file_directory + "/" + file_map[process] + ".root:Events"
    try:
      processed_events = uproot.concatenate([file_string], branches, cut=good_events, library="np")
    except FileNotFoundError:
      print(text_options["yellow"] + "FILE NOT FOUND! " + text_options["reset"], end="")
      print(f"continuing without loading {file_map[process]}...")
      continue
    process_list[process] = {}
    process_list[process]["info"] = processed_events
    if "Generator_weight" not in branches: branches.append("Generator_weight") # only works if Data is first
  return process_list


def get_parent_process(MC_process):
  parent_process = ""
  if   "DY"    in MC_process:
    parent_process = "DY"
  elif "WJets" in MC_process:
    parent_process = "WJ"
  elif "TT"    in MC_process:
    parent_process = "TT"
  elif "ST"    in MC_process:
    parent_process = "ST"
  elif ("WW"   in MC_process or 
        "WZ"   in MC_process or 
        "ZZ"   in MC_process):
    parent_process = "VV"
  return parent_process


def match_objects_to_trigger_bit():
  '''Using the final state object kinematics, check if the filter bit of a used trigger is matched'''
  #FS ditau - two taus, match to ditau
  #FS mutau - one tau, one muon
  # - if not cross-trig, match muon to filter
  # - if cross-trig, use cross-trig filters to match both
  match = False
  # step 1 check fired triggers
  # step 2 ensure correct trigger bit is fired
  # step 3 calculate dR and compare with 0.5
  dR_trig_offline = calculate_dR(trig_eta, trig_phi, off_eta, off_phi)


def calculate_dR(eta1, phi1, eta2, phi2): 
  delta_eta = eta1-eta2
  delta_phi = phi_mpi_pi(lep_phi - MET_phi)
  return np.sqrt(delta_eta*delta_eta + delta_phi*delta_phi)

def phi_mpi_pi(delta_phi):
  '''return phi between a range of negative pi and pi'''
  return 2 * np.pi - delta_phi if delta_phi > np.pi else 2 * np.pi + delta_phi if delta_phi < -1*np.pi else delta_phi


def time_print(*args, **kwargs):
  emoji = np.random.choice([";) ", ":^)", "<.<", ">.>", ":O ", "^.^", "UwU", "owO"])
  time  = datetime.now(timezone.utc).strftime('%H:%M:%S')
  if getlogin() == "ballmond":
    time = emoji + "  " + time
  print(f"{time}", *args, **kwargs)


if __name__ == "__main__":

  lxplus_redirector = "root://cms-xrd-global.cern.ch//"
  eos_user_dir      = "eos/user/b/ballmond/NanoTauAnalysis/analysis/HTauTau_2022_fromstep1_skimmed/"
  lxplus_directory  = lxplus_redirector + eos_user_dir
  # there's no place like home :)
  home_directory    = "/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3SkimmedSamples"
  using_directory   = home_directory
  print(f"CURRENT FILE DIRECTORY : {using_directory}")
  

  # final_state_mode sets dataset, "good_events" filter, and cuts
  #final_state_mode = "ditau"  # 5/1 min full/skim dataset
  final_state_mode = "mutau"  # 10/5 min full/skim dataset
  #final_state_mode = "etau"   # not working yet
  #final_state_mode = "dimuon" # 30 min full dataset ~ 3min testing

  attention(final_state_mode)
  # TODO implement arg parse to handle testing, dataset, interactively showing plots, what to plot
  testing    = True
  show_plots = True
  lumi       = luminosities["2022 G"] if testing else luminosities["2022 F&G"]
  print(f"Testing: {testing}")
  useDeepTauVersion = "2p5"
  print(f"USING DEEP TAU VERSION {useDeepTauVersion}")

  good_events = set_good_events(final_state_mode)

  native_variables = ["MET_pt", "PuppiMET_pt", "nCleanJet", "HTT_dR", "HTT_m_vis",
                      #"HTT_DiJet_MassInv_fromHighestMjj", "HTT_DiJet_dEta_fromHighestMjj",
                      "HTT_H_pt_using_PUPPI_MET"]
  added_mutau_variables  = ["FS_mu_pt", "FS_mu_eta", "FS_tau_pt", "FS_tau_eta", "HTT_mt", "CleanJet_btagWP"]
  added_ditau_variables  = ["FS_t1_pt", "FS_t1_eta", "FS_t2_pt", "FS_t2_eta"]
  added_variables = added_ditau_variables if final_state_mode=="ditau" else added_mutau_variables
  full_variables   = native_variables + added_variables
  vars_to_plot     = ["FS_mu_pt", "FS_mu_eta"] if testing else full_variables

  # TODO: make and store jet branches correctly
  #  i.e. branches above ending in "fromHighestMjj" should only be plotted for events with nJet>2
  #  "nCleanJetGT30" : (8, 0, 8), # GT = Greater Than 
  #
  
  print(f"Plotting {vars_to_plot}!")

  added_by_processing = ["FS_t1_pt", "FS_t2_pt", "FS_t1_eta", "FS_t2_eta",
                         "FS_mu_pt", "FS_mu_eta", "FS_tau_pt", "FS_tau_eta",
                         "HTT_mt"]
  branches = [
              "run", "luminosityBlock", "event",
              "FSLeptons", "Lepton_pt", "Lepton_eta",
              "Lepton_tauIdx", "Lepton_muIdx", 
              "HTT_Lep_pt", "HTT_Tau_pt",
              "TrigObj_filterBits", # for matching...
             ]

  branches += [var for var in native_variables if var not in branches and var not in added_by_processing]
  branches = add_DeepTau_branches(branches, useDeepTauVersion)
  branches = add_trigger_branches(branches, final_state_mode)
  branches = add_final_state_branches(branches, final_state_mode)
  # errors like " 'NoneType' object is not iterable " usually mean you forgot
  # to add some branch relevant to the final state

  # if testing, only use Data from Era G (small, 3 fb^-1), DY, WJ, and VBF
  # else use all Era FG 2022 Data for Tau or Muon dataset, and all available MC samples listed in file_map.py
  process_list = {}
  process_list = fill_process_list(process_list, using_directory, branches, good_events, final_state_mode,
                                   testing=testing)

  # TODO : consider combining the cutting step with processing step.
  # naively i think it would reduce the overall memory consumption,
  # because once a big set of files is loaded, it is immediately reduced,
  # then another big set would be loaded and reduced, as opposed to loading several
  # large data sets and then reducing them all in turn.
  # simply bundle below code into "fill_process_list", with an initialization sequence for containers

  # TODO: move to lxplus, make the jump!
  # UPDATE: 40+ minutes on lxplus, and no output.. difficult situation to debug
  # for now will keep things local and work to scale up and over by the end of October!
  # UPDATE Oct 12: ~25 minutes on lxplus, with working output!
  # idea to immediately cut any loaded data to reduce memory use, allowing quicker plotting
  # of the 25 minutes, 1 was spent plotting, 3 processing, and the rest loading the files (network time)
  # the User Time was 4.5 minutes

  # make and apply cuts to any loaded events, store in new dictionaries for plotting
  protected_processes = ["DataTau", "DataMuon", "ggH", "VBF"]
  data_dictionary, stack_dictionary, signal_dictionary = {}, {}, {}
  for process in process_list: 
    time_print(f"Processing {process}")
    process_events = process_list[process]["info"]
    process_events = append_lepton_indices(process_events)
    process_events = make_final_state_cut(process_events, useDeepTauVersion, final_state_mode)
    process_events = apply_cut(process_events, "pass_cuts")

    if process not in protected_processes:
      stack_dictionary[process] = {
        "PlotEvents": {}, 
        "Generator_weight": process_events["Generator_weight"],
      }
      for var in vars_to_plot:
        stack_dictionary[process]["PlotEvents"][var] = process_events[var]

    elif process == "ggH" or process == "VBF":
      signal_dictionary[process] = {
        "PlotEvents": {},
        "Generator_weight": process_events["Generator_weight"],
        "plot_scaling" : 100 if process == "ggH" else 500 if process == "VBF" else 50,
      }
      for var in vars_to_plot:
        signal_dictionary[process]["PlotEvents"][var] = process_events[var]
  
    elif "Data" in process:
      # overwrites process name of data for homogeneity in later plotting
      data_dictionary["Data"] = {
        "PlotEvents":{}
      }
      for var in vars_to_plot:
        data_dictionary["Data"]["PlotEvents"][var] = process_events[var]
    else:
      print("Process not recognized: {process}")

  time_print("Processing finished!")
  ## end processing loop, begin plotting
  
  # TODO: commit in a sensible place at a sensible configuration
  
  for var in vars_to_plot:
    time_print(f"Plotting {var}")
    # plotting setup
    xbins = make_bins(var)
    hist_ax, hist_ratio = setup_ratio_plot()
    h_data, h_backgrounds, h_signals = [], [], []
  
    # bin Data, MC, signal
    data_variable = data_dictionary["Data"]["PlotEvents"][var]
    data_weights  = np.ones(np.shape(data_variable)) # weights of one for data
    h_data = get_binned_info("Data", data_variable, xbins, data_weights, lumi)

    h_MC_by_process = {}
    for process in stack_dictionary:
      process_variable = stack_dictionary[process]["PlotEvents"][var]
      process_weights  = stack_dictionary[process]["Generator_weight"]
      h_MC_by_process[process] = {}
      h_MC_by_process[process]["BinnedEvents"] = get_binned_info(process, process_variable, xbins, process_weights, lumi)

    # add together subprocesses of each MC family
    h_MC_by_family = {}
    MC_families = ["DY", "WJ", "VV"] if testing else ["DY", "TT", "ST", "WJ", "VV"]
    for family in MC_families:
      h_MC_by_family[family] = {}
      h_MC_by_family[family]["BinnedEvents"] = accumulate_MC_subprocesses(family, h_MC_by_process)
    h_backgrounds = h_MC_by_family

    h_signals = {}
    for signal in signal_dictionary:
      signal_variable = signal_dictionary[signal]["PlotEvents"][var]
      signal_weights  = signal_dictionary[signal]["Generator_weight"]
      h_signals[signal] = {}
      h_signals[signal]["BinnedEvents"] = get_binned_info(signal, signal_variable, xbins, signal_weights, lumi)

    # derive additional quanities (for labels, ratio plot, etc)
    yields = calculate_yields(h_data, h_backgrounds, h_signals)
    
    h_summed_backgrounds = 0
    for background in h_backgrounds:
      h_summed_backgrounds += h_backgrounds[background]["BinnedEvents"]
    
    # plot everything :)
    plot_data(hist_ax, xbins, h_data, lumi)
    plot_MC(hist_ax, xbins, h_MC_by_family, lumi)
    plot_signal(hist_ax, xbins, h_signals, lumi)

    make_ratio_plot(hist_ratio, xbins, h_data, h_summed_backgrounds)
  
    spruce_up_plot(hist_ax, hist_ratio, var, lumi)
    hist_ax.legend()
  
    # TODO: more organized plotting output
    # save and show results
    #if not exist make dir for plots
    #else make dir + "alt" in name
    #  print("WARNING: directory already exists, putting images in alternate: {alternate_dir}")
    
    plt.savefig(str(var) + "_my_happy_plot.png")
  if show_plots: plt.show()

