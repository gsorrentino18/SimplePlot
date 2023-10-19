# Authored by Braden Allmond, Sep 11, 2023

# libraries
import uproot # only used by fill_process_list
import numpy as np
import sys
import matplotlib.pyplot as plt

# explicitly import used functions from user files, grouped roughly by call order and relatedness
from plotting_functions    import setup_ratio_plot, make_ratio_plot, spruce_up_plot
from plotting_functions    import plot_data, plot_MC, plot_signal

from get_and_set_functions import set_good_events, make_bins, get_binned_info
from get_and_set_functions import add_final_state_branches, add_DeepTau_branches, add_trigger_branches
from get_and_set_functions import accumulate_MC_subprocesses

from calculate_functions   import calculate_signal_background_ratio
from utility_functions     import time_print, attention, text_options, make_directory

from cut_and_study_functions import make_final_state_cut, apply_cut, append_lepton_indices

# fill process list is the only function which uses this line
from file_map      import testing_file_map, full_file_map, luminosities

# TODO, change from receiving list to receiving one file
# actually, this receives an empty list, so that's an erroneous argument
#
# would be best to separate out the file_map so you have more control over 
# what is loaded when
#
# would like to load one process at a time, immediately apply cuts, and store that info
# (destroying/freeing resources as early as possible), and move on until finished
# come back to this after TT/binning fixes...
# https://www.codingdeeply.com/delete-object-from-memory-python/
def fill_process_list(process_list, file_directory, branches, good_events, final_state_mode, testing=False):
  '''
  Most important function! Contains the only call to uproot in this library! 
  Loads into memory files relevant to the given 'final_state_mode' by reading
  'file_map' which is a python dictionary maintained in a separate file. 
  uproot.concatenate grabs all files matching the wildcard in 'file_map[process]'
  and loads ONLY the data specified by 'branches' which pass the cut 'good_events'.
  Both 'branches' and 'good_events' are specified in other places and depend on the
  final state mode.
  library='np' loads the data in a numpy array that looks like this
   {{"branch_1" : [event1, event2, event3, ..., eventN]},
    {"branch_2" : [event1, event2, event3, ..., eventN]},
    {"branch_3" : [event1, event2, event3, ..., eventN]}, etc.
    {"branch_N" : [event1, event2, event3, ..., eventN]}}
  This coding library is built using numpy arrays as the default and will not work
  with other types of arrays (although the methods could be copied and rewritten). 
  Note: that a numpy array is generated for each loaded process, which corresponds
  to a set of files. 
  '''
  file_map = testing_file_map if testing else full_file_map
  if final_state_mode == "ditau": del(file_map["DataMuon"])
  elif final_state_mode == "mutau": del(file_map["DataTau"])
  else: 
    print(f"{final_state_mode} isn't a final state that I have a mapping for! Quitting...")
    sys.exit()
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


def match_objects_to_trigger_bit():
  '''
  Current work in progress
  Using the final state object kinematics, check if the filter bit of a used trigger is matched
  '''
  #FS ditau - two taus, match to ditau
  #FS mutau - one tau, one muon
  # - if not cross-trig, match muon to filter
  # - if cross-trig, use cross-trig filters to match both
  match = False
  # step 1 check fired triggers
  # step 2 ensure correct trigger bit is fired
  # step 3 calculate dR and compare with 0.5
  dR_trig_offline = calculate_dR(trig_eta, trig_phi, off_eta, off_phi)


if __name__ == "__main__":
  '''
  Just read the code, it speaks for itself.
  Kidding.
  This is the main block, which calls a bunch of other functions from other files
  and uses local variables and short algorithms to, by final state
  1) load data from files
  2) apply bespoke cuts to reject events
  3) create a lovely plot

  Ideally, if one wants to use this library to make another type of plot, they
  would look at this script and use its format as a template.
  '''


  lxplus_redirector = "root://cms-xrd-global.cern.ch//"
  eos_user_dir      = "eos/user/b/ballmond/NanoTauAnalysis/analysis/HTauTau_2022_fromstep1_skimmed/"
  lxplus_directory  = lxplus_redirector + eos_user_dir
  # there's no place like home :)
  home_directory    = "/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3SkimmedSamples"
  using_directory   = home_directory
  print(f"CURRENT FILE DIRECTORY : {using_directory}")
  
  
  import argparse 
  parser = argparse.ArgumentParser(description='Make a standard Data-MC agreement plot.')
  # store_true : when the argument is supplied, store it's value as true
  # for 'testing' below, the default value is false if the argument is not specified
  parser.add_argument('--testing',     dest='testing',     default=False,   action='store_true')
  parser.add_argument('--hide_plots',  dest='hide_plots',  default=False,   action='store_true')
  parser.add_argument('--hide_yields', dest='hide_yields', default=False,   action='store_true')
  parser.add_argument('--final_state', dest='final_state', default="mutau", action='store')
  parser.add_argument('--plot_dir',    dest='plot_dir',    default="plots", action='store')

  args = parser.parse_args() 
  testing     = args.testing     # False by default, do full dataset unless otherwise specified
  hide_plots  = args.hide_plots  # False by default, show plots unless otherwise specified
  hide_yields = args.hide_yields # False by default, show yields unless otherwise specified
  lumi = luminosities["2022 G"] if testing else luminosities["2022 F&G"]
  useDeepTauVersion = "2p5"

  # final_state_mode affects dataset, 'good_events' filter, and cuts
  final_state_mode = args.final_state # default mutau [possible values ditau, mutau, etau, dimuon]
  plot_dir = make_directory(args.plot_dir, args.final_state, testing=testing) # for output files of plots

  # show info to user
  attention(final_state_mode)
  print(f"Testing: {testing}")
  print(f"USING DEEP TAU VERSION {useDeepTauVersion}")

  good_events = set_good_events(final_state_mode)

  native_variables = ["MET_pt", "PuppiMET_pt", "nCleanJet", "HTT_dR", "HTT_m_vis",
                      #"HTT_DiJet_MassInv_fromHighestMjj", "HTT_DiJet_dEta_fromHighestMjj",
                      "HTT_H_pt_using_PUPPI_MET"]
  #added_mutau_variables  = ["FS_mu_pt", "FS_mu_eta", "FS_tau_pt", "FS_tau_eta", "HTT_mt", "CleanJet_btagWP"]
  added_mutau_variables  = ["FS_mu_pt", "FS_mu_eta", "FS_tau_pt", "FS_tau_eta", "HTT_mt"]
  added_ditau_variables  = ["FS_t1_pt", "FS_t1_eta", "FS_t2_pt", "FS_t2_eta"]
  added_variables  = added_ditau_variables if final_state_mode=="ditau" else added_mutau_variables
  full_variables   = native_variables + added_variables
  vars_to_plot     = ["FS_mu_pt", "FS_mu_eta"] if testing else full_variables

  # TODO: make and store jet branches correctly
  #  i.e. branches above ending in "fromHighestMjj" should only be plotted for events with nJet>2
  #  "nCleanJetGT30" : (8, 0, 8), # GT = Greater Than 
  
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
  # UPDATE Oct 18th: did some memory consumption testing and found that 4GB of memory are used when running 
  # this script in testing mode. Also, files not automatically closed. Would managing this better help?

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
  
    plt.savefig(plot_dir + "/" + str(var) + ".png")

    # calculate and print these quanities only once
    if (var=="FS_mu_pt" or var=="FS_t1_pt"): calculate_signal_background_ratio(h_data, h_backgrounds, h_signals)

  if hide_plots: pass
  else: plt.show()


