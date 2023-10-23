# Authored by Braden Allmond, Sep 11, 2023

# libraries
import numpy as np
import sys
import matplotlib.pyplot as plt
import gc

# explicitly import used functions from user files, grouped roughly by call order and relatedness
from file_functions        import testing_file_map, full_file_map, luminosities
from file_functions        import load_process_from_file, append_to_combined_processes, sort_combined_processes

from cut_and_study_functions import make_final_state_cut, apply_cut, append_lepton_indices

from plotting_functions    import setup_ratio_plot, make_ratio_plot, spruce_up_plot
from plotting_functions    import plot_data, plot_MC, plot_signal

from get_and_set_functions import set_good_events, make_bins, get_binned_info
from get_and_set_functions import add_final_state_branches, add_DeepTau_branches, add_trigger_branches
from get_and_set_functions import accumulate_MC_subprocesses, accumulate_datasets

from calculate_functions   import calculate_signal_background_ratio
from utility_functions     import time_print, attention, make_directory


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
  3) explicitly remove large objects after use
  4) create a lovely plot

  Ideally, if one wants to use this library to make another type of plot, they
  would look at this script and use its format as a template.

  This code sometimes loads very large files, and then makes very large arrays from the data.
  Because of this, I do a bit of memory management, which is atypical of python programs.
  This handling reduces the program burden on lxplus nodes, and subsequently leads to faster results.
  Usually, calling the garbage collector manually like this reduces code efficiency, and if the program
  runs very slowly in the future the memory consumption would be the first thing to check.
  In the main loop below, gc.collect() tells python to remove unused objects and free up resources,
  and del(large_object) lets python know we no longer need an object, and its resources can be
  reacquired at the next gc.collect() call
  '''


  lxplus_redirector = "root://cms-xrd-global.cern.ch//"
  eos_user_dir      = "eos/user/b/ballmond/NanoTauAnalysis/analysis/HTauTau_2022_fromstep1_skimmed/"
  lxplus_directory  = lxplus_redirector + eos_user_dir
  # there's no place like home :)
  home_directory    = "/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3SkimmedSamples"
  home_directory    = "/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3FSSplitSamples/mutau"
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
  #added_mutau_variables  = ["FS_mu_pt", "FS_mu_eta", "FS_tau_pt", "FS_tau_eta", "HTT_mt"]
  added_mutau_variables  = ["FS_mu_pt", "FS_mu_eta", "FS_tau_pt", "FS_tau_eta", "HTT_mt", 
                            "dummy_HTT_Lep_pt", "dummy_HTT_Tau_pt"]
  added_ditau_variables  = ["FS_t1_pt", "FS_t1_eta", "FS_t2_pt", "FS_t2_eta"]
  added_variables  = added_ditau_variables if final_state_mode=="ditau" else added_mutau_variables
  full_variables   = native_variables + added_variables
  vars_to_plot     = ["MET_pt", "HTT_m_vis"] if testing else full_variables

  # TODO: make and store jet branches correctly
  #  i.e. branches above ending in "fromHighestMjj" should only be plotted for events with nJet>2
  #  "nCleanJetGT30" : (8, 0, 8), # GT = Greater Than 
  
  print(f"Plotting {vars_to_plot}!")

  added_by_processing = ["FS_t1_pt", "FS_t2_pt", "FS_t1_eta", "FS_t2_eta",
                         "FS_mu_pt", "FS_mu_eta", "FS_tau_pt", "FS_tau_eta",
                         "HTT_mt"]
  branches = [
              "run", "luminosityBlock", "event", "Generator_weight",
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

  # TODO: move to lxplus, make the jump!
  # UPDATE: 40+ minutes on lxplus, and no output.. difficult situation to debug
  # for now will keep things local and work to scale up and over by the end of October!
  # UPDATE Oct 12: ~25 minutes on lxplus, with working output!
  # idea to immediately cut any loaded data to reduce memory use, allowing quicker plotting
  # of the 25 minutes, 1 was spent plotting, 3 processing, and the rest loading the files (network time)
  # the User Time was 4.5 minutes
  # UPDATE Oct 18th: did some memory consumption testing and found that 4GB of memory are used when running 
  # this script in testing mode. Also, files not automatically closed. Would managing this better help?
  # UPDATE Oct 23rd: managing this better helps a lot, reduces hang time by 2x on machines with fewer resources
  # additionally, skimming samples by FS helps tremendously

  file_map = testing_file_map if testing else full_file_map
  # TODO fix double-counting between datasets 
  #if final_state_mode == "ditau": del(file_map["DataMuon"]) # remove this after FS skimming and combine datasets
  #if final_state_mode == "mutau": del(file_map["DataTau"])

  # make and apply cuts to any loaded events, store in new dictionaries for plotting
  combined_process_dictionary = {}
  for process in file_map: 

    gc.collect()
    new_process_list = load_process_from_file(process, using_directory, 
                                              branches, good_events, final_state_mode,
                                              data=("Data" in process), testing=testing)
    if new_process_list == None: continue

    time_print(f"Processing {process}")
    process_events = new_process_list[process]["info"]
    del(new_process_list)

    process_events = append_lepton_indices(process_events)
    process_events = make_final_state_cut(process_events, useDeepTauVersion, final_state_mode)
    cut_events = apply_cut(process_events, "pass_cuts")
    del(process_events)

    combined_process_dictionary = append_to_combined_processes(process, cut_events, vars_to_plot, 
                                                               combined_process_dictionary)
    del(cut_events)

  # after loop, sort big dictionary into three smaller ones
  data_dictionary, background_dictionary, signal_dictionary = sort_combined_processes(combined_process_dictionary)

  time_print("Processing finished!")
  ## end processing loop, begin plotting
  
  for var in vars_to_plot:
    time_print(f"Plotting {var}")

    xbins = make_bins(var)
    hist_ax, hist_ratio = setup_ratio_plot()

    # accumulate datasets into one flat dictionary called h_data  
    h_data_by_dataset = {}
    for dataset in data_dictionary:
      data_variable = data_dictionary[dataset]["PlotEvents"][var]
      data_weights  = np.ones(np.shape(data_variable)) # weights of one for data
      h_data_by_dataset[dataset] = {}
      h_data_by_dataset[dataset]["BinnedEvents"] = get_binned_info(dataset, data_variable, xbins, data_weights, lumi)
    h_data = accumulate_datasets(h_data_by_dataset)

    # treat each MC process, then group the output by family into flat dictionaries
    # also sum all backgrounds into h_summed_backgrounds to use in ratio plot
    h_MC_by_process = {}
    for process in background_dictionary:
      process_variable = background_dictionary[process]["PlotEvents"][var]
      process_weights  = background_dictionary[process]["Generator_weight"]
      h_MC_by_process[process] = {}
      h_MC_by_process[process]["BinnedEvents"] = get_binned_info(process, process_variable, xbins, process_weights, lumi)
    # add together subprocesses of each MC family
    h_MC_by_family = {}
    MC_families = ["DY", "WJ", "VV"] if testing else ["DY", "TT", "ST", "WJ", "VV"]
    for family in MC_families:
      h_MC_by_family[family] = {}
      h_MC_by_family[family]["BinnedEvents"] = accumulate_MC_subprocesses(family, h_MC_by_process)
    h_backgrounds = h_MC_by_family
    # used for ratio plot
    h_summed_backgrounds = 0
    for background in h_backgrounds:
      h_summed_backgrounds += h_backgrounds[background]["BinnedEvents"]

    # signal is put in a separate dictionary from MC, but they are processed very similarly
    h_signals = {}
    for signal in signal_dictionary:
      signal_variable = signal_dictionary[signal]["PlotEvents"][var]
      signal_weights  = signal_dictionary[signal]["Generator_weight"]
      h_signals[signal] = {}
      h_signals[signal]["BinnedEvents"] = get_binned_info(signal, signal_variable, xbins, signal_weights, lumi)

    # plot everything :)
    plot_data(hist_ax, xbins, h_data, lumi)
    plot_MC(hist_ax, xbins, h_backgrounds, lumi)
    plot_signal(hist_ax, xbins, h_signals, lumi)

    make_ratio_plot(hist_ratio, xbins, h_data, h_summed_backgrounds)
  
    spruce_up_plot(hist_ax, hist_ratio, var, lumi)
    hist_ax.legend()
  
    plt.savefig(plot_dir + "/" + str(var) + ".png")

    # calculate and print these quantities only once
    if (var=="FS_mu_pt" or var=="FS_t1_pt"): calculate_signal_background_ratio(h_data, h_backgrounds, h_signals)

  if hide_plots: pass
  else: plt.show()


