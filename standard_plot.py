# Authored by Braden Allmond, Sep 11, 2023

# libraries
import numpy as np
import sys
import matplotlib.pyplot as plt
import gc

# explicitly import used functions from user files, grouped roughly by call order and relatedness
from file_functions        import testing_file_map, full_file_map, luminosities
from file_functions        import load_process_from_file, append_to_combined_processes, sort_combined_processes

from cut_and_study_functions import set_branches, set_vars_to_plot # set good events should be here
from cut_and_study_functions import apply_final_state_cut, apply_jet_cut, append_lepton_indices

from plotting_functions    import setup_ratio_plot, make_ratio_plot, spruce_up_plot
from plotting_functions    import plot_data, plot_MC, plot_signal

from get_and_set_functions import set_good_events, make_bins, get_binned_info
from get_and_set_functions import accumulate_MC_subprocesses, accumulate_datasets

from calculate_functions   import calculate_signal_background_ratio, yields_for_CSV
from utility_functions     import time_print, make_directory, print_setup_info
 

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

  import argparse 
  parser = argparse.ArgumentParser(description='Make a standard Data-MC agreement plot.')
  # store_true : when the argument is supplied, store it's value as true
  # for 'testing' below, the default value is false if the argument is not specified
  parser.add_argument('--testing',     dest='testing',     default=False,      action='store_true')
  parser.add_argument('--hide_plots',  dest='hide_plots',  default=False,      action='store_true')
  parser.add_argument('--hide_yields', dest='hide_yields', default=False,      action='store_true')
  parser.add_argument('--final_state', dest='final_state', default="mutau",    action='store')
  parser.add_argument('--plot_dir',    dest='plot_dir',    default="plots",    action='store')
  parser.add_argument('--lumi',        dest='lumi',        default="2022 F&G", action='store')

  args = parser.parse_args() 
  testing     = args.testing     # False by default, do full dataset unless otherwise specified
  hide_plots  = args.hide_plots  # False by default, show plots unless otherwise specified
  hide_yields = args.hide_yields # False by default, show yields unless otherwise specified
  lumi = luminosities["2022 G"] if testing else luminosities[args.lumi]
  useDeepTauVersion = "2p5"

  # final_state_mode affects many things automatically, including good_events, datasets, plotting vars, etc.
  final_state_mode = args.final_state # default mutau [possible values ditau, mutau, etau, dimuon]
  jet_mode         = "none"

  #lxplus_redirector = "root://cms-xrd-global.cern.ch//"
  eos_user_dir    = "eos/user/b/ballmond/NanoTauAnalysis/analysis/HTauTau_2022_fromstep1_skimmed/" + final_state_mode
  # there's no place like home :)
  home_dir        = "/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3FSSplitSamples/" + final_state_mode
  using_directory = home_dir # if not <<some env var specific to lxplus>>
 
  good_events  = set_good_events(final_state_mode)
  branches     = set_branches(final_state_mode)
  # jet category define at run time as 0, 1, 2, inclusive (≥0), ≥1, or ≥2
  vars_to_plot = set_vars_to_plot(final_state_mode, jet_mode=jet_mode)
  plot_dir = make_directory(args.plot_dir, args.final_state, testing=testing) # for output files of plots

  # show info to user
  print_setup_info(final_state_mode, lumi, jet_mode, testing, useDeepTauVersion,
                   using_directory, plot_dir,
                   good_events, branches, vars_to_plot)

  file_map = testing_file_map if testing else full_file_map

  # make and apply cuts to any loaded events, store in new dictionaries for plotting
  combined_process_dictionary = {}
  for process in file_map: 

    gc.collect()
    if   final_state_mode == "ditau" and (process=="DataMuon" or process=="DataElectron"): continue
    elif final_state_mode == "mutau" and (process=="DataTau"  or process=="DataElectron"): continue
    elif final_state_mode == "etau"  and (process=="DataTau"  or process=="DataMuon"):     continue

    new_process_list = load_process_from_file(process, using_directory, 
                                              branches, good_events, final_state_mode,
                                              data=("Data" in process), testing=testing)
    if new_process_list == None: continue # skip datasets if nothing is in them

    time_print(f"Processing {process}")
    process_events = new_process_list[process]["info"]
    if len(process_events["run"])==0: continue # skip datasets if nothing is in them
    del(new_process_list)

    process_events = append_lepton_indices(process_events)
    cut_events = apply_final_state_cut(process_events, final_state_mode, useDeepTauVersion)
    if len(cut_events["pass_cuts"])==0: continue # skip datasets if nothing is in them
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

    ### FUNCTION return h_data
    # accumulate datasets into one flat dictionary called h_data  
    h_data_by_dataset = {}
    for dataset in data_dictionary:
      data_variable = data_dictionary[dataset]["PlotEvents"][var]
      data_weights  = np.ones(np.shape(data_variable)) # weights of one for data
      h_data_by_dataset[dataset] = {}
      h_data_by_dataset[dataset]["BinnedEvents"] = get_binned_info(dataset, data_variable, xbins, data_weights, lumi)
    h_data = accumulate_datasets(h_data_by_dataset)
    ###

    # FUNCTION return h_backgrounds and h_summed_backgrounds
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
    MC_families = ["DY", "TT", "ST", "WJ", "VV"]
    for family in MC_families:
      h_MC_by_family[family] = {}
      h_MC_by_family[family]["BinnedEvents"] = accumulate_MC_subprocesses(family, h_MC_by_process)
    h_backgrounds = h_MC_by_family
    # used for ratio plot
    h_summed_backgrounds = 0
    for background in h_backgrounds:
      h_summed_backgrounds += h_backgrounds[background]["BinnedEvents"]
    ###

   
    # FUNCTION return h_signals 
    # signal is put in a separate dictionary from MC, but they are processed very similarly
    h_signals = {}
    for signal in signal_dictionary:
      signal_variable = signal_dictionary[signal]["PlotEvents"][var]
      signal_weights  = signal_dictionary[signal]["Generator_weight"]
      h_signals[signal] = {}
      h_signals[signal]["BinnedEvents"] = get_binned_info(signal, signal_variable, xbins, signal_weights, lumi)
    ###

    # plot everything :)
    plot_data(hist_ax, xbins, h_data, lumi)
    plot_MC(hist_ax, xbins, h_backgrounds, lumi)
    plot_signal(hist_ax, xbins, h_signals, lumi)

    make_ratio_plot(hist_ratio, xbins, h_data, h_summed_backgrounds)
  
    spruce_up_plot(hist_ax, hist_ratio, var, lumi)
    hist_ax.legend()
  
    plt.savefig(plot_dir + "/" + str(var) + ".png")

    # calculate and print these quantities only once
    if (var == "HTT_m_vis"): 
      calculate_signal_background_ratio(h_data, h_backgrounds, h_signals)
      yields_for_CSV(hist_ax, desired_order=["Data", "TT", "WJ", "DY", "VV", "ST", "ggH", "VBF"])

  if hide_plots: pass
  else: plt.show()


