# Authored by Braden Allmond, Sep 11, 2023

# libraries
import numpy as np
import sys
import matplotlib.pyplot as plt
import gc

# explicitly import used functions from user files, grouped roughly by call order and relatedness
from file_functions        import testing_file_map, full_file_map, testing_dimuon_file_map, dimuon_file_map, luminosities
from file_functions        import compare_eras_file_map
from file_functions        import load_process_from_file, append_to_combined_processes, sort_combined_processes

from cut_and_study_functions import set_branches, set_vars_to_plot # TODO set good events should be here
from cut_and_study_functions import apply_cuts_to_process

from plotting_functions    import get_binned_data, get_binned_backgrounds, get_binned_signals
from plotting_functions    import setup_ratio_plot, make_ratio_plot, spruce_up_plot, spruce_up_legend
from plotting_functions    import plot_data, plot_MC, plot_signal

from get_and_set_functions import set_good_events, make_bins

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

def plot_QCD_preview(xbins, h_data, h_summed_backgrounds, h_QCD, h_MC_frac, h_QCD_FF):
  FF_before_after_ax, FF_info_ax = setup_ratio_plot()

  FF_before_after_ax.set_title("QCD Preview")
  FF_before_after_ax.set_ylabel("Events / Bin")
  FF_before_after_ax.minorticks_on()

  FF_before_after_ax.plot(xbins[0:-1], h_data, label="Data",
                          color="black", marker="o", linestyle='none', markersize=3)
  FF_before_after_ax.plot(xbins[0:-1], h_summed_backgrounds, label="MC",
                          color="blue", marker="^", linestyle='none', markersize=3)
  FF_before_after_ax.plot(xbins[0:-1], h_QCD, label="QCD", 
                          color="orange", marker="v", linestyle='none', markersize=4)

  FF_info_ax.plot(xbins[0:-1], h_MC_frac, label="1-MC/Data",
                  color="red", marker="*", linestyle='none', markersize=3)
  FF_info_ax.plot(xbins[0:-1], h_QCD_FF, label="FF from fit",
                  color="green", marker="s", linestyle='none', markersize=3)
  FF_info_ax.axhline(y=1, color='grey', linestyle='--')

  FF_before_after_ax.legend()
  FF_info_ax.legend()


def set_FF_values(final_state_mode, jet_mode):
  '''
  '''
  # should have aiso/iso as well
  FF_values = {
    "ditau" : {
      "0j"  : [7.25594, -0.0289055],
      "1j"  : [6.57284, -0.0214702],
      "2j"  : [5.74548, -0.011847], # Not 2j, ≥2j. need to decide on the convention and propagate. GTE? ≥? two_or_more?
    },
    "mutau" : {
      "0j"  : [1.93213, 0.0160684],
      "1j"  : [1.87358, 0.0166869],
      "2j"  : [2.06906, 0.00789786],
    },
    "etau"  : {#Dummy values
      "0j"  : [1, 1],
      "1j"  : [1, 1],
      "2j"  : [1, 1],
    },
  }
  intercept = FF_values[final_state_mode][jet_mode][0]
  slope     = FF_values[final_state_mode][jet_mode][1]

  return intercept, slope
 

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
  and del(large_object) in related functions lets python know we no longer need an object, and its resources can be
  reacquired at the next gc.collect() call
  '''

  import argparse 
  parser = argparse.ArgumentParser(description='Make a standard Data-MC agreement plot.')
  # store_true : when the argument is supplied, store it's value as true
  # for 'testing' below, the default value is false if the argument is not specified
  parser.add_argument('--testing',     dest='testing',     default=False,       action='store_true')
  parser.add_argument('--hide_plots',  dest='hide_plots',  default=False,       action='store_true')
  parser.add_argument('--hide_yields', dest='hide_yields', default=False,       action='store_true')
  parser.add_argument('--final_state', dest='final_state', default="mutau",     action='store')
  parser.add_argument('--plot_dir',    dest='plot_dir',    default="plots",     action='store')
  parser.add_argument('--lumi',        dest='lumi',        default="2022 F&G",  action='store')
  parser.add_argument('--jet_mode',    dest='jet_mode',    default="Inclusive", action='store')

  args = parser.parse_args() 
  testing     = args.testing     # False by default, do full dataset unless otherwise specified
  hide_plots  = args.hide_plots  # False by default, show plots unless otherwise specified
  hide_yields = args.hide_yields # False by default, show yields unless otherwise specified
  lumi = luminosities["2022 G"] if testing else luminosities[args.lumi]
  useDeepTauVersion = "2p5"

  # final_state_mode affects many things automatically, including good_events, datasets, plotting vars, etc.
  final_state_mode = args.final_state # default mutau [possible values ditau, mutau, etau, dimuon]
  jet_mode         = args.jet_mode # default Inclusive [possible values 0j, 1j, 2j, GTE2j]

  #lxplus_redirector = "root://cms-xrd-global.cern.ch//"
  eos_user_dir    = "eos/user/b/ballmond/NanoTauAnalysis/analysis/HTauTau_2022_fromstep1_skimmed/" + final_state_mode
  # there's no place like home :)
  home_dir        = "/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3PreEEFSSplitSamples/" + final_state_mode
  home_dir        = "/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3FSSplitSamples/" + final_state_mode
  using_directory = home_dir # TODO: add some auto handling for 'if not <<some env var specific to lxplus>>'
 
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
  if final_state_mode == "dimuon": 
    file_map = testing_dimuon_file_map if testing else dimuon_file_map

  file_map = compare_eras_file_map
  print(file_map)

  # make and apply cuts to any loaded events, store in new dictionaries for plotting
  combined_process_dictionary = {}
  for process in file_map: 

    gc.collect()
    if "Data" not in process: continue

    if process=="DataMuonEraF": continue

    if   final_state_mode == "ditau"  and ("Muon" in process or "Electron" in process): continue
    elif final_state_mode == "mutau"  and ("Tau" in process  or "Electron" in process): continue
    elif final_state_mode == "etau"   and ("Tau" in process  or "Muon" in process):     continue
    elif final_state_mode == "dimuon" and not ("Muon" in process or "DY" in process):   continue

    new_process_dictionary = load_process_from_file(process, using_directory, file_map,
                                              branches, good_events, final_state_mode,
                                              data=("Data" in process), testing=testing)
    if new_process_dictionary == None: continue # skip process if empty

    cut_events = apply_cuts_to_process(process, new_process_dictionary, final_state_mode, jet_mode,
                                       DeepTau_version=useDeepTauVersion)
    if cut_events == None: continue

    combined_process_dictionary = append_to_combined_processes(process, cut_events, vars_to_plot, 
                                                               combined_process_dictionary)

  # after loop, sort big dictionary into three smaller ones
  data_dictionary = combined_process_dictionary

  print(data_dictionary)

  #data_dict_eraF = {"DataMuonEraF" : data_dictionary["DataMuonEraF"]}
  data_dict_eraG = {"DataMuonEraG" : data_dictionary["DataMuonEraG"]}

  time_print("Processing finished!")
  ## end processing loop, begin plotting

  for var in vars_to_plot:
    time_print(f"Plotting {var}")

    xbins = make_bins(var)
    hist_ax, hist_ratio = setup_ratio_plot()

    #h_data_eraF = get_binned_data(data_dict_eraF, var, xbins, lumi)
    h_data_eraG = get_binned_data(data_dict_eraG, var, xbins, lumi)

    # plot everything :)
    #plot_data(hist_ax, xbins, h_data_eraF, lumi, label="2022 Era F")
    plot_data(hist_ax, xbins, h_data_eraG, lumi, color="red", label="2022 Era G")

    #make_ratio_plot(hist_ratio, xbins, h_data, h_summed_backgrounds)
  
    spruce_up_plot(hist_ax, hist_ratio, var, lumi)
    spruce_up_legend(hist_ax, final_state_mode="skip_dimuon_handling")

    plt.savefig(plot_dir + "/" + str(var) + ".png")

    # calculate and print these quantities only once
    #if (var == "HTT_m_vis"): 
      #calculate_signal_background_ratio(h_data, h_backgrounds, h_signals)
      #labels, yields = yields_for_CSV(hist_ax, desired_order=["Data", "TT", "WJ", "DY", "VV", "ST", "ggH", "VBF"])
      #print(f"Reordered     Labels: {labels}")
      #print(f"Corresponding Yields: {yields}")

  if hide_plots: pass
  else: plt.show()


