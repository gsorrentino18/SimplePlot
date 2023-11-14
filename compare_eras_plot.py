# Authored by Braden Allmond, Sep 11, 2023

# libraries
import numpy as np
import sys
import matplotlib.pyplot as plt
import gc

# explicitly import used functions from user files, grouped roughly by call order and relatedness
from file_map_dictionary   import testing_file_map, full_file_map, testing_dimuon_file_map, dimuon_file_map
from file_map_dictionary   import compare_eras_file_map
from file_functions        import load_process_from_file, append_to_combined_processes, sort_combined_processes

from luminosity_dictionary import luminosities_with_normtag as luminosities

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
  This script is meant to compare different eras of data in the same plot. 
  It uses the same basic structure and functions as the main plotting script,
  with some additional handling for splitting up dictionaries and passing styles.
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
  using_directory = "multiple_directories"
 
  good_events  = set_good_events(final_state_mode)
  branches     = set_branches(final_state_mode)
  # jet category define at run time as 0, 1, 2, inclusive (≥0), ≥1, or ≥2
  vars_to_plot = set_vars_to_plot(final_state_mode, jet_mode=jet_mode)
  plot_dir = make_directory(args.plot_dir, args.final_state, testing=testing) # for output files of plots

  # show info to user
  print_setup_info(final_state_mode, lumi, jet_mode, testing, useDeepTauVersion,
                   using_directory, plot_dir,
                   good_events, branches, vars_to_plot)

  file_map = compare_eras_file_map
  print(file_map)

  # make and apply cuts to any loaded events, store in new dictionaries for plotting
  combined_process_dictionary = {}
  for process in file_map: 
    print(f"starting process: {process}")
    using_directory = "/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3Unskimmed/"
    gc.collect()
    if "Data" not in process: continue

    # handling for files in different directories
    if (process=="DataMuonEraF") or (process=="DataMuonEraGPrompt"):
      using_directory = "/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3FSSplitSamples/dimuon/DataPromptReco/"

    if (process=="DataMuonEraE") or (process=="DataMuonEraG"):
      using_directory = "/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3FSSplitSamples/dimuon/DataReReco/"

    if (process=="DataMuonEraC") or (process=="DataMuonEraD"):
      using_directory = "/Users/ballmond/LocalDesktop/trigger_gain_plotting/Run3preEEFSSplitSamples/dimuon/"

    print(using_directory)

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

  # only data in this script
  data_dictionary = combined_process_dictionary

  data_dict_eraC = {"DataMuonEraC" : data_dictionary["DataMuonEraC"]}
  data_dict_eraD = {"DataMuonEraD" : data_dictionary["DataMuonEraD"]}
  data_dict_eraE = {"DataMuonEraE" : data_dictionary["DataMuonEraE"]}
  data_dict_eraF = {"DataMuonEraF" : data_dictionary["DataMuonEraF"]}
  data_dict_eraG = {"DataMuonEraG" : data_dictionary["DataMuonEraG"]}
  data_dict_eraGPrompt = {"DataMuonEraGPrompt" : data_dictionary["DataMuonEraGPrompt"]}

  time_print("Processing finished!")
  ## end processing loop, begin plotting

  for var in vars_to_plot:
    time_print(f"Plotting {var}")

    xbins = make_bins(var)
    hist_ax, hist_ratio = setup_ratio_plot()

    h_data_eraC = get_binned_data(data_dict_eraC, var, xbins, lumi)
    h_data_eraD = get_binned_data(data_dict_eraD, var, xbins, lumi)
    h_data_eraE = get_binned_data(data_dict_eraE, var, xbins, lumi)
    h_data_eraF = get_binned_data(data_dict_eraF, var, xbins, lumi)
    h_data_eraG = get_binned_data(data_dict_eraG, var, xbins, lumi)
    h_data_eraGPrompt = get_binned_data(data_dict_eraGPrompt, var, xbins, lumi)

    # plot everything :)
    plot_data(hist_ax, xbins, h_data_eraC, lumi, color="green", label="2022 Era C")
    plot_data(hist_ax, xbins, h_data_eraD, lumi, color="orange", label="2022 Era D")
    plot_data(hist_ax, xbins, h_data_eraE, lumi, color="blue", label="2022 Era E")
    plot_data(hist_ax, xbins, h_data_eraF, lumi, label="2022 Era F (Prompt)")
    plot_data(hist_ax, xbins, h_data_eraG, lumi, color="red", label="2022 Era G", fillstyle="none")
    plot_data(hist_ax, xbins, h_data_eraGPrompt, lumi, color="purple", label="2022 Era G (Prompt)", marker="x")

    #make_ratio_plot(hist_ratio, xbins, h_data_eraC, h_data_eraD) # example
    # reversed dictionary search for era name based on lumi 
    title_era = [key for key in luminosities.items() if key[1] == lumi][0][0]
    title = f"{title_era}, {lumi:.2f}" + r"$fb^{-1}$"
    spruce_up_plot(hist_ax, hist_ratio, var, title)
    spruce_up_legend(hist_ax, final_state_mode="skip_dimuon_handling")

    plt.savefig(plot_dir + "/" + str(var) + ".png")

  if hide_plots: pass
  else: plt.show()


