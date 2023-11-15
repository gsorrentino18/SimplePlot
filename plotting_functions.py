# libraries
import numpy as np
import matplotlib.pyplot as plt

# for plotting errors in MC
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

### README
# this file contains functions to setup plotting interfaces and draw the plots themselves

# user files
from MC_dictionary        import MC_dictionary

from get_and_set_functions import get_midpoints, set_MC_process_info

from get_and_set_functions import get_binned_info, accumulate_MC_subprocesses, accumulate_datasets

from calculate_functions import yields_for_CSV

def plot_data(histogram_axis, xbins, data_info, luminosity, 
              color="black", label="Data", marker="o", fillstyle="full"):
  '''
  Add the data histogram to the existing histogram axis, computing errors in a simple way.
  For data, since points and error bars are used, they are shifted to the center of the bins.
  TODO: The error calculation should be followed up and separated to another function. 
  '''
  sum_of_data = np.sum(data_info)
  stat_error = np.array([np.sqrt(entry) if entry > 0 else 0 for entry in data_info]) #error = √N
  midpoints   = get_midpoints(xbins)
  label = f"Data [{sum_of_data:>.0f}]" if label == "Data" else label
  histogram_axis.errorbar(midpoints, data_info, yerr=stat_error, 
                          color=color, marker=marker, fillstyle=fillstyle, label=label,
                          linestyle='none', markersize=3)
  #below plots without error bars
  #histogram_axis.plot(midpoints, data_info, color="black", marker="o", linestyle='none', markersize=2, label="Data")


def plot_MC(histogram_axis, xbins, stack_dictionary, luminosity,
            custom=False, color="default", label="MC", fill=True):
  '''
  Add background MC histograms to the existing histogram axis. The input 'stack_dictionary'
  contains a list of backgrounds (which should be pre-grouped, normally), the name of which
  determines colors and labels of the stacked output. 

  Since the 'bar' method of matplotlib doesn't necessarily expect histogram data, 
  the final bin edge is omitted so that the size of the xaxis array and the plotted 
  data array are equal. To stack the plots, the 'bottom'  keyword argument is 
  adjusted each iteration of the loop such that it is the top of the previous histogram. 
  '''
  previous_histogram_tops = 0
  weight_per_bin_squared = 0
  for MC_process in stack_dictionary:
    if custom == True:
      pass
    else:
      color, label, _ = set_MC_process_info(MC_process, luminosity)
    current_hist = stack_dictionary[MC_process]["BinnedEvents"]
    # TODO handle variable binning with list of differences
    #print("weight per bin squared")
    #print(weight_per_bin_squared)
    #weight_per_bin_squared += np.array([entry*entry if entry > 0 else 0 for entry in current_hist])
    label += f" [{np.sum(current_hist):>.0f}]"
    bars = histogram_axis.bar(xbins[0:-1], current_hist, width=xbins[1]-xbins[0],
                            color=color, label=label, edgecolor=color if custom else None,
                            bottom=previous_histogram_tops, fill=fill, align='edge')
    #print("current_hist")
    #print(current_hist)
    #print("bars output")
    #print(bars)
    previous_histogram_tops += current_hist # stack

  summed_squared_weights = weight_per_bin_squared
  #print("summed_squared_weights")
  #print(summed_squared_weights)
  # error = √sum(weights_i**2)
  #stat_error = np.array([np.sqrt(squared_weight) if squared_weight > 0 else 0 for squared_weight in summed_squared_weights])
  #print("stat_error")
  #print(stat_error)
  #histogram_axis.errorbar(xbins[0:-1], previous_histogram_tops, yerr=stat_error, fmt="o")
  # error boxes not quite working yet
  #histogram_axis.errorbar(xbins[0:-1], previous_histogram_tops, yerr=5, fmt="o", markersize=3) # faked
  #_ = make_error_boxes(histogram_axis, xbins[0:-1], previous_histogram_tops, # fake
  #                     abs(xbins[1:]-xbins[0:-1])/2, 10*np.ones(previous_histogram_tops.shape),
  #                     facecolor='grey', edgecolor='none', alpha=0.1)

#from https://matplotlib.org/stable/gallery/statistics/errorbars_and_boxes.html#sphx-glr-gallery-statistics-errorbars-and-boxes-py
def make_error_boxes(ax, xdata, ydata, xerror, yerror, facecolor='r',
                     edgecolor='none', alpha=0.5):
    # Loop over data points; create box from errors at each point
    #errorboxes = [Rectangle((x - xe, y - ye), xe.sum(), ye.sum())
    errorboxes = [Rectangle((x, y), xe, ye)
                  for x, y, xe, ye in zip(xdata, ydata, xerror, yerror)]
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor)
    # Add collection to axes
    ax.add_collection(pc)
    # Plot errorbars
    artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
                          fmt='none', ecolor='k')
    return artists
  

def plot_signal(histogram_axis, xbins, signal_dictionary, luminosity):
  '''
  Similar to plot_MC, except signals are not stacked, and the 'stair' method
  of matplotlib DOES expect histogram data, so no adjustment to xbins is necessary.
  '''
  for signal in signal_dictionary:
    color, label, _ = set_MC_process_info(signal, luminosity, scaling=True, signal=True)
    current_hist = signal_dictionary[signal]["BinnedEvents"]
    label += f" [{np.sum(current_hist):>.0f}]"
    stairs = histogram_axis.stairs(current_hist, xbins, color=color, label=label, fill=False)
  #histogram_axis.errorbar(xbins[0:-1], previous_histogram_tops, yerr=5, fmt="o", markersize=3) # faked
  # actually, check if stairs has an error method
  # TODO : it would be easier to store errors when binning events? probably
  # USE fill_between axis method for shaded errors
  # https://stackoverflow.com/questions/60697625/error-bars-as-a-shaded-area-on-matplotlib-pyplot-step


def setup_ratio_plot():
  '''
  Define a standard plot format with a plotting area on top, and a ratio area below.
  The plots share the x-axis, and other functions should handle cosmetic additions/subtractions.
  '''
  gs = gridspec_kw = {'height_ratios': [4, 1], 'hspace': 0}
  fig, (upper_ax, lower_ax) = plt.subplots(nrows=2, sharex=True, gridspec_kw=gridspec_kw)
  return (upper_ax, lower_ax)


def add_CMS_preliminary(axis):
  '''
  Add text to plot following CMS plotting guidelines
  https://twiki.cern.ch/twiki/bin/viewauth/CMS/Internal/FigGuidelines#Example_ROOT_macro_python
  '''
  CMS_text = "CMS"
  axis.text(0.01, 1.02, CMS_text, transform=axis.transAxes, fontsize=16, weight='bold')
  preliminary_text = "Preliminary"
  axis.text(0.12, 1.02, preliminary_text, transform=axis.transAxes, fontsize=16, style='italic')


def add_final_state_and_jet_mode(axis, final_state_mode, jet_mode):
  axis.text(0.09, 0.94, final_state_mode + "-" + jet_mode, transform=axis.transAxes, fontsize=12)


def spruce_up_plot(histogram_axis, ratio_plot_axis, variable_name, title, final_state_mode, jet_mode):
  '''
  Add title and axes labels
  Additionally:
    - hide a zero that overlaps with the upper plot.
    - add a horizontal line at y=1 to the ratio plot
  '''
  add_CMS_preliminary(histogram_axis)
  add_final_state_and_jet_mode(histogram_axis, final_state_mode, jet_mode)
  histogram_axis.set_title(title, loc='right')
  histogram_axis.set_ylabel("Events / Bin")
  histogram_axis.minorticks_on()
  histogram_axis.tick_params(which="both", top=True, right=True, direction="inout")

  yticks = histogram_axis.yaxis.get_major_ticks()
  yticks[0].label1.set_visible(False) # hides a zero that overlaps with the upper plot
  plt.minorticks_on()

  ratio_plot_axis.set_ylim([0.6, 1.4]) # 0.0, 2.0 also make sense
  ratio_plot_axis.set_xlabel(variable_name) # shared axis label
  ratio_plot_axis.set_ylabel("Obs. / Exp.")
  ratio_plot_axis.axhline(y=1, color='grey', linestyle='--')
  ratio_plot_axis.tick_params(bottom=True, right=True, direction="inout")


def spruce_up_legend(histogram_axis, final_state_mode):
  histogram_axis.legend()
  if final_state_mode == "dimuon":
    handles, original_labels = histogram_axis.get_legend_handles_labels()
    labels, yields = yields_for_CSV(histogram_axis)
    save_entry = [i for i,yield_ in enumerate(yields) if yield_ != 0]
    save_handles = [handles[entry] for entry in save_entry]
    save_labels  = [labels[entry] for entry in save_entry]
    histogram_axis.legend(save_handles, save_labels)
    print(f"Legend lables were previously {original_labels}")
    print("Removed samples with yield=0 from legend!")

 
def make_ratio_plot(ratio_axis, xbins, numerator_data, denominator_data):
  '''
  Uses provided numerator and denominator info to make a ratio to add to given plotting axis.
  Errors are also calculated using the same matplotlib function as used in plot_data.
  '''
  ratio = numerator_data/denominator_data
  ratio[np.isnan(ratio)] = 0 # numpy idiom to set "nan" values to 0
  # TODO : technically errors from stack should be individually calculated, not one stack
  statistical_error = [ ratio[i] * np.sqrt( (np.sqrt(numerator_data[i])   / numerator_data[i]) ** 2 +
                                            (np.sqrt(denominator_data[i]) / denominator_data[i]) ** 2 )
                      for i,_ in enumerate(denominator_data)]

  midpoints = get_midpoints(xbins)
  ratio_axis.errorbar(midpoints, ratio, yerr=statistical_error,
                    color="black", marker="o", linestyle='none', markersize=2)


def get_trimmed_Generator_weight_copy(variable, single_background_dictionary, jet_mode):
  '''
  Only to be used with Inclusive jet mode

  When plotting inclusively, we would like to show j1_pt, j2_pt, j3_pt, etc.
  However, these branches are not available for every event, and to plot them we
  also need the Generator_weight branch. So, what we do is make copies of the Generator_weight
  branch for only the events also found in the variable we want to plot. This is
  possible because the jet_mode cut branch is defined and stored in our dictionary
  '''
  Gen_weight = single_background_dictionary["Generator_weight"]
  Cuts       = single_background_dictionary["Cuts"]
  if   ("_1" in variable) and (jet_mode=="Inclusive"):
    # events with ≥1 j
    pass_jet_cut = sorted(np.concatenate([Cuts["pass_1j_cuts"], Cuts["pass_2j_cuts"], Cuts["pass_3j_cuts"]]))
  elif (("_2" in variable) or ("fromHighestMjj" in variable)) and (jet_mode=="Inclusive"):
    # events with ≥2 j
    pass_jet_cut = sorted(np.concatenate([Cuts["pass_2j_cuts"], Cuts["pass_3j_cuts"]]))
  elif ("_3" in variable) and ((jet_mode=="Inclusive") or (jet_mode=="GTE2j")):
    pass_jet_cut = Cuts["pass_3j_cuts"]

    print("gen weight, pass jet cut, and temp weight shapes") # DEBUG
    print(Gen_weight.shape, pass_jet_cut.shape) # DEBUG
  temp_weight = np.take(Gen_weight, pass_jet_cut)
  print(temp_weight.shape) # DEBUG

  return temp_weight


def get_binned_data(data_dictionary, variable, xbins_, lumi_):
  '''
  Standard loop to get only the plotted variable from a dictionary containing data.

  This is written so that it can be extended to use multiple datasets, but the default
  usage is only one dataset at a time. 
  '''
  h_data_by_dataset = {}
  for dataset in data_dictionary:
    data_variable = data_dictionary[dataset]["PlotEvents"][variable]
    data_weights  = np.ones(np.shape(data_variable)) # weights of one for data
    # normalized to Era D, Era E for some reason is still larger than it should be
    if dataset == "DataMuonEraC":
      data_weights = (2.922/4.953)*data_weights
    if dataset == "DataMuonEraD":
      data_weights = (2.922/2.922)*data_weights
    if dataset == "DataMuonEraE":
      data_weights = (2.922/5.672)*data_weights
    if dataset == "DataMuonEraF":
      data_weights = (2.922/17.61)*data_weights
    if (dataset == "DataMuonEraG") or (dataset=="DataMuonEraGPrompt"):
      data_weights = (2.922/3.06)*data_weights
    #print(f"For {dataset} using weights:")
    #print(data_weights)
    h_data_by_dataset[dataset] = {}
    h_data_by_dataset[dataset]["BinnedEvents"] = get_binned_info(dataset, data_variable, 
                                                                 xbins_, data_weights, lumi_)
  h_data = accumulate_datasets(h_data_by_dataset)
  return h_data


def get_binned_backgrounds(background_dictionary, variable, xbins_, lumi_, jet_mode):
  '''
  Treat each MC process, then group the output by family into flat dictionaries.
  Also, sum all backgrounds into h_summed_backgrounds to use in ratio plot.
  '''
  h_MC_by_process = {}
  for process in background_dictionary:
    process_variable = background_dictionary[process]["PlotEvents"][variable]
    if len(process_variable) == 0: continue
    if process == "QCD":  
      process_weights = background_dictionary[process]["FF_weight"]
    #if ( (("JetGT30_" in variable or "fromHighestMjj" in variable) and (jet_mode=="Inclusive")) or 
    #     ((variable == "CleanJetGT30_pt_3" or variable == "CleanJetGT30_eta3") and (jet_mode=="GTE2j")) ):
    ##if ("JetGT30_" in variable) and (jet_mode=="Inclusive"):
    #  process_weights = get_trimmed_Generator_weight_copy(variable, background_dictionary[process], jet_mode)
    else:
      process_weights_gen = background_dictionary[process]["Generator_weight"]
      process_weights_SF  = background_dictionary[process]["SF_weight"]
      process_weights = process_weights_gen*process_weights_SF
    #print("process, variable, variable and weight shapes") # DEBUG 
    #print(process, variable, process_variable.shape, process_weights.shape) # DEBUG
    h_MC_by_process[process] = {}
    h_MC_by_process[process]["BinnedEvents"] = get_binned_info(process, process_variable, 
                                                               xbins_, process_weights, lumi_)
  # add together subprocesses of each MC family
  h_MC_by_family = {}
  # see what processes exist in the dictionary
  all_MC_families  = ["DY", "TT", "ST", "WJ", "VV"]
  used_MC_families = []
  for process in h_MC_by_process:
    for family in all_MC_families:
      if (("WW" in process) or ("WZ" in process) or ("ZZ" in process)) and ("VV" not in used_MC_families):
        used_MC_families.append("VV")
      elif (family in process) and (family not in used_MC_families):
        used_MC_families.append(family)
  for family in used_MC_families:
    h_MC_by_family[family] = {}
    h_MC_by_family[family]["BinnedEvents"] = accumulate_MC_subprocesses(family, h_MC_by_process)
  if "QCD" in background_dictionary.keys():
    h_MC_by_family["QCD"] = {}
    h_MC_by_family["QCD"]["BinnedEvents"] = h_MC_by_process["QCD"]["BinnedEvents"]
  h_backgrounds = h_MC_by_family
  # used for ratio plot
  h_summed_backgrounds = 0
  for background in h_backgrounds:
    h_summed_backgrounds += h_backgrounds[background]["BinnedEvents"]

  return h_backgrounds, h_summed_backgrounds


def get_binned_signals(signal_dictionary, variable, xbins_, lumi_, jet_mode):
  '''
  Signal is put in a separate dictionary from MC, but they are processed very similarly
  '''
  h_signals = {}
  for signal in signal_dictionary:
    signal_variable = signal_dictionary[signal]["PlotEvents"][variable]
    if len(signal_variable) == 0: continue
    #if ("JetGT30_" in variable or "fromHighestMjj" in variable) and (jet_mode=="Inclusive" or jet_mode=="GTE2j"):
    if ( (("JetGT30_" in variable or "fromHighestMjj" in variable) and (jet_mode=="Inclusive")) or 
         ((variable == "CleanJetGT30_pt_3" or variable == "CleanJetGT30_eta3") and (jet_mode=="GTE2j")) ):
    #if ("JetGT30_" in variable) and (jet_mode=="Inclusive"):
      signal_weights = get_trimmed_Generator_weight_copy(variable, signal_dictionary[signal], jet_mode)
    else:
      signal_weights = signal_dictionary[signal]["Generator_weight"]
    h_signals[signal] = {}
    h_signals[signal]["BinnedEvents"] = get_binned_info(signal, signal_variable,
                                                        xbins_, signal_weights, lumi_)
  return h_signals


