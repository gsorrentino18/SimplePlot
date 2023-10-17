# libraries
import numpy as np
import matplotlib.pyplot as plt

# user files
from MC_dictionary        import MC_dictionary

from get_and_set_functions import get_midpoints, set_MC_process_info
from file_map              import luminosities

def plot_data(histogram_axis, xbins, data_info, luminosity):
  '''
  Add the data histogram to the existing histogram axis, computing errors in a simple way.
  For data, since points and error bars are used, they are shifted to the center of the bins.
  TODO: The error calculation should be followed up and separated to another function. 
  '''
  # TODO: understand errors for poisson statistics
  # TODO: understand all of statistics
  #stat_error = np.array([1000/np.sqrt(entry) if entry > 0 else 0 for entry in data_info]) #wrong but scaled up...
  #stat_error  = sum_of_data**2 * np.ones(np.shape(data_info)) #wrong
  sum_of_data = np.sum(data_info)
  stat_error  = np.array([np.sqrt(entry * (1 - entry/sum_of_data)) if entry > 0 else 0 for entry in data_info])
  #stat_error = np.array([np.sqrt(entry) if entry > 0 else 0 for entry in data_info]) # unsure if correct?
  midpoints   = get_midpoints(xbins)
  label = f"Data [{np.sum(data_info):>.0f}]"
  histogram_axis.errorbar(midpoints, data_info, yerr=stat_error, 
                          color="black", marker="o", linestyle='none', markersize=2, label=label)
  #histogram_axis.plot(midpoints, data_info, color="black", marker="o", linestyle='none', markersize=2, label="Data")
  # above plots without error bars


def plot_MC(histogram_axis, xbins, stack_dictionary, luminosity):
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
  for MC_process in stack_dictionary:
    color, label, _ = set_MC_process_info(MC_process, luminosity)
    current_hist = stack_dictionary[MC_process]["BinnedEvents"]
    # TODO handle variable binning with list of differences
    label += f" [{np.sum(current_hist):>.0f}]"
    bars = histogram_axis.bar(xbins[0:-1], current_hist, width=xbins[1]-xbins[0],
                            color=color, label=label, 
                            bottom=previous_histogram_tops, fill=True, align='edge')
    previous_histogram_tops += current_hist # stack
  

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
  #axis.text(0.01, 1.02, CMS_text, transform=hist_ax.transAxes, fontsize=16, weight='bold')
  axis.text(0.01, 1.02, CMS_text, transform=axis.transAxes, fontsize=16, weight='bold')
  preliminary_text = "Preliminary"
  #axis.text(0.12, 1.02, preliminary_text, transform=hist_ax.transAxes, fontsize=16, style='italic')
  axis.text(0.12, 1.02, preliminary_text, transform=axis.transAxes, fontsize=16, style='italic')


def spruce_up_plot(histogram_axis, ratio_plot_axis, variable_name, luminosity):
  '''
  Add title, axes labels, and data and luminosity info.
  Additionally:
    - hide a zero that overlaps with the upper plot.
    - add a horizontal line at y=1 to the ratio plot
  '''
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

 
def make_ratio_plot(ratio_axis, xbins, numerator_data, denominator_data):
  '''
  Uses provided numerator and denominator info to make a ratio to add to given plotting axis.
  Errors are also calculated using the same matplotlib function as used in plot_data.
  '''
  # TODO fix error calculation, should be np.sqrt(entry / sum? )
  numerator_statistical_error   = [1/np.sqrt(entry) if entry > 0 else 0 for entry in numerator_data]
  denominator_statistical_error = [1/np.sqrt(entry) if entry > 0 else 0 for entry in denominator_data]
  combined_statistical_error    = np.add(numerator_statistical_error, denominator_statistical_error)

  midpoints = get_midpoints(xbins)
  ratio_axis.errorbar(midpoints, numerator_data/denominator_data, yerr=combined_statistical_error,
                    color="black", marker="o", linestyle='none', markersize=2)


