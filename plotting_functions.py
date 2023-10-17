# libraries
import numpy as np
import matplotlib.pyplot as plt

# user files
from MC_dictionary        import MC_dictionary

from get_and_set_functions import get_midpoints, set_MC_process_info

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


