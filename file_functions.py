import uproot
import numpy as np

from utility_functions import time_print, text_options
from get_and_set_functions import set_good_events

### README ###
# This file contains the main method to load data from root files
# The wildcarding works for the 'concatenate' function of uproot, and might not in the future.
# This file also contains methods relevant to sorting samples from files.

def load_process_from_file(process, file_directory, file_map, branches, good_events, final_state_mode, data=False, testing=False):
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
  time_print(f"Loading {file_map[process]}")
  file_string = file_directory + "/" + file_map[process] + ".root:Events"
  if data: 
    branches = [branch for branch in branches if branch != "Generator_weight"]
  try:
    processed_events = uproot.concatenate([file_string], branches, cut=good_events, library="np")
  except FileNotFoundError:
    print(text_options["yellow"] + "FILE NOT FOUND! " + text_options["reset"], end="")
    print(f"continuing without loading {file_map[process]}...")
    return None
  process_list = {}
  process_list[process] = {}
  process_list[process]["info"] = processed_events
 
  return process_list


def sort_combined_processes(combined_processes_dictionary):
  data_dictionary, background_dictionary, signal_dictionary = {}, {}, {}
  for process in combined_processes_dictionary:
    if "Data" in process:
      data_dictionary[process]       = combined_processes_dictionary[process]
    elif ("VBF" in process) or ("ggH" in process):
      signal_dictionary[process]     = combined_processes_dictionary[process]
    else:
      background_dictionary[process] = combined_processes_dictionary[process]
  return data_dictionary, background_dictionary, signal_dictionary


def append_to_combined_processes(process, cut_events, vars_to_plot, combined_processes):
  if "Data" not in process:
    combined_processes[process] = {
      "PlotEvents": {}, 
      "Cuts": {},
      "Generator_weight": cut_events["Generator_weight"],
      "SF_weight": np.ones(cut_events["Generator_weight"].shape)
    }
  elif "Data" in process:
    combined_processes[process] = { 
      "PlotEvents": {},
      "Cuts": {},
    }
  for var in vars_to_plot:
    combined_processes[process]["PlotEvents"][var] = cut_events[var]

  for cut in ["pass_cuts", 
              "pass_0j_cuts", "pass_1j_cuts", "pass_2j_cuts", "pass_3j_cuts",
              "pass_GTE2j_cuts"]:
    if cut in cut_events.keys():
      combined_processes[process]["Cuts"][cut] = cut_events[cut]

  return combined_processes

# not used by default
# only to be used if final states are plotted at the same time and datasets are simultaneously loaded
def reject_duplicate_events(final_state_mode, dataset):
  '''
   disable any triggers in set_good_events
   then manually adding them based on dataset here

   below is confusing to read
   BUT, i put the dataset at the top of each that the final_state_mode
   gets the most data from. Then subleading datasets, vetoing each used trigger
   along the way to remove the overlap in datasets.
  '''
  good_events = set_good_events(final_state_mode, disable_triggers=True)
  
  if final_state_mode == "ditau":
    if dataset == "DataTau":
      good_events += " & (Trigger_ditau)" # implied OR with mutau/etau by convention
    elif dataset == "DataMuon":
      good_events += " & (Trigger_mutau) & (Trigger_ditau==0)" # implied OR with etau
    elif dataset == "DataElectron":
      good_events += " & (Trigger_etau)  & (Trigger_ditau==0) & (Trigger_mutau==0)"

  elif final_state_mode == "mutau":
    if dataset == "DataMuon":
      good_events += " & (Trigger_mutau)" # implied OR with ditau/etau by convention
    elif dataset == "DataTau":
      good_events += " & (Trigger_ditau) & (Trigger_mutau==0)" # implied OR with etau
    elif dataset == "DataElectron":
      good_events += " & (Trigger_etau)  & (Trigger_mutau==0) & (Trigger_ditau==0)"

  elif final_state_mode == "etau":
    if dataset == "DataElectron":
      good_events += " & (Trigger_etau)" # implied OR with ditau/mutau by convention
    elif dataset == "DataTau":
      good_events += " & (Trigger_ditau) & (Trigger_etau==0)" # implied OR with mutau
    elif dataset == "DataMuon":
      good_events += " & (Trigger_mutau) & (Trigger_etau==0) & (Trigger_ditau==0)"

  print(f"Rejecting duplicate events in Data sample with temporary selection \n {good_events}")

  return good_events


