import uproot
import numpy as np

from utility_functions import time_print, text_options
from get_and_set_functions import set_good_events

### README ###
# This file contains mappings of process names (shared with XSec.py) to wildcards for related samples.
# The wildcarding works for the 'concatenate' function of uproot, and might not in the future.
# The testing_file_map is a reduced number of files for faster processing times.
# Additionally, this file contains methods relevant to loading and sorting samples from files.

# fb ^ -1
# useful webpage : https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#Year_2022
# without normtag
luminosities = {
  "2022 C&D" : 7.875,
  "2022 C"   : 4.953, # 2022 pre EE
  "2022 D"   : 2.922, # 2022 pre EE
  "2022 E"   : 5.672, # 2022 EE
  "2022 F"   : 17.61, # 2022 EE
  "2022 G"   :  3.06, # 2022 EE
  "2022 F&G" : 20.67,
  "2022 EFG" : 26.342,
}
# with normtag (default in future)
# TODO: update this and study plots, 2% effect everywhere
luminosities_w_normtag = {
  "2022 C&D" : 8.077,
  "2022 C"   : 5.0707,
  "2022 D"   : 3.0063,
  "2022 E"   : 5.8783,
  "2022 F"   : 18.0070,
  "2022 G"   : 3.1219,
  "2022 F&G" : 21.1289,
  "2022 EFG" : 27.0072,
}

pre2022_file_map = {
  "DataMuon" : "Data/Muon*",

  "DYInc"    : "DY/DYJets*",
}

testing_file_map = {
  "DataTau"  : "Data/Tau*G*",
  "DataMuon" : "Data/Muon*G*",
  "DataElectron" : "Data/EGamma*G*",

  "DYInc" : "DY/DYJets*part1",
  "WJetsInc" : "WJ/WJ*part1",

  "WWTo2L2Nu" : "VV/WWTo2L2Nu*",
  "WWTo4Q"    : "VV/WWTo4Q*",
  "WWToLNu2Q" : "VV/WWToLNu2Q*",
  "WZTo3LNu"  : "VV/WZTo3LNu*",
  "WZTo2L2Q"  : "VV/WZTo2L2Q*",
  "WZToLNu2Q" : "VV/WZToLNu2Q*",
  "ZZTo2L2Nu" : "VV/ZZTo2L2Nu*",
  "ZZTo2L2Q"  : "VV/ZZTo2L2Q*", 
  "ZZTo2Nu2Q" : "VV/ZZTo2Nu2Q*",
  "ZZTo4L"    : "VV/ZZTo4L*", 

  "VBF"   : "Signal/VBF*",
}

compare_eras_file_map = {
  "DataMuonEraC" : "Data*/Muon*C*",
  "DataMuonEraD" : "Data*/Muon*D*",
  "DataMuonEraE" : "Data/Muon*E*",
  "DataMuonEraF" : "Data/Muon*F*",
  "DataMuonEraG" : "Data/Muon*G*",
  "DataMuonEraGPrompt" : "Data/Muon*G*",
  "DataTauEraF" : "Data/Tau*F*",
  "DataTauEraG" : "Data/Tau*G*",
  "DataElectronEraF" : "Data/Electron*F*",
  "DataElectronEraG" : "Data/Electron*G*",
}

dimuon_file_map = {
  "DataMuon" : "Data/Muon*",
  "DYInc"    : "DY/DY*part*",
}

testing_dimuon_file_map = {
  "DataMuon" : "Data/Muon*G*",
  "DYInc"    : "DY/DYJets*part1",
}

full_file_map = {
  "DataTau"  : "Data/Tau*",
  "DataMuon" : "Data/Muon*",
  "DataElectron" : "Data/EGamma*",

  "DYInc" : "DY/DY*part*",
  #"DYJetsToLL_M-50_1J" : "DY/DYJetsToLL_M-50_1J*",
  #"DYJetsToLL_M-50_2J" : "DY/DYJetsToLL_M-50_2J*",
  #"DYJetsToLL_M-50_3J" : "DY/DYJetsToLL_M-50_3J*",
  #"DYJetsToLL_M-50_4J" : "DY/DYJetsToLL_M-50_4J*",

  "TTTo2L2Nu"         : "TT/TTTo2L2Nu*",
  "TTToFullyHadronic" : "TT/TTToFullyHadronic*",
  "TTToSemiLeptonic"  : "TT/TTToSemiLeptonic*",

  "ST_TWminus_2L2Nu"   : "ST/ST_TWminus_2L2Nu*",
  "ST_TbarWplus_2L2Nu" : "ST/ST_TbarWplus_2L2Nu*",
  "ST_TWminus_4Q"      : "ST/ST_TWminus_4Q*",
  "ST_TbarWplus_4Q"    : "ST/ST_TbarWplus_4Q*",
  "ST_TWminus_LNu2Q"   : "ST/ST_TWminus_LNu2Q",
  "ST_TbarWplus_LNu2Q" : "ST/ST_TbarWplus_LNu2Q",

  "WJetsInc" : "WJ/WJ*",
  #"WJetsToLNu_1J" : "WJ/W1JetsToLNu*",
  #"WJetsToLNu_2J" : "WJ/W2JetsToLNu*",
  #"WJetsToLNu_3J" : "WJ/W3JetsToLNu*",
  #"WJetsToLNu_4J" : "WJ/W4JetsToLNu*",

  "WWTo2L2Nu" : "VV/WWTo2L2Nu*",
  "WWTo4Q"    : "VV/WWTo4Q*",
  "WWToLNu2Q" : "VV/WWToLNu2Q*",
  "WZTo3LNu"  : "VV/WZTo3LNu*",
  "WZTo2L2Q"  : "VV/WZTo2L2Q*",
  "WZToLNu2Q" : "VV/WZToLNu2Q*",
  "ZZTo2L2Nu" : "VV/ZZTo2L2Nu*",
  "ZZTo2L2Q"  : "VV/ZZTo2L2Q*", 
  "ZZTo2Nu2Q" : "VV/ZZTo2Nu2Q*",
  "ZZTo4L"    : "VV/ZZTo4L*", 

  "VBF"   : "Signal/VBF*",
  "ggH"   : "Signal/ggH*",
}

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


