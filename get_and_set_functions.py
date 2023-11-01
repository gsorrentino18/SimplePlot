# libraries
import numpy as np

### README
# this file contains simple functions to get or set and return simple objects
# examples include specific logic for reading and updating a dictionary,
# combining dictionary entries of the same type, 
# returning a list based on values defined in a dictionary

# user files
from binning_dictionary   import binning_dictionary
from MC_dictionary        import MC_dictionary
from triggers_dictionary  import triggers_dictionary

from calculate_functions  import calculate_underoverflow


def make_bins(variable_name):
  '''
  Create a linear numpy array to use for histogram binning.
  Information for binning is referenced from a python dictionary in a separate file.
  A check is made on bin edges to see if they end in 1 or 0.1, which generally 
  result in better plots with edges that align with axes tickmarks.
  
  This method returns only linearly spaced bins
  '''
  nbins, xmin, xmax = binning_dictionary[variable_name]
  check_uniformity = (xmax-xmin)/nbins
  if (check_uniformity % 1 != 0 and check_uniformity % 0.1 != 0):
    print(f"nbins, xmin, xmax : {nbins}, {xmin}, {xmax}")
    print(f"(xmax-xmin)/nbins = {check_uniformity}, results in bad bin edges and centers")
  xbins = np.linspace(xmin, xmax, nbins+1)
  return xbins


def get_midpoints(input_bins):
  '''
  From an input array of increasing values, return the values halfway between each value.
  The input array is size N, and the output array is size N-1
  '''
  midpoints = []
  for i, ibin in enumerate(input_bins):
    if (i+1 != len(input_bins)):
      midpoints.append( ibin + (input_bins[i+1] - ibin)/2 )
  midpoints = np.array(midpoints)
  return midpoints


def set_MC_process_info(process, luminosity, scaling=False, signal=False):
  '''
  Obtain process-specific styling and scaling information.
  MC_dictionary is maintained in a separate file.
  '''
  color = MC_dictionary[process]["color"]
  label = MC_dictionary[process]["label"]
  if scaling:
  # factor of 1000 comes from lumi and XSec units of fb^-1 = 10E15 b^-1 and pb = 10E-12 b respectively
    plot_scaling = MC_dictionary[process]["plot_scaling"] # 1 for all non-signal processes by default
    scaling = 1000. * plot_scaling * luminosity * MC_dictionary[process]["XSec"] / MC_dictionary[process]["NWevents"]
    if process=="QCD": scaling = 1
  if signal:
    label += " x" + str(plot_scaling)
  return (color, label, scaling)


def get_binned_info(process_name, process_variable, xbins, process_weights, luminosity):
  '''
  Take in a list of events and produce a histogram (values binned in a numpy array).
  'scaling' is either set to 1 for data (no scaling) or retrieved from the MC_dictionary.
  Underflows and overflows are included in the first and final bins of the output histogram by default.
  Note: 'process_variable' is a list of events
  '''
  scaling = 1 if "Data" in process_name else set_MC_process_info(process_name, luminosity, scaling=True)[2]
  weights = scaling*process_weights
  underflow, overflow = calculate_underoverflow(process_variable, xbins, weights)
  binned_values, _    = np.histogram(process_variable, xbins, weights=weights)
  binned_values[0]   += underflow
  binned_values[-1]  += overflow
  return binned_values


def accumulate_MC_subprocesses(parent_process, process_dictionary):
  '''
  Add up separate MC histograms for processes belonging to the same family.
  For example, with three given inputs of the same family, the output is the final line:
    WWToLNu2Q = [0.0, 1.0, 5.5, 0.5]
    WZTo2L2Nu = [0.0, 2.0, 7.5, 0.2]
    ZZTo4L    = [0.0, 3.0, 4.5, 0.1]
    --------------------------------
    VV        = [0.0, 6.0, 17.5, 0.8]
  Inputs not belonging to the specified 'parent_process' are ignored,
  therefore, this function is called once for each parent process
  '''
  accumulated_values = 0
  for MC_process in process_dictionary:
    if get_parent_process(MC_process) == parent_process:
      accumulated_values += process_dictionary[MC_process]["BinnedEvents"]
  return accumulated_values

def accumulate_datasets(dataset_dictionary):
  '''
  Very similar to accumulate_MC_subproceses
  Written to add datasets (Muon, Tau, EGamma, MuonEG) together
  '''
  accumulated_values = 0
  for dataset in dataset_dictionary:
    accumulated_values += dataset_dictionary[dataset]["BinnedEvents"]

  return accumulated_values


def get_parent_process(MC_process):
  '''
  Given some process, return a corresponding parent_process, effectively grouping
  related processes (i.e. DYInclusive, DY1, DY2, DY3, and DY4 all belong to DY).
  TODO: simplify this code, it is currently written in a brain-dead way
  '''
  parent_process = ""
  if   "DY"    in MC_process:  parent_process = "DY"
  elif "WJets" in MC_process:  parent_process = "WJ"
  elif "TT"    in MC_process:  parent_process = "TT"
  elif "ST"    in MC_process:  parent_process = "ST"
  elif ("WW"   in MC_process or 
        "WZ"   in MC_process or 
        "ZZ"   in MC_process): parent_process = "VV"
  else:
    print("No matching parent process for {MC_process}, continuing as individual process...")
  return parent_process


def set_good_events(final_state_mode, disable_triggers=False):
  '''
  Return a string defining a 'good_events' flag used by uproot to preskim input events
  to only those passing these simple requirements. 'good_events' changes based on
  final_state_mode, and the trigger condition is removed if a trigger study is 
  being conducted (since requiring the trigger biases the study).
  '''
  good_events = ""
  if disable_triggers: print("*"*20 + " removed trigger requirement " + "*"*20)

  # relevant definitions from NanoTauAnalysis /// modules/TauPairSelector.py
  # HTT_SRevent and HTT_ARevent require opposite sign objects
  # HTT_SRevent = ((pdgIdPair < 0) 
  #            and ( ((LeptonIso < 0.2) and (abs(pdgIdPair)==11*13)) or (LeptonIso < 0.15)) 
  #            and TauPassVsJet and (self.leptons[finalpair[1]].pt > 15))
  # HTT_ARevent = ((pdgIdPair < 0) 
  #            and ( ((LeptonIso < 0.2) and (abs(pdgIdPair)==11*13)) or (LeptonIso < 0.15)) 
  #            and (not TauPassVsJet) and (self.leptons[finalpair[1]].pt > 15))
  #     # All SR requirements besides TauPassVsJet
  # HTT_SSevent = ((pdgIdPair > 0) 
  #            and ( ((LeptonIso < 0.2) and (abs(pdgIdPair)==11*13)) or (LeptonIso < 0.15)) 
  #            and TauPassVsJet and (self.leptons[finalpair[1]].pt > 15)) 
  #     # All SR requirements besides opposite sign
  
  # apply FS cut separately so it can be used with reject_duplicate_events
  if final_state_mode == "ditau":
    good_events = "(HTT_SRevent) & (METfilters) & (LeptonVeto==0) & (abs(HTT_pdgId)==15*15) & (Trigger_ditau)"
    if disable_triggers: good_events = good_events.replace(" & (Trigger_ditau)", "")

  elif final_state_mode == "mutau":
    good_events = "(HTT_SRevent) & (METfilters) & (LeptonVeto==0) & (abs(HTT_pdgId)==13*15) & (Trigger_mutau)"
    if disable_triggers: good_events = good_events.replace(" & (Trigger_mutau)", "")

  elif final_state_mode == "etau":
    good_events = "(HTT_SRevent) & (METfilters) & (LeptonVeto==0) & (abs(HTT_pdgId)==11*15) & (Trigger_etau)"
    if disable_triggers: good_events = good_events.replace(" & (Trigger_etau)", "")

  elif final_state_mode == "dimuon":
    # lepton veto must be applied manually for this final state
    good_events = "(METfilters) & (HTT_pdgId==-13*13) & (HLT_IsoMu24)"
    if disable_triggers: good_events = good_events.replace(" & (HLT_IsoMu24)", "")

  return good_events


