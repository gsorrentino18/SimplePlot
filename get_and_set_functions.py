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


