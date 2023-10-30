### README
# This is a mapping of final_state_mode to HLT triggers included in branches "Trigger_[final_state_mode]"
# Very similar in design to "TriggerList.py" in NanoTauAnalysis.

triggers_dictionary = {
  "ditau" : [
             "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
             "HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1", # Run2, only for Era F study
             "HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1",  # Run2, only for Era F study
          ],
  "mutau" : [
             "HLT_IsoMu24", 
             "HLT_IsoMu27",
             "HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1", 
             "HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1", # Run2, only for Era F study
          ],
  "etau"  : [
             #"HLT_Ele30_WPTight_Gsf",
             "HLT_Ele32_WPTight_Gsf",
             "HLT_Ele35_WPTight_Gsf",
             "HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1",
          ],
  "dimuon": [
             "HLT_IsoMu24",
             "HLT_IsoMu27",
          ],
}
