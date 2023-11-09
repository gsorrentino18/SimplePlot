import uproot

branches = ["MET_pt"]

skim_subset = "(METfilters) & (abs(HTT_pdgId)==13*15)"
this_string = "(HTT_SRevent) & (METfilters) & (LeptonVeto==0) & (abs(HTT_pdgId)==13*15) & (Trigger_mutau)"
good_events = this_string

file = uproot.concatenate("../Run3FSSplitSamples/mutau/Data/Muon_*.root:Events",
                          branches, cut=good_events, library="np")

#nEvents = len(file["HTT_SRevent"])
nEvents = len(file["MET_pt"])
print(nEvents)
