### README ###
# binning for different variables are defined below and are separated by use-case
# All variables are assumed to be linearly binned.
 
binning_dictionary = {
#  var  : (nBins, xmin, xmax),

#  ditau
  "FS_t1_pt"   : (36, 0, 180),
  "FS_t1_eta"  : (30, -3, 3),
  "FS_t2_pt"   : (24, 0, 120),
  "FS_t2_eta"  : (30, -3, 3),

#  mutau/etau
  "FS_mu_pt"   : (40, 0, 120),
  "FS_mu_eta"  : (30, -3, 3),
  "FS_el_pt"   : (40, 0, 120),
  "FS_el_eta"  : (30, -3, 3),
  "FS_tau_pt"  : (36, 0, 180),
  "FS_tau_eta" : (30, -3, 3),
  #"dummy_HTT_Lep_pt" : (40, 0, 120),
  #"dummy_HTT_Tau_pt" : (30, 0, 180),

# common, calculated on the fly
  "FS_mt"        : (20, 0, 200),
  "nCleanJetGT30" : (8, 0, 8), # GT = Greater Than 

# common, from branches
  "MET_pt"      : (30, 0, 150),
  "PuppiMET_pt" : (30, 0, 150),
  "nCleanJet"   : (8, 0, 8),
  "CleanJet_pt" : (30, 20, 200),
  "CleanJet_eta": (50, -5, 5),
  "HTT_DiJet_MassInv_fromHighestMjj" : (30, 0, 1500),
  "HTT_DiJet_dEta_fromHighestMjj"    : (35, 0, 7),
  "HTT_H_pt_using_PUPPI_MET"         : (30, 0, 300),
  "HTT_dR"     : (60, 0, 6),
  "HTT_m_vis" : (30, 0, 300), # this is m_lep
}
