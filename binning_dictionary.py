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

# dimuon
  "FS_m1_pt"   : (40, 0, 120),
  "FS_m2_pt"   : (40, 0, 120),
  "FS_m1_eta"  : (30, -3, 3),
  "FS_m2_eta"  : (30, -3, 3),

# common, calculated on the fly
  "FS_mt"        : (20, 0, 200),
  "nCleanJetGT30" : (8, 0, 8), # GT(E) = Greater Than (Equal to)
  #TODO: delete these
  "pass_0j_cuts"  : (100, 0, 100), 
  "pass_1j_cuts"  : (100, 0, 100),
  "pass_2j_cuts"  : (100, 0, 100), 
  "pass_3j_cuts"  : (100, 0, 100),
  "pass_GTE2j_cuts"    : (100, 0, 100),
  ##
  "CleanJetGT30_pt_1"  : (60, 0, 300),
  "CleanJetGT30_pt_2"  : (60, 0, 300),
  "CleanJetGT30_pt_3"  : (60, 0, 300),
  "CleanJetGT30_eta_1" : (50, -5, 5),
  "CleanJetGT30_eta_2" : (50, -5, 5),
  "CleanJetGT30_eta_3" : (50, -5, 5),

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
