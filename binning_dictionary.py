import numpy as np

### README ###
# binning for different variables are defined below and are separated by use-case
# All variables are assumed to be linearly binned.

label_dictionary = {
  "FS_t1_pt" : r'Leading Tau $p_T$ [GeV]',
  "FS_t1_eta" : r'Leading Tau $\eta$',
  "FS_t1_phi" : r'Leading Tau $\phi$',
  "FS_t1_DeepTauVSjet" : r'Leading Tau DeepTau Vs Jet',
  "FS_t1_DeepTauVSmu"  : r'Leading Tau DeepTau Vs Muon',
  "FS_t1_DeepTauVSe"   : r'Leading Tau DeepTau Vs Electron',
  "FS_t1_dxy"  : r'Leading Tau $D_{xy}$',
  "FS_t1_dz"   : r'Leading Tau $D_Z$',
  "FS_t1_chg"  : r'Leading Tau Charge',
  "FS_t1_DM"   : r'Leading Tau Decay Mode',

  "FS_t2_pt" : r'Sub-leading Tau $p_T$ [GeV]',
  "FS_t2_eta" : r'Sub-leading Tau $\eta$',
  "FS_t2_phi" : r'Sub-leading Tau $\phi$',
  "FS_t2_DeepTauVSjet" : r'Sub-leading Tau DeepTau Vs Jet',
  "FS_t2_DeepTauVSmu"  : r'Sub-leading Tau DeepTau Vs Muon',
  "FS_t2_DeepTauVSe"   : r'Sub-leading Tau DeepTau Vs Electron',
  "FS_t2_dxy"  : r'Sub-leading Tau $D_{xy}$',
  "FS_t2_dz"   : r'Sub-leading Tau $D_Z$',
  "FS_t2_chg"  : r'Sub-leading Tau Charge',
  "FS_t2_DM"   : r'Sub-leading Tau Decay Mode',
   

  "MET_pt"      : r'MET [GeV]',
  "PuppiMET_pt" : r'PUPPI MET [GeV]',
  "nCleanJetGT30"   : r'Number of Jets',
  #"CleanJet_pt" : '',
  #"CleanJet_eta": '',
  "HTT_DiJet_MassInv_fromHighestMjj" : r'Dijet Mass [GeV]',
  "HTT_DiJet_dEta_fromHighestMjj"    : r'|$\Delta\eta$|',
  "HTT_H_pt_using_PUPPI_MET"         : r'Higgs $p_T$ [GeV]',
  "HTT_dR"     : r'$\Delta$R',
  "HTT_m_vis-KSUbinning" : r'$m_{vis}$ [GeV]',
  "HTT_m_vis-SFbinning" : r'$m_{vis}$ [GeV]',
  "HTT_pT_l1l2" : r'$p_T^{ll}$',
  "FastMTT_PUPPIMET_mT" : r'Fast MTT Transverse Mass [GeV]',
  "FastMTT_PUPPIMET_mass" : r'Fast MTT Mass [GeV]',
  "PV_npvs"    : r'Number of Primary Vertices',

  "FS_el_pt"   : r'Ele $p_T$ [GeV]',
  "FS_el_eta"  : r'Ele $\eta$',
  "FS_el_phi"  : r'Ele $\phi$',
  "FS_el_iso"  : r'Ele Isolation',
  "FS_el_dxy"  : r'Ele $D_{xy}$',
  "FS_el_dz"   : r'Ele $D_{Z}$',
  "FS_el_chg"  : r'Electron charge',

  "FS_el_hoe"      : r'Ele H/E',
  "FS_el_lostHits" : r'Ele Lost Hits',
  "FS_el_mass"     : r'Elecron mass',
  "FS_el_r9"      : r'Electron R9',
  "FS_el_sieie"    : r'Ele $\sigma_{i \eta i \eta}$',
  "FS_el_isochg"   : r'Ele Charged Isolation',
  "FS_el_dr03TkSumPt"         : r'Ele dr03TkSumPt',
  "FS_el_dr03EcalRecHitSumEt" : r'Ele dr03EcalRecHitSumEt',
  "FS_el_dr03HcalDepth1TowerSumEt" : r'Ele dr03HcalDepth1TowerSumEt',

  "FS_tau_pt"  : r'Tau $p_T$ [GeV]',
  "FS_tau_eta" : r'Tau $\eta$',
  "FS_tau_phi" : r'Tau $\phi$',
  "FS_tau_dxy" : r'Tau $D_{xy}$',
  "FS_tau_dz"  : r'Tau $D_{Z}$',
  "FS_tau_chg" : r'Tau charge',
  "FS_tau_DM"  : r'Tau DM',

  "FS_mt"         : r'Transverse mass',
  "HTT_dR"     : r'$D_{R}$',
  "HTT_m_vis-KSUbinning" : r'$m_{vis}$ [GeV]',
  "HTT_m_vis-SFbinning" :  r'$m_{vis}$ [GeV]' ,
  "HTT_Lep_iso" : r'Lepton isolation'

}
 
binning_dictionary = {
#  var  : (nBins, xmin, xmax),

#  ditau
  "FS_t1_pt"   : (36, 0, 180),
  "FS_t1_eta"  : (30, -3, 3),
  "FS_t1_phi"  : (32, -3.2, 3.2),
  "FS_t1_DeepTauVSjet" : (8, 1, 9),
  "FS_t1_DeepTauVSmu"  : (4, 1, 5),
  "FS_t1_DeepTauVSe"   : (8, 1, 9),
  "FS_t1_dxy"  : (50, 0, 0.20),
  "FS_t1_dz"   : (50, 0, 0.25),
  "FS_t1_chg"  : (5, -2, 2),
  "FS_t1_DM"   : (20, 0, 19),

  "FS_t2_pt"   : (24, 0, 120),
  "FS_t2_eta"  : (30, -3, 3),
  "FS_t2_phi"  : (32, -3.2, 3.2),
  "FS_t2_DeepTauVSjet" : (8, 1, 9),
  "FS_t2_DeepTauVSmu"  : (4, 1, 5),
  "FS_t2_DeepTauVSe"   : (8, 1, 9),
  "FS_t2_dxy"  : (50, 0, 0.20),
  "FS_t2_dz"   : (50, 0, 0.25),
  "FS_t2_chg"  : (5, -2, 2),
  "FS_t2_DM"   : (20, 0, 19),

#  mutau/etau
  "FS_mu_pt"   : (40, 0, 120),
  "FS_mu_eta"  : (30, -3, 3),
  "FS_mu_phi"  : (32, -3.2, 3.2),
  "FS_mu_iso"  : (25, 0, 1),
  "FS_mu_dxy"  : (50, 0, 0.05),
  "FS_mu_dz"   : (50, 0, 0.25),
  "FS_mu_chg"  : (5, -2, 2),

  "FS_el_pt_newbins"   : (80, 26, 46), #40, 0, 120
  "FS_el_pt"   : (40, 20, 100), #75, 25, 120
  "FS_el_eta"  : (30, -3, 3),
  "FS_el_phi"  : (32, -3.2, 3.2),
  "FS_el_iso"  : (30, 0, 0.3),
  "FS_el_dxy"  : (30, 0, 0.15),
  "FS_el_dz"   : (15, 0, 0.15),
  "FS_el_chg"  : (5, -2, 2),

  "FS_el_hoe"      : (40,0,0.02),
  "FS_el_lostHits" : (9,0,8), 
  "FS_el_mass"     : (40, -0.002, 0.002),
  "FS_el_r9"      : (25, 0.2, 1.2),
  "FS_el_sieie"    : (30, 0, 0.05),
  "FS_el_isochg"   : (30, 0, 0.3),
  "FS_el_dr03TkSumPt"         : (40, 0,0.4),
  "FS_el_dr03EcalRecHitSumEt" : (30, 0, 0.3),
  "FS_el_dr03HcalDepth1TowerSumEt" : (20, 0, 1.0),

  "FS_tau_pt"  : (50, 20, 100),#36, 0, 180
  "FS_tau_pt_newbins"  : (60, 25, 40), #36, 0, 180
  "FS_tau_eta" : (30, -3, 3),
  "FS_tau_phi" : (32, -3.2, 3.2),
  "FS_tau_dxy" : (50, 0, 0.20),
  "FS_tau_dz"  : (50, 0, 0.25),
  "FS_tau_chg" : (5, -2, 2),
  "FS_tau_DM"  : (20, 0, 19),

  "FS_tau_rawPNetVSjet" : (50, 0, 1),
  "FS_tau_rawPNetVSmu"  : np.array([0, 0.95, 0.96, 0.97, 0.98, 0.99, 1]), # plot logx
  "FS_tau_rawPNetVSe"   : np.array([0, 0.95, 0.96, 0.97, 0.98, 0.99, 1]),

  "FS_t1_rawPNetVSjet" : (50, 0, 1),
  "FS_t1_rawPNetVSmu"  : np.array([0, 0.95, 0.96, 0.97, 0.98, 0.99, 1]), # plot logx
  "FS_t1_rawPNetVSe"   : np.array([0, 0.95, 0.96, 0.97, 0.98, 0.99, 1]),

  "FS_t2_rawPNetVSjet" : (50, 0, 1),
  "FS_t2_rawPNetVSmu"  : np.array([0, 0.95, 0.96, 0.97, 0.98, 0.99, 1]), # plot logx
  "FS_t2_rawPNetVSe"   : np.array([0, 0.95, 0.96, 0.97, 0.98, 0.99, 1]),

  "mutau_TnP" : {
    "FS_mu_pt"   : (40, 0, 120),
    "FS_mu_eta"  : (7, -2.7, 2.7),
    "FS_mu_phi"  : (16, -3.2, 3.2),
 
    "FS_tau_pt" : np.array([0, 20, 25, 30, 35, 40, 50, 70, 150]),
    "FS_tau_eta"  : (7, -2.7, 2.7),
    "FS_tau_phi" : (16, -3.2, 3.2),
  },
  "pass_tag"   : (2, 0, 2),
  "pass_probe" : (2, 0, 2),

# dimuon
  "FS_m1_pt"   : (60, 0, 300),
  "FS_m1_eta"  : (99, -2.5, 2.5),
  "FS_m1_phi"  : (64, -3.2, 3.2),
  "FS_m1_iso"  : (25, 0, 1),
  "FS_m1_dxy"  : (50, 0, 0.05),
  "FS_m1_dz"   : (50, 0, 0.25),
  "FS_m1_chg"  : (5, -2, 2),

  "FS_m2_pt"   : (60, 0, 300),
  "FS_m2_eta"  : (99, -2.5, 2.5),
  "FS_m2_phi"  : (64, -3.2, 3.2),
  "FS_m2_iso"  : (25, 0, 1),
  "FS_m2_dxy"  : (50, 0, 0.05),
  "FS_m2_dz"   : (50, 0, 0.25),
  "FS_m2_chg"  : (5, -2, 2),
  "FS_m_vis_tight" : (80, 70, 110), # TODO add this

# common, calculated on the fly
  "FS_mt"         : (40, 0, 200),
  "nCleanJetGT30" : (8, 0, 8), # GT(E) = Greater Than (Equal to)
  "CleanJetGT30_pt_1"  : (60, 0, 300),
  "CleanJetGT30_pt_2"  : (60, 0, 300),
  "CleanJetGT30_pt_3"  : (60, 0, 300),
  "CleanJetGT30_eta_1" : (50, -5, 5),
  "CleanJetGT30_eta_2" : (50, -5, 5),
  "CleanJetGT30_eta_3" : (50, -5, 5),
  "CleanJetGT30_phi_1" : (16, -3.2, 3.2),
  "CleanJetGT30_phi_2" : (16, -3.2, 3.2),
  "CleanJetGT30_phi_3" : (16, -3.2, 3.2),
  "FS_mjj"    : (30, 0, 1500),
  "FS_detajj" : (31, 0, 7),

# common, from branches
  "MET_pt"      : (30, 0, 150),
  "PuppiMET_pt" : (50, 0, 150),
  "nCleanJet"   : (8, 0, 8),
  "CleanJet_pt" : (30, 20, 200),
  "CleanJet_eta": (50, -5, 5),
  "HTT_DiJet_MassInv_fromHighestMjj" : (30, 0, 1500),
  "HTT_DiJet_dEta_fromHighestMjj"    : (35, 0, 7),
  "HTT_H_pt_using_PUPPI_MET"         : (30, 0, 300),
  "HTT_dR"     : (60, 0, 6),
  "HTT_m_vis-KSUbinning" : (40, 0, 200),
  "HTT_m_vis-SFbinning"  : (12, 60, 120), #40,0,200
  "HTT_Lep_iso" : (25, 0, 1),
  "HTT_pT_l1l2" : (30, 0, 150),
  "FastMTT_PUPPIMET_mT" : (40, 0, 400),
  "FastMTT_PUPPIMET_mass" : (20, 0, 400),
  "FS_t1_flav" : (11, 0, 11),
  "FS_t2_flav" : (11, 0, 11),
  "PV_npvs"    : (30, 0, 90),
}
