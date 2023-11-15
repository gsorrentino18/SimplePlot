### README ###
# This file contains mappings of process names (shared with XSec.py) to wildcards for related samples.
# The :testing" file maps are subsets of full filelists for faster processing times.

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
  "DYInc"    : "DY/DYJets*part*",
}

testing_dimuon_file_map = {
  "DataMuon" : "Data/Muon*G*",
  #"DYInc"    : "DY/DYJets*part1",
  "DYInc1"    : "DY/DYJets*part1",
  "DYInc2"    : "DY/DYJets*part2",
  "DYInc3"    : "DY/DYJets*part3",
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
