### README ###
# This is a mapping of process name (shared with XSec.py) to wildcards for related samples.
# The wildcarding works for the 'concatenate' of uproot, and might not in the future.
# The testing_file_map is a reduced number of files for faster processing times.

testing_file_map = {
  "DataTau"  : "Data/Tau*G*",
  "DataMuon" : "Data/Muon*G*",

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

full_file_map = {
  "DataTau"  : "Data/Tau*",
  "DataMuon" : "Data/Muon*",

  "DYInc" : "DY/DY*part*",
  "DYJetsToLL_M-50_1J" : "DY/DYJetsToLL_M-50_1J*",
  "DYJetsToLL_M-50_2J" : "DY/DYJetsToLL_M-50_2J*",
  "DYJetsToLL_M-50_3J" : "DY/DYJetsToLL_M-50_3J*",
  "DYJetsToLL_M-50_4J" : "DY/DYJetsToLL_M-50_4J*",

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
  "WJetsToLNu_1J" : "WJ/WJetsToLNu_1J*",
  "WJetsToLNu_2J" : "WJ/WJetsToLNu_2J*",
  "WJetsToLNu_3J" : "WJ/WJetsToLNu_3J*",
  "WJetsToLNu_4J" : "WJ/WJetsToLNu_4J*",

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
 
