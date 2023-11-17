from XSec import XSecRun3 as XSec

### README ###
# this file contains small dictionaries combining information for MC processes 
# MC_dictionary SHOULD contain one entry for each type of sample, it should    
# it references the smaller dictionaries color_dictionary and label_dictionary 
# for additionaly information. Cross-sections are read in from XSec.py, and
# NWEvents are hand-coded (to be removed later). Plot scaling for all 
# non-signal processes are set to 1 by default.

color_dictionary = {
  "ggH": "#0000ff", # blue
  "VBF": "#ff0000", # red
  "DY" : "#ffc942", # DY yellow
  "TT" : "#9999cc", # light purple 
  "ST" : "#660099", # dark purple
  "WJ" : "#e44e4e", # dark red
  "VV" : "#808080", # grey
  "QCD": "#ffc9ff", # pink
}
label_dictionary = {
 "ggH" : "ggH",
 "VBF" : "VBF",
 "DY"  : "DY",
 "TT"  : "TT",
 "ST"  : "ST",
 "WJ"  : "WJ",
 "VV"  : "VV",
 "QCD" : "QCD",
}
 
MC_dictionary = {
  "DY" : {"label": "DY", "color": color_dictionary["DY"]},
  "TT" : {"label": "TT", "color": color_dictionary["TT"]},
  "ST" : {"label": "ST", "color": color_dictionary["ST"]},
  "WJ" : {"label": "WJ", "color": color_dictionary["WJ"]},
  "VV" : {"label": "VV", "color": color_dictionary["VV"]},

  #########################################################
  # above are dummy dictionaries for grouped subprocesses #
  # below are real dictionaries for subprocesses          #
  #########################################################

  "ggH" : {"XSec": XSec["ggH_TauTau"], "NWevents": 19289464.362, 
           "label": "ggH", "color": color_dictionary["ggH"],
           "plot_scaling" : 100},

  "VBF" : {"XSec": XSec["VBF_TauTau"], "NWevents": 2402853.147599998, 
           "label": "VBF", "color": color_dictionary["VBF"],
           "plot_scaling" : 500},

  #TODO: This is for 2022 pre EE dimuon! Somehow store this info and make it switchable!
  #"DYInc" : {"XSec": XSec["DYJetsToLL_M-50"], "NWevents": 1925298812420.0,
  #"DYInc" : {"XSec": XSec["DYJetsToLL_M-50"], "NWevents": 6534910623420.0, # new DYInc MC from Dennis (???)
  "DYInc" : {"XSec": XSec["DYJetsToLL_M-50"], "NWevents": 3222224508811.3296, # old August samples
             "label": "DY", "color": color_dictionary["DY"],
             "plot_scaling" : 1.12}, # k-factor to NNLO

  "TTTo2L2Nu"         : {"XSec": XSec["TTTo2L2Nu"], "NWevents": 6750783685.336004,
                         "label": "TT", "color": color_dictionary["TT"],
                         "plot_scaling" : 1},
  "TTToFullyHadronic" : {"XSec": XSec["TTToFullyHadronic"], "NWevents": 62787618087.51002,
                         "label": "TT", "color": color_dictionary["TT"],
                         "plot_scaling" : 1},
  "TTToSemiLeptonic"  : {"XSec": XSec["TTToSemiLeptonic"], "NWevents": 89738778580.06801,
                         "label": "TT", "color": color_dictionary["TT"],
                         "plot_scaling" : 1},

  "ST_TWminus_2L2Nu"   : {"XSec": XSec["ST_TWminus_2L2Nu"], "NWevents": 31630228.619640004,
                          "label": "ST", "color": color_dictionary["ST"],
                          "plot_scaling" : 1},
  "ST_TbarWplus_2L2Nu" : {"XSec": XSec["ST_TbarWplus_2L2Nu"], "NWevents": 32536952.371500004,
                          "label": "ST", "color": color_dictionary["ST"],
                          "plot_scaling" : 1},
 
  "ST_TWminus_4Q"      : {"XSec": XSec["ST_TWminus_4Q"], "NWevents": 221087105.7696,
                          "label": "ST", "color": color_dictionary["ST"],
                          "plot_scaling" : 1},
  "ST_TbarWplus_4Q"    : {"XSec": XSec["ST_TbarWplus_4Q"], "NWevents": 219993701.15240002,
                          "label": "ST", "color": color_dictionary["ST"],
                          "plot_scaling" : 1},
 
  "ST_TWminus_LNu2Q"   : {"XSec": XSec["ST_TWminus_LNu2Q"], "NWevents": 612917819.6544,
                          "label": "ST", "color": color_dictionary["ST"],
                          "plot_scaling" : 1},
  "ST_TbarWplus_LNu2Q" : {"XSec": XSec["ST_TbarWplus_LNu2Q"], "NWevents": 605855578.9104003,
                          "label": "ST", "color": color_dictionary["ST"],
                          "plot_scaling" : 1},

  "WJetsInc"      : {"XSec": XSec["WJetsToLNu"], "NWevents": 43876130657994.01,
                     "label": "WJ", "color": color_dictionary["WJ"],
                     "plot_scaling" : 1},
  "WJetsToLNu_1J" : {"XSec": XSec["WJetsToLNu"], "NWevents": 1302140190027.9631,
                     "label": "WJ", "color": color_dictionary["WJ"],
                     "plot_scaling" : 1},
  "WJetsToLNu_2J" : {"XSec": XSec["WJetsToLNu"], "NWevents": -999999.99999, # sample DNE
                     "label": "WJ", "color": color_dictionary["WJ"],
                     "plot_scaling" : 1},
  "WJetsToLNu_3J" : {"XSec": XSec["WJetsToLNu"], "NWevents": 320401283585.0479,
                     "label": "WJ", "color": color_dictionary["WJ"],
                     "plot_scaling" : 1},
  "WJetsToLNu_4J" : {"XSec": XSec["WJetsToLNu"], "NWevents": 35153418064.277626,
                     "label": "WJ", "color": color_dictionary["WJ"],
                     "plot_scaling" : 1},

  "WWTo2L2Nu"    : {"XSec": XSec["WWTo2L2Nu"], "NWevents": 264682597.7669993,
                    "label": "VV", "color": color_dictionary["VV"],
                    "plot_scaling" : 1},
  "WWTo4Q"       : {"XSec": XSec["WWTo4Q"], "NWevents": 4983316073.657991,
                    "label": "VV", "color": color_dictionary["VV"],
                    "plot_scaling" : 1},
  # not in XSec yet
  #"WWToLNu2Q"    : {"XSec": XSec["WWToLNu2Q"], "NWevents": 0,
  #                  "label": "VV", "color": color_dictionary["VV"],
  #                  "plot_scaling" : 1},

  "WZTo3L1Nu"    : {"XSec": XSec["WZTo3L1Nu"], "NWevents": 48195008.27371989,
                    "label": "VV", "color": color_dictionary["VV"],
                    "plot_scaling" : 1},
  "WZTo2L2Q"     : {"XSec": XSec["WZTo2L2Q"], "NWevents": 113225054.20599996, 
                    "label": "VV", "color": color_dictionary["VV"],
                    "plot_scaling" : 1},
  "WZToLNu2Q"    : {"XSec": XSec["WZToLNu2Q"], "NWevents": 478386870.3699988,
                    "label": "VV", "color": color_dictionary["VV"],
                    "plot_scaling" : 1},

  "ZZTo2L2Nu"    : {"XSec": XSec["ZZTo2L2Nu"], "NWevents": 51651062.90956998,
                    "label": "VV", "color": color_dictionary["VV"],
                    "plot_scaling" : 1},
  "ZZTo2L2Q"     : {"XSec": XSec["ZZTo2L2Q"], "NWevents": 343533147.0025,
                    "label": "VV", "color": color_dictionary["VV"],
                    "plot_scaling" : 1},
  "ZZTo2Nu2Q"    : {"XSec": XSec["ZZTo2Nu2Q"], "NWevents": 0,
                    "label": "VV", "color": color_dictionary["VV"],
                    "plot_scaling" : 1},
  "ZZTo4L"       : {"XSec": XSec["ZZTo4L"], "NWevents": 69812228.61073996,
                    "label": "VV", "color": color_dictionary["VV"],
                    "plot_scaling" : 1},

  "QCD"   : {"XSec": 1, "NWevents": 1,
             "label": "QCD", "color": color_dictionary["QCD"],
             "plot_scaling" : 1}, 

}
