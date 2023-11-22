from triggers_dictionary import triggers_dictionary

def add_trigger_branches(branches_, final_state_mode):
  '''
  Helper function to add HLT branches used by a given final state
  '''
  for trigger in triggers_dictionary[final_state_mode]:
    branches_.append(trigger)
  return branches_


def add_DeepTau_branches(branches_, DeepTauVersion):
  '''
  Helper function to add DeepTauID branches
  '''
  if DeepTauVersion == "2p1":
    for DeepTau_v2p1_branch in ["Tau_idDeepTau2017v2p1VSjet", "Tau_idDeepTau2017v2p1VSmu", "Tau_idDeepTau2017v2p1VSe"]:
      branches_.append(DeepTau_v2p1_branch)

  elif DeepTauVersion == "2p5":
    for DeepTau_v2p5_branch in ["Tau_idDeepTau2018v2p5VSjet", "Tau_idDeepTau2018v2p5VSmu", "Tau_idDeepTau2018v2p5VSe"]:
      branches_.append(DeepTau_v2p5_branch)

  else:
    print(f"no branches added with argument {DeepTauVersion}. Try 2p1 or 2p5.")

  return branches_


