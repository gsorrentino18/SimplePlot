import numpy as np

### README
# this file contains functions to perform simple calculations and return or print the result


def calculate_underoverflow(events, xbins, weights):
  '''
  Count the number of events falling outside (below and above) the specified bins. 
  For data, an array of ones is passed to the 'weights' variable.
  For MC, event weights must be passed correctly when the function is called.
  '''
  count_bin_values = [-999999., xbins[0], xbins[-1], 999999.]
  values, bins = np.histogram(events, count_bin_values, weights=weights)
  underflow_value, overflow_value = values[0], values[-1]
  return underflow_value, overflow_value


def calculate_signal_background_ratio(data, backgrounds, signals):
  '''
  Calculate and display signal-to-background ratios.
  '''
  yields = [np.sum(data)]
  total_background, total_signal = 0, 0
  for process in backgrounds: 
    process_yield = np.sum(backgrounds[process]["BinnedEvents"])
    yields.append(process_yield)
    total_background += process_yield
  for signal in signals:
    signal_yield = np.sum(backgrounds[process]["BinnedEvents"])
    yields.append(signal_yield)
    total_signal += signal_yield

  print("signal-to-background information")
  print(f"S/B      : {total_signal/total_background:.3f}")
  print(f"S/(S+B)  : {total_signal/(total_signal+total_background):.3f}")
  print(f"S/√(B)   : {total_signal/np.sqrt(total_background):.3f}")
  print(f"S/√(S+B) : {total_signal/np.sqrt(total_signal+total_background):.3f}")


def calculate_mt(lep_pt, lep_phi, MET_pt, MET_phi):
  '''
  Calculates the experimental paricle physicist variable of "transverse mass"
  which is a measure a two-particle system's mass when known parts (neutrinos)
  are missing. 
  Notably, there is another variable called "transverse mass" which is what
  ROOT.Mt() calculates. This is not the variable we are interested in and instead
  calculate the correct transverse mass by hand. Either form below is equivalenetly valid.
  '''
  # used in mutau, etau, emu
  #delta_phi = phi_mpi_pi(lep_phi - MET_phi)
  #mt = np.sqrt(2 * lep_pt * MET_pt * (1 - np.cos(delta_phi) ) ) 
  #alternate calculation, same quantity up to 4 decimal places
  lep_x = lep_pt*np.cos(lep_phi)
  lep_y = lep_pt*np.sin(lep_phi)
  MET_x = MET_pt*np.cos(MET_phi)
  MET_y = MET_pt*np.sin(MET_phi)
  sum_pt_2  = (lep_pt + MET_pt) * (lep_pt + MET_pt)
  sum_ptx_2 = (lep_x + MET_x) * (lep_x + MET_x)
  sum_pty_2 = (lep_y + MET_y) * (lep_y + MET_y)
  mt_2      = (sum_pt_2 - sum_ptx_2 - sum_pty_2)
  if mt_2 < 0:
    # this is floating-point calculation error, all values with 0.01 of 0
    #print(mt_2, sum_pt_2, sum_ptx_2, sum_pty_2)
    mt = 0
  else: 
    mt = np.sqrt(mt_2) 
  return mt

#def ROOT_mt(lep_pt, lep_eta, lep_phi, lep_mass, MET_pt, MET_phi): #import ROOT to use this function
#  Lep_vec = ROOT.TLorentzVector()
#  Lep_vec.SetPtEtaPhiE(lep_pt, lep_eta, lep_phi, lep_mass)
#  MET_vec = ROOT.TLorentzVector(MET_pt, MET_phi)
#  ROOT_mt = (Lep_vec + MET_vec).Mt()
#  return ROOT_mt
  

def calculate_dR(eta1, phi1, eta2, phi2): 
  '''return value of delta R cone defined by two objects'''
  delta_eta = eta1-eta2
  delta_phi = phi_mpi_pi(lep_phi - MET_phi)
  return np.sqrt(delta_eta*delta_eta + delta_phi*delta_phi)


def phi_mpi_pi(delta_phi):
  '''return phi between a range of negative pi and pi'''
  return 2 * np.pi - delta_phi if delta_phi > np.pi else 2 * np.pi + delta_phi if delta_phi < -1*np.pi else delta_phi


def yields_for_CSV(histogram_axis, desired_order=[]):
    handles, labels = histogram_axis.get_legend_handles_labels()
    desired_order    = labels if desired_order == [] else desired_order
    reordered_labels = []
    corresponding_yields = []
    for compare_label in desired_order:
      for original_label in labels:
        if compare_label in original_label:
          reordered_labels.append(original_label)
          label_yield_start = original_label.find("[")
          label_yield_end   = original_label.find("]")
          label_yield       = original_label[label_yield_start+1:label_yield_end]
          corresponding_yields.append(int(label_yield))
    print(f"Legend Labels: {reordered_labels}")
    print(f"Yields: {corresponding_yields}")


