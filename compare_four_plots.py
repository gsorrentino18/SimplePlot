import numpy as np
from PIL import ImageTk, Image
import matplotlib.pyplot as plt

variables = ["FS_mu_eta", "FS_tau_pt", "PuppiMET_pt", "HTT_m_vis",
             "FS_mu_pt",  "HTT_H_pt_using_PUPPI_MET", "HTT_mt",    
             "nCleanJet", "FS_tau_eta", "HTT_dR", "MET_pt"]

for variable in variables:
  print(f"setting up plots for {variable}")
  no_cut_img       = np.array(Image.open("mutau_plots_from_19-10_at_1327/"+variable+".png"))
  one_trig_cut_img = np.array(Image.open("mutau_plots_from_19-10_at_1335/"+variable+".png"))
  all_trig_cut_img = np.array(Image.open("mutau_plots_from_19-10_at_1342/"+variable+".png"))
  mt_cut_img       = np.array(Image.open("mutau_plots_from_19-10_at_1350/"+variable+".png"))

  images = [no_cut_img, one_trig_cut_img, all_trig_cut_img, mt_cut_img]
  image_titles = ["pass good events, no cuts", 
                  "tau plus single muon from trigger", 
                  "plus single muon from cross trigger", 
                  "plus mt cut"]

  fig = plt.figure(figsize=(9,9))
  for i in range(0,len(images)):
    ax = fig.add_subplot(2, 2, i+1)
    ax.title.set_text(image_titles[i])
    plt.imshow(images[i])
    plt.axis('off')
    plt.tight_layout()
  filename = "compare-four-"+variable
  print(f"saving as {filename}")
  plt.savefig(filename+".png")

print("finished")
