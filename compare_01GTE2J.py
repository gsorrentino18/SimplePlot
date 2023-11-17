import numpy as np
from PIL import ImageTk, Image
import matplotlib.pyplot as plt

variables = ["FS_t1_dxy", "FS_t1_eta", "FS_t1_pt", "FS_t2_dz",
             "FS_t2_phi", "HTT_dR", "nCleanJetGT30", "FS_t1_dz", "FS_t1_phi",
             "FS_t2_dxy", "FS_t2_eta", "FS_t2_pt", "HTT_m_vis"]

common_name = "FS_plots/2p1_best_QCD_correct_ratio_err_ditau_"

for variable in variables:
  print(f"setting up plots for {variable}")
  set1 = np.array(Image.open(common_name + "0j/"+variable+".png"))
  set2 = np.array(Image.open(common_name + "1j/"+variable+".png"))
  set3 = np.array(Image.open(common_name + "GTE2j/"+variable+".png"))

  images = [set1, set2, set3]
  image_titles = ["0j", "1j", "GTE2j"]

  fig = plt.figure(figsize=(15,5))
  for i in range(0,len(images)):
    ax = fig.add_subplot(1, 3, i+1)
    ax.title.set_text(image_titles[i])
    plt.imshow(images[i])
    plt.axis('off')
    plt.tight_layout()
  filename = "FS_plots/compare-jet-modes-"+variable
  print(f"saving as {filename}")
  plt.savefig(filename+".png")

print("finished")
