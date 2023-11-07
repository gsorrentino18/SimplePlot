import numpy as np
from PIL import ImageTk, Image
import matplotlib.pyplot as plt

variables = ["QCD_FF_0j", "QCD_FF_1j", "QCD_FF_2j",
             "QCD_FF_0j_aiso", "QCD_FF_1j_aiso", "QCD_FF_2j_aiso"]

for variable in variables:
  print(f"setting up plots for {variable}")
  set1 = np.array(Image.open("QCD2018/"+variable+".png"))
  set2 = np.array(Image.open("QCD2022postEE_DT2p1/"+variable+".png"))
  set3 = np.array(Image.open("QCD2022postEE_DT2p5/"+variable+".png"))

  images = [set1, set2, set3]
  image_titles = ["2018, DeepTau 2p1", 
                  "2022, DeepTau 2p1", 
                  "2022, DeepTau 2p5"]

  fig = plt.figure(figsize=(9,9))
  for i in range(0,len(images)):
    ax = fig.add_subplot(1, 3, i+1)
    ax.title.set_text(image_titles[i])
    plt.imshow(images[i])
    plt.axis('off')
    plt.tight_layout()
  filename = "compare-three-"+variable
  print(f"saving as {filename}")
  plt.savefig(filename+".png")

print("finished")
