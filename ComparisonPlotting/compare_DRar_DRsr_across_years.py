import numpy as np
from PIL import ImageTk, Image
import matplotlib.pyplot as plt


variablesSR = ["QCD_ditau_region_DRsr_0j_set1_tau1pt",
               "QCD_ditau_region_DRsr_1j_set1_tau1pt",
               "QCD_ditau_region_DRsr_2j_set1_tau1pt"]

variablesAR = ["QCD_ditau_region_DRar_0j_set1_tau1pt", 
               "QCD_ditau_region_DRar_1j_set1_tau1pt", 
               "QCD_ditau_region_DRar_2j_set1_tau1pt"]

for variableSR, variableAR in zip(variablesSR, variablesAR):
  print(f"setting up plots for {variableSR} and {variableAR}")
  set1 = np.array(Image.open("QCD2018_DRsr/"+variableSR+".png"))
  set2 = np.array(Image.open("QCD2022postEE_DT2p1_DRsr/"+variableSR+".png"))
  set3 = np.array(Image.open("QCD2022postEE_DT2p5_DRsr/"+variableSR+".png"))

  set4 = np.array(Image.open("QCD2018_DRar/"+variableAR+".png"))
  set5 = np.array(Image.open("QCD2022postEE_DT2p1_DRar/"+variableAR+".png"))
  set6 = np.array(Image.open("QCD2022postEE_DT2p5_DRar/"+variableAR+".png"))

  images = [set1, set2, set3, set4, set5, set6]
  image_titles = ["2018, DeepTau 2p1, Numerator", 
                  "2022, DeepTau 2p1, Numerator", 
                  "2022, DeepTau 2p5, Numerator",
                  "2018, DeepTau 2p1, Denominator",
                  "2022, DeepTau 2p1, Denominator",
                  "2022, DeepTau 2p5, Denominator",
                 ]

  fig = plt.figure(figsize=(9,9))
  for i in range(0,len(images)):
    ax = fig.add_subplot(2, 3, i+1)
    ax.title.set_text(image_titles[i])
    plt.imshow(images[i])
    plt.axis('off')
    plt.tight_layout()
  variable = variableSR.replace("QCD_ditau_region_DRsr_","")
  variable = variableSR.replace("_set1_tau1pt","")
  filename = "compare-six-regions"+variable
  print(f"saving as {filename}")
  plt.savefig(filename+".png")

print("finished")
