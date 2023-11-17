import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

width = 8
height = 16
angle1 = 55
angle2 = 125
pos = 1

name_top_L = "ditau"
name_top_R = "ditau+jet"
name_bot_L = "VBF ditau" # headache, how to delineate
name_bot_R = "tight VBF ditau"

#ellipse_(top/bottom)_(Left/Right)
ell_top_L = Ellipse(xy=[-pos,  pos], width=width, height=height, angle=angle1, linestyle="solid")
ell_top_R = Ellipse(xy=[pos,   pos], width=width, height=height, angle=angle2, linestyle="dotted")
ell_bot_L = Ellipse(xy=[-pos-2.4, -pos-1], width=width, height=height, angle=angle1+10, linestyle="dashed")
ell_bot_R = Ellipse(xy=[pos+2.4,  -pos-1], width=width, height=height, angle=angle2-10, linestyle=(0, (1,10)))

ellipse_list = [ell_top_L,   ell_top_R,
                ell_bot_L,   ell_bot_R]
color_list   = [[255, 0, 0], [0, 255, 0],
                [0, 0, 255], [120, 120, 0]]

fig = plt.figure(0)
ax = fig.add_subplot(111, aspect='equal')
ax.axis('off')

ax.text(-pos-6, pos+6, name_top_L)
ax.text(pos+3, pos+6, name_top_R)
ax.text(-pos-14, -pos+3, name_bot_L)
ax.text(pos+10, -pos+3, name_bot_R)

# centers
#ax.text(-pos, pos, "x")
#ax.text(pos, pos, "x")
#ax.text(-pos-2.4, -pos-1, "x")
#ax.text(pos+2.4, -pos-1, "x")

# 13 regions
ax.text(-4.25, 4, "A", horizontalalignment="center")
ax.text(4.25, 4, "B", horizontalalignment="center")
ax.text(-7.5, -1, "C", horizontalalignment="center")
ax.text(7.5, -1, "D", horizontalalignment="center")

ax.text(0, 2.5, "AB", horizontalalignment="center")
ax.text(-5.5, 1.5, "AC", horizontalalignment="center")
ax.text(5.5, 1.5, "BD", horizontalalignment="center")
ax.text(0, -5.75, "CD", horizontalalignment="center")
ax.text(-4.5, -2.6, "BC", horizontalalignment="center")
ax.text(4.5, -2.6, "AD", horizontalalignment="center")

ax.text(-3, 0, "ABC", horizontalalignment="center")
ax.text(3, 0, "ABD", horizontalalignment="center")
ax.text(-2.5, -3.9, "BCD", horizontalalignment="center")
ax.text(2.5, -3.9, "ACD", horizontalalignment="center")

ax.text(0, -2, "ABCD", horizontalalignment="center")

ax.text(-9, -10, "ex A", horizontalalignment="center")
ax.text(-9, -11, "val", horizontalalignment="center")
ax.text(-3, -10, "ex B", horizontalalignment="center")
ax.text(-3, -11, "val", horizontalalignment="center")
ax.text(3, -10, "ex C", horizontalalignment="center")
ax.text(3, -11, "val", horizontalalignment="center")
ax.text(9, -10, "ex D", horizontalalignment="center")
ax.text(9, -11, "val", horizontalalignment="center")


for ellipse, color in zip(ellipse_list, color_list):
  ax.add_artist(ellipse)
  ellipse.set_alpha(0.3)
  ellipse.set_clip_box(ax.bbox)
  ellipse.set_edgecolor("black")
  ellipse.set_facecolor( [rgb_val/255 for rgb_val in color] )

square_lim = 12
ax.set_xlim(-square_lim, square_lim)
ax.set_ylim(-square_lim, square_lim)

plt.show()





