import matplotlib
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tkinter import filedialog	 #library to allow interactive file selection
from tkinter import Tk
import os						 #operating system library to allow filename manipulation

#This script allows you to inspect previous data and rotate the figure of interest. The resulting plot shows the finalised cell classification (outside versus inside).
#To add if wanted- channel data to allow inspection of cell protein content- but I am hesitant to do this because the colour schemes may lead to mis-intepretation of the data

#open numpy data file
root = Tk()
root.filename = filedialog.askopenfilename()   #open dialogue box to allow use to locate and open data file
print(root.filename)
fname = root.filename
root.destroy()
#unpack data from numpy file
data = np.load(fname, allow_pickle=True)
cell_id = data[:,0]
points = data[:,1:4]
outside = data[:,4]
s = np.shape(cell_id)
num_cells =s[0]
print(num_cells)


#plot figure as before in finalfig
fig = plt.figure()
ax = fig.gca(projection='3d')
p = ax.scatter(points[:, 0], points[:, 1], points[:, 2], facecolors=['lightgrey'] * num_cells, edgecolors=['white'] * num_cells, picker=5, s=[150] * num_cells, alpha=1)  # 5 points tolerance
#for i in range(num_cells):
 #   ax.text(points[i,0],points[i,1],points[i,2],str(int(cell_id[i])))                     #for each cell, label cell point with corresponding cell ID

hull = ConvexHull(points)
for i in hull.simplices:
    plt.plot(points[i, 0], points[i, 1], points[i, 2], '-', color='darkgrey')  #again plot a rough mesh for visualisation of embryo shape

ec_array=['white'] * num_cells
for ii in range(num_cells):
    if outside[ii]==1:
        ec_array[ii]='black'

p._edgecolor3d= ec_array
fig.canvas.draw()
#ax.axis('equal')
ax.set_title('Final classification')
ax.set_xticklabels(['']*10)
ax.set_yticklabels(['']*10)
ax.set_zticklabels(['']*10)

plt.show()