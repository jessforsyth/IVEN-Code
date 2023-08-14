import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import ConvexHull
from tkinter import filedialog	 #library to allow interactive file selection
matplotlib.use('TkAgg')
from matplotlib.widgets import Button
from matplotlib.widgets import RadioButtons

# change plotting parameters
s = 150                             # size of markers in scatter plot
cmap = matplotlib.cm.spring         # colour map of points, see others at https://matplotlib.org/stable/tutorials/colors/colormaps.html
cell_ID_labels = 'on'               # change to 'off' if you want to remove the cell ID labels on the plot


########################################################################################################################
# Create open file dialog
print('\nDisplay final IVEN figure -')
print('Select the IVEN2out_.... file you would like to display: ')
filepath = filedialog.askopenfilename()

# create structures
class data():  # structure to store all imported data
    pass
class plotting():  # structure to store IVEN results/outputs
    pass

# load in data from excel file
data.all = pd.read_csv(filepath, header=[0])    #read in data as array, my data table has row 1 as headings/titles,
data.all = data.all.drop(columns=['nbr_ids'])   # we can't display this so ignore it
data.headings = data.all.keys()
data.xyz = np.asarray(data.all[data.headings[1:4]])
data.num_cells = len(data.all['ID'])

plotting.edge_col = np.multiply(np.ones([data.num_cells, 4]), [1,1,1,1])
if 'outside_bool' in data.all:
    for cell in range(data.num_cells):
        if data.all.loc[cell, 'outside_bool'] == 1:
            plotting.edge_col[cell, :] = [0,0,0,0]

plotting.sizes = np.ones(data.num_cells) * s

position = [-0.22, -0.12, 1.1, 1.1]


###### Here is the code for plotting the figure ##############
matplotlib.rcParams['font.family'] = 'Arial'

fig = plt.figure('IVEN - Manual Correction Stage', figsize=(6, 5))
ax = fig.add_subplot(111) #, projection='3d')
fig.subplots_adjust(right=1, left=0, bottom=0, top=1, wspace=0, hspace=0)



em_ax = fig.add_axes(position, projection='3d')
sct_plt = em_ax.scatter(data.xyz[:, 0], data.xyz[:, 1], data.xyz[:, 2],
                        edgecolors=plotting.edge_col, picker=True, s=plotting.sizes, linewidth=2, alpha=1)
em_ax.set_xlim(np.min(data.xyz[:, 0]), np.max(data.xyz[:, 0]))
em_ax.set_ylim(np.min(data.xyz[:, 1]), np.max(data.xyz[:, 1]))
em_ax.set_zlim(np.min(data.xyz[:, 2]), np.max(data.xyz[:, 2]))
em_ax.set_box_aspect([1.0, 1.0, 1.0])
em_ax.patch.set_alpha(0.0)

hull = ConvexHull(data.xyz)  # generate ConvexHull of points to show ConvexHull meshwork
for i in hull.simplices:
    i = np.append(i, i[0])
    plt.plot(data.xyz[i, 0], data.xyz[i, 1], data.xyz[i, 2], '-', color='darkgrey',
             linewidth=1.5)  # plot meshwork

if cell_ID_labels=='on':
    for cell in range(data.num_cells):
        em_ax.text(data.all.loc[cell, data.headings[1]], data.all.loc[cell, data.headings[2]],
                   data.all.loc[cell, data.headings[3]], str(data.all.loc[cell, 'ID']), fontsize=10)

em_ax.grid(False)
plt.axis('off')

add_done_but = fig.add_axes([0.88, 0.01, 0.1, 0.075])
done_but = Button(add_done_but, 'Finish')

add_rad = fig.add_axes([0.64, 0.1, 0.32, 0.03 * (len(data.headings)-4)])
rad = RadioButtons(
    ax=add_rad,
    labels=data.headings[4:],
    radio_props={'facecolor':['magenta'], 's':[40]})

sm = matplotlib.cm.ScalarMappable(cmap=cmap)
cb_ax = fig.add_axes([0.48, 0.85, 0.47, 0.04])
cbar = fig.colorbar(sm, cax=cb_ax, ticks=[0, 1], orientation='horizontal')# , label='Property value')
plt.annotate('Property value', xy=(0.656,0.902), xycoords='figure fraction')

def radio_click(prop_event):
    min_prop_val = np.min(data.all[prop_event])
    max_prop_val = np.max(data.all[prop_event])
    trans_data = np.divide(data.all[prop_event] - min_prop_val, max_prop_val - min_prop_val)
    plotting.face_col = cmap(trans_data)
    sct_plt.set_facecolors(plotting.face_col)
    cbar.ax.set_xticklabels([str(round(min_prop_val, 4)), str(round(max_prop_val, 4))])
    fig.canvas.draw()

def done_click(event):
    plt.close()

radio_click(data.headings[4])
done_but.on_clicked(done_click)
rad.on_clicked(radio_click)

fig.canvas.draw()
plt.show()