# Internal Versus External Neighbourhood quantification v2.0
# August 2023
# Jessica E. Forsyth

import numpy as np
from tkinter import filedialog	 #library to allow interactive file selection
#from tkinter import Tk
import pandas as pd				 #library for data handling
import os						 #operating system library to allow filename manipulation
import funcs_IVEN2 as f2                #IVEN specific function library



# Create open file dialog
filepath = filedialog.askopenfilename()
#filepath = "C:/Users/Admin/Documents/IVEN-2/em1_7June19E3.5PDxCD110amEcadSOX2CDX2.xls"
print('\nIVEN analysis started -')
print('Analysing file : ', filepath)

# create structure
class data():  # structure to store all imported data
    pass
class results():  # structure to store IVEN results/outputs
    pass
class plotting():  # structure to store IVEN results/outputs
    pass

# load in data from excel file
data.all = pd.read_excel(filepath, header=[0])    #read in data as array, my data table has row 1 as headings/titles,
data.headings = data.all.keys()
data.id = data.all[data.headings[0]]
data.xyz = np.asarray(data.all[data.headings[1:4]])
data.properties = np.asarray(data.all[data.headings[4:-1]])
data.num_cells = len(data.id)

results.all = pd.DataFrame.copy(data.all)
results.dist_all = pd.DataFrame()


# find convex hull of points and define outside cells vers1
[results.outside_bool1, results.outside_ids1] = f2. classify_outside(data.xyz, data.num_cells)
# interactive plot to select and correct outside versus inside cells
results = f2.manual_correct_outside(data, results, plotting)

# ask what thresholding method you want to use? (none, value, automatic (cell type specific), automatic (all))
# ask if you would like to select cavity points within embryo (later implement automatic selection of these)
results = f2.analysis_prefs(data, results, plotting)

# generate raw nbr matrix
results.nbr_matrix1 = f2.nbr_matrix(data)
# calculate threshold
results = f2.eval_threshold(data, results)
print('...Neighbour threshold method- %s \n...Threshold value(s) -' % results.thresh_method)
print('-----',results.thresh_val)

# neighbourhood evaluation (count number of nbrs)
results = f2.check_nbrs(data, results)
print('...Removed %3d neighbours.' % int((np.sum(results.nbr_matrix1) - np.sum(results.nbr_matrix2))/2))

# generate new DT including cavity points
# any cells with connection to cavity points are highlighted as potentially cavity adjacent cells
if len(results.cavity_pts)>0:
    results = f2.dt_withcav(data, results)

# count number of nbrs which are outside/inside
# output cell neighbour IDs
# number of nbrs
# number of nbrs that are inside, outside, cavity adjacent

print('...Compiling data to .xls/.csv files.')
results = f2.compile_results(data, results)
results = f2.compile_distances(data, results)

fname_only = os.path.basename(filepath)
directory = os.path.dirname(filepath)
output_fname = directory + '/IVEN2out_' + fname_only[0:-4] + '.csv'
output_fname_dist = directory + '/IVEN2dists_' + fname_only[0:-4] + '.csv'

results.all.to_csv(output_fname, index=False)
results.dist_all.to_csv(output_fname_dist, index=False)

#plot embryo with cell classifications
f2.final_figure(directory, fname_only, data, results, plotting)

print('Analysis complete. View final figure and cell properties by running the display_final_fig.py file. ')
print('...Results summary saved to -', output_fname)
print('...Distances summary saved to -', output_fname_dist)