#Internal Versus External Neighbourhood quantification- IVEN
#January 2020
#Jessica E. Forsyth

########################################################################################################################

#first load modules into python, other modules needed are included in the funcs.py file
import numpy as np				 #standard numpy array library for array manipulation
from tkinter import filedialog	 #library to allow interactive file selection
from tkinter import Tk
import pandas as pd				 #library for data handling
import os						 #operating system library to allow filename manipulation
import funcs                     #IVEN specific function library
import matplotlib
#Here we quickly check that the users are using the correct version of matplotlib. 
#In later releases of matplotlib, the way the interactive plot extracts data appears to be different
vers=matplotlib.__version__
if vers!='3.1.2':
	print('Current version of matplotlib librry:')
	print(vers)
	print('Please change your matplotlib version to 3.1.2 using pip install matplotlib==3.1.2 ')
	exit()

root = Tk()
root.filename = filedialog.askopenfilename()   #open dialogue box to allow use to locate and open data file
print(root.filename)
fname = root.filename
root.destroy()

data = pd.read_excel(fname, header=[0,1])    #read in data as array, my data table has row 1 and row 2 as headings/titles, so I wrote 0,1 , if you only have one title row, just write 0
data = np.asarray(data)
dshape = np.shape(data)

#unpack data into variables
num_cells = dshape[0]  #infer number of cells from size of data table
cell_id = data[:, 0]		 #cell_ID output from segmentation or allocated during data import
x = data[:, 1]     	 #x coordinate of cell position
y = data[:, 2]  		 #y coordinate of cell position
z = data[:, 3]			 #z coordinate of cell position
ch1 = data[:, 4]        #channel 1 intensity
ch2 = data[:, 5]		 #channel 2 intensity
ch3 = data[:, 6]		 #channel 3 intensity
ch4 = data[:, 7]		 #channel 4 intensity
ch2_a = np.array(data[:, 8], dtype=np.float32)    #adjusted channel 2 intensity (corrected for depth attenuation)
ch3_a = np.array(data[:, 9], dtype=np.float32)    #adjusted channel 3 intensity (corrected for depth attenuation)
ch4_a = np.array(data[:, 10], dtype=np.float32)   #adjusted channel 4 intensity (corrected for depth attenuation)


points = np.array([x, y, z])   #combine x,y,z coordinates
points = points.T

val1 = np.reshape(ch2_a/np.max(ch2_a), [num_cells, 1])  #normalise channel 2 data
val2 = np.reshape(ch3_a/np.max(ch3_a), [num_cells, 1])  #normalise channel 3 data
val3 = np.reshape(ch4_a/np.max(ch4_a), [num_cells, 1])  #normalise channel 4 data


#now the data is unpacked, we can start running the IVEN pipeline
[outside] = funcs.inside_outside(points, num_cells)								#first calculate which cells are inside and outside using the convex hull
[outside] = funcs.manual_correction(points, num_cells, outside, val1, val2, val3)		#load and display the manual correction window
[nbrs] = funcs.num_nbrs(points, num_cells)									#calculate the number of neighbours of each cell using the DT
#threshold = funcs.threshold_f(points, num_cells, nbrs,outside)					#calculate threshold based on distances between nbrs
threshold = funcs.threshold_f_cellspec(points, num_cells, nbrs,outside)					#calculate threshold based on distances between nbrs inside versus outside specific
[nbrs2] = funcs.nbr_check(points, num_cells, nbrs, threshold,outside)     			#threshold the neighbourhood counts based on distance (correction for the cavity)
neighbour_count = np.sum(nbrs2, axis=1)										#count the total number of neighbours for each cell along row (this is important if using cell spec thresholding) 
                                                                            #if you don't count along the rows, then you count the number of neighbours using the WRONG threshold- i..e case where you have TE and ICM cell, maybe come back to ? 

#now calculate the number of neighbours of cells which are outside (this was used as an extension to the outside analysis (to show mural and polar TE cells)
nbr_comp = np.zeros([num_cells, 1])   #initate array to store data
for c1 in range(num_cells):
    for c2 in range(num_cells):
        if nbrs2[c1, c2] == 1:        #check if cellj is a neighbour of celli
            if outside[c2, 0] == 1:        #if it is a neighbour, check if the cell is TE or INS
                nbr_comp[c1, 0] += 1  #if it is TE, add one to the count


#now prepare data for saving
fname_only = os.path.basename(fname)
dire = os.path.dirname(fname)
fname_out = dire + '/nbrs_' + fname_only  #file name that includes the directory and filename of the excel document

#pepare data for excel table ({'Name':variable,.....}) if you add more to the table, this is where you have to change the output
data_out = pd.DataFrame({'Cell ID': cell_id, 'x': x, 'y': y, 'z': z, 'CH1': ch1, 'CH2': ch2, 'CH3': ch3, 'CH4': ch4, 'CH2_adj': ch2_a, 'CH3_adj': ch3_a, 'CH4_adj': ch4_a, 'outside': outside[:, 0], '# nbrs': neighbour_count, '#nbrs that are outside': nbr_comp[:, 0]})
data_out = data_out[['Cell ID', 'x', 'y', 'z', 'CH1', 'CH2', 'CH3', 'CH4', 'CH2_adj', 'CH3_adj', 'CH4_adj', 'outside', '# nbrs', '#nbrs that are outside']]
data_out = data_out.sort_values(by='outside')  #sort data by outside vs INS, choose whatever you like!
data_out.to_excel(fname_out, sheet_name='test', index=False)

funcs.show_fig_final(points, num_cells, outside, cell_id)  #now display the final figure with final allocations and CELL_ID labels for reference

funcs.save_data(points, num_cells, outside, cell_id,dire,fname_only)

print('Finished!')
