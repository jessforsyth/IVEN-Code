#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import matplotlib
import numpy as np
import scipy.stats
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math


#This file contains all of the functions needed to run IVEN, make sure this file is in the same directory as the run_IVEN.py file.


def inside_outside(points, num_cells):
    #this function first categorises cells as inside or outside using the standard Convex Hull algorithm
 
    hull = ConvexHull(points)
    outside_points = np.unique(hull.simplices)                #identify all unique outside points
    inside_points = []
    outside = np.ones([num_cells, 1])
    for i in range(0, num_cells):                            #create boolean lookup array for cell identity (outside or inside)
        if i not in outside_points:
            inside_points.append(i)
            outside[i, 0] = 0
    print(np.sum(outside), 'outside cells identified')
    return[outside]

def manual_correction(points,num_cells, outside,val1,val2,val3):
    #this function displays the points along with the Convex Hull mesh and cell identities
    #cells that are outside are outlined in black, and cells outlined in white are inside cells
    #if cell classification need to be changed, simply click on the cell
    #to aid manual checking/re-classification, channel intensities can be imported and displayed using 'g','r' and 'i'
    print(' ', end='\n')
    print('Please manually correct the outside classification.....')
    print('Click on cells to change their classification, black edge= outside, white edge=inside.')
    print('Press "g","r", or "i" to visualise the green, red and far-red channels respectively.')
    print('Press enter when finished.')
    
    #green colour gradient : this creates a normalised colour scheme for the green channel data
    fc = np.zeros([100, 4])
    c1 = np.array([1, 1, 1])
    c2 = np.array([0, 0.5, 0])
    for i in range(0, 100, 1):
        c = c1 + ((i / 100) * (c2 - c1))
        fc[i, 0:3] = c[:]
    order1t = np.round(val1, 2) * 100
    order1 = order1t.astype(int)
    fc1 = fc[order1-1, :]
    fc1 = np.reshape(fc1, [num_cells, 4])

    #red colour gradient : this creates a normalised colour scheme for the red channel data
    fc = np.zeros([100, 4])
    c1 = np.array([1, 1, 1])
    c2 = np.array([1, 0, 0])
    for i in range(0, 100, 1):
        c = c1 + ((i / 100) * (c2 - c1))
        fc[i, 0:3] = c[:]
    order2t = np.round(val2, 2) * 100
    order2 = order2t.astype(int)
    fc2 = fc[order2 - 1, :]
    fc2 = np.reshape(fc2, [num_cells, 4])

    #far-red colour gradient : this creates a normalised colour scheme for the far-red channel data
    fc = np.zeros([100, 4])
    c1 = np.array([0, 0, 0])
    c2 = np.array([1, 1, 1])
    for i in range(0, 100, 1):
        c = c1 + ((i / 100) * (c2 - c1))
        fc[i, 0:3] = c[:]
    order3t = np.round(val3, 2) * 100
    order3 = order3t.astype(int)
    fc3 = fc[order3-1, :]
    fc3 = np.reshape(fc3, [num_cells, 4])

    fig = plt.figure()
    ax = fig.gca(projection='3d')  #coordinates are in three dimensions so set up axes for 3D

    p = ax.scatter(points[:, 0], points[:, 1], points[:, 2], facecolors=['lightgrey'] * num_cells, edgecolors=['white'] * num_cells,picker=5, s=[50] * num_cells, alpha=1)  #plot cells as 3D scatter plot.
    hull = ConvexHull(points)                                                #generate ConvexHull of points to show ConvexHull meshwork
    for i in hull.simplices:
        plt.plot(points[i, 0], points[i, 1], points[i, 2], '-', color='darkgrey')  #plot meshwork

    ec_array = ['white'] * num_cells
    for ii in range(num_cells):
        if outside[ii, 0] == 1:
            ec_array[ii] = 'black'

    p._edgecolor3d = ec_array

    #p._edgecolor3d = (np.ones([num_cells, 4]) * np.reshape((1 - outside[:, 0]), [num_cells, 1]))  #specify edgecolours of markers
    fig.canvas.draw()
    ec = p.get_edgecolors()                                                       #find edgecolours of scatter plot points
    #ax.axis('equal')                                                              #set axes equal (otherwise embryos may look skewed due to misleading axes)
    ax.set_title('Correct the cell classification if necessary')
    ax.set_xticklabels(['']*10)
    ax.set_yticklabels(['']*10)
    ax.set_zticklabels(['']*10)

    def onpick(event):                                                            # this function is the response for when cells are clicked on the plot
        ich = event.ind[0]
        if outside[ich, 0] == 1:
            ec[ich, :] = (1, 1, 1, 1)                                                   #if the selected cell is already a outside cell, reassign as an inside cell
            outside[ich, 0] = 0                                                            #change cell classification and marker edge colour
            print('Cell assigned to inside group ', end='\r')
        else:
            ec[ich, :] = (0, 0, 0, 1)                                                   #if the selected cell is already a inside cell, reassign as an TE cell
            outside[ich, 0] = 1                                                            #change cell classification and marker edge colour
            print('Cell assigned to outside group', end='\r')

        p._edgecolor3d = ec                                                       #update edgecolour list
        fig.canvas.draw()
        return outside

    def press(event):                                                        #this function responds to the pressing of keyboard keys- used to display different channel or close the figure
        if event.key == 'enter':                                             #if you press 'enter' this will close the figure
            plt.close(fig)
            plt.ioff()  # turn interactive off
            print(' ', end='\n')
            print('Manual correction complete.')

        if event.key == 'g':   #green channel                               #press 'g' on the keyboard to display the green-channel intensities by changing the facecolour of the cell markers
            print('Green channel displayed [low(white)--->high(green)]')
            p._facecolor3d = fc1
            fig.canvas.draw()

        if event.key == 'r':   #red channel                                 #press 'r' on the keyboard to display the red-channel intensities by changing the facecolour of the cell markers
            print('Red channel displayed [low(white)--->high(red)]')
            p._facecolor3d = fc2
            fig.canvas.draw()

        if event.key == 'i':   #far red channel, cant use f
            print('Far-red channel displayed [low(black)--->high(white)]')  #press 'i' on the keyboard to display the farred-channel intensities by changing the facecolour of the cell markers
            p._facecolor3d = fc3
            fig.canvas.draw()
        return()

    fig.canvas.mpl_connect('pick_event', onpick)
    fig.canvas.mpl_connect('key_press_event', press)

    plt.show()    
    return[outside]

def num_nbrs(points, num_cells):
    #this function now calculates the original number of neighbours for each cell (without thresholding to account for potential cavities)
    #uses the DT and infers from the tetrahedrons which cells are neighbours to each other.

    tri = Delaunay(points)                                        #generate the Delaunay Triangulation
    len(tri.simplices)

    nbrs = np.zeros((num_cells, num_cells))                        #initate array to store neighbourhood matrix

    for i in range(0, num_cells):
        for j in range(0, len(tri.simplices)):
            if i in tri.simplices[j, :]:
                nbrs[i, tri.simplices[j, ]] = 1                   #find corresponding term in neighbourhood matrix and change to 1, to indicate neighbour
                nbrs[i, i] = 0                                    #but make sure a cell is not counted as a neighbour of itself

    return[nbrs]

def threshold_f(points, num_cells, nbrs,outside):
    #this function calculates the desired threshold to account for the cavity
    #this is based on cell diameters and an imposed threshold.
    #due to the cells reducing in size over the preimplantation period, we select the threshold based on cell number (i.e. embryo stage)
    
    D=[]
    for i in range(0, num_cells):
        for j in range(0, num_cells):
            if nbrs[i, j] == 1:                       #if cells are neighbours then compare distances
                sqdx = (points[j, 0]-points[i, 0]) ** 2
                sqdy = (points[j, 1]-points[i, 1]) ** 2
                sqdz = (points[j, 2]-points[i, 2]) ** 2
                d = math.sqrt(sqdx+sqdy+sqdz)        #calculating the distance between cells
                D.append(d)

    p75=np.percentile(D,75)
    k=0.5
    iqr=scipy.stats.iqr(D)	

    thresholdval = p75 + (k* iqr)
    print('Distance threshold: ', thresholdval)

    threshold=[thresholdval,thresholdval]
 
    return threshold

def threshold_f_cellspec(points, num_cells, nbrs,outside):
    #this function calculates the desired threshold to account for the cavity
    #this calculates a threshold for inside cells and then for outside cells 
    #this is based on cell diameters and an imposed threshold.
    #due to the cells reducing in size over the preimplantation period, we select the threshold based on cell number (i.e. embryo stage)
    
    Doutside=[]
    Dinside=[]
    for i in range(0, num_cells):
        for j in range(0, num_cells):
            if nbrs[i, j] == 1:                       #if cells are neighbours then compare distances
                sqdx = (points[j, 0]-points[i, 0]) ** 2
                sqdy = (points[j, 1]-points[i, 1]) ** 2
                sqdz = (points[j, 2]-points[i, 2]) ** 2
                d = math.sqrt(sqdx+sqdy+sqdz)        #calculating the distance between cells
                if outside[i,0]==1:
                    Doutside.append(d)
                else:
                    Dinside.append(d)

    p75outside=np.percentile(Doutside,75)
    k=0.5
    iqroutside=scipy.stats.iqr(Doutside)	
    thresholdoutside = p75outside + (k* iqroutside)

    p75inside=np.percentile(Dinside,75)
    iqrinside=scipy.stats.iqr(Dinside)	
    thresholdinside = p75inside + (k* iqrinside)

    print('Distance threshold outside: ', thresholdoutside)
    print('Distance threshold inside ', thresholdinside)

    threshold=[thresholdoutside,thresholdinside]
    
    return threshold

def nbr_check(points, num_cells, nbrs, threshold,outside):
    #we now check our orignal neighbourhood numbers by using the threshold caluclated above on line 187.
    #we calculate the edistance between every cell and every neighbour it has, if this distance is greater than the threshold, the cell is discounted as a neighbour
    #this stage could be ignored if not cavities are present

    nbrs2 = np.zeros([num_cells, num_cells])
    for i in range(0, num_cells):
        for j in range(0, num_cells):
            if nbrs[i, j] == 1:                       #if cells are neighbours then compare distances
                sqdx = (points[j, 0]-points[i, 0]) ** 2
                sqdy = (points[j, 1]-points[i, 1]) ** 2
                sqdz = (points[j, 2]-points[i, 2]) ** 2
                d = math.sqrt(sqdx+sqdy+sqdz)        #calculating the distance between cells
                if outside[i,0]==1 and outside[j,0]==1:
                    if d > threshold[0]: #i.e. untrue neighbours use outside threshold
                        #print(threshold[0])
                        #print(i)
                        #print(j)
                        nbrs2[i, j] = 0
                    else:                              #i.e. true neighbours
                        nbrs2[i, j] = 1
                else:
                    if d> threshold[1]:   #use inside cell threshold
                        #print(threshold[1])
                        #print(i)
                        #print(j)
                        nbrs2[i, j] = 0
                    else:                              #i.e. true neighbours
                        nbrs2[i, j] = 1               

    return[nbrs2]                                  #return true/corrected neighbour matrix

def show_fig_final(points, num_cells, outside, cell_id):
    #this function displays the final classificaion of inside and outside cells, plus labels each point with the cell id
    #unfortunately matplotlib does not yet have a function to export the interactive figure, so any final copies of this plot should be caputured as png files at this stage
    #use user GUI to save in the appropriate rotations/fields of view

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    p = ax.scatter(points[:, 0], points[:, 1], points[:, 2], facecolors=['lightgrey'] * num_cells, edgecolors=['white'] * num_cells, picker=5, s=[50] * num_cells, alpha=1)  # 5 points tolerance
    #for i in range(num_cells):
     #   ax.text(points[i,0],points[i,1],points[i,2],str(int(cell_id[i])))                     #for each cell, label cell point with corresponding cell ID

    hull = ConvexHull(points)
    for i in hull.simplices:
        plt.plot(points[i, 0], points[i, 1], points[i, 2], '-', color='darkgrey')  #again plot the mesh of the convex hull for visualisation of embryo shape

    ec_array=['white'] * num_cells
    for ii in range(num_cells):
        if outside[ii,0]==1:
            ec_array[ii]='black'

    p._edgecolor3d= ec_array
    fig.canvas.draw()
    #ax.axis('equal')
    ax.set_title('Final classification')
    ax.set_xticklabels(['']*10)
    ax.set_yticklabels(['']*10)
    ax.set_zticklabels(['']*10)

    def press(event):
        if event.key == 'enter':                 #to close figure press enter as before
            plt.close(fig)
            plt.ioff()                           # turn interactive off
            print(' ', end='\n')
            print('Figure closed.')
        return()

    print('Press enter to close figure and save data. \n ')
    fig.canvas.mpl_connect('key_press_event', press)
    plt.show()
    return

def save_data(points, num_cells, outside, cell_id, dire, fname_only):
    data = np.concatenate((np.reshape(cell_id, [num_cells, 1]), np.reshape(points[:, 0], [num_cells, 1]), np.reshape(points[:, 1], [num_cells, 1]), np.reshape(points[:, 2], [num_cells, 1]), np.reshape(outside, [num_cells, 1])),axis=1)
    fname= dire + '/figdata_' + fname_only
    np.save(fname, data,allow_pickle=True)
    return

