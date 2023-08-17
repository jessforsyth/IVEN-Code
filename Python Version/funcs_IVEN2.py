import numpy as np
import matplotlib
import scipy
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
from matplotlib.widgets import Button
from matplotlib.widgets import RadioButtons
from matplotlib.widgets import TextBox

def classify_outside(pts, n):
    hull = ConvexHull(pts)
    outside_ids = np.unique(hull.simplices)
    outside_bool = np.zeros(n)
    outside_bool[outside_ids] = 1
    return [outside_bool, outside_ids]

def plot_embryo(fig, data, plotting, position):
    em_ax = fig.add_axes(position, projection='3d')
    sct_plt = em_ax.scatter(data.xyz[:, 0], data.xyz[:, 1], data.xyz[:, 2], facecolors=plotting.face_col,
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

    em_ax.grid(False)
    plt.axis('off')

    return em_ax, sct_plt

def manual_correct_outside(data, results, plotting):
    results.outside_bool2 = np.copy(results.outside_bool1)

    matplotlib.rcParams['font.family'] = 'Arial'
    fig = plt.figure('IVEN - Manual Correction Stage', figsize=(8, 5))
    ax = fig.add_subplot(111) #, projection='3d')
    fig.subplots_adjust(right=1, left=0, bottom=0, top=1, wspace=0, hspace=0)


    ax.annotate('Check and correct outside cell classification -', xy=(0.02, 0.95), xycoords='figure fraction', fontsize=14)
    ax.annotate('Select cells to toggle cell classification:\n - Black edge = outside\n - White edge= inside',
                xy=(0.02, 0.84), xycoords='figure fraction', fontsize=11)
    ax.annotate('Use buttons to visualise \nimported cell properties.',
                xy=(0.56, 0.1), xycoords='figure fraction', fontsize=11)

    plotting.face_col = np.ones([data.num_cells, 4])
    plotting.edge_col = np.ones([data.num_cells, 4])
    plotting.edge_col[results.outside_ids1, 0:3] = [0, 0, 0]
    s = 150
    plotting.sizes = s * np.ones(data.num_cells)

    embryo_position = [-0.22, -0.12, 0.95, 1.1]
    em_ax, sct_plt = plot_embryo(fig, data, plotting, embryo_position)

    add_done_but = fig.add_axes([0.88, 0.01, 0.1, 0.075])
    done_but = Button(add_done_but, 'Finish')

    add_rad = fig.add_axes([0.78, 0.1, 0.2, 0.03* (len(data.headings)-4)])
    rad = RadioButtons(
        ax=add_rad,
        labels=data.headings[4:],
        radio_props={'facecolor':['magenta'], 's':[40]})


    def onpick(event):
        sel_id = event.ind[0]
        plotting.edge_col[sel_id, 0:3] = 1 - plotting.edge_col[sel_id, 0:3]
        sct_plt.set_edgecolors(plotting.edge_col)
        #update outside indices
        results.outside_bool2[sel_id] = 1 - results.outside_bool2[sel_id]
        fig.canvas.draw()

    cmap = matplotlib.cm.spring
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
    fig.canvas.mpl_connect('pick_event', onpick)
    done_but.on_clicked(done_click)
    rad.on_clicked(radio_click)

    fig.canvas.draw()
    plt.show()

    # update list of outside cell ids
    results.outside_ids2 = np.where(results.outside_bool2 == 1)[0]
    results.inside_ids2 = np.where(results.outside_bool2 == 0)[0]

    results.all['outside_bool'] = results.outside_bool2
    diff = sum(abs(results.outside_bool1 - results.outside_bool2))
    print('...Number of cell allocations manually changed : %3d out of %3d' % (int(diff), int(data.num_cells)))

    return results

def analysis_prefs(data, results, plotting):
    results.attempted_cav_num =0

    matplotlib.rcParams['font.family'] = 'Arial'
    results.cavity_pts =[]
    results.guess_cavity_adj = np.zeros(data.num_cells)
    plotting.all_col = [255 / 256, 102 / 256, 204 / 256, 1]
    plotting.cav_col = [255/256, 204/256, 0, 1]

    # generate canvas
    fig = plt.figure('IVEN - Analysis Preferences', figsize=(8, 5))
    ax = fig.add_subplot(111)# , projection='3d')
    fig.subplots_adjust(right=1, left=0, bottom=0, top=1, wspace=0, hspace=0)

    # add buttons and options
    ax.annotate('Select cavity points -', xy=(0.02, 0.95), xycoords='figure fraction', fontsize=14)
    ax.annotate('Select cells around a single cavity. Press the \nnext button to store the cavity \npoint.',
                xy=(0.02, 0.84), xycoords='figure fraction', fontsize=11)

    ax.annotate('Select thresholding method -', xy=(0.5, 0.95), xycoords='figure fraction', fontsize=14)
    add_rad2 = fig.add_axes([0.5, 0.7, 0.45, 0.2])
    rad2 = RadioButtons(
        ax=add_rad2,
        active=0,
        labels=['Automatic (cell position dependent)', 'Automatic (cell position independent)', 'Manual', 'None'],
        label_props={'fontsize':[11]},
        radio_props={'facecolor': ['orange'], 's': [40]})

    add_done_but = fig.add_axes([0.88, 0.01, 0.1, 0.075])
    done_but = Button(add_done_but, 'Finish')
    axthresh_box = fig.add_axes([0.64, 0.53, 0.1, 0.045])
    thresh_box = TextBox(axthresh_box, label='Parameters', textalignment='center')


    embryo_position = [-0.07, -0.05, 0.6, 1.1]
    plotting.face_col = np.multiply(np.ones([data.num_cells, 4]), plotting.all_col)
    em_ax, sct_plt = plot_embryo(fig, data, plotting, embryo_position)

    add_next_but = fig.add_axes([0.25, 0.01, 0.16, 0.075])
    next_but = Button(add_next_but, 'Done/Next cavity')
    add_remove_but = fig.add_axes([0.05, 0.01, 0.16, 0.075])
    remove_but = Button(add_remove_but, 'Remove last cavity')

    def remove_cav_pt(event):
        if len(results.cavity_pts) > 0:
            print('...Cavity %2d centre coordinates removed' % int(len(results.cavity_pts)))

            em_ax.scatter(results.cavity_pts[int(len(results.cavity_pts)) - 1][0],
                          results.cavity_pts[int(len(results.cavity_pts)) - 1][1],
                          results.cavity_pts[int(len(results.cavity_pts)) - 1][2],
                          facecolors='r', s=plotting.sizes[0] * 1.25, alpha=1, linewidth=2.5, marker='x')

            # output cavity points on screen
            txt1 = ('-----------------------------------------------------------------------')
            ypos = 0.44 - (0.05 * (len(results.cavity_pts) - 1))
            ax.annotate(txt1, xy=(0.5, ypos), xycoords='figure fraction', fontsize=11, color='r')
            fig.canvas.draw()

            #delete cavity point
            results.cavity_pts.pop(len(results.cavity_pts)-1)


    def calc_cav_pt(event):
        results.attempted_cav_num = results.attempted_cav_num + 1
        cav_bool = results.guess_cavity_adj.astype(bool)
        results.cavity_pts.append(np.mean(data.xyz[cav_bool, :], axis=0))
        plotting.face_col = np.multiply(np.ones([data.num_cells, 4]), plotting.all_col)
        sct_plt.set_facecolors(plotting.face_col)
        results.guess_cavity_adj = np.zeros(data.num_cells)

        print('...Cavity %2d centre coordinates stored' % int(len(results.cavity_pts)))
        print('-----', results.cavity_pts[int(len(results.cavity_pts)) - 1])

        em_ax.scatter(results.cavity_pts[int(len(results.cavity_pts)) - 1][0],
                               results.cavity_pts[int(len(results.cavity_pts)) - 1][1],
                               results.cavity_pts[int(len(results.cavity_pts)) - 1][2],
                               facecolors='k', s=plotting.sizes[0]*1.25, alpha=1, linewidth=2.5, marker='+')

        # output cavity points on screen
        txt1 = ('Cavity %2d - [%.3f, %.3f, %.3f]' % (int(len(results.cavity_pts)),
                results.cavity_pts[int(len(results.cavity_pts)) - 1][0],
                results.cavity_pts[int(len(results.cavity_pts)) - 1][1],
                results.cavity_pts[int(len(results.cavity_pts)) - 1][2]))
        ypos = 0.44 - (0.05 * (results.attempted_cav_num - 1))
        ax.annotate(txt1, xy=(0.5, ypos), xycoords='figure fraction', fontsize=11)
        fig.canvas.draw()

    def done_click(event):
        plt.close()

    def onpick(event):
        sel_id = event.ind[0]
        if list(plotting.face_col[sel_id, 0:4]) == plotting.all_col:
            plotting.face_col[sel_id, 0:4] = plotting.cav_col
        else:
            plotting.face_col[sel_id, 0:4] = plotting.all_col
        sct_plt.set_facecolors(plotting.face_col)
        #update guesses of cavity adjacent cells
        results.guess_cavity_adj[sel_id] = 1- results.guess_cavity_adj[sel_id]

        # plot inferred cavity point?

        fig.canvas.draw()

    def manual_input(val):
        results.thresh_val = val

    def threshold_sel(event):
        results.thresh_method = event
        if event == 'Automatic (cell position dependent)':
            plotting.thresh_annotation.set_text('Threshold automatically calculated for inside/\noutside cells. Enter k-value and press enter.')
            thresh_box.set_val('0.5')
            thresh_box.label.set_text('k-value   ')
        if event == 'Automatic (cell position independent)':
            plotting.thresh_annotation.set_text('Threshold automatically calculated for all cells.\nEnter k-value and press enter.')
            thresh_box.set_val('0.5')
            thresh_box.label.set_text('k-value   ')
        if event == 'None':
            plotting.thresh_annotation.set_text('No threshold set for neighbourhood evaluation.\n')
            thresh_box.set_val('-')
            thresh_box.label.set_text(' ')
        if event == 'Manual':
            plotting.thresh_annotation.set_text('Enter a threshold value in the box below and \npress enter.')
            thresh_box.label.set_text('Threshold    ')
            thresh_box.set_val('10')


    next_but.on_clicked(calc_cav_pt)
    remove_but.on_clicked(remove_cav_pt)
    thresh_box.on_submit(manual_input)
    fig.canvas.mpl_connect('pick_event', onpick)
    done_but.on_clicked(done_click)
    rad2.on_clicked(threshold_sel)

    plotting.thresh_annotation = plt.annotate('', xy=(0.5, 0.61), xycoords='figure fraction', fontsize=11)
    threshold_sel('Automatic (cell position dependent)')

    plt.show()

    return results

def nbr_matrix(data):
    # generate DT
    dt = Delaunay(data.xyz)

    nbr_matrix = np.zeros((data.num_cells, data.num_cells))
    for cell in range(data.num_cells):
        for sim in range(0, len(dt.simplices)):
            if cell in dt.simplices[sim,:]:
                nbr_matrix[cell, dt.simplices[sim,]] = 1
        nbr_matrix[cell, cell] = 0 # stop double counting

    return nbr_matrix

def eval_threshold(data, results):

    results.dist_matrix1 = np.zeros((data.num_cells, data.num_cells))
    results.dist_matrix1[:] = 'nan'

    for cell1 in range(data.num_cells):
        for cell2 in range(data.num_cells):
            if results.nbr_matrix1[cell1, cell2] == 1:
                # calc distance
                results.dist_matrix1[cell1, cell2] = np.linalg.norm(data.xyz[cell1,:] - data.xyz[cell2,:])

    #check if we have any outside/inside cells classified
    if np.sum(results.outside_bool2) == 0 and results.thresh_method=='Automatic (cell position dependent)':
        print('!!! No outside cells classified - thresholding method changed to "Automatic (cell position independent)" ')
        results.thresh_method = 'Automatic (cell position independent)'
    if np.sum(results.outside_bool2) == data.num_cells and results.thresh_method=='Automatic (cell position dependent)':
        print('!!! No inside cells classified - thresholding method changed to "Automatic (cell position independent)" ')
        results.thresh_method = 'Automatic (cell position independent)'


    if results.thresh_method=='None':
        results.thresh_val = np.nanmax(results.dist_matrix1) * 100  #big number so this doesn't actually do anything
        results.all['threshold'] = ['-'] * data.num_cells

    if results.thresh_method=='Manual':
        results.thresh_val = float(results.thresh_val)
        results.all['threshold'] = np.ones(data.num_cells) * results.thresh_val

    if results.thresh_method=='Automatic (cell position dependent)':
        dists_out = results.dist_matrix1[results.outside_ids2, : ]
        dists_in = results.dist_matrix1[results.inside_ids2, : ]

        k = float(results.thresh_val)
        p75_out = np.nanpercentile(dists_out, 75)
        iqr_out = scipy.stats.iqr(dists_out, nan_policy='omit')
        p75_in = np.nanpercentile(dists_in, 75)
        iqr_in = scipy.stats.iqr(dists_in, nan_policy='omit')

        results.thresh_val = [p75_out + (k * iqr_out), p75_in + (k * iqr_in)]

        results.all['threshold'] = np.ones(data.num_cells) * results.thresh_val[1]
        for cell in range(data.num_cells):
            if results.outside_bool2[cell] == 1:
                results.all.loc[cell, 'threshold'] = results.thresh_val[0]


    if results.thresh_method=='Automatic (cell position independent)':
        k = float(results.thresh_val)
        p75 = np.nanpercentile(results.dist_matrix1, 75)
        iqr = scipy.stats.iqr(results.dist_matrix1, nan_policy='omit')

        results.thresh_val = p75 + (k * iqr)
        results.all['threshold'] = np.ones(data.num_cells) * results.thresh_val

    return results

def check_nbrs(data, results):
    results.nbr_matrix2 = np.copy(results.nbr_matrix1)
    results.dist_matrix2 = np.copy(results.dist_matrix1)

    for cell1 in range(data.num_cells):
        for cell2 in range(data.num_cells):
            if results.thresh_method == 'Automatic (cell position dependent)':
                if results.outside_bool2[cell1] == 1 or results.outside_bool2[cell2] == 1:
                    thresh = results.thresh_val[0]
                else:
                    thresh = results.thresh_val[1]
            else:
                thresh = results.thresh_val

            if results.nbr_matrix1[cell1, cell2] == 1 and results.dist_matrix1[cell1, cell2] > thresh:
                # remove nbr from nbr_matrix1
                # remove dist from dist_matrix1
                results.nbr_matrix2[cell1, cell2] = 0
                results.dist_matrix2[cell1, cell2] = 'nan'

    return results

def dt_withcav(data, results):
    print('...Identifying cavity adjacent cells.')
    pts = np.concatenate([data.xyz, results.cavity_pts])
    n = len(pts)
    n_cav = len(results.cavity_pts)

    dt = Delaunay(pts)

    nbr_matrix = np.zeros([n, n])
    for cell in range(n):
        for sim in range(0, len(dt.simplices)):
            if cell in dt.simplices[sim,:]:
                nbr_matrix[cell, dt.simplices[sim,]] = 1
        nbr_matrix[cell, cell] = 0 # stop double counting

    cav_adj_ids=[]
    for i in range(n_cav):
        i_cav = data.num_cells + i
        ids = np.where(nbr_matrix[i_cav, :] == 1)
        cav_adj_ids.extend(ids[0])

    results.cav_adj_ids = np.unique(cav_adj_ids)
    results.cav_adj_ids = results.cav_adj_ids[results.cav_adj_ids < data.num_cells]
    results.cav_adj_bool = np.zeros(data.num_cells)
    results.cav_adj_bool[results.cav_adj_ids] = 1
    results.all['cavity_adj_bool'] = results.cav_adj_bool
    print('----- %2d cavity adjacent cells identified.' % len(results.cav_adj_ids))

    return results

def compile_results(data, results):
    # number of neighbours
    results.all['num_nbrs'] = np.sum(results.nbr_matrix2, axis=1)

    # number of neighbours of each type
    results.all['num_nbrs_outside'] = np.zeros(data.num_cells)
    results.all['num_nbrs_inside'] = np.zeros(data.num_cells)
    for cell1 in range(data.num_cells):
        for cell2 in range(data.num_cells):
            if results.nbr_matrix2[cell1, cell2] == 1:
                if results.outside_bool2[cell2] == 1:
                    results.all.loc[cell1, 'num_nbrs_outside'] = results.all.loc[cell1, 'num_nbrs_outside'] + 1
                else:
                    results.all.loc[cell1, 'num_nbrs_inside'] = results.all.loc[cell1, 'num_nbrs_inside'] + 1

    # neighbour list
    results.all['nbr_ids'] = [''] * data.num_cells
    temp = []

    for cell1 in range(data.num_cells):
        nbrs = np.where(results.nbr_matrix2[cell1, :] == 1)
        results.all.loc[cell1,'nbr_ids'] = str(nbrs[0])

    # median and range of cell-nbr distances
    results.all['nbr_dist_mean'] = np.zeros(data.num_cells)
    results.all['nbr_dist_range'] = np.zeros(data.num_cells)

    for cell1 in range(data.num_cells):
        results.all.loc[cell1, 'nbr_dist_mean'] = np.nanmean(results.dist_matrix2[cell1, :])
        results.all.loc[cell1, 'nbr_dist_range'] = np.nanmax(results.dist_matrix2[cell1,:]) - np.nanmin(results.dist_matrix2[cell1,:])


    return results

def compile_distances(data, results):
    tot_num_nbrs = int(np.sum(results.nbr_matrix2) / 2)

    results.dist_all['cell_id1'] = np.zeros(tot_num_nbrs)
    results.dist_all['cell_id2'] = np.zeros(tot_num_nbrs)

    # num_props = np.shape(data.properties)[1]
    # prop_heads = data.headings[4:-1]
    # for i in range(num_props):
    #
    #     results.dist_all[prop_heads[i]] = np.zeros(int(tot_num_nbrs))

    results.dist_all['outside_bool1'] = np.zeros(tot_num_nbrs)
    results.dist_all['outside_bool2'] = np.zeros(tot_num_nbrs)
    if len(results.cavity_pts)>0:
        results.dist_all['cavity_adj_bool1'] = np.zeros(tot_num_nbrs)
        results.dist_all['cavity_adj_bool2'] = np.zeros(tot_num_nbrs)

    results.dist_all['nbr_dist'] = np.zeros(tot_num_nbrs)

    counter = 0
    for cell1 in range(data.num_cells):
        for cell2 in range(cell1+1, data.num_cells):
            if results.nbr_matrix2[cell1, cell2] == 1 :
                results.dist_all.loc[counter, 'cell_id1'] = results.all.loc[cell1, 'ID']
                results.dist_all.loc[counter, 'cell_id2'] = results.all.loc[cell2, 'ID']

                results.dist_all.loc[counter, 'outside_bool1'] = results.all.loc[cell1, 'outside_bool']
                results.dist_all.loc[counter, 'outside_bool2'] = results.all.loc[cell2, 'outside_bool']

                if len(results.cavity_pts) > 0:
                    results.dist_all.loc[counter, 'cavity_adj_bool1'] = results.all.loc[cell1, 'cavity_adj_bool']
                    results.dist_all.loc[counter, 'cavity_adj_bool2'] = results.all.loc[cell2, 'cavity_adj_bool']

                results.dist_all.loc[counter, 'nbr_dist'] = results.dist_matrix2[cell1, cell2]

                counter += 1

    return results

def final_figure(directory, fname_only, data, results, plotting):
    print('...Press enter to close the figure.')
    plotting.face_col = np.multiply(np.ones([data.num_cells, 4]), plotting.all_col)

    if len(results.cavity_pts) > 0:
        for cell in range(data.num_cells):
            if results.all.loc[cell, 'cavity_adj_bool'] == 1:
                plotting.face_col[cell, :] = plotting.cav_col

    embryo_position = [0, 0, 1, 1]
    fig = plt.figure('IVEN - final classification', figsize=(5, 5))
    em_ax, sct_plt = plot_embryo(fig, data, plotting, embryo_position)

    for cell in range(data.num_cells):
        em_ax.text(data.all.loc[cell, data.headings[1]], data.all.loc[cell, data.headings[2]],
                   data.all.loc[cell, data.headings[3]], str(data.all.loc[cell, 'ID']), fontsize=10)

    # add wait to close figure with 'Enter button'
    def press(event):
        if event.key == 'enter':  # to close figure press enter as before
            plt.close(fig)
        return ()

    output_fname_img = directory + '/IVEN2img_' + fname_only[0:-4] + '.png'
    fig.canvas.draw()
    fig.canvas.mpl_connect('key_press_event', press)
    plt.savefig(output_fname_img, dpi=500)
    plt.show()


    return plotting



