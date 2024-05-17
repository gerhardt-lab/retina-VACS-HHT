import seaborn as sns
from pylab import plt
import os
import scipy
from scipy import stats
import shutil
import numpy as np
import math
from skimage import io
import scipy.misc
import matplotlib.colors as colors
import pandas as pd
from statannotations.Annotator import Annotator


sns.set_context('poster')
sns.set_style('white')
sns.set_style("ticks", {"xtick.major.size": 1, "ytick.major.size": 1})
"""
sns.jointplot(rearranged_data['GFP']['CDH5CreERT2 YFP']['P3_P7']['distances_AV_0to1200'], rearranged_data['GFP']['CDH5CreERT2 YFP']['P3_P7']['distances_R_0to1200'], kind='kde', color='#ff7f0e', space=0, ylim=[0,1200], xlim=[0,1]).set_axis_labels('$\phi_{v-a}$', '$r \; [\mu m]$')
plt.savefig(dir_plots + 'R0to1200_ctr_kdePlot_' + str(tp) + '.png', format="png", bbox_inches = "tight", dpi=150)
plt.savefig(dir_plots + 'R0to1200_ctr_kdePlot_' + str(tp) + '.pdf', format="pdf", bbox_inches = "tight", dpi=150)
plt.savefig(dir_plots + 'R0to1200_ctr_kdePlot_' + str(tp) + '.eps', format="eps", bbox_inches = "tight", dpi=150)

sns.jointplot(rearranged_data['GFP']['CDH5CreERT2;YFP;YesFloxFlox']['P3_P7']['distances_AV_0to1200'], rearranged_data['GFP']['CDH5CreERT2;YFP;YesFloxFlox']['P3_P7']['distances_R_0to1200'], kind='kde', space=0, ylim=[0,1200], xlim=[0,1]).set_axis_labels('$\phi_{v-a}$', '$r \; [\mu m]$')
plt.savefig(dir_plots + 'R0to1200_YesKO_kdePlot_' + str(tp) + '.png', format="png", bbox_inches = "tight", dpi=150)
plt.savefig(dir_plots + 'R0to1200_YesKO_kdePlot_' + str(tp) + '.pdf', format="pdf", bbox_inches = "tight", dpi=150)
plt.savefig(dir_plots + 'R0to1200_YesKO_kdePlot_' + str(tp) + '.eps', format="eps", bbox_inches = "tight", dpi=150)

plt.xlim(0, 1e-52)
sns.distplot(p_values[tp]['GFP']['distances_AD_0to1200'], bins=100)
"""
#======================================================================================================================
#=                                                 Plotting histograms                                                =
#======================================================================================================================

def plot_all_hists(parameters, all_data, df, exp_cond, dir_plots):
    print('Plotting combined histograms...')
    plot_combined_hists(parameters, all_data, dir_plots)
    #compare_normalization_sloppily_histograms(parameters, all_data, dir_append, dir_plots)
    
    if parameters['radial_bins'] != []:
        print('Plotting histograms for radial bins...')
        plot_combined_hists_for_radial_bins(parameters, all_data, dir_plots)
    
    if parameters['plot_kde']:
        print('Plotting kde plots...')
        kde_plots(parameters, all_data, dir_plots)
        
    # plot individual histograms if specified in parameter file       
    if parameters['plot_individual_hist']:
        print('Plotting individual histograms...')
        individual_hists_allPixels(parameters, df[df[exp_cond]==1], dir_plots)

def plot_combined_hists(parameters, all_data, dir_plots):
    hfont = {'fontname':'Helvetica'}
    time_points = parameters['time_points']
    dir_name = dir_plots
    
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    
    types = ['AV', 'AD', 'VD', 'R', 'AD_outSF', 'AV_outSF', 'AV_inSF']
    xlabels = ['$\phi_{v-a}$', '$distance \; to \; artery \; [\mu m]$' , '$distance \; to \; vein \; [\mu m]$', '$r \; [\mu m]$', '$distance \; to \; artery \; [\mu m]$', '$\phi_{v-a}$', '$\phi_{v-a}$']
    
    for typ in range(len(types)):
        plt.figure()
        for i in range(len(time_points)):
            if time_points[i] in all_data.keys():
                thisLabel = time_points[i].split("_")[1]
                median = all_data[time_points[i]]['median_' + types[typ]]
                sns.distplot(all_data[time_points[i]]['distances_' + types[typ]], bins=30, color=all_data[time_points[i]]['color'], label=thisLabel)
                plt.axvline(x=median, color = all_data[time_points[i]]['color'])
            else:
                print('Function plot_combined_hists: There is no data for '+str(time_points[i]) + ' in the dataframe.')
        plt.xlabel(xlabels[typ], **hfont)
        plt.title(types[typ])
        plt.legend(loc = 'upper right', prop={'size': 16})
        plt.tick_params(axis='both', which='major', labelsize=20)

        plt.savefig(dir_name + 'combined_' + types[typ] + '_hist_' + str(time_points) + '.png', format="png", bbox_inches = "tight", dpi=150)
        plt.savefig(dir_name + 'combined_' + types[typ] + '_hist_' + str(time_points) + '.pdf', format="pdf", bbox_inches = "tight", dpi=150)
        plt.close()

def plot_combined_hists_for_radial_bins(parameters, all_data, dir_plots):
    hfont = {'fontname':'Helvetica'}
    time_points = parameters['time_points']
    dir_name = dir_plots + '/histograms_radial_bins/'
    
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    
    types = ['AV', 'AD', 'VD']
    xlabels = ['$\phi_{v-a}$', '$distance \; to \; artery \; [\mu m]$' , '$distance \; to \; vein \; [\mu m]$']
    
    for typ in range(len(types)):
        for l in parameters['radial_bins']:
            plt.figure()
            for i in range(len(time_points)):
                if time_points[i] in all_data.keys():
                    thisLabel = time_points[i].split("_")[1]
                    median = all_data[time_points[i]]['median_' + types[typ]+ '_' + str(l[0]) + 'to' + str(l[1])]
                    sns.distplot(all_data[time_points[i]]['distances_' + types[typ] + '_' + str(l[0]) + 'to' + str(l[1])], bins=30, color=all_data[time_points[i]]['color'], label=thisLabel)
                    plt.axvline(x=median, color = all_data[time_points[i]]['color'])
                else:
                    print('Function plot_combined_hists: There is no data for '+str(time_points[i]) + ' in the dataframe.')
            plt.xlabel(xlabels[typ], **hfont)
            plt.legend(loc = 'upper right', prop={'size': 16})
            plt.tick_params(axis='both', which='major', labelsize=20)
            plt.title('(' + str(l[0]) + ' - ' + str(l[1]) + ') ' + '$\mu m$', fontsize = 18)
    
            plt.savefig(dir_name + 'combined_' + types[typ] + '_hist_' + str(l[0]) + 'to' + str(l[1]) + '_' + str(time_points) + '.png', format="png", bbox_inches = "tight", dpi=150)
            plt.savefig(dir_name + 'combined_' + types[typ] + '_hist_' + str(l[0]) + 'to' + str(l[1]) + '_' + str(time_points) + '.pdf', format="pdf", bbox_inches = "tight", dpi=150)
            plt.close()            
        
        for i in range(len(time_points)):
            plt.figure()
            ax = plt.gca()
            if time_points[i] in all_data.keys():
                for l in parameters['radial_bins']:
                    color = next(ax._get_lines.prop_cycler)['color']
                    thisLabel = '(' + str(l[0]) + ' - ' + str(l[1]) + ') ' + '$\mu m$'
                    median = all_data[time_points[i]]['median_' + types[typ]+ '_' + str(l[0]) + 'to' + str(l[1])]
                    #percentile25 = np.percentile(all_data[time_points[i]]['distances_' + types[typ]+ '_' + str(l[0]) + 'to' + str(l[1])], 25)
                    #percentile75 = np.percentile(all_data[exp_conds[i]][tp]['distances_' + types[typ]+ '_' + str(l[0]) + 'to' + str(l[1])], 75)
                    sns.distplot(all_data[time_points[i]]['distances_' + types[typ] + '_' + str(l[0]) + 'to' + str(l[1])], bins=30, color=color, label=thisLabel)
                    plt.axvline(x=median, color=color)
                    #plt.axvline(x=percentile25, ls='--', color=color)
                    #plt.axvline(x=percentile75, ls=':', color=color)
                plt.xlabel(xlabels[typ], **hfont)
                plt.legend(loc = 'upper right', prop={'size': 16})
                plt.tick_params(axis='both', which='major', labelsize=20)
                plt.title(time_points[i], fontsize = 18)
    
                plt.savefig(dir_name + str(time_points[i]) + '_' + types[typ] + '_hist_radial_bins_' + str(parameters['radial_bins']) + '.png', format="png", bbox_inches = "tight", dpi=150)
                plt.savefig(dir_name + str(time_points[i]) + '_' + types[typ] + '_hist_radial_bins_' + str(parameters['radial_bins']) + '.pdf', format="pdf", bbox_inches = "tight", dpi=150)
                plt.close()  

def plot_compare_exp_conditions_histograms(parameters, all_data, tp, dir_plots):
    hfont = {'fontname': 'Helvetica'}
    exp_conds = parameters['load_conditions']
    dir_name = dir_plots

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    types = ['AV', 'AD', 'VD', 'R', 'AV_0to1200']#, 'AD_outSF', 'AV_outSF', 'AV_inSF'
    xlabels = ['$\phi_{v-a}$', '$distance \; to \; artery \; [\mu m]$', '$distance \; to \; vein \; [\mu m]$',
               '$r \; [\mu m]$', '$\phi_{v-a}$']#'$distance \; to \; artery \; [\mu m]$', '$\phi_{v-a}$', , '$\phi_{v-a}$'

    for typ in range(len(types)):
        plt.figure()
        ax = plt.gca()
        for i in range(len(exp_conds)):
            thisLabel = exp_conds[i]
            color = next(ax._get_lines.prop_cycler)['color']
            median = all_data[exp_conds[i]][tp]['median_' + types[typ]]
            percentile25 = np.percentile(all_data[exp_conds[i]][tp]['distances_' + types[typ]], 25)
            percentile75 = np.percentile(all_data[exp_conds[i]][tp]['distances_' + types[typ]], 75)
            sns.distplot(all_data[exp_conds[i]][tp]['distances_' + types[typ]], bins=30,
                             color=color, label=thisLabel)
            plt.axvline(x=median, color=color)
            plt.axvline(x=percentile25, ls='--', color=color)
            plt.axvline(x=percentile75, ls=':', color=color)

        plt.xlabel(xlabels[typ], **hfont)
        plt.title(tp + ', ' + types[typ])
        plt.legend(loc='upper right', prop={'size': 16})
        plt.tick_params(axis='both', which='major', labelsize=20)

        plt.savefig(dir_name + 'combined_' + types[typ] + '_hist_' + str(exp_conds) + '_' + tp + '.png', format="png",
                    bbox_inches="tight", dpi=150)
        plt.savefig(dir_name + 'combined_' + types[typ] + '_hist_' + str(exp_conds) + '_' + tp + '.pdf', format="pdf",
                    bbox_inches="tight", dpi=150)
        plt.close()


def kde_plots(parameters, all_data, dir_plots):
    time_points = parameters['time_points']
    dir_name = dir_plots + '/kde_plots/'
    
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    types = ['AV', 'AD', 'VD']
    xlabels = ['$\phi_{v-a}$', '$distance \; to \; artery \; [\mu m]$' , '$distance \; to \; vein \; [\mu m]$']

    for tp in time_points:
        if tp in all_data.keys():
            for art in range(len(types)):
                plt.figure()
                g = sns.jointplot(all_data[tp]['distances_' + types[art]], all_data[tp]['distances_R'], kind='kde', space=0, color=all_data[tp]['color'], ylim=[0,2200], fill='True').set_axis_labels(xlabels[art], '$r \; [\mu m]$')
                for i in range(len(time_points)):
                    g.ax_joint.axvline(all_data[time_points[i]]['median_' + types[art]],ls='--', lw= 1.5, color=all_data[time_points[i]]['color'])
                    g.ax_joint.axhline(all_data[time_points[i]]['median_R'],ls='--', lw= 1.5, color=all_data[time_points[i]]['color'])
                    g.ax_marg_x.axvline(all_data[time_points[i]]['median_' + types[art]],ls='--', lw= 1.5, color=all_data[time_points[i]]['color'])
                    g.ax_marg_y.axhline(all_data[time_points[i]]['median_R'],ls='--', lw= 1.5, color=all_data[time_points[i]]['color'])
                plt.title(tp, fontsize = 18)
                g.savefig(dir_name + 'kde_' + types[art] + '_R_hist_' + str(tp) + '_' + str(time_points) + '.png', format="png", bbox_inches = "tight", dpi=150)
                g.savefig(dir_name + 'kde_' + types[art] + '_R_hist_' + str(tp) + '_' + str(time_points) + '.pdf', format="pdf", bbox_inches = "tight", dpi=150)
                plt.close()
        else:
            print('Function kde_AV_R: There is no data for '+ str(tp) + ' in the dataframe.')

def cond_kde_plots(parameters, all_data, dir_plots):
    time_points = parameters['time_points']
    conds = parameters['load_conditions']
    
    types = ['AV', 'AV_0to1200']#, 'AV_outSF'
    radii = ['R', 'R_0to1200']#, 'R_outSF'
    xlabels = ['$\phi_{v-a}$', '$\phi_{v-a}$']#, '$\phi_{v-a}$']
    ylims = [[0,2200], [0,1200]]#, [0,1200]]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    dir_name = dir_plots + '/cond_kde_plots/'

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    
    for tp in time_points: 
        for typ in range(len(types)):
            for cond in range(len(conds)):
                plt.figure()      
                sns.jointplot(all_data[conds[cond]][tp]['distances_' + types[typ]], all_data[conds[cond]][tp]['distances_' + radii[typ]], kind='kde', space=0, color=colors[cond], ylim=ylims[typ], fill='True').set_axis_labels(xlabels[typ], '$r \; [\mu m]$')
                #plt.title(conds[cond] + types[typ])
                plt.savefig(dir_name + types[typ] + '_' + conds[cond] + '_kdePlot_' + str(tp) + '.png', format="png", bbox_inches = "tight", dpi=150)
                plt.savefig(dir_name + types[typ] + '_' + conds[cond] + '_kdePlot_' + str(tp) + '.pdf', format="pdf", bbox_inches = "tight", dpi=150)
                plt.close()
    
def individual_hists_allPixels(parameters, df, dir_plots):
    time_points = np.unique(df['TP'].values)
    dir_name = dir_plots + 'individual_hists/'

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    
    types = ['AV', 'AD', 'VD']
    
    for typ in types:
        for tp in range(len(time_points)):
            numRows = math.ceil(df[df['TP'] == time_points[tp]].shape[0]/3)
            fig, axes = plt.subplots(numRows,3, sharex=True, sharey=True)
            count = 0
            for i, retina in df[df['TP'] == time_points[tp]].iterrows():
                sns.distplot(retina['distances_' + typ], bins=30, ax=axes.flat[count])
                count += 1
            plt.suptitle(str(time_points[tp]))
            plt.savefig(dir_name + 'individual_' + typ +'_hists_allPixels_' + str(time_points[tp]) + '.png', format="png", bbox_inches = "tight", dpi=150)
            plt.savefig(dir_name + 'individual_' + typ +'_hists_allPixels_' + str(time_points[tp]) + '.pdf', format="pdf", bbox_inches = "tight", dpi=150)
            plt.close()

#-----------------
def compare_normalization_sloppily_histograms(parameters, all_data, dir_plots):
    time_points = all_data.keys()
    dir_name = dir_plots
    
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    
    types = ['AV', 'R', 'AD', 'VD']
    
    numRows = len(time_points)
    fig, axes = plt.subplots(numRows,4, figsize=(30,20))#, sharex=True, sharey=True)
    plt.subplots_adjust(hspace = 0.6)
    count = 0
    for tp in time_points:
        for typ in range(len(types)):
            v_mps = all_data[tp]['v_' + types[typ] + '_mps']
            hist = all_data[tp]['mean_' + types[typ] + '_GFP']
            
            distances = []
            
            for i in range(len(hist)):
                number = int(100000*hist[i])
                for j in range(number):
                    distances.append(v_mps[i])
            
            g = sns.distplot(all_data[tp]['distances_' + types[typ]], bins=len(v_mps), color='orange', ax=axes[count][typ], label='Our approach')
            sns.distplot(distances, bins=int(len(v_mps)*0.6), color='green', ax=axes[count][typ], label='w/ np.hist fct')
            g.set_title(types[typ]+', '+tp, fontsize = 18)
            g.legend(loc = 'best', prop={'size': 15})
            
        count += 1            
    plt.savefig(dir_name + '/compare_normalizations.png', format="png", bbox_inches = "tight", dpi=150)  
    plt.savefig(dir_name + '/compare_normalizations.pdf', format="pdf", bbox_inches = "tight", dpi=150)
    plt.close()          

#======================================================================================================================
#=                                              Plotting coordinate system                                            =
#======================================================================================================================
def plot_coordinate_system(parameters, key_file, experimentID):
    if experimentID in key_file["ExperimentID"].unique():
        experiment_df  = key_file[key_file["ExperimentID"] == experimentID]
        row_drawn = experiment_df[experiment_df["Drawn"] == 1].iloc[0]
        scale_factor = key_file[key_file["ExperimentID"] == experimentID].iloc[0][28]
        
        out_dir = parameters['out_dir'] + 'processed/plots/'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        drawn_filename = parameters["data_dir"] + row_drawn['filename']
    
        img_map = np.array(io.imread(drawn_filename)).astype('bool')
        print(img_map.shape)
        # problem that accoured with Yi's retinas when they were RGB
        if img_map.ndim != 3:
            print('WARNING! Something was wrong with the draw file dimensions...')
            img_map = np.moveaxis(img_map[:,:,:,0], 0, -1)
        # masks that included the sprounting front showed a different dimension
        if img_map.shape[0] < img_map.shape[2]:
            print('WARNING! Something was wrong with the draw file dimensions.')
            print('Old dimensions: ' + str(img_map.shape))
            img_map = np.moveaxis(img_map[:,:,:], 0, -1)
            print('New dimensions: ' + str(img_map.shape))
        if img_map.shape[2] > 4:
            if np.count_nonzero(img_map[:,:,4] == True) > (img_map.shape[0]*img_map.shape[1]/2):
                img_map[:,:,4] = np.invert(img_map[:,:,4])
            img_map[:,:,4] = np.invert(scipy.ndimage.morphology.binary_fill_holes(img_map[:,:,4]))
        
        pos_on = 0 	
        pos_artery = 1
        pos_vein = 2
        pos_mask = 3

        img_map[:, :, pos_mask] = scipy.ndimage.morphology.binary_fill_holes(img_map[:, :, pos_mask])

        artery_distance = scale_factor*scipy.ndimage.morphology.distance_transform_edt(np.invert(img_map[:, :, pos_artery]))
        vein_distance = scale_factor*scipy.ndimage.morphology.distance_transform_edt(np.invert(img_map[:, :, pos_vein]))
        optical_distance = scale_factor*scipy.ndimage.morphology.distance_transform_edt(np.invert(img_map[:, :, pos_on]))
        va_pos = vein_distance / (artery_distance + vein_distance)
    
        #####################################################
        ########### plot AV coordinate system ###############
        #####################################################
    
        fig = plt.figure()
		# mask image
        va_pos_masked = np.ma.masked_where(img_map[:, :, pos_mask]==0, va_pos)
        # plot image with defined colormap palette_R
        palette = plt.cm.coolwarm
        palette.set_bad('w', 1.0)
        im = plt.imshow(va_pos_masked, interpolation='bilinear', cmap=palette,
		                norm=colors.Normalize(vmin=0.0, vmax=1.0),
		                aspect='auto', origin='lower')
		
        # plot contour lines
        CS = plt.contour(va_pos_masked, [0.2, 0.5, 0.8], colors=('b', 'gray', 'r'), linewidths=1)

    	# label contour lines
        plt.clabel(CS, [0.2, 0.5, 0.8], inline=1, fmt='%1.1f', fontsize=8)
        		
        # add color bar
        cbar = fig.colorbar(im)
        # remove axis labels
        im.axes.get_xaxis().set_visible(False)
        im.axes.get_yaxis().set_visible(False)
		
        plt.axis('off')
        plt.savefig(out_dir + 'coordinate_system_ExpID_' + str(experimentID) + '_AV.png', format="png", bbox_inches = "tight", dpi=150)
        plt.savefig(out_dir + 'coordinate_system_ExpID_' + str(experimentID) + '_AV.pdf', format="pdf", bbox_inches = "tight", dpi=150)
        plt.close()
    
        #####################################################
        ########### plot radial coordinate system ###########
        #####################################################    

        fig = plt.figure()
		
        # mask image
        r_pos_masked = np.ma.masked_where(img_map[:, :, pos_mask]==0, optical_distance)
		
        # plot image with defined colormap palette_R
        palette_R = plt.cm.RdPu
        palette_R.set_bad('white', 1.0)
        im = plt.imshow(r_pos_masked, interpolation='bilinear',
		                cmap=palette_R, norm=colors.Normalize(vmin=0.0, vmax=np.max(r_pos_masked.flatten())),
		                aspect='auto', origin='lower')
		
        # plot contour lines
        CS = plt.contour(r_pos_masked, [1000, 2000, 3000], linewidths=2, cmap=palette_R)
		
        # label contour lines
        plt.clabel(CS, [1000, 2000], inline=1, fmt='%1.1f', fontsize=12)
		
        # add color boar
        cbar = fig.colorbar(im)
		
        plt.axis('off')
		
        # remove axis labels
        im.axes.get_xaxis().set_visible(False)
        im.axes.get_yaxis().set_visible(False)
        plt.savefig(out_dir + 'coordinate_system_ExpID_' + str(experimentID) + '_radial.png', format="png", bbox_inches = "tight", dpi=150)
        plt.savefig(out_dir + 'coordinate_system_ExpID_' + str(experimentID) + '_radial.pdf', format="pdf", bbox_inches = "tight", dpi=150)
        plt.close()
    else:
        print('Experiment ID %s is not listed in the key file.' % experimentID)        

def single_retina_mean_box_plot(parameters, key_file, dir_append):

        #inputDir = "data/processed/python_results/"
    inputDir = parameters["out_dir"] + "processed/python_results/" 
    outputDir = parameters["out_dir"] + "processed/plots/" 

    ax_labels = {'AV': '$\phi_{v-a}$', 'R': '$r \; [\mu m]$'}
    
    #folderList = dict()
    #for path, subdirs, files in os.walk(inputDir):
    #    if 'ExperimentID' in path and dir_append in path:
    #        folderLista.append(path)
 

    #folderList = []
    #for path, subdirs, files in os.walk(inputDir):
    #    if 'ExperimentID' in path and dir_append in path:
    #        folderList.append(path)
    #folder_list = np.unique(folderList)

    #res_list = []

    ######
    ## Compute plot features
    ###### 


    plot_features = pd.DataFrame()
    counter = 0
    
    for index, row in key_file.iterrows():
        #res_folder in folder_list:
        res_folder_base = inputDir + "ExperimentID_" + str(row['ExperimentID']) + "_"
 
        if not row[dir_append]:
            continue            
       
        for time_points in parameters['time_points']:
            res_folder = res_folder_base + time_points + "/" + dir_append  + "/"
            if os.path.exists(res_folder):
                # print("Folder exists:")
                # print(res_folder)
                distances_AV = np.load(res_folder + '/distances_AV.npy')
                # print('This is what the AV numpy array looks like:')
                # print(distances_AV)
                # print(np.isnan(distances_AV).any())
                # print('This is what the R numpy array looks like:')
                distances_R = np.load(res_folder + '/distances_R.npy')
                # print(distances_R)
                # print(np.isnan(distances_R).any())
                #
                # print("Shape of distances_AV")
                # print(distances_AV.shape)
                #
                # print("Shape of distances_R")
                # print(distances_R.shape)
                # print("Mean of distances_R")
                # print(np.mean(distances_R))

                plot_features.at[counter, "ExperimentID"] = row['ExperimentID']
                plot_features.at[counter, "time_interval"] = time_points

                if row['Collection Point'] == 15:
                    plot_features.at[counter, "time_point"] = "P15"
                if row['Collection Point'] == 7:
                    plot_features.at[counter, "time_point"] = "P7"
                if row['Collection Point'] == 8:
                    plot_features.at[counter, "time_point"] = "P8"
                
                for condition in parameters["load_conditions"]:
                    if row[condition]:
                        plot_features.at[counter, "condition"] = condition
                
                plot_features.at[counter, "mean_AV"] = np.mean(distances_AV)
                plot_features.at[counter, "mean_R"] = np.mean(distances_R)
                
                plot_features.at[counter, "median_AV"] = np.median(distances_AV)
                plot_features.at[counter, "median_R"] = np.median(distances_R)
 

                counter += 1
            else:
                print("Folder does not exist:")
                print(res_folder)
 
    print(str(plot_features["mean_AV"]))
    print(plot_features.keys())

    plot_dir = outputDir + "/box_plots/"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    time_points = ['P8', 'P15']
    time_intervals = parameters['time_points']
    conditions = parameters["load_conditions"]
 
    width = 10
    length = 9
   
    stats_df = pd.DataFrame()
    
    data_P7 = plot_features[plot_features["time_point"]=="P7"]
    data_P15 = plot_features[plot_features["time_point"]=="P15"]
    print(data_P15['ExperimentID'])
    test_smad4 = stats.ttest_ind(data_P7['mean_AV'], data_P15['mean_AV'], axis=0, equal_var=False,
                                nan_policy='propagate', permutations=None, random_state=None, alternative='two-sided')
    print("Welch's test Smad4")
    print(test_smad4)

    try:
        data_P7_yesfloxflox = data_P7[data_P7["condition"] == conditions[1]]
        data_P7_ctr= data_P7[data_P7["condition"] == conditions[0]]

        test_P7 = stats.ttest_ind(data_P7_ctr['mean_AV'], data_P7_yesfloxflox['mean_AV'], axis=0, equal_var=False, nan_policy='propagate', permutations=None, random_state=None, alternative='two-sided')
    
        print("Welch's test P7")
        print(test_P7)
    except:
        print("Welch's t_test to compare conditions was not performed for P7 as there was something wrong (E.g. not enough conditions to compare).")

    data_P15 = plot_features[plot_features["time_point"]=="P15"]

    try:
        data_P15_yesfloxflox = data_P15[data_P15["condition"] == conditions[1]]
        data_P15_ctr= data_P15[data_P15["condition"] == conditions[0]]


        test_P15 = stats.ttest_ind(data_P15_ctr['mean_AV'], data_P15_yesfloxflox['mean_AV'], axis=0, equal_var=False, nan_policy='propagate', permutations=None, random_state=None, alternative='two-sided')

        print("Welch's test P15")
        print(test_P15)
    except:
        print("Welch's t_test was not performed for P15 as there was something wrong (E.g. not enough conditions to compare).")

    # plot of the mean (per retina sample) AV values

    fig, ax = plt.subplots(figsize=(width, length))
    
    sns.boxplot(y="time_point", x="mean_AV", data=plot_features, hue='condition', orient = "h")
    sns.swarmplot(y="time_point", x="mean_AV", data=plot_features, hue='condition', orient = "h", size=15.0, color="k", dodge = True)
    #sns.boxplot(y="time_point", x="mean_AV", data=plot_features, orient="h")
    #sns.swarmplot(y="time_point", x="mean_AV", data=plot_features, orient = "h", size=15.0, color="k", dodge = True)

    # stats
    pairs=[((time_points[0],conditions[0]),(time_points[0],conditions[1])),
            ((time_points[1],conditions[0]),(time_points[1],conditions[1]))]
    # TODO Read time points to plot from parameters file instead of hardcoding here
    #pairs = [('P8', 'P15')]
    annotator = Annotator(ax, pairs, data=plot_features, y="time_point", x="mean_AV", hue="condition", orient = "h")
                                        #order = time_intervals_order, hue = "condition", hue_order = [condition_other,"control"])
    annotator.configure( test="t-test_welch", text_format="star", loc="inside")
    annotator.apply_and_annotate()
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], loc = 'upper right', fontsize = 20)
    
    ax.set_xlabel(ax_labels['AV']) #+ " median per sample")
    ax.set_ylabel("time")
    plt.tight_layout()   


    plot_path = plot_dir + "mean_distances_AV.pdf"
    plot_path_png = plot_dir + "mean_distances_AV.png"
    plt.savefig(plot_path)
    plt.savefig(plot_path_png)

    # TRIAL


        # plot of the mean (per retina sample) AV values
    
    fig, ax = plt.subplots(figsize=(10, 9))
    #sns.boxplot(x="time_point", y="mean_AV", data=plot_features, hue = "condition", palette=color_palette)
    #sns.boxplot(y="time_point", x="median_AV", data=plot_features, hue='condition', orient = "h")
    g_boxplot =sns.boxplot(y="time_point", x="median_AV", data=plot_features, hue = "condition", orient = "h")
    g_swarm = sns.swarmplot(y="time_point", x="median_AV", data=plot_features, hue = "condition", orient = "h", size=15.0, color="k", dodge = True)
    
    pairs=[((time_points[0],conditions[0]),(time_points[0],conditions[1])),
            ((time_points[1],conditions[0]),(time_points[1],conditions[1]))]
    #pairs = [('P8', 'P15')]

    annotator = Annotator(ax, pairs, data=plot_features, y="time_point", x="median_AV", hue='condition', orient = "h")
    #annotator = Annotator(ax, pairs, data=plot_features, y="time_point", x="median_AV", orient="h")
    annotator.configure( test="t-test_welch", text_format="star", loc="inside")
    annotator.apply_and_annotate()

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], loc = 'upper right', fontsize = 20)
    
    ax.set_xlabel(ax_labels['AV']) #+ " median per sample")
    ax.set_ylabel("time")
    plt.tight_layout()   

    plot_path = plot_dir + "median_distances_AV.pdf"
    plot_path_png = plot_dir + "median_distances_AV.png"
    plt.savefig(plot_path)
    plt.savefig(plot_path_png)

    # plot of the mean (per retina sample) R values
    
    fig, ax = plt.subplots(figsize=(9, 10))
    #sns.boxplot(x="time_point", y="mean_AV", data=data, hue = "condition", palette=color_palette)
    sns.boxplot(x="time_point", y="mean_R", data=plot_features, hue = "condition")
    #sns.boxplot(x="time_point", y="mean_R", data=plot_features)
    sns.swarmplot(x="time_point", y="mean_R", data=plot_features, hue = "condition", size=15.0, color="k", dodge = True)
    #sns.swarmplot(x="time_point", y="mean_R", data=plot_features, size=15.0, color="k", dodge=True)

    pairs=[((time_points[0],conditions[0]),(time_points[0],conditions[1])),
            ((time_points[1],conditions[0]),(time_points[1],conditions[1]))]
    #pairs = [('P8', 'P15')]

    annotator = Annotator(ax, pairs, data=plot_features, x="time_point", y="mean_R", hue="condition")
    # annotator = Annotator(ax, pairs, data=plot_features, x="time_point", y="mean_R")
    annotator.configure( test="t-test_welch", text_format="star", loc="inside")
    annotator.apply_and_annotate()

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], loc = 'lower left', fontsize = 15)
    
    ax.set_xlabel(ax_labels['R']) #+ " median per sample")
    ax.set_ylabel("time")
    plt.tight_layout()   


    plot_path = plot_dir + "distances_R.pdf"
    plot_path_png = plot_dir + "distances_R.png"
    plt.savefig(plot_path)
    plt.savefig(plot_path_png)


    print(plot_features)
    plot_features.to_csv(plot_dir + "plot_features.csv", index = False)


    return 0
    
    

#======================================================================================================================
#=                                            Plots for accessing velocities                                          =
#======================================================================================================================


#======================================================================================================================
#=                                            Plots for accessing velocities                                          =
#======================================================================================================================



"""
def velocity_plot(parameters, all_data):
    
"""
