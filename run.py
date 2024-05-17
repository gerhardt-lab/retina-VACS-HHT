import sys
from skimage import io
import pandas as pd
import time

sys.path.append("src/")
import functions
import plotResults

parameters = functions.get_parameters()
key_file = functions.get_key_file(parameters)

key_file.to_csv(parameters["out_dir"] + "key_file.csv")

#======================================================================================================================
#=                                            STEP 1: COMPUTE ALL RESULTS                                             =
#======================================================================================================================

print('Computing results...')
functions.compute_res(parameters, key_file)
print('Finished computing.')


#======================================================================================================================
#=                                              STEP 2: LOAD RESULTS                                                  =
#======================================================================================================================

# dataframe gets too big when loading all results together. Splitting it into 'GFP' and 'GFP X ERG' instead
# initializing a dataframe for that
all_data = {}
rearranged_data = {}
for image_type in parameters['image_types']:
    rearranged_data[image_type] = {}
    try: 
        print()
        print('++++++++++++++++++++++++ %s ++++++++++++++++++++++++' %image_type)
        #print('Loading results...')
        #print("image type column")
        ##print(key_file[image_type])
        #print("key_file extract")
        #print(key_file[key_file[image_type] == 1])
        all_data[image_type] = pd.merge(key_file[key_file[image_type] == 1],functions.load_res(parameters, image_type))
    except:
        print('Error: %s is not a valid image type.' %image_type)
        break

#======================================================================================================================
#=                              STEP 3: PLOT RESULTS FOR DIFFERENT EXPERIMENTAL CONDITIONS                            =
#======================================================================================================================

    print("key file columns")
    print(key_file.columns)

    plotResults.single_retina_mean_box_plot(parameters, key_file, image_type)

    # plot results depending on what experimental condition(s) is/are specified in parameters file ()
    for i in parameters['load_conditions']:
        print()
        print('EXPERIMENTAL CONDITION: ' + i)
        print()
    
        folder = i
            
        #try:
        print('Rearranging results for plotting...')
        rearranged_data[image_type][i] = functions.rearrange_results_for_plotting(parameters,all_data[image_type][all_data[image_type][i]==1])
        dir_plots = parameters['out_dir'] + 'processed/plots/pixels_per_retina_' + str(parameters['pixel_no_per_retina']) + '/' + image_type + '/' + folder + '/'
        print()
        # Plot histograms and kde plots
        plotResults.plot_all_hists(parameters, rearranged_data[image_type][i], all_data[image_type], i, dir_plots)

    # plot
    dir_plots = parameters['out_dir'] + 'processed/plots/pixels_per_retina_' + str(parameters['pixel_no_per_retina']) + '/' + image_type + '/'
    print('Plotting histograms for different experimental conditions...')
    for tp in parameters['time_points']:
        plotResults.plot_compare_exp_conditions_histograms(parameters, rearranged_data[image_type], tp, dir_plots)
    
    if parameters['plot_cond_kde']:
        print('Plotting kde plots for different conditions...')
        plotResults.cond_kde_plots(parameters, rearranged_data[image_type], dir_plots)

# 
#print('KS test...')
#p_values = {tp : {} for tp in parameters['time_points']}
#for tp in parameters['time_points']:
#    p_values[tp] = functions.confidence_KS_test(parameters, all_data, tp, parameters['load_conditions'])
    
functions.job_finished() 
