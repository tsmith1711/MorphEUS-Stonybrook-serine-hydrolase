morpheus-serine-hydrolase

The contents of this git include:

The workspaces generated from applying on thl, aa691, aa692, aa701, and aa702 onto the 32 drug morphEUS workspace.

The workspaces for the high dose,low dose, and joint profiles are found in their own corresponding folders. The file path for these spaces is:

analysis>bayes>stonybrook

To recreate the figures used in the paper run recreate_figures.m. This script needs to be ran in an iterative manner according to the order desired for each drug in the figure. For each run change the name on line 50 to the desired drug name then run the script. Once the script has been ran to include all desired data columns copy and paste lines 152 - 154 in the command window. 

names of drugs in this study when entering in line 50:

d691, d692, d701, d702, thl_3x

To plot the broad category classification for the compounds desired, paste the following in the command line:

make_connectivity_heatmap(blank_table,false);
xtickangle(90)
colormap(parula_plus)

To plot the nearest neighbor drug profiles for the compounds desired, paste the following in the command line:

make_connectivity_heatmap(blank_table_drug,false);
xtickangle(90)
colormap(parula_plus)

The lines of code above are lines 152-154 of the recreate_figures.m script.

The raw data workspace that contains all of the feature data can be found in the workspaces folder titled "stonybrook". This file contains all of the data from the aa691,692,701,702, and thl treatments. Each of these profiles were applied onto the morphEUS space previously published (using the classification_trials.m script) to generate the workspaces (in the bayes folder) used to generate the figures for this paper. 

