%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MorphEUS classification trials %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script to perform classification trials and MorphEUS classifications
%%% Description: 
% prompts user if doing joint or individual profile, asks for other
% settings, asks user to select drugs to include in calculations, asks user
% to select drug(s) to apply, then performs MorphEUS classification trials
% saves result at end of maxrun number runs
% loads joint profile, takes 80 untreated from each 025x and 3x workspace,
% runs backwards select, picks out the best result (decided by highest
% percent success, if there is a tie, then with the lowest variable number)
% records results from what was pointing to what, then repeat the process
% from the top with a new random set of untreated
% will save result in bayes and a log of all outputs in logs


%% clear workspace and set up paths 
clear

add_all_paths

rng('default')

% turn off TeX interpreter 
set(groot, 'defaultAxesTickLabelInterpreter', 'none')

% get user input for name of variable to save at end 
prompt = {'Enter a name for the variable information to be stored'};
dlgtitle = 'Save bayes_struct as:';
definput = {'bayes_struct'};
user_save_title = inputdlg(prompt,dlgtitle,[1 60],definput);


%% easy way to get and change the workspaces
%%% updated so that you only need joint workspace name and the other two
%%% will automatically take and separate the workspace approproiately 

%x3_workspace_name = "full_3x_workspace_with_stonybrook";
% if you don't have most updated for 3x, will default to full_3x_with_redos_no_A22_drug_name_only
%x025_workspace_name = "full_025x_workspace_with_stonybrook";
% if you don't have most updated for 025x, will default to full_025x_with_redos_v3_drug_names_only

%%%% BASE WORKSPACES %%%%
%%%% workspaces must be located in subfolder workspaces 
% requried spaces -- main data workspace 
morpheus_workspace = "morpheus_workspace_03232020.mat";
% second morpheus workspace
morpheus2_workspace = "morpheus2_baseworkspace_11132020.mat";
% third morphues workspace
morpheus3_workspace = "morpheus3_workspace_img_01062021.mat";
% morphEUS3 standard workspace
morpheus3_standard_workspace = "morphEUS3_standard_basespace_05142021.mat";

base_workspaces = [morpheus_workspace morpheus2_workspace morpheus3_workspace morpheus3_standard_workspace];

% allow user to select base workspace 
[indexes_chosen, any_chosen_tf] = listdlg('ListString',base_workspaces,...
            'Name',"Select base workspace",'ListSize',[220 300],'PromptString',"base_workspaces");

joint_workspace_name = base_workspaces(indexes_chosen);

% all applied workspaces 
condition_workspace_name = "all_condition_data_03232020.mat";
timecourse_workspace_name = "timecourse_data_03232020.mat";
clinical_workspace_name = "clinical_data_formatted_03232020.mat";  
stonybrook2h_workspace_name = "stonybrook_2h_03232020.mat";
stonybrook24h_workspace_name = "stonybrook_24h_03232020.mat";
applied_table_low_workspace_name = "low_dose_applied_07212020.mat";
applied_table_high_workspace_name = "high_dose_applied_07212020.mat";
batch_workspace_name = "Batch_experiment_workspace_10132020.mat";
newdrug_workspace_name = "newdrug_experiment_workspace_img_11192020.mat";
morpheus2_workspace_name = "morpheus2_11052020.mat";
morpheus3_standard_workspace_name = "morpheus3_standard_applied_05142021.mat";


applied_workspace_list = {timecourse_workspace_name, ...
    clinical_workspace_name,stonybrook24h_workspace_name,stonybrook2h_workspace_name,...
    condition_workspace_name,applied_table_low_workspace_name,applied_table_high_workspace_name,...
    batch_workspace_name,newdrug_workspace_name,morpheus2_workspace_name,morpheus3_standard_workspace_name}; 

all_workspaces = {joint_workspace_name, timecourse_workspace_name, ...
    clinical_workspace_name,stonybrook24h_workspace_name,stonybrook2h_workspace_name,...
    condition_workspace_name,batch_workspace_name,newdrug_workspace_name,morpheus2_workspace_name,morpheus3_standard_workspace_name}; 

%% create a log
% get time string for a unique log name/variable name at the end 
time_string = strrep(datestr(clock),":","-");

logname = strcat("./logs/bayes log ",user_save_title, time_string,".txt");

diary(logname)
% for clarity of reading diary 
warning off 

%% some other settings 
with_randomized_labels = false; %Different than random_comapre (just for this script)

% we are doing backwards select, so do feature reduction is true
do_feature_reduction = true;

 % set to on to create, off to not create (will be invisible with off)
create_cm = 'off'; 
show_pca = false; 

% settings for knn_helper
tit_add = "";
knn_type = 'median';
knn_with_apply = false;

% yes we are doing bayesian stuff
bayesian = true;

% for new52 joint profile
use_default = false;

% joint or individual profile question
joint_profile_question = questdlg('Joint or individual profiles', ...
	'joint profiles?', ...
	'joint','individual','joint');
switch joint_profile_question
    case 'joint'
        do_joint_profile = true;
    case 'individual'
        do_joint_profile = false; 
end 

% if not joint profile, 3x or 025x
if ~do_joint_profile
    dose_question = questdlg('3x or 025x', ...
	'dose?', ...
	'3x','025x','3x');
    switch dose_question
        case '3x'
            x3_dose = true;
        case '025x'
            x3_dose = false; 
    end 
    
else % if do joint profiles, do flattened or separate?    
    flatten_question = questdlg('Flattened or separate?', ...
        '???',...
        'flattened','separate','flattened');
    switch flatten_question
        case 'flattened'
            flattened_joint = true;
        case 'separate'
            flattened_joint = false; 
    end
end 

% the usual 
variables_to_create = {'disc',  'include_DMSO', 'large_group','nn2',  'no_TVN','apply_a_workspace','replicates_random','random_at_end','no_FM', 'no_syto','remove_feature_count', 'remove_heterogeneity_features','debug_mode'};
default_values =      [true,      false,          true,        false, false,     false,              false,                 false,    false,       false,          false,       false,   false];
disp("include_DMSO = include DMSO in TVN transformation")
disp("large group = categorize by larger groups instead of finer groups")
disp("nn2 = use second nearest neighbors as well")
disp("no_TVN = do not TVN transform")
disp("apply_a_workspace =  get prompted to apply a workspace")
disp("replicates_random = run this with random labels for each replicate")
disp("random_at_end = randomize the labels at the end, after doing backwards select")
disp("no_FM = do not use any FM features (both feature count and foci features)")
disp("no_syto = do not use any syto features (both feature count and foci features)")
disp("remove_feature_count = do not use feature count features")
disp("remove_heterogeneity_features = remove Q1, Q3 and IQR features") 
disp("debug_mode = run in debug mode")

target_values = checkboxList("Select Settings", variables_to_create, default_values);

initialize_variables


% we've been using essentially this set so let's have it get
% suggested
suggestion = {'Mer','Amp','Ctax','INH','EMB','ETA','IMI','Van','Cyc','Del',...
'Lev','Mox','Clz','MIT','Olf','Kan','Amk','Cam','Cla','Dox','Gent',...
'Strep','Tet','Lin','Pre','CCCP','Cer','Mon','Nig','Thi','RifT','BDQ',... 
'RIF','THL','water','Untreated'};

%% do first run separate to not mess with parfor

% do an initial load of workspace to check number of untreated that can max
% be used 
chosen_workspace = joint_workspace_name;
chosen_workspaces = {joint_workspace_name};

load_workspaces


num_untreated_in_wksp = sum(strcmp(final_data_table.DRUG,"Untreated"));
if do_joint_profile
    % divide by 2 bc after flattening will be half
    num_untreated_in_wksp = num_untreated_in_wksp/2;
end

% if debugging, only do 5 runs 
% num rows = number of untreated rows to keep 
% if debugging, do a small amount so it can run faster 
if debug_mode
    maxrun = 5; 
    num_rows = 32;
    number_untreated_used = 80;
else
    maxrun = 70;
    num_rows = 80;
    number_untreated_used = 80;
    if num_untreated_in_wksp < num_rows 
        num_rows = num_untreated_in_wksp;
        number_untreated_used = num_untreated_in_wksp;
    end
end 
runct = 1;

% don't want to just copy paste code again, this is the same content as the
% loop just over separately 
% make something to keep track of the variables
all_variable_sets = cell(maxrun,1); 
% make something to keep track of the results_tables
all_results_tables = cell(maxrun,1);

tic 

bayes_loop_content

%% run the rest of them -> use that cluster bby
% ah oh well we tried parfor is not working bby
for runct = 2:maxrun

    % since the exact same as before, just copied into a separate document 
    bayes_loop_content
   
end


%% auto save because this will take a while

if do_joint_profile
   disp("JOINT PROFILE")
   save_title = strcat("./bayes/",user_save_title," joint ", time_string);
else
    disp("INDIVIDUAL PROFILE")
    if x3_dose
        disp("3x dose")
         save_title = strcat("./bayes/",user_save_title," 3x ", time_string);
    else
        disp("025x dose")
         save_title = strcat("./bayes/",user_save_title," 025x ", time_string);
    end
end

% want to indicate in the title if we had random replciates 
if replicates_random
   save_title = strcat(save_title, " replicates randomized"); 
end

if random_at_end
    save_title = strcat(save_title, " random at end");   
end 

% let's save the git info as well
try
    git_info = getGitInfo();
    git_hash = git_info.hash;
catch
    git_hash = "unable to save git information--likely run using a previous commit";
    disp("unable to save git information") 
    
end

% if there was an applied drug, save that
if apply 
    save(save_title,'bayes_struct','all_variable_sets','all_results_tables','git_hash','applied_drug','final_data_table')

else
    % save a sample final data table (even if it's just from last run)
    save(save_title,'bayes_struct','all_variable_sets','all_results_tables','git_hash','final_data_table')

end 


disp(strcat("Saved bayes_struct as ", save_title))

diary off 
toc 