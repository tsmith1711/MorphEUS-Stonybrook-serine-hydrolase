%% this works for combining any columns from a drug-drug or drug-moa table
% instead of pulling one column from multiple figures, just pull one column
% from each data set and make the figure all at once. better for saving
% space and makes things feel less weird 

% if it's the first run, generate 'firstrun' var so that a table can be
% created 
if ~exist('firstrun','var') 
    firstrun = true;
    
end 
% first thing to do is run create_connectivity_map
create_connectivity_map

%% Recreate figure S8 B
% just do this part the first time
%blank_table = var_to_add;

% comment out if not doing a drug table
%blank_table_drug = var_to_add_drug; 

%% fig 8 
% streamlined it past this but want to keep this and just have it be out of
% the way for archival purposess
while false 
% filename is bayews_file
valid_name = genvarname(bayes_file); 


%two_way_table_bi_moa_pct

% find the variable with 6h
%target_var = contains(two_way_table_bi_moa_pct.Properties.VariableNames,"6h");
%target_var = contains(two_way_table_bi_moa_pct.Properties.VariableNames,"mox");
%target_var = contains(two_way_table_bi_moa_pct.Properties.VariableNames,"beda");
%target_var = contains(two_way_table_bi_moa_pct.Properties.VariableNames,"rifamp");
target_var = contains(two_way_table_bi_moa_pct.Properties.VariableNames,"d6");


% grab that var
var_to_add = two_way_table_bi_moa_pct(:,target_var);

drug_name = var_to_add.Properties.VariableNames{1};
end 

%% now need to do this for fig 4 -- this has both MoA and drug so need 
% this is also the section for the stonybrook part 
% to also get symmetric_table 

search_for = "d702";

is_low_dose = contains(bayes_file,"025x");
is_high_dose = contains(bayes_file,"3x");
is_joint = contains(bayes_file,"joint");

if is_low_dose
    suffix = "_LD";
elseif is_high_dose
    suffix = "_HD";
elseif is_joint
    suffix = "_JP";
else
    suffix = "_NA";
end 

% filename is bayes_file
valid_name = genvarname(bayes_file); 
target_var = contains(two_way_table_bi_moa_pct.Properties.VariableNames,search_for);
target_var_drug = contains(symmetric_table.Properties.VariableNames,search_for);

% grab that var
var_to_add = two_way_table_bi_moa_pct(:,target_var);

% remove dummy final var that was used to keep the applied drug on the
% bottom (if it exists)
if any(contains(var_to_add.Properties.RowNames,"zzz"))
    var_to_add('zzz',:) = [];

end 

% need to rename with JP or HD or LD
drug_name = var_to_add.Properties.VariableNames{1};

var_to_add.Properties.VariableNames{drug_name} = char(strcat(drug_name,suffix));

% for drug, remove the row value from the table
symmetric_table(drug_name,:) = [];

var_to_add_drug = symmetric_table(:,target_var_drug);

var_to_add_drug.Properties.VariableNames{drug_name} = char(strcat(drug_name,suffix));


%% add it to the list 
if firstrun
    % if it's the first run, create the table
    
    % just do this part the first time
    blank_table = var_to_add;

    % comment out if not doing a drug table
    blank_table_drug = var_to_add_drug; 
    
    % future runs should go to the false sections
    firstrun = false;

    
else 
    % all subsequent runs want to add to the table 
    blank_table = horzcat(blank_table,var_to_add);

    % comment out if not doing a drug table
    blank_table_drug = horzcat(blank_table_drug,var_to_add_drug);
    
end 

%% save the workspsace
%save('workspace_for_fig_s8b.mat','blank_table')
%save('workspace_for_fig_4.mat','blank_table')
save('workspace_for_fig_S7.mat','blank_table')

%save('workspace_for_fig_4_drug.mat','blank_table_drug')
save('workspace_for_fig_S7.mat','blank_table_drug','blank_table')

%% when ready to put into plots 
% have a while false there just so that these don't auto-run


% new map -> parula plus some purple 
parula_plus = parula;

parula_plus = parula_plus(6:end,:);

% need to add purple color to beginning -> fade from first blue color to a
% purple color
first_blue = parula_plus(1,:);
purple_hue = [59/255 14/255 101/255];
vec = [0;100];
raw = [first_blue; purple_hue];
N = 14; 
blue_to_purple = interp1(vec,raw,linspace(100,0,N),'pchip');

parula_plus = [blue_to_purple; parula_plus];
while false 
make_connectivity_heatmap(blank_table,false);
xtickangle(90)
colormap(parula_plus)


%% drug connectivity heatmap

make_connectivity_heatmap(blank_table_drug,false);
xtickangle(90)
colormap(parula_plus)

end 