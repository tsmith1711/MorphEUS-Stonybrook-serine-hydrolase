%% knn helper script 
% designed to be used by a bigger knn script. cuts down on repetitive code
% need to set up knn_data and knn_drugs ahead of time

%User settings 

% distance metric to usee
d_metric = 'cosine';

% boolean whether to create graph
create_graph = true;

% boolean whether or not to create confusion matrix
create_cm = true;

% which settings to use -- more general settings or more specific settings
large_group = true;
% true: use only the large 5 groups
% false: break into smaller groups 

%% Set up workspace 
% set up colormap
vec = [100; 83.3; 66.6; 50; 33.3; 16.6; 0;];
hex = ['#006687'; '#00b0e9';'#87e2ff';'#feffff';'#ff7575';'#b00000';'#812020'];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = 128;
map = interp1(vec,raw,linspace(100,0,N),'pchip');

% categories for all drugs 
proteins = ["AZT","Amk","Cam","Cla","Dox","Gent","Kan","Lin","Strep","Tet","Tig",...
    "Lin_025x","Lin_05x","Lin_075x","Lin_1x","Lin_3x"];
lipids = ["Tri","Nis","Nig","Mon","Gra","Cer","CCCP","THL"];
dnas = ["Dau","Lev","MIT","Mox","Nal","Olf","Sulf","Nit",...
    "MOX_025x","MOX_05x","MOX_075x","MOX_1x","MOX_3x"];
controls = ["water","Untreated","NaOH","MeOH","EtOH","DMSO"];
rnaps = ["RIF","RifT","RifT_025x","RifT_05x","RifT_075x","RifT_1x","RifT_3x",...
    "RIF_025x","RIF_05x","RIF_075x","RIF_1x","RIF_3x",];
cell_walls = ["A22","Amp","Carb","Ctax","Cyc","Del","EMB","ETA",...
    "IMI","INH","Mer","Oxa","Pip","Pre","Van",...
    "EMB_025x","EMB_05x","EMB_075x","EMB_1x","EMB_3x",...
    "INH_025x","INH_05x","INH_075x","INH_1x","INH_3x",...
    "Pre_025x","Pre_05x","Pre_075x","Pre_1x","Pre_3x"];
%others = ["PZA", "Ver", "Thi", "Clz"];

% smaller settings
peptidoglycan = ["Mer","Amp","Ctax","IMI","Pip","Oxa","Carb","Van","Cyc","A22"];
mycolic_acid = ["INH","ETA","Del","EMB","Pre",...
    "EMB_025x","EMB_05x","EMB_075x","EMB_1x","EMB_3x",...
    "INH_025x","INH_05x","INH_075x","INH_1x","INH_3x",...
    "Pre_025x","Pre_05x","Pre_075x","Pre_1x","Pre_3x"];
s50_subunit = ["AZT","Cla","Cam","Lin",...
    "Lin_025x","Lin_05x","Lin_075x","Lin_1x","Lin_3x"];
s30_subunit = ["Gent","Kan","Amk","Strep","Tet","Dox","Tig"];
fas_inhibitor = ["Cer","THL","Tri"];
atp_synthesis = ["BDQ","CCCP","Clz","Gra","Mon","Nig","Nis","PZA",...
    "BDQ_1x","BDQ_025x","BDQ_05x","BDQ_075x","BDQ_3x",...
    "Clz_1x","Clz_025x","Clz_05x","Clz_075x","Clz_3x",...
    "PZA_1x","PZA_025x","PZA_05x","PZA_075x","PZA_3x"];

efflux = ["Ver","Thi"];

% colors for different groups 
proteins_col = [0 200 0]/255;
lipids_col = [2 0 200]/255;
dnas_col = [220 0 0]/255;
controls_col = [0 0 0]/255;
rnaps_col = [240 240 0]/255;
cell_walls_col = [255 51 255]/255;
other = [0 153 153]/255;

% smaller settings color
peptidoglycan_col = [127 26 127]/255;
mycolic_acid_col = [255 51 255]/255;
s50_subunit_col = [0 200 0]/255;
s30_subunit_col = [0 110 0]/255;
fas_inhibitor_col = [0 179 179]/255;
atp_synthesis_col = [0 255 255]/255;
efflux_col = [0 0 255]/255; 

%% do knn 
% will rewrite this with the neighbors 
nearest_neighbors = knn_drugs;

% preallocate space to store distances
knn_distances = zeros(length(knn_drugs),1);

for i = 1:size(knn_data,1)
    % get individual drug
    individual_row = knn_data(i,:);

    % get remainder
    remaining_rows = knn_data;
    remaining_rows(i,:) = [];

    % get the drugs for 'remaining_rows'
    remaining_drugs = knn_drugs;
    remaining_drugs(i,:) = [];

    % do knn
    [idx,distances] = knnsearch(remaining_rows,individual_row,'Distance',d_metric);
    knn_distances(i) = distances; 
    nearest_neighbors(i) = remaining_drugs(idx);

end 

% store data in results table 
results_table = table(knn_drugs,nearest_neighbors,knn_distances,'VariableNames',{'DRUG','NEIGHBOR','DISTANCE'});


%% Create confusion matrix 
% will overwrite with information in the loop
drug_moas = results_table.DRUG;
neighbor_moas = results_table.DRUG;

% go through and find moa of each drug and neighbor in order to create
% confusion matrix (don't neeed MoA of others and unknowns)
for k = 1:height(results_table)
    drug1 = results_table.DRUG(k);
    drug1 = drug1{1};
    drug2 = results_table.NEIGHBOR(k);
    drug2 = drug2{1};

    % one setting for large group 
    if large_group
        if ismember(drug1,dnas)
            drug_moas(k) = {'dna'};
        % protein 
        elseif ismember(drug1,proteins)
            drug_moas(k) = {'protein'};
        % lipid
        elseif ismember(drug1,lipids)
            drug_moas(k) = {'lipid'};
        % cell wall
        elseif ismember(drug1,cell_walls)
            drug_moas(k) = {'cell_wall'};
        % rnap
        elseif ismember(drug1,rnaps)
            drug_moas(k) = {'rnap'};
        % control
        elseif ismember(drug1,controls)
            drug_moas(k) = {'control'};
        % else if not one with a "known" moa 
        else
            drug_moas(k) = {'other'};
        end

        % neighbors 
        % dna 
        if ismember(drug2,dnas)
            neighbor_moas(k) = {'dna'};
        % protein 
        elseif ismember(drug2,proteins)
            neighbor_moas(k) = {'protein'};
        % lipid
        elseif ismember(drug2,lipids)
            neighbor_moas(k) = {'lipid'};
        % cell wall
        elseif ismember(drug2,cell_walls)
            neighbor_moas(k) = {'cell_wall'};
        % rnap
        elseif ismember(drug2,rnaps)
            neighbor_moas(k) = {'rnap'};
        % control
        elseif ismember(drug2,controls)
            neighbor_moas(k) = {'control'};
        % else if not one with a "known" moa 
        else
            neighbor_moas(k) = {'other'};
        end
    else 
        
        % one setting for small group
        if ismember(drug1,peptidoglycan)
            drug_moas(k) = {'peptidoglycan'};
        elseif ismember(drug1,mycolic_acid)
            drug_moas(k) = {'mycolic_acid'};
        elseif ismember(drug1,dnas)
            drug_moas(k) = {'dna'};
        elseif ismember(drug1,s50_subunit)
            drug_moas(k) = {'s50_subunit'};
        elseif ismember(drug1,s30_subunit)
            drug_moas(k) = {'s30_subunit'};
        elseif ismember(drug1,controls)
            drug_moas(k) = {'control'};
        elseif ismember(drug1,atp_synthesis)
            drug_moas(k) = {'atp_synthesis'};
        elseif ismember(drug1,fas_inhibitor)
            drug_moas(k) = {'fas_inhibitor'};
        elseif ismember(drug1,rnaps)
            drug_moas(k) = {'rnaps'};
        elseif ismember(drug1,efflux)
            drug_moas(k) = {'efflux'};
        % else if not one with a "known" moa 
        else
            drug_moas(k) = {'other'};
        end
        
        % neighbors 
        if ismember(drug2,peptidoglycan)
            neighbor_moas(k) = {'peptidoglycan'};
        elseif ismember(drug2,mycolic_acid)
            neighbor_moas(k) = {'mycolic_acid'};
        elseif ismember(drug2,dnas)
            neighbor_moas(k) = {'dna'};
        elseif ismember(drug2,s50_subunit)
            neighbor_moas(k) = {'s50_subunit'};
        elseif ismember(drug2,s30_subunit)
            neighbor_moas(k) = {'s30_subunit'};
        elseif ismember(drug2,controls)
            neighbor_moas(k) = {'control'};
        elseif ismember(drug2,atp_synthesis)
            neighbor_moas(k) = {'atp_synthesis'};
        elseif ismember(drug2,fas_inhibitor)
            neighbor_moas(k) = {'fas_inhibitor'};
        elseif ismember(drug2,rnaps)
            neighbor_moas(k) = {'rnaps'};
        elseif ismember(drug2,efflux)
            neighbor_moas(k) = {'efflux'};
        % else if not one with a "known" moa 
        else
            neighbor_moas(k) = {'other'};
        end
        
         
    end

end

% list of classes
if large_group
    class_list = {'cell_wall','dna','protein','lipid','rnap','control','other'};
else
    class_list = {'peptidoglycan','mycolic_acid','dna','s50_subunit','s30_subunit','atp_synthesis',...
        'fas_inhibitor','rnaps','efflux','control'};
end 

% if create confusion matrix
if create_cm
    figure
    % generate confusion matrix with row summaries 
    cm = confusionchart(drug_moas,neighbor_moas , ...
    'RowSummary','row-normalized');
    % organize so that confusion matrix is always in the same ordere 
    try
        sortClasses(cm,class_list)
    catch
        disp("can't sort labels")
    end
    cm_tit = strcat("confusion matrix -- ", d_metric, " ", num2str(numvar), " variables ",tit_add);
    title(cm_tit)
    
    % get success percent
    % get confusion matrixfrom the confusion chart object
    cmat = cm.NormalizedValues;
    % find the sum of the rows to get true classes
    true_class = sum(cmat,2);
    % find the sum of the colimns to get predicted class
    predicted_class = sum(cmat,1);
    if large_group
        % get all of the correct 'other' category
        correct_others = cmat(end,end);
        % add up along diagonal to get all correct
        all_correct = trace(cmat);
        % good correct does not includee correct from others category
        good_correct = all_correct - correct_others;
        % get all of the drugs not including the 'other' category
        all_with_no_others = sum(true_class(1:end-1));

        pct_correct = round((good_correct/all_with_no_others)*100,2); 
        disp(strcat(num2str(pct_correct), "% for ", num2str(numvar), " variables"))
    else
        % no other category
        
        % add up along diagonal to get all correct
        all_correct = trace(cmat);

        all_drugs = sum(true_class);
        
        pct_correct = round((all_correct/all_drugs)*100,2);
        disp(strcat(num2str(pct_correct), "% for ", num2str(numvar), " variables"))
    end
end 

%% Create graph

if create_graph 
    knn_digraph =  digraph(results_table.DRUG,results_table.NEIGHBOR,results_table.DISTANCE);

    % get colors for each node 
    node_colors = [];
    for row = 1:height(knn_digraph.Nodes)
        name = table2array(knn_digraph.Nodes(row,1));
        % use the correct color array to change the node
        
        if large_group
            if ismember(name,proteins) 
                node_colors = [node_colors; proteins_col];

            elseif ismember(name, lipids)           
                node_colors = [node_colors; lipids_col];

            elseif ismember(name, dnas)
                node_colors = [node_colors; dnas_col];

            elseif ismember(name, controls)
                node_colors = [node_colors; controls_col];

            elseif ismember(name, rnaps)
                node_colors = [node_colors; rnaps_col];

            elseif ismember(name, cell_walls)
                node_colors = [node_colors; cell_walls_col];

            else 
                node_colors = [node_colors; other];
            end
        else 
            if ismember(name,peptidoglycan) 
                node_colors = [node_colors; peptidoglycan_col];

            elseif ismember(name, mycolic_acid)           
                node_colors = [node_colors; mycolic_acid_col];

            elseif ismember(name, dnas)
                node_colors = [node_colors; dnas_col];

            elseif ismember(name, controls)
                node_colors = [node_colors; controls_col];

            elseif ismember(name, s50_subunit)
                node_colors = [node_colors; s50_subunit_col];

            elseif ismember(name, s30_subunit)
                node_colors = [node_colors; s30_subunit_col];
                
            elseif ismember(name, fas_inhibitor)
                node_colors = [node_colors; fas_inhibitor_col];

            elseif ismember(name, atp_synthesis)
                node_colors = [node_colors; atp_synthesis_col];

            elseif ismember(name, efflux)
                node_colors = [node_colors; efflux_col];
                
            elseif ismember(name, rnaps)
                node_colors = [node_colors; rnaps_col];
                
            else 
                node_colors = [node_colors; other];
            end
        end 
    end

    % create figure
    h2= figure;
  %  subplot(1,2,2)
   knn_graph = plot(knn_digraph,'Layout','force','LineWidth',2.5,'NodeColor',node_colors);

    % color based on weight
    knn_graph.EdgeCData= knn_digraph.Edges.Weight;

    % attempt to change size based on weight
   % layout(knn_graph,'force','WeightEffect','direct')

    % change background color
    set(gca,'Color',[.8 .8 .8])
    if strcmp(knn_type,'median')
        if knn_with_apply
            tit = strcat("KNN after PCA - medians ", addedDrugStr," applied ",d_metric);
        else  
            % if we made this one by randomly swapping labels, want to show
            % that in the title
            try
                tit = strcat("KNN after PCA - medians ",d_metric, " ", num2str(numvar)," variables ",tit_add);
            catch
                tit = strcat("KNN after PCA - medians ",d_metric, " ", num2str(numvar)," variables ");
            end 
        end
    else
        if knn_with_apply
            tit = strcat("KNN after PCA ", addedDrugStr," applied ",d_metric);
        else
            tit = strcat("KNN after PCA ",d_metric);
        end
    end
    
    title(tit,'Interpreter','none')
    % add colorbar
    colorbar
    colormap(map)
    % resize image
    set(h2, 'Position', [100, 100, 800, 700])
end 