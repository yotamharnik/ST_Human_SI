%% Info
% a Master script to call all functions needed to align diffrent layer of
% data and produce a visium structer for downstream analysis.
% this script is dependent on the input files in "\Analysis\Data\RAW\"
% which are the relveant files from the spaceranger output folders

%% add the path to all the necessary functions and defince working directory
current_dir = cd;
addpath([current_dir,'\fun\']);


%% run the import_visium_data_yotam function for each patient
% define some fixed parameters for all samples (same as function defult):
Z_THRESH_MITO = 4;
Z_THRESH_SUM_UMI = -4;
HEMO_FRAC = inf;
BOUNDARY_PRCTILE = 100;
DIST_THRESH = inf;


% set the path for the files:
input_path = [current_dir,'\Data\RAW\'];

pati = {'P1','P5','P8','P9','P10','P11','P12','P13'};
pPub = sprintfc('P%d', 1:length(pati));
metadata = readtable([input_path,'metadata.csv']);


%% loop over all patients to import the data

v=cell(1,(length(pati)));
for i=1:length(v)
    v{i} = import_visium_data_funcion(input_path,pati{i},...
                                    Z_THRESH_MITO,Z_THRESH_SUM_UMI,HEMO_FRAC,...
                                    BOUNDARY_PRCTILE,DIST_THRESH);
    
    v{i}.patient = pati{i};
    v{i}.patient_publication = pPub{i};
    v{i}.metadata = metadata(i,:);
    v{i}.coor = v{i}.coor .* v{i}.json.tissue_hires_scalef;     % asjust coordianted to the slide image
end

disp(upper(['***** finished importing visium data for ',num2str(length(v)),' patient *****']));

%% import pixel classifier annotations

% set path to files
pxcl_dir = [current_dir,'\Data\RAW\Pixel_classifier_output.txt'];
pxcl_tab = readtable(pxcl_dir);
for i = 1:length(v)
    v{i} = parse_pixel_classifier(v{i},pxcl_tab);
end

disp(upper(['***** finished incorporating pixel classification for ',num2str(length(v)),' patient *****']));


%% Preform digital annotation along the submocusa - villus axis

QTH=0.25;
nBINS = 7;
VILLI_LEN_THRESH = 3;
for i=1:length(v)
    if strcmpi(v{i}.patient,'P10')
        v{i} = compute_villi_spots_height_p10(v{i},nBINS,VILLI_LEN_THRESH,current_dir);
        v{i}.dist_zon_struct.villus_base_spots = ~ismember(1:length(v{5}.spot_name),v{5}.dist_zon_struct.spot_index);
    else
        v{i} = compute_villi_spots_height(v{i}...
            ,[input_path,'manual_annotation\crypt_base_',v{i}.patient,'.csv'],...
            nBINS,VILLI_LEN_THRESH);
    end
    zontxt = ceil(length(find(v{i}.dist_zon_struct.qval<QTH))/length(v{i}.gene_name)*100);
    disp([v{i}.patient,': ',num2str(zontxt),'% of the genes are zonated [qval<',num2str(QTH),']',]);
    
    % add a villus_bin field 
    v{i}.villus_bin = nan(length(v{i}.spot_name),1);
    v{i}.villus_bin(v{i}.dist_zon_struct.spot_index) = v{i}.dist_zon_struct.prctile_bin_ind-1;
    v{i}.villus_bin(v{i}.dist_zon_struct.zon_struct.ind_muscularis_mucosa) = -1;
    v{i}.villus_bin(v{i}.dist_zon_struct.zon_struct.ind_submucosa) = -2;
end

% plot the digital annotation for all patient
close all hidden
vs = sprintfc('V%d', 1:nBINS-1);
EXPORT2R = 0;
bin_lbls = {'SubMucosa','MusMucosa','Crypt',vs{:}};
figure;
for i=1:length(v)
    v{i}.villus_bin_lbl = bin_lbls;
    t = v{i};

    nexttile;
    scatter(t.coor(:,1),t.coor(:,2),15,t.villus_bin+3,'filled','MarkerFaceAlpha',.75);
    c = colorbar;
    c.Ticks = unique(t.villus_bin(~isnan(t.villus_bin)))+3;
    c.TickLabels = t.villus_bin_lbl;
    colormap(tab20(length(t.villus_bin_lbl)));
    title(t.patient)

    % export for R
    if EXPORT2R
        tab  = table(t.spot_name',t.villus_bin);
        writetable(tab,[current_dir,'\Data\processed\',t.patient,'_villus_bin.csv']);
    end
end
sgtitle('Digital annotation');

%% Compute the apical-basel eta for each patient
BIN_TH = [0 6];
edges = 5;

for i=1:length(v)
    v{i}.apicome = compute_apical_basal_bias(v{i},BIN_TH,edges);
end

%% Compute the mean and sem by zones for each paitent

for i=1:length(v)
    U = unique(v{i}.villus_bin(~isnan(v{i}.villus_bin)));
    v{i}.zonal_mat = NaN(length(v{i}.gene_name),length(U));
    v{i}.zonal_sem = NaN(length(v{i}.gene_name),length(U));    
    for j=1:length(U)
        zon_mat = v{i}.mat_norm(:,v{i}.villus_bin==U(j));
        v{i}.zonal_mat(:,j) = mean(zon_mat,2);
        v{i}.zonal_sem(:,j) = std(zon_mat,[],2)./sqrt(size(zon_mat,2));
    end
end

%% Compute the zonation median over all patient
EXP_THRESH = 0;
vall= compute_zonation_all_patients(v,EXP_THRESH);
vall.apicome = compute_apicome_all_patients(v,0);

%% Import the signature matrix and reconstruct zonation for each cell-type in the signature matrix
processed_data_path = [current_dir,'\Data\processed\'];

cell_type_sig_csv = 'burc_teich_mean_exp_cell_types.csv';
cell_type_size_csv = 'burc_teich_size_anno.csv';
cell_type_struct = import_cell_type_signature_matrix(processed_data_path,cell_type_sig_csv,cell_type_size_csv);

lineage_sig_csv = 'burc_teich_integrated_mean_exp_lineages.csv';
lineage_size_csv = 'burc_teich_integrated_size_lineages.csv';
lineage_struct = import_lineage_signature_matrix(processed_data_path,lineage_sig_csv,lineage_size_csv);

for i = 1:length(v)
    v{i}.cell_type_struct = category_specific_reconstruction(v{i},cell_type_struct,2,10^-5);
    v{i}.lineage_struct = category_specific_reconstruction(v{i},lineage_struct,2,10^-5);
end

disp(upper(['***** finished zonation reconstruction for ',num2str(length(v)),' patient *****']));

% add the cell type and linage data to vall
vall.cell_type_struct = v{1}.cell_type_struct;
vall.lineage_struct = v{1}.lineage_struct;

%% save **uncoment to save - takes time and space :)
% save X:\Common\Lab_Papers\zonation_human_villus\Analysis\Data\processed\v_0323.mat v -v7.3
% save X:\Common\Lab_Papers\zonation_human_villus\Analysis\Data\processed\vall_0323.mat vall -v7.3

