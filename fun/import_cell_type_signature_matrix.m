function category_struct = import_cell_type_signature_matrix(folder_path,categoty_sig_csv,categoty_size_csv)

% load and prase category signature matrix
T = readtable([folder_path,categoty_sig_csv],'ReadVariableNames',1);
T2 = readtable([folder_path,categoty_size_csv]);
category_struct.gene_name = T.Var1;
category_struct.mat_sig = table2array(T(:,2:end));
category_struct.categories = T.Properties.VariableNames(2:end);
category_struct.categories_cell_count = T2.Freq;

% rearrange columns (cell type order)
col_ord = 1:length(category_struct.categories);
category_struct.mat_sig = category_struct.mat_sig(:,col_ord);
category_struct.categories = category_struct.categories(col_ord);
category_struct.categories_cell_count = category_struct.categories_cell_count(col_ord);

% Refine categories resolution
% 1. rename endothelial cells
category_struct.categories(1:4) = {'Capillary EC_Endo','Lymphatic EC_Endo','Arterial EC_Endo','Venous EC_Endo'};
% 2. rename enterocytes
category_struct.categories(contains(category_struct.categories,'_Epi')) = ...
    {'Best4_Epi','Enterocytes [Early]_Epi','Enteroendocrine_Epi','Goblet_Epi','Enterocytes[ Interm]_Epi','StemCell_Epi',...
    'Enterocytes [Mature]_Epi','Paneth_Epi','TransAmplifying_Epi','Tuft_Epi'};
% 3. rename menschymal cells
category_struct.categories(contains(category_struct.categories,'_Msnc')) = ...
    {'Fibroblasts_Msnc','Pericytes_Msnc'};
% 4. reoder cells and lineages
catorder = [10,13,6,9,11,5,12,8,7,14,15:19,22:24,1:4,20:21,25];
category_struct.mat_sig = category_struct.mat_sig(:,catorder);
category_struct.categories = category_struct.categories(catorder);
category_struct.categories_cell_count = category_struct.categories_cell_count(catorder);     


end
