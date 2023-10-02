function category_struct = import_lineage_signature_matrix(folder_path,categoty_sig_csv,categoty_size_csv)


% load and prase category signature matrix
T = readtable([folder_path,categoty_sig_csv],'ReadVariableNames',1);
T2 = readtable([folder_path,categoty_size_csv]);
category_struct.gene_name = T.Var1;
category_struct.mat_sig = table2array(T(:,2:end));
category_struct.categories = T.Properties.VariableNames(2:end);
category_struct.categories_cell_count = T2.Var3;

% rearrange columns (cell type order)
col_ord = [2 3 1 7 6 4 5];
category_struct.mat_sig = category_struct.mat_sig(:,col_ord);
category_struct.categories = category_struct.categories(col_ord);
category_struct.categories_cell_count = category_struct.categories_cell_count(col_ord);

end
