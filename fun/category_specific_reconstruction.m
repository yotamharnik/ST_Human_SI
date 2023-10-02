function categoires_singature = category_specific_reconstruction(t,category_mat,THratio,THexp)
% a master function, that calls two functions and loops over all categories
%
% fun1:   find_category_specific_genes
% fun2:   reconstruct_zonation_intestine
%
% input:  is a categoires_singature structer including cellanneal results and the category_mat
%         (same that was used for ca). THratio and THexp are optional to deterimane
%         how stringet to ve when selscting specific genes.
% output: is a zo struncter with cell for each category. the cell is a
%         structer that holds all the zonation data for that category.


%% determine the subset of category specific genes
spcf_genes_names = cell(1,length(category_mat.categories));
spcf_genes_inds = cell(1,length(category_mat.categories));
top5 = cell(1,length(category_mat.categories));

% find specific genes for each category
if nargin < 3
    THratio = 1;
    THexp = 5*10^-6;
end

cel_type_signature_gene_indices=cell(1,length(category_mat.categories));
for i = 1:length(category_mat.categories)
    [spcf_genes_names{i},spcf_genes_inds{i},top5{i}] = ...
        find_category_specific_genes(category_mat.categories(i),...
        category_mat.categories,category_mat.gene_name,category_mat.mat_sig,THratio,THexp);
    cel_type_signature_gene_indices{i}=[];
    for j=1:length(spcf_genes_names{i})
        indd=find(strcmpi(t.gene_name,spcf_genes_names{i}{j}));
        if ~isempty(indd)
            cel_type_signature_gene_indices{i}=[cel_type_signature_gene_indices{i} indd'];
        end
    end
end

% save the category_mat 
categoires_singature = category_mat;
categoires_singature.cel_type_signature_gene_indices = cel_type_signature_gene_indices;

% display results summary 
for i=1:length(category_mat.categories)
    disptab_coln{i} = [category_mat.categories{i},'_TOP_5_genes'];
    disptab_spcn{i} = [num2str(length(spcf_genes_inds{i})),' specific genes total'];
end

end