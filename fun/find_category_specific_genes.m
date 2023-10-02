function [spcf_names,spcf_inds,top5] = find_category_specific_genes(category,categories,genes,mat,THratio,THexp,PLOT_FLAG,DISP_FLAG)
% this function gets a catergory, a signature matrix mat and gene names,
% and gives back the names and indices of genes which are category specific
% 'genes' are expected to be rownames of 'mat'
% 'categrories' are expected to be colnames of 'mat'
% category is the category to find specific genes for

% eg. run 'import_signature_matrix) and then
% find_category_specific_genes(category_mat.categories(3),...
% category_mat.categories,category_mat.gene_name,category_mat.mat_sig,5,10^-4,1);

mat=mat./sum(mat);
if nargin < 7
    PLOT_FLAG = 0;
end

if nargin < 8
    DISP_FLAG = 0;
end
CAT_FLAG = find(strcmpi(categories,category));

m1 = mat(:,CAT_FLAG);
m2 = max(mat(:,setdiff(1:length(categories),CAT_FLAG)),[],2);
spcf_inds = find(m1>THexp & m1>THratio*m2);
spcf_names = genes(spcf_inds);
[~,ord] = sort(mat(spcf_inds,CAT_FLAG),'descend');
if length(ord)<5
    top5 = genes(spcf_inds(ord));
else
    top5 = genes(spcf_inds(ord(1:5)));
end

if PLOT_FLAG
    figure;
    loglog(m1,m2,'.')
    hold on;
    loglog(m1(spcf_inds),m2(spcf_inds),'r.');
    xlabel(categories(CAT_FLAG));
    ylabel(['All beside ',categories(CAT_FLAG)]);
    legend('all genes','Category specific');
end

if DISP_FLAG
    disp(['Found ',num2str(length(spcf_inds)),' specific genes for the ',cell2mat(category),' category.']);
    disp(cell2table(top5))
end
end