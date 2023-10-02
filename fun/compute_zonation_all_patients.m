 function zon_struct_combined=compute_zonation_all_patients(v,EXP_THRESH)

if nargin<2
   EXP_THRESH=10^-5;
end

t=v{1}; %same gene_name for all patients
xtl = t.villus_bin_lbl;

mat_all=[];
mat_norm_all=[];
spots_all=[];
indicator_all_patient=[];
indicator_all_vill_bin=[];
indicator_all_eta_bin=[];
indicator_all_pxcl_max=[];
for i=1:length(v)
    % find relevant spots
    indin = find(~isnan(v{i}.villus_bin));
    indin_eta = find(~isnan(v{i}.apicome.eta_bin));

    % hold the relvant features
    spots = strcat([v{i}.patient,'_'],v{i}.spot_name(indin));
    mat            = v{i}.mat(:,indin);
    mat_norm       = v{i}.mat_norm(:,indin);
    indicator_patient = i*ones(1,length(indin));
    indicator_vill_bin = v{i}.villus_bin(indin)';
    indicator_all_eta_bin = nan(length(indin),1);
    indicator_all_eta_bin(indin_eta) = v{i}.apicome.eta_bin(indin_eta);

    mat_all     =[mat_all mat];
    mat_norm_all=[mat_norm_all mat_norm];
    spots_all = [spots_all spots];
    indicator_all_patient=[indicator_all_patient indicator_patient];
    indicator_all_vill_bin =[indicator_all_vill_bin indicator_vill_bin];
    indicator_all_eta_bin =[indicator_all_eta_bin' indicator_all_eta_bin'];
    indicator_all_pxcl_max = [ ]
end

indin = find(max(mat_norm_all,[],2)>=EXP_THRESH);

% compute the median amd MAD of all paitents
imamat = NaN(length(t.gene_name),length(xtl));
imasem = NaN(length(t.gene_name),length(xtl));
% also compute the mean and sem
imamat_mn = NaN(length(t.gene_name),length(xtl));
imasem_se = NaN(length(t.gene_name),length(xtl));


pval_spots = NaN(length(t.gene_name),1);
qval_spots = NaN(length(t.gene_name),1);
pval_combined = NaN(length(t.gene_name),1);
qval_combined = NaN(length(t.gene_name),1);
for i=1:length(indin)
    if mod(i,500)==0
        display([num2str(i) ' of ' num2str(length(indin))]);
    end
    vec = zeros(length(v),length(xtl));
    for j = 1:length(v)
        vec(j,:) = v{j}.zonal_mat(indin(i),:);
    end
    imamat(indin(i),:) = median(vec);
    imasem(indin(i),:) = median(abs(vec-median(vec)));

    imamat_mn(indin(i),:) = mean(vec);
    imasem_se(indin(i),:) = median(abs(vec-median(vec)));


    pval_spots(indin(i))=kruskalwallis(mat_norm_all(indin(i),~isnan(indicator_all_patient)),indicator_all_patient(~isnan(indicator_all_patient)),'off');
    pval_combined(indin(i))=pfast([v{1}.dist_zon_struct.pval(indin(i)),...
                                   v{2}.dist_zon_struct.pval(indin(i)),...
                                   v{3}.dist_zon_struct.pval(indin(i)),...
                                   v{4}.dist_zon_struct.pval(indin(i)),...
                                   v{5}.dist_zon_struct.pval(indin(i)),...
                                   v{6}.dist_zon_struct.pval(indin(i)),...
                                   v{7}.dist_zon_struct.pval(indin(i)),...
                                   v{8}.dist_zon_struct.pval(indin(i))]);
end
qval_spots=mafdr(pval_spots,'bhfdr','true');
qval_combined=mafdr(pval_combined,'bhfdr','true');


% sort by com for nice presentation
com=calculate_com_mat(imamat);
[~,ord]=sort(com);

zon_struct_combined.gene_name=v{1}.gene_name;
zon_struct_combined.zone_ttls=xtl;
zon_struct_combined.md=imamat;
zon_struct_combined.mad=imasem;
zon_struct_combined.com=com;
zon_struct_combined.pval_spots=pval_spots;
zon_struct_combined.qval_spots=qval_spots;
zon_struct_combined.pval_combined=pval_combined;
zon_struct_combined.qval_combined=qval_combined;
zon_struct_combined.mat_all=mat_all;
zon_struct_combined.mat__norm_all=mat_norm_all;
zon_struct_combined.spot_all=spots_all;
zon_struct_combined.indicator_all_patient=indicator_all_patient;
zon_struct_combined.indicator_all_vill_bin=indicator_all_vill_bin;
zon_struct_combined.indicator_all_eta_bin=indicator_all_eta_bin;
zon_struct_combined.EXP_THRESH=EXP_THRESH;
