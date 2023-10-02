function api_struct_combined=compute_apicome_all_patients(v,EXP_THRESH)

if nargin<2
   EXP_THRESH=10^-5;
end

t=v{1}; %same gene_name for all patients
xtl = {'LP','','Basal','','Apical'};

mat_all=[];
spots_all=[];
indicator_all_patient=[];
indicator_all_eta_bin=[];
for i=1:length(v)
    % find relevant spots
    indin = find(~isnan(v{i}.apicome.eta_bin));

    % hold the relvant features
    spots = strcat([v{i}.patient,'_'],v{i}.spot_name(indin));
    mat       = v{i}.mat_norm(:,indin);
    indicator_patient = i*ones(1,length(indin));
    indicator_eta_bin = v{i}.apicome.eta_bin(indin)';

    mat_all=[mat_all mat];
    spots_all = [spots_all spots];
    indicator_all_patient=[indicator_all_patient indicator_patient];
    indicator_all_eta_bin =[indicator_all_eta_bin indicator_eta_bin];
end

indin = find(max(mat_all,[],2)>=EXP_THRESH);

% compute the mean of all paitents
imamat = NaN(length(t.gene_name),length(xtl));
imasem = NaN(length(t.gene_name),length(xtl));
pval_spots = NaN(length(t.gene_name),1);
qval_spots = NaN(length(t.gene_name),1);
for i=1:length(indin)
    if mod(i,500)==0
        display([num2str(i) ' of ' num2str(length(indin))]);
    end
    vec = zeros(length(v),length(xtl));
    for j = 1:length(v)
        vec(j,:) = v{j}.apicome.mn_eta_bin(indin(i),:);
    end
    imamat(indin(i),:) = mean(vec);
    imasem(indin(i),:) = std(vec)/sqrt(length(vec));

    pval_spots(indin(i))=kruskalwallis(mat_all(indin(i),~isnan(indicator_all_patient)),indicator_all_patient(~isnan(indicator_all_patient)),'off');
end
qval_spots=mafdr(pval_spots,'bhfdr','true');


% sort by com for nice presentation
com=calculate_com_mat(imamat);
[~,ord]=sort(com);

api_struct_combined.gene_name=v{1}.gene_name;
api_struct_combined.zone_ttls=xtl;
api_struct_combined.mn=imamat;
api_struct_combined.sem=imasem;
api_struct_combined.com=com;
api_struct_combined.pval_spots=pval_spots;
api_struct_combined.qval_spots=qval_spots;
api_struct_combined.mat_all=mat_all;
api_struct_combined.spot_all=spots_all;
api_struct_combined.indicator_all_patient=indicator_all_patient;
api_struct_combined.indicator_all_eta_bin=indicator_all_eta_bin;
api_struct_combined.EXP_THRESH=EXP_THRESH;
