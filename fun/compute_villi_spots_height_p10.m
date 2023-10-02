function t = compute_villi_spots_height_p10(t,nBINS,NUM_PER_VILLI,current_dir)

addpath ([current_dir,'\fun\'])
if nargin <2
    nBINS=5;
end

if nargin <3
    NUM_PER_VILLI=0;
end

%% read spots annotations and muscualris serosa spots
tab = readtable([current_dir,'\Data\RAW\manual_annotation\p10_dist2submuc.csv']);
tab_serosa = readtable([current_dir,'\Data\RAW\manual_annotation\P10_Graph-Based.csv']);
tab_serosa.Var1 = replace(tab_serosa.Var1,',Cluster','');

t.dist2sm = cell(length(t.spot_name),1);
ind_serosa = [];
for j=1:length(t.spot_name)
    % annotation
    str=t.spot_name{j};
    str(findstr(str,'_'))='-';
    indd=find(strcmpi(tab.Barcode,str));
    if ~isempty(indd)
        t.dist2sm(j)=tab.dist2submuc(indd);
    else
        t.dist2sm(j)={'non_epi_out'};
    end

    %serosa 
    indd_s=find(strcmpi(tab_serosa.Var1,str));
    if ~isempty(indd_s)
        ind_serosa = [ind_serosa; j];
    end
end
%% assign spots to villi - loop over inner and outer

strvec = {'inn';'out'};
ind_submuc = cell(length(strvec),1);
ind_epi = cell(length(strvec),1);
D = cell(length(strvec),1);
d_closest_sub = cell(length(strvec),1);
closest_sub = cell(length(strvec),1);

f2 = figure;
for i=1:length(strvec)    
    str2cmp = unique(t.dist2sm(contains(t.dist2sm,strvec{i})));

    % compute the distances from each spot to the closest sub-mucosal spot
    ind_submuc{i}=find(strcmpi(t.dist2sm,str2cmp{2}));
    ind_epi{i}=find(strcmpi(t.dist2sm,str2cmp{1}));

    % remove from the analysis the muscularis serosa spots
    ind_submuc{i} = setdiff(ind_submuc{i},ind_serosa);

    D{i}=pdist2(t.coor(ind_epi{i},:),t.coor(ind_submuc{i},:));
    % for each spot find the minimal dist
    [d_closest_sub{i},closest_sub{i}]=min(D{i},[],2);
    closest_sub{i}=ind_submuc{i}(closest_sub{i});

    % also determine the submucosal spots closest to the epithel
    Dsub=pdist2(t.coor(ind_submuc{i},:),t.coor(ind_epi{i},:));
    [dist_epi,ind_epi_closest]=min(Dsub,[],2);
    zon_struct.ind_muscularis_mucosa{i}=ind_submuc{i}(find(dist_epi<50));
    zon_struct.ind_submucosa{i}=ind_submuc{i}(find(dist_epi>=50)); %median(dist_epi)

    subplot(2,3,1*i)
    scatter(t.coor(ind_epi{i},1),t.coor(ind_epi{i},2),25,closest_sub{i},'filled');
    c = colorbar;
    c.Label.String = 'Closest Submuc spot idx';
    title('2 - Assign spots to a villus');
    box on;
    set(gcf,'units','normalized','position',[0.3026    0.2630    0.2901    0.4463])
    hold on;
    scatter(t.coor(zon_struct.ind_muscularis_mucosa{i},1),t.coor(zon_struct.ind_muscularis_mucosa{i},2),'k.');
    scatter(t.coor(zon_struct.ind_submucosa{i},1),t.coor(zon_struct.ind_submucosa{i},2),'r.');
    axis square

    subplot(2,3,2*i);
    scatter(t.coor(ind_epi{i},1),t.coor(ind_epi{i},2),25,d_closest_sub{i},'filled');
    c = colorbar;
    c.Label.String = 'Distance from SM spot';
    title('3 - Compute villus hight per spot');
    box on;
    set(gcf,'units','normalized','position',[0.3026    0.2630    0.2901    0.4463])
    hold on;
    scatter(t.coor(ind_submuc{i},1),t.coor(ind_submuc{i},2),'w.');
    axis square

    subplot(2,3,3*i)
    imshow(t.slide_image);
    set(gca,'YDir','normal')

    % compute all pairwise distances in microns
    DD=squareform(pdist(t.coor));
    for j=1:length(DD)
        DD(j,j)=NaN;
    end
    vec=nanmin(DD,[],2);
    mn=min(vec(vec>0));
    DD=DD*100/mn;
    d_closest_sub{i}=d_closest_sub{i}*100/mn;
end

%% merge inner and outer

ind_submuc_all = [ind_submuc{1};ind_submuc{2}];
ind_epi_all = [ind_epi{1};ind_epi{2}];
d_closest_sub_all = [d_closest_sub{1};d_closest_sub{2}];
closest_sub_all = [closest_sub{1};closest_sub{2}];
zon_struct.ind_muscularis_mucosa = [zon_struct.ind_muscularis_mucosa{1};zon_struct.ind_muscularis_mucosa{2}];
zon_struct.ind_submucosa = [zon_struct.ind_submucosa{1};zon_struct.ind_submucosa{2}];

%% divide into villi
indices=closest_sub_all;
villi_unique=unique(indices);

% for each unique villus maintain the distance to base vectors
dist2base=cell(length(villi_unique),1);
villi_indices=cell(length(villi_unique),1);
for i=1:length(villi_unique)
    ind=find(closest_sub_all==villi_unique(i));
    dist2base{i}=d_closest_sub_all(ind);
    villi_indices{i}=ind;
end

% remove villi with less than NUM_PER_VILLI spots
clear L;
for i=1:length(villi_indices),
    L(i)=length(villi_indices{i});
end
ind2include=find(L>=NUM_PER_VILLI);
villi_unique=villi_unique(ind2include);
dist2base=dist2base(ind2include);
villi_indices=villi_indices(ind2include);

% show the remaining villi
figure;
scatter(t.coor(ind_epi_all,1),t.coor(ind_epi_all,2),25,closest_sub_all,'filled');
c = colorbar;
c.Label.String = 'Closest Submuc spot idx';
title('show the remaining villi');
box on;
set(gcf,'units','normalized','position',[0.3026    0.2630    0.2901    0.4463])
hold on;
scatter(t.coor(ind_submuc_all,1),t.coor(ind_submuc_all,2),'k.');

hold on;
for i=1:length(villi_indices)
   scatter(t.coor(villi_unique(i),1),t.coor(villi_unique(i),2),'go')
    scatter(t.coor(ind_epi_all(villi_indices{i}),1),t.coor(ind_epi_all(villi_indices{i}),2),'r.')
end


%% assign each spot the length of the villi it belongs to and save data
clear spot_struct;
spot_struct.ind=[];
spot_struct.dist2base=[];
spot_struct.vil_length=[];
for i=1:length(villi_indices)
    spot_struct.ind=[spot_struct.ind;ind_epi_all(villi_indices{i})];
    spot_struct.dist2base=[spot_struct.dist2base;d_closest_sub_all(villi_indices{i})];
    spot_struct.vil_length=[spot_struct.vil_length;repmat(max(d_closest_sub_all(villi_indices{i})),length(villi_indices{i}),1)];
end

VIL_LEN_BINS=[0 400;0 inf];
prctile_bin_vec=[0:(100/nBINS):100];

zon_struct.mat=cell(1,size(VIL_LEN_BINS,1));
spot_struct.prctile_bin_ind = zeros(length(spot_struct.ind),1);
for i=1:size(VIL_LEN_BINS,1)
    indin=find(spot_struct.vil_length>=VIL_LEN_BINS(i,1) & spot_struct.vil_length<=VIL_LEN_BINS(i,2));
    ind_spots=spot_struct.ind(indin);
    dist_bins=prctile(spot_struct.dist2base(indin),prctile_bin_vec);
    zon_struct.mat{i}=NaN(length(t.gene_name),length(dist_bins)-1);
    zon_struct.mat_se{i}=NaN(length(t.gene_name),length(dist_bins)-1);
    for j=1:length(dist_bins)-1
        ind2include=indin(find(spot_struct.dist2base(indin)>=dist_bins(j) & spot_struct.dist2base(indin)<=dist_bins(j+1)));
        if isempty(ind2include)
            continue;
        end
        ind2average=spot_struct.ind(ind2include);
        zon_struct.mat{i}(:,j)=mean(t.mat_norm(:,ind2average),2);
        zon_struct.mat_se{i}(:,j)=std(t.mat_norm(:,ind2average),[],2)/sqrt(length(ind2average));
        spot_struct.prctile_bin_ind(ind2include) = j;
    end
end

% compute p-values
indicator = spot_struct.prctile_bin_ind';
for i=1:size(t.mat_norm,1)
        vec = t.mat_norm(i,spot_struct.ind);
        zon_struct.pval(i) = kruskalwallis(vec,indicator,'off');
end

zon_struct.qval=mafdr(zon_struct.pval,'bhfdr','true');

if isfield(t,'dist_zon_struct')
    t = rmfield(t,'dist_zon_struct');
end

t.dist_zon_struct.spot_struct = spot_struct;
t.dist_zon_struct.spot_index = spot_struct.ind;
t.dist_zon_struct.spot_dist2base = spot_struct.dist2base;
t.dist_zon_struct.spot_vil_length = spot_struct.vil_length;
t.dist_zon_struct.prctile_bin_ind = spot_struct.prctile_bin_ind;
t.dist_zon_struct.mat = zon_struct.mat;
t.dist_zon_struct.mat_se = zon_struct.mat_se;
t.dist_zon_struct.zon_struct = zon_struct;
t.dist_zon_struct.pval = zon_struct.pval;
t.dist_zon_struct.qval = zon_struct.qval;
t.dist_zon_struct.zones_lables = {'Crypt'  'V1'  'V2'  'V3'  'V4'};
t.dist_zon_struct.closest_sub = closest_sub;

end
