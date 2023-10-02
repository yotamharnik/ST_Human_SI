function t = compute_villi_spots_height(t,SM_spots_input,nBINS,VILLI_LEN_THRESH)

if nargin <3
    nBINS=5;
end
if nargin <4
    VILLI_LEN_THRESH=0;
end

%% add the mannualy SM spots annotation
SM_spots = zeros(length(t.spot_name),1);
SM_spots_barcodes = read_spots2exclude(SM_spots_input);
ind_SM = find(ismember(t.spot_name',SM_spots_barcodes));
SM_spots(ind_SM) = 1;

%% assign spots to villi

f1 = figure;
scatter(t.coor(:,1),t.coor(:,2),25,SM_spots,'filled');
c = colorbar;
c.Ticks = unique(SM_spots);
c.TickLabels = {'MUCOSA','SUBMUCOSA'};
title('1 - identify mocusa');
box on;
set(gcf,'units','normalized','position',[0.3026    0.2630    0.2901    0.4463])


% compute the distances from each spot to the closest sub-mucosal spot
ind_submuc=find(SM_spots);
ind_epi=setdiff(1:size(t.coor,1),ind_submuc);
D=pdist2(t.coor,t.coor(ind_submuc,:));

% for each spot find the minimal dist
[d_closest_sub,closest_sub]=min(D,[],2);
closest_sub=ind_submuc(closest_sub);

% also determine the submucosal spots closest to the epithel
Dsub=pdist2(t.coor(ind_submuc,:),t.coor(ind_epi,:));
[dist_epi,ind_epi_closest]=min(Dsub,[],2);
zon_struct.ind_muscularis_mucosa=ind_submuc(find(dist_epi<median(dist_epi)));
zon_struct.ind_submucosa=ind_submuc(find(dist_epi>=median(dist_epi)));


f2 = figure;
scatter(t.coor(:,1),t.coor(:,2),25,closest_sub,'filled');
c = colorbar;
c.Label.String = 'Closest Submuc spot idx';
title('2 - Assign spots to a villus');
box on;
set(gcf,'units','normalized','position',[0.3026    0.2630    0.2901    0.4463])
hold on;
scatter(t.coor(zon_struct.ind_muscularis_mucosa,1),t.coor(zon_struct.ind_muscularis_mucosa,2),'k.');
scatter(t.coor(zon_struct.ind_submucosa,1),t.coor(zon_struct.ind_submucosa,2),'r.');

f3 = figure;
scatter(t.coor(:,1),t.coor(:,2),25,d_closest_sub,'filled');
c = colorbar;
c.Label.String = 'Distance from SM spot';
title('3 - Compute villus hight per spot');
box on;
set(gcf,'units','normalized','position',[0.3026    0.2630    0.2901    0.4463])
hold on;
scatter(t.coor(ind_submuc,1),t.coor(ind_submuc,2),'w.');

f4 = figure;
scatter(t.coor(:,1),t.coor(:,2),25,t.pxcl.classes_mat(:,4),'filled');
c = colorbar;
c.Label.String = "SM pixels spot fraction";
title('4 - Show pixel calssifier fraction');
box on;
set(gcf,'units','normalized','position',[0.3026    0.2630    0.2901    0.4463])
hold on;
scatter(t.coor(ind_submuc,1),t.coor(ind_submuc,2),'k.');

% distFig('Pos','W','Rows',2,'Cols',2,'Extra','Ignore')

% compute all pairwise distances in microns
D=squareform(pdist(t.coor));
for i=1:length(D)
    D(i,i)=NaN;
end
vec=nanmin(D,[],2);
mn=min(vec(vec>0));
D=D*100/mn;
d_closest_sub=d_closest_sub*100/mn;

%% divide into villi
indin=ind_epi;
indices=closest_sub(indin);
villi_unique=unique(indices);

% postprocess - for each villi_unique entry find all spots that belong to
% it and impose that every spot there has another spot within 100um from it
include=zeros(length(villi_unique),1);
for i=1:length(villi_unique)
    ind=find(closest_sub==villi_unique(i));
    ind=intersect(ind,indin);
    d=D(ind,ind);
    mn=nanmin(d,[],2);
    if all(mn<=110)
        include(i)=1;
    end
end

% show the remaining villi
f2 = figure;
scatter(t.coor(:,1),t.coor(:,2),25,closest_sub,'filled');
c = colorbar;
c.Label.String = 'Closest Submuc spot idx';
title('show the remaining villi');
box on;
set(gcf,'units','normalized','position',[0.3026    0.2630    0.2901    0.4463])
hold on;
% scatter(t.coor(ind_submuc,1),t.coor(ind_submuc,2),'k.');
scatter(t.coor(ind_submuc,1),t.coor(ind_submuc,2),'ko');

ind2plot=[];
vil=villi_unique(find(include));
for i=1:length(vil),
    ind2plot=[ind2plot ;find(closest_sub==vil(i))];
end

hold on;scatter(t.coor(ind2plot,1),t.coor(ind2plot,2),'r.')

villi_unique=vil;

% for each unique villus maintain the distance to base vectors
ind_epithel=ind_epi;

dist2base=cell(length(villi_unique),1);
villi_indices=cell(length(villi_unique),1);
for i=1:length(villi_unique)
    ind=find(closest_sub==villi_unique(i));
    ind=intersect(ind,ind_epithel);
    dist2base{i}=d_closest_sub(ind);
    villi_indices{i}=ind;
end

% remove villi with less than NUM_PER_VILLI spots
clear L;
for i=1:length(villi_indices),
    L(i)=length(villi_indices{i});
end
ind2include=find(L>=VILLI_LEN_THRESH);
villi_unique=villi_unique(ind2include);
dist2base=dist2base(ind2include);
villi_indices=villi_indices(ind2include);

% assign each spot the length of the villi it belongs to
clear spot_struct;
spot_struct.ind=[];
spot_struct.dist2base=[];
spot_struct.vil_length=[];
for i=1:length(villi_indices)
    spot_struct.ind=[spot_struct.ind;villi_indices{i}];
    spot_struct.dist2base=[spot_struct.dist2base;d_closest_sub(villi_indices{i})];
    spot_struct.vil_length=[spot_struct.vil_length;repmat(max(d_closest_sub(villi_indices{i})),length(villi_indices{i}),1)];
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


t.dist_zon_struct.spot_struct = spot_struct;
t.dist_zon_struct.spot_index = spot_struct.ind;
t.dist_zon_struct.spot_dist2base = spot_struct.dist2base;
t.dist_zon_struct.spot_vil_length = spot_struct.vil_length;
t.dist_zon_struct.prctile_bin_ind = spot_struct.prctile_bin_ind;

t.dist_zon_struct.mat = zon_struct.mat;
t.dist_zon_struct.mat_se = zon_struct.mat_se;
t.dist_zon_struct.pval = zon_struct.pval;
t.dist_zon_struct.qval = zon_struct.qval;
t.dist_zon_struct.zon_struct = zon_struct;
t.dist_zon_struct.zones_lables = {'Crypt'  'V1'  'V2'  'V3'  'V4'};
t.dist_zon_struct.closest_sub = closest_sub;

end
