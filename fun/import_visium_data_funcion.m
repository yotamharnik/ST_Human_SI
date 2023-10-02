function t=import_visium_data_funcion(input_path,patient,Z_THRESH_MITO,Z_THRESH_SUM_UMI,HEMO_FRAC,BOUNDARY_PRCTILE,DIST_THRESH)

if nargin<7
    DIST_THRESH=Inf;
end
if nargin<6
    BOUNDARY_PRCTILE=100; % the spot is removed if it is beyond this distance from the spot center x or y 
end
if nargin<5
    %HEMO_FRAC=0.01;
    HEMO_FRAC=inf ; % we don't want to remove spots with high fraction of hemoglobin genes since they often mark vessels
end
if nargin<4
    Z_THRESH_SUM_UMI=-4;
end

if nargin<3
    Z_THRESH_MITO=4;
end

% read spot position
display('Reading spot positions');
T=readtable([input_path,patient,'_tissue_positions_list.csv']);
barcode_space=table2cell(T(:,1));
coor=table2array(T(:,2:end));

% import the microscopy image
im=imread([input_path,patient,'_tissue_hires_image.png']);
display('Importing UMI counts');
T=readtable([input_path,patient,'_counts_UTT.csv']);

t.slide_image=im;
t.spot_name=T.Properties.VariableNames(2:end);
t.gene_name=table2cell(T(:,1));
t.coor=NaN*ones(length(t.spot_name),2);
t.mat=table2array(T(:,2:end));
t.mat_orig = t.mat;

% prefrom bg substruction 
bg_vector = compute_background_counts(input_path,patient,'mean');
if ~isempty(bg_vector)
    t.mat = t.mat-bg_vector;
    t.mat(find(t.mat<0)) = 0;
end
t.bg_vector = bg_vector;

t.mat_norm=t.mat./sum(t.mat);

% assign the coordinates of each spot according to the barcode_space
for i=1:length(t.spot_name)
    str=t.spot_name{i};
    str(findstr(str,'_'))='-';
    indd=find(strcmpi(barcode_space,str));
    if ~isempty(indd)
        t.coor(i,:)=coor(indd,end:-1:end-1);
    end
end


%% remove problematic spots
% 1. low UMI counts
s=log10(sum(t.mat));
nexttile; subplot(2,2,1:2); hist(s,100);
xlabel('log10(#UMIs)');
Z_log_sum_umi=(s-mean(s))/std(s);
indin=find(Z_log_sum_umi>Z_THRESH_SUM_UMI);

% 2. high mitochondrial content
ind=find(contains(lower(t.gene_name),'mt-'));
frac_mito=sum(t.mat_norm(ind,:));
Z_mito=(frac_mito-mean(frac_mito))./std(frac_mito);
indin=intersect(indin,find(Z_mito<Z_THRESH_MITO));

% 3. high hemoglobin expression (blood vessels)
ind=find(contains((t.gene_name),'Hba-'));
ind=union(ind,find(contains((t.gene_name),'Hbb-')));
frac_hemo=sum(t.mat_norm(ind,:));
indin=intersect(indin,find(frac_hemo<HEMO_FRAC));

% Further remove boundary spots
mn=min(t.coor)+(max(t.coor)-min(t.coor))/2;
dist2centerX=max(abs(t.coor(:,1)-mn(1)),[],2);
dist2centerY=max(abs(t.coor(:,2)-mn(2)),[],2);
ind_boundary_spots=union(find(dist2centerX>prctile(dist2centerX,BOUNDARY_PRCTILE)),find(dist2centerY>prctile(dist2centerY,BOUNDARY_PRCTILE)));
indin=setdiff(indin,ind_boundary_spots);

% remove spots selected to be excluded if specified in function argumaemts
% in cluster annotation
spots_to_exclude = read_spots2exclude([input_path,'\manual_annotation\spots2exclude_',lower(patient),'.csv']);
ind_exc = find(ismember(t.spot_name,spots_to_exclude'));
indin = setdiff(indin,ind_exc);  


t.spot_name=t.spot_name(indin);
t.coor=t.coor(indin,:);

t.mat_orig = t.mat_orig(:,indin);
t.mat=t.mat(:,indin);
t.mat_norm=t.mat./sum(t.mat);
t.json = readjson([input_path,patient,'_scalefactors_json.json']);