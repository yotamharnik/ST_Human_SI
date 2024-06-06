function imamat=plot_zonation_carpet_vall_panel(vall,gns,zones,exp_filter,dynamic_range_filter,PLOT_FLAG,SORT_FLAG,NEXTTILE)

if nargin<8
    NEXTTILE=1;
end

if nargin<7
    SORT_FLAG=1;
end

if nargin<6
    PLOT_FLAG=1;
end

if nargin<5
    dynamic_range_filter=0;
end

if nargin<4
    exp_filter=0;
end
if nargin<3
    zones=1:9;
end

xtl = vall.zone_ttls(zones);
[ia,ib] = ismember(upper(gns),vall.gene_name);
indin = ib(ia);


% set the subplotting size
n = ceil(sqrt(length(indin)));


% compute the mean of all paitents
imamat = vall.md(indin,zones);
imasem = vall.mad(indin,zones);

if exp_filter>0
    indin2 = find(max(imamat,[],2)>exp_filter);
    dr=max(imamat,[],2)./min(imamat,[],2);
    indin2=intersect(indin2,find(dr>=dynamic_range_filter));
    imamat = imamat(indin2,:);
    imasem = imasem(indin2,:);
    indin = indin(indin2);
end

% sort by com for nice presentation
[~,mxind]=max(imamat,[],2);
mxind=(mxind-min(mxind))./(max(mxind)-min(mxind));
if SORT_FLAG
    com=calculate_com_mat(imamat);
    [~,ord]=sort(com+mxind);
else
    ord = 1:size(imamat,1);
end
if PLOT_FLAG
    if NEXTTILE
        nexttile;
    end
    imagesc(imamat(ord,:)./max(imamat(ord,:),[],2));
    set(gca,'xtick',1:length(xtl),'xticklabel',xtl);
    set(gca,'ytick',1:length(ord),'yticklabel',vall.gene_name(indin(ord)));
    set(gca,'fontsize',7);
    title('Median of all paitents');
    c = colorbar;
    c.Label.String = 'Max norm';
end
for i=1:size(vall.md(:,zones),2)+1
    xline(i-0.5,'-k','LineWidth',0.5);
end
for i=1:length(gns)+1
    yline(i-0.5,'-k','LineWidth',0.5);
end
