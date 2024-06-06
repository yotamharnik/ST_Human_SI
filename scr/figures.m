%% add the path to all the necessary functions and defince working directory
current_dir = cd;
addpath([current_dir,'\fun\']);

% load processed data
load([current_dir,'\Data\processed\v_0323.mat']);
load([current_dir,'\Data\processed\vall_0323.mat']);

%% F1 A-E: plot the zonation carpet + profiles

% set the colors and thresholds
clrs = inferno(4);
CARPET_TH = 10^-5;
DR_THRESH = 0.5;
zones2plot = 3:9;

m =  vall.md(:,zones2plot);
se = vall.mad(:,zones2plot);
dr = (max(m,[],2)-min(m,[],2))./mean(m,2);
[mmax] = nanmax(m,[],2);

indin = find(mmax>CARPET_TH & dr > DR_THRESH);
m=m(indin,:)./mmax(indin,:);
se=se(indin,:)./mmax(indin,:);
[~,zone_max] = max(m,[],2);
zone_max =(zone_max -min(zone_max ))./(max(zone_max )-min(zone_max ));
[y,ord] = sort(calculate_com_mat(m)+zone_max);
GNS = vall.gene_name(indin(ord));

gns = {'REG1A','PIGR','DGAT1','SLC5A1','APOA4',... % enterocytes
       'DMBT1','GSTA1','RBP2','SLC2A5','AQP10',...     % enterocytes
       'TPH1','GCG','SCT','CCK','MLN',...          % EEC
       'SPINK4','MUC2','FCGBP','CLCA1','ZG16',...  % Goblet
       'TRPM5','SH2D6','PLCG2','GABRP','DEFB1'};   % tuft
ind_gns = find_indices_in_mat(GNS,gns);

% 5 rows by 6 columns
sp1 = sort([1:6:30]);

figure('Position',[488         230        1051         595]);
% carpet
subplot(5,6,sp1)
imagesc(m(ord,:));
set(gca,'XTick',1:size(m,2),'xticklabel',vall.zone_ttls(zones2plot));
[~,ord_ytk] = sort(ind_gns);
set(gca,'YTick',ind_gns(ord_ytk),'yticklabel',GNS(ind_gns(ord_ytk)));
qq = length(find(vall.pval_combined(indin)<0.05))/length(indin)
hold on;
xline([1:length(zones2plot)-1]+0.5);
title([num2str(length(indin)),' genes']);
colormap(viridis);
c = colorbar;
c.Position = [0.2370 0.9109 0.0114 0.0538];


% get the current tick labeks
ticklabels = get(gca,'YTickLabel');
% prepend a color for each tick label
ticklabels_new = cell(size(ticklabels));
for i = 1:length(ticklabels)
    if ord_ytk(i)<=10
        ticklabels_new{i} = ['\color[rgb]{',regexprep(num2str(clrs(1,:)),'\s+',','),'}',ticklabels{i}];
    elseif ord_ytk(i)>10 & ord_ytk(i)<=15
        ticklabels_new{i} = ['\color[rgb]{',regexprep(num2str(clrs(2,:)),'\s+',','),'}',ticklabels{i}];
    elseif ord_ytk(i)>15 & ord_ytk(i)<=18
        ticklabels_new{i} = ['\color[rgb]{',regexprep(num2str(clrs(3,:)),'\s+',','),'}',ticklabels{i}];
    elseif ord_ytk(i)>18
        ticklabels_new{i} = ['\color[rgb]{',regexprep(num2str(clrs(4,:)),'\s+',','),'}',ticklabels{i}];
    end
end
% set the tick labels
set(gca, 'YTickLabel', ticklabels_new,'TickLabelInterpreter','tex');


sp2 = sort([2:6:30, 3:6:30, 4:6:30 5:6:30 6:6:30]);



% individual proflies
ax = NaN
for i =1:length(gns)
    ax(i) = subplot(5,6,sp2(i));
    ind_gg = find(strcmpi(vall.gene_name,gns{i}));
    vec = vall.md(ind_gg,zones2plot);
    sem = vall.mad(ind_gg,zones2plot);
    if i<=10
        plot_patch(zones2plot,vec/mean(vec),sem/mean(vec),clrs(1,:),0,1);
    elseif i>10 & i<=15
        plot_patch(zones2plot,vec/mean(vec),sem/mean(vec),clrs(2,:),0,1);
    elseif i>15 & i<=20
        plot_patch(zones2plot,vec/mean(vec),sem/mean(vec),clrs(3,:),0,1);
    elseif i>20
        plot_patch(zones2plot,vec/mean(vec),sem/mean(vec),clrs(4,:),0,1);
    end
    ylim([0 max(ylim)]);
    xlim([min(zones2plot) max(zones2plot)]);
    if sp2(i)>5*5
        set(gca,'XTick',zones2plot,'XTickLabel',vall.zone_ttls(zones2plot),'FontSize',8,'XTickLabelRotation',45);
    else
        set(gca,'XTick',zones2plot,'XTickLabel',{},'FontSize',8);
    end
    box on;
    grid on;
    title(gns{i},'FontSize',12,'FontWeight','normal');
    %     plot_gene_zonation_all_patients(v,gns{i},'mean');
    if intersect(sp2(i),2:6:30)
        ylabel('Expression');
    end
end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.523, 0.96,'Enterocytyes','Color',clrs(1,:),'FontSize',12,'FontWeight','bold');
text(0.521, 0.608,'Enteroendocrine','Color',clrs(2,:),'FontSize',12,'FontWeight','bold');
text(0.552, 0.438,'Golbet','Color',clrs(3,:),'FontSize',12,'FontWeight','bold');
text(0.556, 0.266,'Tuft','Color',clrs(4,:),'FontSize',12,'FontWeight','bold');

% Export Source Data files:
% SD1 = array2table(m,'VariableNames',vall.zone_ttls(zones2plot));
% SD1.gene_name = GNS;
% writetable(SD1,[current_dir,'\SD\Source_Data_Fig.1.xlsx'],'FileType','spreadsheet');

%% panel G+H plot codex panels zonation proflies
NORM = 'max';
zone2plot = 3:9;
fig = figure;
tt = tiledlayout(1,5);
nexttile(1,[1 2]);
gns = {'CA9','NDRG1','SLC2A1'};

clm = {'b','g','r'};
l=[];
for i = 1:length(gns)
    gg = gns{i};
    ind = find(strcmpi(vall.gene_name,gg));
    if ~isempty(ind)
        if strcmpi(NORM,'mean')
            plot_patch(1:size(vall.md(:,zone2plot),2),vall.md(ind,zone2plot)./mean(vall.md(ind,zone2plot)),vall.mad(ind,zone2plot)./mean(vall.md(ind,zone2plot)),clm{i});
        elseif strcmpi(NORM,'max')
            plot_patch(1:size(vall.md(:,zone2plot),2),vall.md(ind,zone2plot)./max(vall.md(ind,zone2plot)),vall.mad(ind,zone2plot)./max(vall.md(ind,zone2plot)),clm{i});
        else
            plot_patch(1:size(vall.md(:,zone2plot),2),vall.md(ind,zone2plot),vall.mad(ind,zone2plot),clm{i});
        end
        ylim([0 max(ylim)]);
    end
    set(gca,'XTick',1:size(vall.md(:,zone2plot),2),'XTickLabel',vall.zone_ttls(zone2plot),'FontSize',10)
    box on;
    grid on;
    ylim([0 max(ylim)]);
    l = [l;{gg};{''}];
end
legend(l{:},'location','best');
ylabel(['Expression [',NORM,' Norm]'],'FontSize',12);'['
xlabel('Zone','FontSize',12);
title('mRNA zonation proflies');
set(gcf, 'position',1.0e+03 *[0.4266    0.4202    1.0088    0.3784]);

%% I – zonation of metabolic processes (carpet) + EDF10
% gene sets based on burclaff et al 2022

Zt = 'Fatty-acid transport and processing';
Z = {'CD36','FABP2','FABP1','APOA4','APOA1','PCK1','MTTP',...
    'GK','PCK2','APOB','PDK4','DGAT1','APOC3','ACSL5','LPGAT1',...
    'CREB3L3','MFSD2A','SCARB1','SLC27A4','DGAT2'};

At = 'Glucose transport';
A = {'SLC5A1','SLC2A5','SLC2A2','SLC5A9','SLC23A1','SLC2A1'};

Bt = 'Choelsterol transport';
B = {'NPC1L1','SCARB1','ABCG5','ABCG8'};

Ct = 'Amino-acid transport';
C = {'SLC15A1','SLC1A1','SLC36A1','SLC3A1','SLC6A19','SLC7A7','SLC7A8','SLC7A9'};

Dt = 'Ions-acid transport';
D = {'SLC26A3','SLC26A2','SLC4A4','SLC20A2','SLC9A3R1','SCNN1D','SCNN1A'};

Et = 'Digestive enzymes';
E = {'ANPEP','SI','LCT','MGAM','ENPEP','DPP4','PEPD','LIPA','TMPRSS15'};

Ft= 'Aquaporins';
F = {'AQP1','AQP3','AQP7','AQP10','AQP11'};

Gt = 'Xenobiotic metabolism';
G = {'CYP3A4','CYP3A5','CPY1A1','CYP2C9','CYP2C18','CBR1','ABCC2','ABCG2','FMO5'};

Ht = 'Heavy-metals metabolism';
H = [vall.gene_name(startsWith(vall.gene_name,'MT1'));{'MT2A';'STEAP2'}];
H = H([2:4,6:end]);

It = 'Butanoate metabolism';
I = {'SLC16A1','OXCT2','HMGCS2','OXCT1','ALDH7A1','ALDH5A1','ACAT2','ACSM3'};

Jt = 'Glutathione metabolism';
J = {'RRM2','ODC1','GSTA2','GSTA1','GSTA4','GGCT','GSTO2','GPX2','IDH2','GSS'};

Kt = 'Antimicrobial response';
K = {'DEFA5','DEFA6','REG3A','REG1A','LYZ','ITLN2','DMBT1','SLC12A2','GP2','RETNLB','GPX2','ITLN1','ITGB8','GSTA2',...
    'SNHG1','GAS5','AQP1'};

Lt = 'Sulfotransferases';
L = v{1}.gene_name(startsWith(v{1}.gene_name,'SULT'));

% %GOCC_MHC_CLASS_I_PROTEIN_COMPLEX
M = {'HFE','HLA-A','HLA-B','HLA-C','HLA-E','HLA-F','HLA-G','HLA-H','MR1','B2M'} 
Mt= 'MHC Class I';

% %GOCC_MHC_CLASS_II_PROTEIN_COMPLEX
N = {'HLA-DMA','HLA-DMB','HLA-DOA','HLA-DOB','HLA-DPA1','HLA-DPB1','HLA-DQA1',...
    'HLA-DQA2','HLA-DQB1','HLA-DQB2','HLA-DRA','HLA-DRB1','HLA-DRB3','HLA-DRB4','HLA-DRB5','B2M','CD74'};
Nt = 'MHC Class II';

EXP_TH = 10^-5;
DR_TH = 0.5;
SORT_FLAG = 1;

GNS = {Z,A,B,C,D,E,F,G,H,I,J,K,L,M,N};
TTLS = {Zt,At,Bt,Ct,Dt,Et,Ft,Gt,Ht,It,Jt,Kt,Lt,Mt,Nt};

zones = 3:9;
figure('Position',[0    0.0463    1.0000    0.8646]);
md=[];
for i =1:length(GNS)
    imamat=plot_zonation_carpet_vall_panel(vall,GNS{i},zones,EXP_TH,DR_TH,1,SORT_FLAG);
    title(TTLS{i});
    md=[md;median(imamat./max(imamat,[],2))];
    c = get( ancestor(gca, 'axes'), 'Colorbar');
    if i==1
        %         c.Position = [0.93 0.168 0.022 0.7];
        c.Position = [0.9262 0.1119 0.0381 0.1119];

    else
        c.Visible = 'off';
    end
end
set(gcf, 'units','normalized','outerposition',[0.2344    0.0741    0.4760    0.8287]);
colormap(viridis);
sgtitle('Zonation patterns of metabolic patterns');

md_com=calculate_com_mat(md(:,1:end));
[y,ord]=sort(md_com);
[~,mxind]=max(md,[],2);
[~,ord]=sort(mxind,'descend')
ord=ord(end:-1:1)


figure;
imagesc(md(ord,1:end)./max(md(ord,1:end),[],2));
set(gca,'ytick',1:length(ord),'yticklabel',TTLS(ord));
set(gca,'xtick',1:length(zones),'xticklabel',vall.zone_ttls(zones));

for i=1:size(md,2)-1
    xline(i+0.5,'-k');
end
for i=1:size(md,1)-1
    yline(i+0.5,'-k');
end
colormap(viridis);
c = colorbar;
c.Position = [0.9152 0.1097 0.0381 0.1006];
% sgtitle('Zonation patterns of metabolic patterns');
set(gcf,'Position',[680   531   454   447]);

% Export Source Data files:
SD1i = array2table(md(ord,1:end)./max(md(ord,1:end),[],2),'VariableNames',vall.zone_ttls(zones2plot));
SD1i.pathway_name = TTLS(ord)';
writetable(SD1i,[current_dir,'\SD\Source_Data_Fig.1i.xlsx'],'FileType','spreadsheet');


GNS = GNS(ord);
TTLS = TTLS(ord);

figure('Position',[0    0.0463    1.0000    0.8646]);
md=[];
for i =1:length(GNS)
    imamat=plot_zonation_carpet_vall_panel(vall,GNS{i},zones,EXP_TH,1,SORT_FLAG);
    title(TTLS{i});
    md=[md;median(imamat./max(imamat,[],2))];
    c = get( ancestor(gca, 'axes'), 'Colorbar');
    if i==1
        %         c.Position = [0.93 0.168 0.022 0.7];
        c.Position = [0.9262 0.1119 0.0381 0.1119];

    else
        c.Visible = 'off';
    end
end
set(gcf, 'units','normalized','outerposition',[0.2344    0.0380    0.3542    0.8833]);
colormap(viridis);
sgtitle('Zonation patterns of metabolic patterns');


%% F2 A+B : plot the CM biosyntatic pathway 

GNSA = {'BSCL2','TMEM159'};
GNSB = {'MOGAT3','MOGAT2','DGAT1'};
GNSC = {'CIDEB','MTTP','CD36','FABP2','DGAT2','APOB'};

zones = 3:9;
figure;
clrs = plasma(10);
plot_zonation_gene_set_avrege(vall,GNSA,zones,0,clrs(1,:))
plot_zonation_gene_set_avrege(vall,GNSB,zones,0,clrs(5,:))
plot_zonation_gene_set_avrege(vall,GNSC,zones,0,clrs(10,:))
set(gcf,'Position',[488.0000  460.2000  347.4000  301.8000]);
box on;
grid on;
ylabel('Expression','FontSize',12);

%
figure;
SDF2a = plot_zonation_carpet_vall_panel(vall,GNSA,3:9,0,1,1);
set(gca,'DataAspectRatio',[1 1 1]);
colormap(viridis);
set(gcf,'Position',[488.0000  460.2000  347.4000  301.8000]);

figure;
SDF2b = plot_zonation_carpet_vall_panel(vall,GNSB,3:9,0,1,1,0);
set(gca,'DataAspectRatio',[1 1 1]);
colormap(viridis);
set(gcf,'Position',[488.0000  460.2000  347.4000  301.8000]);

figure;
SDF2c = plot_zonation_carpet_vall_panel(vall,GNSC,3:9,0,1,1,0);
set(gca,'DataAspectRatio',[1 1 1]);
colormap(viridis);
set(gcf,'Position',[488.0000  460.2000  347.4000  301.8000]);

% % Export Source Data files:
% SDF2a = array2table(SDF2a,'VariableNames',vall.zone_ttls(zones));
% SDF2a.gene_name = GNSA';
% SDF2b = array2table(SDF2b,'VariableNames',vall.zone_ttls(zones));
% SDF2b.gene_name = GNSB';
% SDF2c = array2table(SDF2c,'VariableNames',vall.zone_ttls(zones));
% SDF2c.gene_name = GNSC';
% 
% SDF2 = [ SDF2a; SDF2b; SDF2c];
% writetable(SDF2,[current_dir,'\SD\Source_Data_Fig.2.xlsx'],'FileType','spreadsheet','Sheet','Fig2b');

%% panel C - dgat1 and dgat2
NORM = 'mean';
fig = figure;
GNS = {'DGAT1','DGAT2'};
clrs = plasma(length(GNS));

zone2plot = 3:9;
mat = [];
l=[];
for i = 1:length(GNS)
    gg = GNS{i};
    ind = find(strcmpi(vall.gene_name,gg));
    if ~isempty(ind)
        mat=[mat ;vall.md(ind,zone2plot)];
        if strcmpi(NORM,'mean')
            plot_patch(1:size(vall.md(:,zone2plot),2),vall.md(ind,zone2plot)./mean(vall.md(ind,zone2plot)),vall.mad(ind,zone2plot)./mean(vall.md(ind,zone2plot)),rgb2hex(clrs(i,:)));
        elseif strcmpi(NORM,'max')
            plot_patch(1:size(vall.md(:,zone2plot),2),vall.md(ind,zone2plot)./max(vall.md(ind,zone2plot)),vall.mad(ind,zone2plot)./max(vall.md(ind,zone2plot)),rgb2hex(clrs(i,:)));
        else
            plot_patch(1:size(vall.md(:,zone2plot),2),vall.md(ind,zone2plot),vall.mad(ind,zone2plot),rgb2hex(clrs(i,:)));
        end
        ylim([0 max(ylim)]);
    end
    set(gca,'XTick',1:size(vall.md(:,zone2plot),2),'XTickLabel',vall.zone_ttls(zone2plot),'FontSize',8)
    box on;
    grid on;
    axis tight
    ylim([0 max(ylim)]);
    l = [l;{gg};{''}];
end
legend(l{:},'location','northwest');
han=axes(fig,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Expression [Mean Norm]','FontSize',12);
xlabel(han,'Zone','FontSize',12);
title(han,'');

set(gcf,'Position',[ 680   633   426   345]);

%% Panels E + F - Iron uptake and release
NORM = 'MEAN';
zone2plot = 3:9;

A = {'CYBRD1','SLC11A2','TFRC'};... % luminal transport
    B = {'SLC40A1','HEPH','FTH1','FTL'};... %iron storage and secration
GNS = {A,B};
TTLS = {'Iron Uptake','Iron storage and release'};
clm = [{plasma(length(A))},{viridis(length(B))}];
zones = 3:9;

fig = figure;
tt = tiledlayout(2,1);
for j =1:length(GNS)
    nexttile;
    gns = GNS{j};
    l=[];
    for i = 1:length(gns)
        gg = gns{i};
        ind = find(strcmpi(vall.gene_name,gg));
        if ~isempty(ind)
            if strcmpi(NORM,'mean')
                plot_patch(1:size(vall.md(:,zone2plot),2),vall.md(ind,zone2plot)./mean(vall.md(ind,zone2plot)),vall.mad(ind,zone2plot)./mean(vall.md(ind,zone2plot)),rgb2hex(clm{j}(i,:)));
            elseif strcmpi(NORM,'max')
                plot_patch(1:size(vall.md(:,zone2plot),2),vall.md(ind,zone2plot)./max(vall.md(ind,zone2plot)),vall.mad(ind,zone2plot)./max(vall.md(ind,zone2plot)),rgb2hex(clm{j}(i,:)));
            else
                plot_patch(1:size(vall.md(:,zone2plot),2),vall.md(ind,zone2plot),vall.mad(ind,zone2plot),rgb2hex(clm{j}(i,:)));
            end
            ylim([0 max(ylim)]);
        end
        set(gca,'XTick',1:size(vall.md(:,zone2plot),2),'XTickLabel',vall.zone_ttls(zone2plot),'FontSize',8)
        box on;
        grid on;
        ylim([0 max(ylim)]);
        l = [l;{gg};{''}];
    end
    L = legend(l{:},'location','best');
    L.Box = 'off';
    title(TTLS{j});
end
han=axes(fig,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Expression [Mean Norm]','FontSize',12);
xlabel(han,'Zone','FontSize',12);
title(han,'');
set(gcf, 'position',[632   426   196   337]);

% % Export Source Data files:
% SDF2e = plot_zonation_carpet_vall_panel(vall,A,3:9,0,1,1,0);
% SDF2e = array2table(SDF2e,'VariableNames',vall.zone_ttls(zone2plot));
% SDF2e.gene_name = A';
% 
% SDF2f = plot_zonation_carpet_vall_panel(vall,B,3:9,0,1,1,0);
% SDF2f = array2table(SDF2f,'VariableNames',vall.zone_ttls(zone2plot));
% SDF2f.gene_name = B';
% 
% writetable(SDF2e,[current_dir,'\SD\Source_Data_Fig.2.xlsx'],'FileType','spreadsheet','Sheet','Fig2e');
% writetable(SDF2f,[current_dir,'\SD\Source_Data_Fig.2.xlsx'],'FileType','spreadsheet','Sheet','Fig2f');

