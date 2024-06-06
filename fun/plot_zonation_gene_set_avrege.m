function plot_zonation_gene_set_avrege_vall(vall,gns,zones,EXP_THRESH,clr)

addpath 'X:\Common\Lab_Papers\zonation_human_villus\Analysis\Code\functions'

if nargin<4
    clr='b';
end

if nargin<3
    EXP_THRESH=5*10^-6;
end
xtl = vall.zone_ttls(zones);

indin = find(ismember(vall.gene_name,gns));
mx=max(vall.md(:,zones),[],2);
indin=intersect(indin,find(mx>EXP_THRESH));
m =vall.md(indin,zones);

mnorm = m./sum(m,2);
vec = median(mnorm,1);
mad = median(abs(mnorm-median(mnorm)));

plot_patch(1:size(vec,2),vec./max(vec),mad./max(vec),clr);
set(gca,'XTick',1:size(vec,2),'xticklabels',xtl,'FontSize',8);
xtickangle(-45)
% xlim([1 max(xlim)]);
ylim([0 max(ylim)])
% ylabel({'Madien of all genes in pathway';'[Max norm]'});


