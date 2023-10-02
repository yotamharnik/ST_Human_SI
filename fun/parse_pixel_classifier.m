function t = parse_pixel_classifier(t,pxcl_tab)
% this function get a visium struct and incorperates pixel classifier
% results into the structer


%% read classifier results

class_names = {'Apical','Basel','Musmucosa','Lumen','Stroma','Submucosa'};
t.pxcl.class_names = class_names;

% find revlevant rows
ind_spots = find(startsWith(pxcl_tab.barcodes,[t.patient,'_']));
tab = pxcl_tab(ind_spots,:);
tab.barcodes = replace(tab.barcodes,[t.patient,'_'],'');
tab.barcodes = replace(tab.barcodes,'-','_');

indin = zeros(length(t.spot_name),1);
for i=1:length(t.spot_name)
    ind = find(strcmpi(tab.barcodes,t.spot_name{i}));
    if ~isempty(ind)
        indin(i) = ind;
    end
end

mat = table2array(tab(:,1:length(class_names)));
t.pxcl.classes_mat = NaN(length(t.spot_name),6);
t.pxcl.classes_mat(find(indin~=0),:) = mat(indin(indin~=0),:);

% make some QC plots

figure('Units','normalized','Position',[0    0.05    1.0000    0.6]);
for i=1:length(class_names)
    subplot(2,6,i);
    scatter(t.coor(:,1),t.coor(:,2),15,t.pxcl.classes_mat(:,i),'filled');
    title(t.pxcl.class_names{i});
    subplot(2,6,i+6);
    hist(t.pxcl.classes_mat(:,i),100);
    title(t.pxcl.class_names{i});
    xlabel('Mean pixel probability');
end

[mx,mxind]=max(t.pxcl.classes_mat,[],2);
figure('Units','normalized','Position',[0.3438    0.6398    0.3620    0.2593]);
subplot(1,2,1);
hist(mx,100);
xlabel('max probability');
ylabel('# spots');
subplot(1,2,2);
hist(mxind,1:length(t.pxcl.class_names));
set(gca,'xtick',1:length(t.pxcl.class_names),'xticklabel',t.pxcl.class_names,'XTickLabelRotation',45);
xlim([.5 length(t.pxcl.class_names)+.5]);
ylabel('# spots');

end