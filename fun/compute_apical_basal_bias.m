function apicome =  compute_apical_basal_bias(t,BIN_TH,edges)

% this function (1) computes the apical-lamina_propria eta for all spots in a
% given zone in bin_th
% and (2) computes the zonation mean and sem for all spots binned to
% number of edges

if nargin <2
    BIN_TH = [2 6];
end

if nargin <3
    edges = 5;
end
%% compute eta

% eta = (apical+lumen/apical+lumen+basel+lp)
api = t.pxcl.classes_mat(:,1);
bas = t.pxcl.classes_mat(:,2);
lp  = t.pxcl.classes_mat(:,5);
lum = t.pxcl.classes_mat(:,4);
eta = (api+lum)./(api+lum+bas+lp);

indin = find(t.villus_bin>=BIN_TH(1) & t.villus_bin<=BIN_TH(2));
indout = setdiff(1:length(t.villus_bin),indin);

figure('Units','normalized','Position',[0 0 1 1]);

% Full image
SZ = 15;
nexttile;
imshow(t.slide_image); hold on;
scatter(t.coor(indin,1),t.coor(indin,2),SZ,eta(indin),'filled','MarkerFaceAlpha',.5,'MarkerEdgeColor','flat');
c = colorbar;;
c.Ticks = [min(eta(indin)) max(eta(indin))];
c.TickLabels = {'LP','Apical'};
title([t.patient,' - Villus bin >=',num2str(BIN_TH)]);

%% bin the eta and plot

Y = NaN(length(t.spot_name),1);
[Y(indin),TH] = discretize(eta(indin),edges);

% Blow up
xx = mean(t.coor(:,1));
yy = mean(t.coor(:,1));
ind1 = find(inpolygon(t.coor(:,1),t.coor(:,2),[xx-20 xx+30],[yy-20 yy+30]));

WIN = 100;
SZ = 1000;
nexttile;
imshow(t.slide_image); hold on;
scatter(t.coor(indin,1),t.coor(indin,2),SZ,Y(indin),'filled','MarkerFaceAlpha',.5,'MarkerEdgeColor','flat');
c = colorbar;;
c.Ticks = [min(eta(indin)) max(eta(indin))];
c.TickLabels = {'LP','Apical'};
axis([min(t.coor(ind1,1))-WIN max(t.coor(ind1,1))+WIN min(t.coor(ind1,2))-WIN max(t.coor(ind1,2))+WIN]);

colormap(plasma);



% plot the distrabution of etas per villus bin
figure;
bin2inc = unique(t.villus_bin(indin));
leg = [];
for i=1:length(bin2inc)
    histogram(eta(t.villus_bin==bin2inc(i)),50);
    hold on;
    ylabel('eta');
    leg = [leg {t.villus_bin_lbl{bin2inc(i)+3}}];
end
legend(leg);

xline(TH,'--r');
sgtitle('Apicome eta zonal dependent distrabution')

%% compute the apicome "zonation" table

apicome.gene_name = t.gene_name;
apicome.eta_spots = eta;
apicome.eta_bin = Y;
apicome.mn_eta_bin = zeros(length(t.gene_name),edges);
apicome.se_eta_bin = zeros(length(t.gene_name),edges);
apicome.BIN_TH = BIN_TH;

for i=1:length(apicome.gene_name)
        for j=1:edges
            apicome.mn_eta_bin(i,j) = mean(t.mat_norm(i,apicome.eta_bin==j));
            apicome.se_eta_bin(i,j) = std(t.mat_norm(i,apicome.eta_bin==j))/sqrt(sum(apicome.eta_bin==j));
        end
end

apicome.com = calculate_com_mat(apicome.mn_eta_bin);
end