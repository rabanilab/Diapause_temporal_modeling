function plot_data_two(Xid,logE1,tm1,M0_1,M1_1,M2_1,DE_1,logE2,tm2,M0_2,M1_2,M2_2,DE_2,resdir,minT)

if (nargin < 15)
    minT = 5;
end

kmin1 = (tm1 < minT);
kmin2 = (tm2 < minT);

utm = min([tm1 tm2]):0.25:max([tm1 tm2]);
umin = (utm < minT);


% plot heatmaps: ALL data
[~,i] = cluster_sort([logE1 logE2],2);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
subplot(1,6,1);
x = logE1 - nanmean(logE1(:,kmin1),2);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
colormap(gene_colormap(1));
set(gca,'xtick',1:size(tm1,2),'xticklabel',tm1,'ytick',[]);
xtickangle(90);
set(gca,'color',[0.8 0.8 0.8]);
subplot(1,6,2);
x = logE2 - nanmean(logE2(:,kmin2),2);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
colormap(gene_colormap(1));
set(gca,'xtick',1:size(tm2,2),'xticklabel',tm2,'ytick',[]);
xtickangle(90);
set(gca,'color',[0.8 0.8 0.8]);
subplot(1,6,6);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
set(gca,'xtick',[],'ytick',1:100:size(Xid,1));
colorbar;
ylabel(sprintf('n=%d',size(Xid,1)));
saveas(h, [resdir '/expr2.all.jpg'],'jpg');


% plot heatmaps: ALL data, by model
M1(DE_1==0,:) = repmat(M0_1(DE_1==0),1,size(utm,2));
M1(DE_1==1,:) = logistic_eval(M1_1(DE_1==1,:),utm);
M1(DE_1==2,:) = sigmo_eval(M2_1(DE_1==2,:),utm);
M2(DE_2==0,:) = repmat(M0_2(DE_2==0),1,size(utm,2));
M2(DE_2==1,:) = logistic_eval(M1_2(DE_2==1,:),utm);
M2(DE_2==2,:) = sigmo_eval(M2_2(DE_2==2,:),utm);
[~,i] = cluster_sort([M1 M2],2);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
subplot(1,6,1);
x = M1 - nanmean(M1(:,umin),2);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
colormap(gene_colormap(1));
set(gca,'xtick',1:size(utm,2),'xticklabel',utm,'ytick',[]);
xtickangle(90);
set(gca,'color',[0.8 0.8 0.8]);
subplot(1,6,2);
x = M2 - nanmean(M2(:,umin),2);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
colormap(gene_colormap(1));
set(gca,'xtick',1:size(utm,2),'xticklabel',utm,'ytick',[]);
xtickangle(90);
set(gca,'color',[0.8 0.8 0.8]);
subplot(1,6,3);
x = logE1 - nanmean(logE1(:,kmin1),2);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
colormap(gene_colormap(1));
set(gca,'xtick',1:size(tm1,2),'xticklabel',tm1,'ytick',[]);
xtickangle(90);
set(gca,'color',[0.8 0.8 0.8]);
subplot(1,6,4);
x = logE2 - nanmean(logE2(:,kmin2),2);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
colormap(gene_colormap(1));
set(gca,'xtick',1:size(tm2,2),'xticklabel',tm2,'ytick',[]);
xtickangle(90);
set(gca,'color',[0.8 0.8 0.8]);
subplot(1,6,6);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
set(gca,'xtick',[],'ytick',1:100:size(Xid,1));
colorbar;
ylabel(sprintf('n=%d',size(Xid,1)));
saveas(h, [resdir '/expr2.all.bymodel.jpg'],'jpg');

