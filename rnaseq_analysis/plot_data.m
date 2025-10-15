function plot_data(Xid,logE,tm,M0,M1,M2,DE,Err,Rsq,resdir,minT)

if (nargin < 11)
    minT = 5;
end

% heatmaps
utm = min(tm):0.25:max(tm);
umin = (utm < minT);
kmin = (tm < minT);

% plot heatmaps: ALL data
[~,i] = cluster_sort(logE,2);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
subplot(1,6,1);
x = logE - nanmean(logE(:,kmin),2);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
colormap(gene_colormap(1));
set(gca,'xtick',1:size(tm,2),'xticklabel',tm,'ytick',[]);
xtickangle(90);
set(gca,'color',[0.8 0.8 0.8]);
subplot(1,6,6);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
set(gca,'xtick',[],'ytick',1:100:size(Xid,1));
colorbar;
ylabel(sprintf('n=%d',size(Xid,1)));
saveas(h, [resdir '/expr.all.jpg'],'jpg');

% plot heatmaps: DE data
[~,i] = cluster_sort(logE(DE>0,:),2);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
subplot(1,6,1);
x = logE(DE>0,:) - nanmean(logE(DE>0,kmin),2);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
colormap(gene_colormap(1));
set(gca,'xtick',1:size(tm,2),'xticklabel',tm,'ytick',[]);
xtickangle(90);
set(gca,'color',[0.8 0.8 0.8]);
subplot(1,6,6);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
set(gca,'xtick',[],'ytick',1:100:size(Xid,1));
colorbar;
ylabel(sprintf('n=%d',sum(DE>0)));
saveas(h, [resdir '/expr.de.jpg'],'jpg');

% plot heatmaps: DE data, by model
M(DE==0,:) = repmat(M0(DE==0),1,size(utm,2));
M(DE==1,:) = logistic_eval(M1(DE==1,:),utm);
M(DE==2,:) = sigmo_eval(M2(DE==2,:),utm);
[~,i] = cluster_sort(M(DE>0,:),2);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
subplot(1,6,1);
x = M(DE>0,:) - nanmean(M(DE>0,umin),2);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
colormap(gene_colormap(1));
set(gca,'xtick',1:size(utm,2),'xticklabel',utm,'ytick',[]);
xtickangle(90);
set(gca,'color',[0.8 0.8 0.8]);
subplot(1,6,2);
x = logE(DE>0,:) - nanmean(logE(DE>0,kmin),2);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
colormap(gene_colormap(1));
set(gca,'xtick',1:size(tm,2),'xticklabel',tm,'ytick',[]);
xtickangle(90);
set(gca,'color',[0.8 0.8 0.8]);
subplot(1,6,6);
im = imagesc(x(i,:),[-4 4]);
im.AlphaData = ~isnan(x(i,:));
set(gca,'xtick',[],'ytick',1:100:size(Xid,1));
colorbar;
ylabel(sprintf('n=%d',sum(DE>0)));
saveas(h, [resdir '/expr.de.bymodel.jpg'],'jpg');



% plot examples
mkdir([resdir '/sigmoid_fit_examples']);

f1 = (Rsq >= 0.9) + (Err <= 2) > 0;
k = f1.*(max(logE,[],2) >= 9).*(DE==0)==1;
sum(k)
k = find(k);
if (max(size(k))>25)
    k = k(1:25);
end
plot_genes(logE,tm,Xid,M0,M1,M2,DE,Rsq,Err,k,[resdir '/sigmoid_fit_examples/fit0_']);
k = f1.*(max(logE,[],2) >= 8).*(DE==1)==1;
sum(k)
k = find(k);
if (max(size(k))>25)
    k = k(1:25);
end
plot_genes(logE,tm,Xid,M0,M1,M2,DE,Rsq,Err,k,[resdir '/sigmoid_fit_examples/fit1_']);
k = f1.*(max(logE,[],2) >= 6).*(DE==2)==1;
sum(k)
k = find(k);
if (max(size(k))>25)
    k = k(1:25);
end
plot_genes(logE,tm,Xid,M0,M1,M2,DE,Rsq,Err,k,[resdir '/sigmoid_fit_examples/fit2_']);

f0 = (Rsq <= 0.5) + (Err >= 10) > 0;
k = f0.*(f1==0).*(max(logE,[],2) >= 9).*(DE==0)==1;
sum(k)
k = find(k);
if (max(size(k))>25)
    k = k(1:25);
end
plot_genes(logE,tm,Xid,M0,M1,M2,DE,Rsq,Err,k,[resdir '/sigmoid_fit_examples/nofit0_']);
k = f0.*(f1==0).*(max(logE,[],2) >= 7).*(DE==1)==1;
sum(k)
k = find(k);
if (max(size(k))>25)
    k = k(1:25);
end
plot_genes(logE,tm,Xid,M0,M1,M2,DE,Rsq,Err,k,[resdir '/sigmoid_fit_examples/nofit1_']);
k = f0.*(f1==0).*(max(logE,[],2) >= 5).*(DE==2)==1;
sum(k)
k = find(k);
if (max(size(k))>25)
    k = k(1:25);
end
plot_genes(logE,tm,Xid,M0,M1,M2,DE,Rsq,Err,k,[resdir '/sigmoid_fit_examples/nofit2_']);
