result_dir = 'Results';
mkdir(result_dir);
rng('default');

% -----------------------------------------------------------------------
% qc
% -----------------------------------------------------------------------
X = importdata('data/VitaminD_FPKM_RNAseq_star2206.csv');
rid = X.textdata(2:end,1);
cid = regexprep(X.textdata(1,2:end),'"','');
tid = str2double(regexprep(strtok(cid,'_'),'X',''));
[~,j] = sort(tid);
tid = tid(j);
cid = cid(j);
FPKM = X.data(:,j);

% reorder samples 1&2
tid = tid([2 1 3:end]);
cid = cid([2 1 3:end]);
FPKM = FPKM(:,[2 1 3:end]);

% % remove sample X11, X8
tid = tid([1:8 10:11 13:end]);
cid = cid([1:8 10:11 13:end]);
FPKM = FPKM(:,[1:8 10:11 13:end]);

n = size(FPKM,2);
minE = prctile(FPKM(:),50);
FPKM(FPKM<minE) = minE;
logFPKM = log2(FPKM);
[cid; num2cell([min(logFPKM);max(logFPKM);sum(logFPKM>log2(minE),1)])]'

k = (sum(FPKM,2) >= n*minE).*(max(FPKM,[],2) >= 4*minE) == 1;
fprintf('Expressed genes: %d (%.1f%%)\n',sum(k),100*sum(k)/size(k,1))

% histograms
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

x = -1:0.5:12;
clear y;
for i = 1:n
    y(i,:) = hist(logFPKM(k,i),x);
    y(i,:) = y(i,:)./sum(y(i,:));
end
plot(x,y,'-','linewidth',2,'marker','.','markersize',15);
xlabel('log FPKM');
ylabel('fraction');
legend(cid,'location', 'bestOutside');
title(sprintf('expressed genes = %d (%.1f%%)',sum(k),100*sum(k)./size(k,1)));
saveas(h, [result_dir '/hist.jpg'],'jpg');

clf;
for i = 1:n-1
    subplot(5,5,i);
    dscatter(logFPKM(k,i),logFPKM(k,i+1));
    axis square;
    set(gca,'xlim',[-1 12],'ylim',[-1 12]);
    xlabel(cid(i));
    ylabel(cid(i+1));
end
saveas(h, [result_dir '/corr.xy.jpg'],'jpg');

clf;
C = corr(logFPKM(k,:));
imagesc(C,[0.75 1]);
colormap(intensity_colormap);
colorbar;
axis square;
set(gca,'xtick',1:n,'xticklabel',cid);
set(gca,'ytick',1:n,'yticklabel',cid);
saveas(h, [result_dir '/corr.hitmap.jpg'],'jpg');

close all;


% -----------------------------------------------------------------------
% temporal order
% -----------------------------------------------------------------------
% PCA
[pc,~,pcv] = pca(logFPKM(k,:));

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for j = 1:4
    subplot(2,2,j);
    hold on;
    plot(pc(:,j),pc(:,j+1),'.','markersize',25);
    for i = 1:n
        text(pc(i,j),pc(i,j+1),cid{i},'FontSize', 14);
    end
    hold off;
    xlabel(sprintf('PC %d (var = %.1f%%)',j,pcv(j)));
    ylabel(sprintf('PC %d (var = %.1f%%)',j+1,pcv(j+1)));
    axis square;
    title(sprintf('n = %d', sum(k)));
end
saveas(h, [result_dir '/pca.jpg'],'jpg');
close all;

% PCA temporal order
[Si,Sd,h] = traveling_salesman(pc(:,2),pc(:,3));
tid2 = 1+cumsum(Sd)*12;
saveas(h, [result_dir '/tso.jpg'],'jpg');

C = corr(logFPKM(k,:));
[~,tid1] = cluster_sort(C);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
hold on;
plot(tid,tid1,'.k','markersize',20);
plot(tid,tid2,'.r','markersize',20);
xlabel('original time');
ylabel('inferred time');
legend({'by correlation' 'by PCA'},'location','bestOutside');
set(gca,'ylim',[0 25],'xlim',[0 25]);
axis square;
saveas(h, [result_dir '/times.jpg'],'jpg');

clf;
imagesc(C(tid1,tid1),[0.75 1]);
colormap(intensity_colormap);
colorbar;
axis square;
set(gca,'xtick',1:n,'xticklabel',cid(tid1));
set(gca,'ytick',1:n,'yticklabel',cid(tid1));
saveas(h, [result_dir '/corr.hitmap.o1.jpg'],'jpg');

clf;
imagesc(C(Si,Si),[0.75 1]);
colormap(intensity_colormap);
colorbar;
axis square;
set(gca,'xtick',1:n,'xticklabel',cid(Si));
set(gca,'ytick',1:n,'yticklabel',cid(Si));
saveas(h, [result_dir '/corr.hitmap.o2.jpg'],'jpg');

close all;

% Entropy
E = entropy(logFPKM(k,:));

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
hold on;
plot(tid,E,'.','markersize',25);
P = robustfit(tid,E);
t = min(tid):max(tid);
y = t*P(2) + P(1);
plot(t,y,'-k','linewidth',2);
hold off;
axis tight;
set(gca,'ytick',-10:0.5:10,'ylim',[-8.5 -6.5],'fontsize',20);
xlabel('sample order');
ylabel('entropy');
saveas(h, [result_dir '/entropy.jpg'],'jpg');

clf;
hold on;
plot(1:n,E(tid1),'.','markersize',25);
P = robustfit(1:n,E(tid1));
t = min(tid):max(tid);
y = t*P(2) + P(1);
plot(t,y,'-k','linewidth',2);
hold off;
axis tight;
set(gca,'ytick',-10:0.5:10,'ylim',[-8.5 -6.5],'fontsize',20);
xlabel('sample order');
ylabel('entropy');
saveas(h, [result_dir '/entropy.o1.jpg'],'jpg');

clf;
hold on;
plot(tid2,E(Si),'.','markersize',25);
P = robustfit(tid2,E(Si));
t = min(tid):max(tid);
y = t*P(2) + P(1);
plot(t,y,'-k','linewidth',2);
hold off;
axis tight;
set(gca,'ytick',-10:0.5:10,'ylim',[-8.5 -6.5],'fontsize',20);
xlabel('sample order');
ylabel('entropy');
saveas(h, [result_dir '/entropy.o2.jpg'],'jpg');

close all;

% -----------------------------------------------------------------------
% clustering
% -----------------------------------------------------------------------

% k-means
kn = 25;
S = logFPKM(k,:) - mean(logFPKM(k,:),2);

[cidx, ctrs] = kmeans(S,kn,'dist','corr','rep',5,'disp','final','MaxIter',250,'Replicates',25);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for c = 1:kn
    subplot(5,5,c);
    plot(1:n,S((cidx == c),:)','-b');
    axis tight;
    set(gca,'ylim',[-4 4]);
    title(sprintf('n=%d',sum(cidx == c)));
end
sgtitle('K-Means Clustering of Profiles');
saveas(h, [result_dir '/kmeans.jpg'],'jpg');
write_text_file([result_dir '/kmeans.txt'],[rid(k) num2cell(cidx)]);
close all;

[cidx, ctrs] = kmeans(S(:,tid1),kn,'dist','corr','rep',5,'disp','final','MaxIter',250,'Replicates',25);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for c = 1:kn
    subplot(5,5,c);
    plot(1:n,S((cidx == c),tid1)','-b');
    axis tight;
    set(gca,'ylim',[-4 4]);
    title(sprintf('n=%d',sum(cidx == c)));
end
sgtitle('K-Means Clustering of Profiles');
saveas(h, [result_dir '/kmeans.o1.jpg'],'jpg');
write_text_file([result_dir '/kmeans.o1.txt'],[rid(k) num2cell(cidx)]);
close all;

[cidx, ctrs] = kmeans(S(:,Si),kn,'dist','corr','rep',5,'disp','final','MaxIter',250,'Replicates',25);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for c = 1:kn
    subplot(5,5,c);
    plot(tid2,S((cidx == c),Si)','-b');
    axis tight;
    set(gca,'ylim',[-4 4]);
    title(sprintf('n=%d',sum(cidx == c)));
end
sgtitle('K-Means Clustering of Profiles');
saveas(h, [result_dir '/kmeans.o2.jpg'],'jpg');
write_text_file([result_dir '/kmeans.o2.txt'],[rid(k) num2cell(cidx)]);
close all;
