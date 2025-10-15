mkdir('results');
minE = log2(0.1);

% ---------------------------------------------------------------------
% load data
% ---------------------------------------------------------------------
load_norm_data = 0;

if (load_norm_data)

    X = importdata('data/RNAseq_data.csv');
    M = 0.005*X.data(:,2:end)';
    Xid = X.textdata(1,2:end)';
    tm = X.data(:,1)';
    cid = regexprep(cellstr(num2str(tm')),' ','')';
    S = repmat({'diap'},size(tm));
    [Xid,~,logE] = normalize_fpkm(Xid,Xid,M,zeros(size(M)),S,cid,'results',1,1,2.^minE);

    save('results/data.norm.mat','Xid','logE','cid','tm');
end

% -----------------------------------------------------------------------
% temporal order (pseudo-time)
% -----------------------------------------------------------------------
tso_run = 1;

if (tso_run)
    load('results/data.norm.mat','Xid','logE','cid','tm');
    n = size(logE,2);

    % PCA
    [pc,~,pcv] = pca(logE);

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
        title(sprintf('n = %d', size(logE,1)));
    end
    saveas(h, 'results/tso.pca.jpg','jpg');
    saveas(h, 'results/tso.pca.eps','epsc');

    % PCA temporal order
    [Si,Sd,h] = traveling_salesman(pc(:,2),pc(:,3));
    tid2 = 1+cumsum(Sd)*12;
    saveas(h, 'results/tso.jpg','jpg');
    saveas(h, 'results/tso.svg','svg');

    save('results/tso.mat','Si');

    C = corr(logE);
    [~,tid1] = cluster_sort(C);

    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    hold on;
    plot(tm,tid1,'.k','markersize',20);
    plot(tm,tid2,'.r','markersize',20);
    xlabel('original time');
    ylabel('inferred time');
    legend({'by correlation' 'by PCA'},'location','bestOutside');
    set(gca,'ylim',[0 25],'xlim',[0 25]);
    axis square;
    saveas(h, 'results/tso.times.jpg','jpg');

    clf;
    imagesc(C(tid1,tid1),[0.75 1]);
    colormap(intensity_colormap);
    colorbar;
    axis square;
    set(gca,'xtick',1:n,'xticklabel',cid(tid1));
    set(gca,'ytick',1:n,'yticklabel',cid(tid1));
    saveas(h, 'results/tso.corr.hitmap.o1.jpg','jpg');

    clf;
    imagesc(C(Si,Si),[0.75 1]);
    colormap(intensity_colormap);
    colorbar;
    axis square;
    set(gca,'xtick',1:n,'xticklabel',cid(Si));
    set(gca,'ytick',1:n,'yticklabel',cid(Si));
    saveas(h, 'results/tso.corr.hitmap.o2.jpg','jpg');

    % Entropy
    E = entropy(logE);

    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    hold on;
    plot(tm,E,'.','markersize',25);
    P = robustfit(tm,E);
    t = min(tm):max(tm);
    y = t*P(2) + P(1);
    plot(t,y,'-k','linewidth',2);
    hold off;
    axis tight;
    set(gca,'ytick',-10:0.5:10,'ylim',[-8.5 -6.5],'fontsize',20);
    xlabel('sample order');
    ylabel('entropy');
    saveas(h, 'results/tso.entropy.jpg','jpg');

    clf;
    hold on;
    plot(1:n,E(tid1),'.','markersize',25);
    P = robustfit(1:n,E(tid1));
    t = min(tm):max(tm);
    y = t*P(2) + P(1);
    plot(t,y,'-k','linewidth',2);
    hold off;
    axis tight;
    set(gca,'ytick',-10:0.5:10,'ylim',[-8.5 -6.5],'fontsize',20);
    xlabel('sample order');
    ylabel('entropy');
    saveas(h, 'results/tso.entropy.o1.jpg','jpg');

    clf;
    hold on;
    plot(tid2,E(Si),'.','markersize',25);
    P = robustfit(tid2,E(Si));
    t = min(tm):max(tm);
    y = t*P(2) + P(1);
    plot(t,y,'-k','linewidth',2);
    hold off;
    axis tight;
    set(gca,'ytick',-10:0.5:10,'ylim',[-8.5 -6.5],'fontsize',20);
    xlabel('sample order');
    ylabel('entropy');
    saveas(h, 'results/tso.entropy.o2.jpg','jpg');

    close all;
end

% ---------------------------------------------------------------------
% Sigmoid Model Fit
% ---------------------------------------------------------------------
sigmoid_model_fit = 1;

if (sigmoid_model_fit)
    load('results/data.norm.mat','Xid','logE','tm');
    load('results/tso.mat','Si');
    logE = logE(:,Si);
    fit_temporal_models(Xid,logE,tm,'results');
end
