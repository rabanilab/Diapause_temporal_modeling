function plot_genes(logE,tm,Xid,M0,M1,M2,DE,Rsq,Err,k,output_pref)

tx = min(tm):0.1:max(tm);
minE = floor(min(logE(:)));
maxE = ceil(max(logE(:)));

h = figure;
for i = k'
    clf;
    hold on;
    plot(tm,logE(i,:),'.k','markersize',25);
    m = repmat(M0(i,:),size(tx));
    plot(tx,m,'-m','linewidth',1);
    m = logistic_eval(M1(i,:),tx);
    plot(tx,m,'-r','linewidth',1);
    m = sigmo_eval(M2(i,:),tx);
    plot(tx,m,'-b','linewidth',1);
    hold off;
    xlabel('pseudo-time');
    ylabel('log2(FPKM)');
    axis tight;
    set(gca,'ylim',[minE maxE],'fontsize',14);
    title(sprintf('%s M=%d, r2=%.2f, e=%.2f',Xid{i},DE(i,:),Rsq(i),Err(i)));
    legend({'data' 'M0' 'M1' 'M2'},'location','bestOutside','box','off');
    saveas(h, [output_pref Xid{i} '.jpg'],'jpg');
end
close all;
