function fit_temporal_models(Xid,logE,tm,resdir,minT)

if (nargin < 5)
    minT = 5;
end

sig = prctile(std(logE,1,2),70);
llr_alpha = 0.01;
mkdir(resdir);

% ---------------------------------------------------------------------
% Model fit
% ---------------------------------------------------------------------
fprintf('Model fitting: data %d x %d\n', size(logE));

if (exist([resdir '/model_fit.mat'],'file'))
    load([resdir '/model_fit.mat'],'M0','X0','R0','E0','M1','X1','R1','E1','M2','X2','R2','E2');
else

    % fit a constant expression level per gene
    M0 = mean(logE,2);
    X0 = repmat(M0,1,size(logE,2));
    R0 = zeros(size(M0));
    E0 = sum((X0-logE).^2,2);

    % fit a single sigmoid function per gene
    % M1 = <b,h0,h1,t1> (parameters)
    [M1,E1] = logistic_fit(logE,tm,10);
    X1 = logistic_eval(M1,tm);
    R1 = diag(corr(X1',logE')).^2;
    sum(R1 > 0.8)

    % fit a double sigmoid function per gene
    % M2 = <b,h0,h1,h2,t1,t2> (parameters)
    [M2,E2] = sigmo_fit(logE,tm,10);
    X2 = sigmo_eval(M2,tm);
    R2 = diag(corr(X2',logE')).^2;
    sum(R2 > 0.8)

    save([resdir '/model_fit.mat'],'M0','X0','R0','E0','M1','X1','R1','E1','M2','X2','R2','E2');
end

% ---------------------------------------------------------------------
% Model selection
% ---------------------------------------------------------------------
fprintf('Model selection: sigma=%.2f, alpha=%.2f\n', sig, llr_alpha);

if (exist([resdir '/model_fit_select.mat'],'file'))
    load([resdir '/model_fit_select.mat'],'P','DE','Rsq','Err');
else

    % likelihood ratio test
    L0 = sum(log2(normpdf(logE,X0,sig)),2);
    L1 = sum(log2(normpdf(logE,X1,sig)),2);
    L2 = sum(log2(normpdf(logE,X2,sig)),2);

    fprintf('rows with L1<L0: %d\n', sum(L1<L0))
    [~,Q(:,1)] = lratiotest(L1,min(L0,L1),size(M1,2)-size(M0,2));
    fprintf('rows reject L0 for L1: %d\n', sum(Q(:,1)<llr_alpha));
    fprintf('rows with L2<L1: %d\n', sum(L2<L1))
    [~,Q(:,2)] = lratiotest(L2,min(L1,L2),size(M2,2)-size(M1,2));
    fprintf('rows reject L1 for L2: %d\n', sum(Q(:,2)<llr_alpha));
    fprintf('rows with L2<L0: %d\n', sum(L2<L0))
    [~,Q(:,3)] = lratiotest(L2,min(L0,L2),size(M2,2)-size(M0,2));
    fprintf('rows reject L0 for L2: %d\n', sum(Q(:,3)<llr_alpha));

    DE = zeros(size(Q,1),1);
    DE(Q(:,1)<llr_alpha) = 1;
    DE((Q(:,1)<llr_alpha).*(Q(:,2)<llr_alpha)==1) = 2;
    DE((Q(:,1)>=llr_alpha).*(Q(:,3)<llr_alpha)==1) = 2;
    Rsq = zeros(size(Q,1),1);
    Rsq(DE==1) = R1(DE==1);
    Rsq(DE==2) = R2(DE==2);
    Err = E0;
    Err(DE==1) = E1(DE==1);
    Err(DE==2) = E2(DE==2);

    P = M2; % <b,h0,h1,h2,t1,t2>
    P(DE==1,:) = [M1(DE==1,1) M1(DE==1,2) M1(DE==1,3) zeros(sum(DE==1),1) M1(DE==1,4) inf(sum(DE==1),1)];
    P(DE==0,:) = [zeros(sum(DE==0),1) M0(DE==0,1) zeros(sum(DE==0),1) zeros(sum(DE==0),1) inf(sum(DE==0),1) inf(sum(DE==0),1)];
    save([resdir '/model_fit_select.mat'],'P','DE','Rsq','Err');
end

% model selection
f = (Rsq >= 0.8) + (Err <= 5) > 0;
num2cell([sum(f) sum(f)./size(DE,1)])

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
subplot(1,2,1);
y = hist(DE,0:2);
bar(y./sum(y));
set(gca,'xtick',1:3,'xticklabel',{'M0' 'M1' 'M2'},'ylim',[0 1],'fontsize',15);
ylabel('fraction of genes');
title(sprintf('all genes (n=%d)', sum(y)));
for i = 1:max(size(y))
    text(i,0.05+y(i)./sum(y),sprintf('n=%d (%.0f %%)',y(i),100*y(i)/sum(y)),'fontsize',15);
end
subplot(1,2,2);
y1 = hist(DE(f),0:2);
bar(y1./sum(y1));
set(gca,'xtick',1:3,'xticklabel',{'M0' 'M1' 'M2'},'ylim',[0 1],'fontsize',15);
ylabel('fraction of genes');
title(sprintf('selected genes with a good fit (n=%d)', sum(y1)));
for i = 1:max(size(y1))
    text(i,0.05+y1(i)./sum(y1),sprintf('n=%d (%.0f %%)',y1(i),100*y1(i)/sum(y1)),'fontsize',15);
end
saveas(h,[resdir '/model_fit_select.jpg'],'jpg');

% model parameters
h = plot_model_param(M0,M1,M2,DE);
saveas(h,[resdir '/model_fit_param.all.jpg'],'jpg');
h = plot_model_param(M0(f,:),M1(f,:),M2(f,:),DE(f,:));
saveas(h,[resdir '/model_fit_param.good_fit.jpg'],'jpg');

% plot model fits
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
subplot(2,2,1);
x = 0:0.05:1;
hold on;
y = hist(Rsq,x);
plot(x,y./sum(y),'-k','linewidth',2);
y = hist(Rsq(DE==1),x);
plot(x,y./sum(y),'-b','linewidth',2);
y = hist(Rsq(DE==2),x);
plot(x,y./sum(y),'-r','linewidth',2);
hold off;
ylabel('fraction of genes');
xlabel('r-squared');
set(gca,'fontsize',14);
legend({'all' 'M1' 'M2'},'box','off');
subplot(2,2,2);
dscatter(Rsq,log2(Err),'MSIZE',20);
axis square;
xlabel('rsq');
ylabel('mse; log2');
set(gca,'fontsize',14);
subplot(2,2,3);
dscatter(Rsq(DE==1),log2(Err(DE==1)),'MSIZE',20);
axis square;
xlabel('rsq');
ylabel('mse; log2');
set(gca,'fontsize',14);
title(sprintf('n=%d (Model = 1)',sum(DE==1)));
subplot(2,2,4);
dscatter(Rsq(DE==2),log2(Err(DE==2)),'MSIZE',20);
axis square;
xlabel('rsq');
ylabel('mse; log2');
set(gca,'fontsize',14);
title(sprintf('n=%d (Model = 2)',sum(DE==2)));
saveas(h,[resdir '/model_fit_sigmoid.jpg'],'jpg');

% param values
write_text_file([resdir '/model_fit_param.txt'], ...
    [{'id' 'model' 'rsq' 'mse' 'b' 'h0' 'h1' 'h2' 't1' 't2'};...
    [Xid num2cell([DE Rsq Err P])]]);

write_text_file([resdir '/model_fit_param.DE.txt'], ...
    [{'id' 'model' 'rsq' 'mse' 'b' 'h0' 'h1' 'h2' 't1' 't2'};...
    [Xid(DE>0) num2cell([DE(DE>0,:) Rsq(DE>0,:) Err(DE>0,:) P(DE>0,:)])]]);

% plot model
plot_data(Xid,logE,tm,M0,M1,M2,DE,Err,Rsq,resdir,minT);

close all;
