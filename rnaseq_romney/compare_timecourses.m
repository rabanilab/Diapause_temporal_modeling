function compare_timecourses(Xid,logE,tm,u,u2,resdir,minT)

if (nargin < 7)
    minT = 5;
end

llr_alpha = 0.01;
sig = prctile(std(logE,1,2),70);
fprintf('Compare timecourses: sigma=%.2f, alpha=%.2f\n', sig, llr_alpha);


% fit a model to both conditions
load(['results.' u2 '/model_fit.mat'],'M0','X0','M1','X1','M2','X2');
load(['results.' u2 '/model_fit_select.mat'],'P','DE','Rsq','Err');
sig = prctile(std(logE,1,2),70);
Lx(DE==0,1) = sum(log2(normpdf(logE(DE==0,:),X0(DE==0,:),sig)),2);
Lx(DE==1,1) = sum(log2(normpdf(logE(DE==1,:),X1(DE==1,:),sig)),2);
Lx(DE==2,1) = sum(log2(normpdf(logE(DE==2,:),X2(DE==2,:),sig)),2);
Nx = ones(size(DE));
Nx(DE==1) = 4;
Nx(DE==2) = 6;

% fit a model to each condition (two models)
[logE1,tm1,~,M0_1,X0_1,M1_1,X1_1,M2_1,X2_1,P1,DE1,R1,E1] = load_tcourse(u{1}); % escape
L1(DE1==0,1) = sum(log2(normpdf(logE1(DE1==0,:),X0_1(DE1==0,:),sig)),2);
L1(DE1==1,1) = sum(log2(normpdf(logE1(DE1==1,:),X1_1(DE1==1,:),sig)),2);
L1(DE1==2,1) = sum(log2(normpdf(logE1(DE1==2,:),X2_1(DE1==2,:),sig)),2);
N1 = ones(size(DE1));
N1(DE1==1) = 4;
N1(DE1==2) = 6;

[logE2,tm2,~,M0_2,X0_2,M1_2,X1_2,M2_2,X2_2,P2,DE2,R2,E2] = load_tcourse(u{2}); % diapause
L2(DE2==0,1) = sum(log2(normpdf(logE2(DE2==0,:),X0_2(DE2==0,:),sig)),2);
L2(DE2==1,1) = sum(log2(normpdf(logE2(DE2==1,:),X1_2(DE2==1,:),sig)),2);
L2(DE2==2,1) = sum(log2(normpdf(logE2(DE2==2,:),X2_2(DE2==2,:),sig)),2);
N2 = ones(size(DE2));
N2(DE2==1) = 4;
N2(DE2==2) = 6;

Ly = L1+L2;
Ny = N1+N2;

% likelihood ratio
[~,Q(:,1)] = lratiotest(Ly,min(Lx,Ly),max(Ny-Nx,1));
fprintf('rows reject Lx (=single model) for Ly (=two models): %d\n', sum(Q(:,1)<llr_alpha));
DExy = zeros(size(Q,1),1);
DExy(Q(:,1)<llr_alpha) = 1;

% plot single model fit
plot_data(Xid(DExy==0),logE(DExy==0,:),tm,M0(DExy==0,:),M1(DExy==0,:),...
    M2(DExy==0,:),DE(DExy==0,:),Err(DExy==0,:),Rsq(DExy==0,:),resdir,minT);

write_text_file([resdir '/model_fit_param.one.txt'], ...
    [{'id' 'model' 'rsq' 'mse' 'b' 'h0' 'h1' 'h2' 't1' 't2'}; ...
    [Xid(DExy==0,:) num2cell([DE(DExy==0,:) Rsq(DExy==0,:) Err(DExy==0,:) P(DExy==0,:)])]]);

% plot two models fit
plot_data_two(Xid(DExy==1,:),logE1(DExy==1,:),tm1,M0_1(DExy==1,:),M1_1(DExy==1,:),M2_1(DExy==1,:),DE1(DExy==1,:),...
    logE2(DExy==1,:),tm2,M0_2(DExy==1,:),M1_2(DExy==1,:),M2_2(DExy==1,:),DE2(DExy==1,:),resdir,minT);

write_text_file([resdir '/model_fit_param.two.txt'], ...
    [{'id' 'model' 'rsq' 'mse' 'b' 'h0' 'h1' 'h2' 't1' 't2' 'rsq' 'mse' 'b' 'h0' 'h1' 'h2' 't1' 't2'}; ...
    [Xid(DExy==1,:) num2cell([DExy(DExy==1,:) ...
    R1(DExy==1,:) E1(DExy==1,:) P1(DExy==1,:) ...
    R2(DExy==1,:) E2(DExy==1,:) P2(DExy==1,:)])]]);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
subplot(1,2,1);
hold on;
dscatter(Ly,Lx,'MSIZE',25);
line([-600 0],[-600 0],'LineStyle','-','color','k','linewidth',1.2);
hold off;
axis square;
set(gca,'xlim',[-600 0],'ylim',[-600 0],'fontsize',18);
xlabel('two models');
ylabel('one model');
title(sprintf('log likelihood; #(Ly<=Lx)=%d\n', sum(Ly<=Lx)));
subplot(1,2,2);
hold on;
dscatter(Ny,Nx,'MSIZE',100);
line([0 15],[0 15],'LineStyle','-','color','k','linewidth',1.2);
hold off;
axis square;
set(gca,'xlim',[0 15],'ylim',[0 15],'fontsize',18);
xlabel('two models');
ylabel('one model');
title(sprintf('number of parameters;#(Ny<=Nx)=%d\n', sum(Ny<=Nx)));
saveas(h, [resdir '/tcourses.fit.jpg'],'jpg');

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

N = [DE1 DE2];
[s,~,t] = unique(N(DExy==1,:),'rows');
c = accumarray(t,1);
[[{'cnt' 'pct'} u];num2cell([c round(100*c./sum(c)) s])]

k = (DExy==1).*(DE1==0).*(DE2==0) == 1;
hold on;
dscatter(M0_1(k),M0_2(k),'MSIZE',20);
line([0 14],[0 14],'LineStyle','-','color','k','linewidth',1.2);
hold off;
xlabel(u{1});
ylabel(u{2});
axis square;
set(gca,'xlim',[0 14],'ylim',[0 14],'fontsize',18);
saveas(h, [resdir '/tcourses.0.0.jpg'],'jpg');

close all;


function [logEi,tmi,Xid,M0,X0,M1,X1,M2,X2,P,DE,Rsq,Err] = load_tcourse(u)
    
load('results/data.norm.mat','Xid','logE','tm','cid');
S = strtok(cid,'_');

j = strcmp(S,u);
logEi = logE(:,j);
tmi = tm(:,j);
load(['results.' u '/model_fit.mat'],'M0','X0','M1','X1','M2','X2');
load(['results.' u '/model_fit_select.mat'],'P','DE','Rsq','Err');

