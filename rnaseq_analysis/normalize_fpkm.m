function [G1,G2,logM,logP] = normalize_fpkm(G1,G2,M,P,S,T,resdir,norm1,norm2,minE,minFPKM)

if (nargin < 8)
    norm1 = 1;
end
if (nargin < 9)
    norm2 = 1;
end
if (nargin < 10)
    minE = 0.2; % not in logscale
end
if (nargin < 11)
    minFPKM = 2;
end
NORMF = 60;

u = unique(S);
nc = max(size(u));

% first step: local normalization between samples within a dataset
% using control genes
MN = M;
PN = P;
if (norm1)
    for i = 1:nc
        j = strcmp(S,u{i});

        h = figure;
        scrsz = get(0,'ScreenSize');
        set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
        n1 = -1;
        n2 = 0;
        r = 1;
        while ((n1<n2) && (r<=15))
            x = MN(:,j) + PN(:,j);
            THR = max(prctile(x(x>minE),NORMF),2.^minFPKM);
            mnX = mean(x,2);
            stdX = std(x,1,2);
            x1 = log2(mnX);
            x2 = log2(stdX);
            k = (~isnan(x1)).*(~isnan(x2)).*(~isinf(x1)).*(~isinf(x2)) == 1;
            p = polyfit(x1(k.*(x1>2)==1),x2(k.*(x1>2)==1),1);
            c = (x1 >= log2(THR)).*(x2 <= polyval(p,x1)-1)==1;
            n1 = n2;
            n2 = sum(c);
            fprintf('norm %s [step %d]: controls=%d\n',regexprep(u{i},'_',' '),r,n2);

            if (n2==0)
                N = 1;
            else
                N = mean(x(c,:));
                N = N./mean(N);
            end
            MN(:,j) = (MN(:,j)./N);
            PN(:,j) = (PN(:,j)./N);

            subplot(3,5,r);
            if (sum(k)>0)
                dscatter(x1(k),x2(k));
                hold on;
                plot(x1,polyval(p,x1),'-k','linewidth',2);
                plot(x1,polyval(p,x1)-1,'-r','linewidth',2);
            end
            hold off;
            axis square;
            axis tight;
            xlabel('mean');
            ylabel('stdv');
            title(sprintf('norm %s, %d (n=%d, c=%d)',regexprep(u{i},'_',' '),r,sum(k),sum(c)));
            r = r+1;
        end
        saveas(h, [resdir '/norm.' u{i} '.jpg'],'jpg');
    end
end

% second step: global normalization between dataset
% by median expression
PN2 = PN;
MN2 = MN;
if (norm2)
    Mall = [];
    for i = 1:nc
        j = strcmp(S,u{i});
        RW = PN2(:,j) + MN2(:,j);
        Mall(1,i) = mean(RW(RW>minE));
    end

    for i = 1:nc
        j = strcmp(S,u{i});
        N = Mall(1,i)./mean(Mall(1,:));
        PN2(:,j) = PN2(:,j)./N;
        MN2(:,j) = MN2(:,j)./N;
        RW = PN2(:,j) + MN2(:,j);
        Mall(2,i) = mean(RW(RW>minE));
    end
    Mall
end
PN2(PN2<minE) = minE;
MN2(MN2<minE) = minE;
logP = log2(PN2);
logM = log2(MN2);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for i = 1:nc
    subplot(1,nc,i);
    j = strcmp(S,u{i});
    n = sum(j);
    hold on;
    RW = M(:,j);
    RW(RW<minE) = minE;
    logRW = log2(RW);
    boxplot(logM(:,j),'widths',1,'positions',(1:n) + 0.35,'Colors','r', ...
        'plotstyle','compact','MedianStyle','target','Symbol','.k','OutlierSize',4,'Jitter',0.1,'whisker',2);
    boxplot(logRW,'labels',T(j),'widths',1,'positions',1:n,'Colors','b',...
        'plotstyle','compact','MedianStyle','target','Symbol','.k','OutlierSize',4,'Jitter',0.1,'whisker',2);
    hold off;
    xlabel('samples');
    ylabel('FPKM; log2');
    set(gca,'ylim',[log2(minE) 15],'fontsize',16);
    title(regexprep(u{i},'_','-'));
end
saveas(h, [resdir '/norm.' num2str(NORMF) '.jpg'],'jpg');

% select expressed genes by maximal expression
maxE = [];
for i = 1:nc
    k = strcmp(S,u{i});
    maxE(:,i) = max([logP(:,k) logM(:,k)],[],2);
end

k = (sum(maxE>=minFPKM,2) >= floor(nc/2)+1);
fprintf('select by maximal expression: %d genes\n', sum(k));
logP = logP(k,:);
logM = logM(k,:);
G1 = G1(k,:);
G2 = G2(k,:);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for i = 1:nc
    subplot(1,nc,i);
    j = strcmp(S,u{i});
    hold on;
    n = sum(j);
    boxplot(logM(:,j),'widths',1,'positions',(1:n) + 0.35,'Colors','b', ...
        'plotstyle','compact','MedianStyle','target','Symbol','.k','OutlierSize',4,'Jitter',0.1,'whisker',2);
    boxplot(logP(:,j),'labels',T(j),'widths',1,'positions',1:n,'Colors','r',...
        'plotstyle','compact','MedianStyle','target','Symbol','.k','OutlierSize',4,'Jitter',0.1,'whisker',2);
    hold off;
    axis tight;
    xlabel('samples');
    ylabel(sprintf('FPKM; log2 (n=%d)',sum(k)));
    set(gca,'ylim',[log2(minE) 15],'fontsize',16);
    title(regexprep(u{i},'_',' '));
end
saveas(h, [resdir '/fpkm_hist_box.jpg'],'jpg');

% text files
for i = 1:nc
    j = strcmp(S,u{i});
    L = [{'ensid' 'gid'} strcat('M_',T(j)) strcat('P_',T(j))];
    write_text_file([resdir '/data.normFPKM.' u{i} '.txt'],[L;[G1 G2 num2cell([logM(:,j) logP(:,j)])]]);
end

close all;
