function h = plot_model_param(M0,M1,M2,DE)
% M0 = <h0>
% M1 = <b,h0,h1,t1>
% M2 = <b,h0,h1,h2,t1,t2>

s = hist(DE,[0 1 2]);
minT = min([M1(:,4);M2(:,5);M2(:,6)]);
maxT = max([M1(:,4);M2(:,5);M2(:,6)]);
minB = min([M1(:,1);M2(:,1)]);
maxB = max([M1(:,1);M2(:,1)]);
minX = min([M0;M1(:,2);M1(:,3);M2(:,2);M2(:,3);M2(:,4)]);
maxX = max([M0;M1(:,2);M1(:,3);M2(:,2);M2(:,3);M2(:,4)]);



h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

subplot(2,2,1);
x = minT:maxT;
hold on;
y = hist(M1(DE==1,4),x);
plot(x,y./sum(y),'-r','LineWidth',1.5);
y = hist(M2(DE==2,4),x);
plot(x,y./sum(y),'-b','LineWidth',1.5);
y = hist(M2(DE==2,5),x);
plot(x,y./sum(y),'-c','LineWidth',1.5);
hold off;
set(gca,'fontsize',18);
L = strcat({'M1' 'M2-1' 'M2-2'}', regexprep(cellstr(strcat(';n=',num2str([s(2) s(3) s(3)]'))),' ',''));
legend(L,'Box','off');
xlabel('fitted onset time')
ylabel('fraction of genes')
title('onset time');

subplot(2,2,2);
x = minB:0.1:maxB;
hold on;
y = hist(M1(DE==1,1),x);
plot(x,y./sum(y),'-r','LineWidth',1.5);
y = hist(M2(DE==2,1),x);
plot(x,y./sum(y),'-b','LineWidth',1.5);
hold off;
set(gca,'fontsize',18);
L = strcat({'M1' 'M2'}', regexprep(cellstr(strcat(';n=',num2str([s(2) s(3)]'))),' ',''));
legend(L,'Box','off');
xlabel('fitted slope; log2')
ylabel('fraction of genes')
title('slope');

subplot(2,2,3);
x = minX:maxX;
hold on;
y = hist(M0(DE==0,1),x);
plot(x,y./sum(y),'-k','LineWidth',1.5);
y = hist(M1(DE==1,2),x);
plot(x,y./sum(y),'-r','LineWidth',1.5);
y = hist(M1(DE==1,3),x);
plot(x,y./sum(y),'-m','LineWidth',1.5);
hold off;
set(gca,'fontsize',18);
L = strcat({'M0' 'M1-1' 'M1-2'}', regexprep(cellstr(strcat(';n=',num2str([s(1) s(2) s(2)]'))),' ',''));
legend(L,'Box','off');
xlabel('expression level; log2')
ylabel('fraction of genes')
title('expression level');

subplot(2,2,4);
hold on;
y = hist(M2(DE==2,2),x);
plot(x,y./sum(y),'-b','LineWidth',1.5);
y = hist(M2(DE==2,3),x);
plot(x,y./sum(y),'-c','LineWidth',1.5);
y = hist(M2(DE==2,4),x);
plot(x,y./sum(y),'-g','LineWidth',1.5);
hold off;
set(gca,'fontsize',18);
L = strcat({'M2-1' 'M2-2' 'M2-3'}', regexprep(cellstr(strcat(';n=',num2str([s(3) s(3) s(3)]'))),' ',''));
legend(L,'Box','off');
xlabel('expression level; log2')
ylabel('fraction of genes')
title('expression level');
