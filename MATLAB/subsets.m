function S = subsets(n)

S = [];
for k = n:-1:0
    %fprintf('>>%d,%d\n',n, k);
    p = all_subsets(n,k);
    S = [S; p];
end



function p = all_subsets(n,k)

%fprintf('%d,%d\n',n, k);

if (k > n)
    p = [];
elseif (k == 0)
    p = zeros(1,n);
elseif (n==1 && k==1)
    p = 1;
elseif (n==1 && k==0)
    p = 0;
else
    p1 = all_subsets(n-1, k-1);
    p2 = all_subsets(n-1, k);
    p = [[ones(size(p1,1),1) p1]; [zeros(size(p2,1),1) p2]];
end

