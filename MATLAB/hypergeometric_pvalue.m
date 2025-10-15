function p = hypergeometric_pvalue(X,M,K,N,lower)
% p = hypergeometric_pvalue(X,M,K,N,lower)
%
% X = number of items with the desired characteristic in sample
% M = size of the population
% K = number of items with the desired characteristic in population
% N = size of sample drawn
% lower = calculate the lower tail (by default calculates the upper tail)
%
% Upper tail = p is the probability of drawing X or MORE items with the 
%              desired characteristic when drawing N out of M items in 
%              which K items has the desired characteristic
% Lower tail = p is the probability of drawing X or LESS items with the 
%              desired characteristic when drawing N out of M items in 
%              which K items has the desired characteristic

if (nargin < 5)
    lower = 0;
end

if (lower) 
    p = hygecdf(X,M,K,N);
else
    p = hygecdf(X,M,K,N,'upper') + hygepdf(X,M,K,N);
end
