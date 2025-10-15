function p = binomial_pvalue(X,N,F,lower)
% p = binomial_pvalue(X,N,P)
%
% X = number of items with the desired characteristic in sample
% N = size of sample drawn
% F = frequency of items with the desired characteristic in population
% lower = calculate the lower tail (by default calculates the upper tail)
%
% Upper tail = p is the probability of observing X or MORE items with the 
%              desired characteristic when drawing N items from a population
%              in which the frequency of the desired characteristic is P
% Lower tail = p is the probability of observing X or LESS items with the 
%              desired characteristic when drawing N items from a population
%              in which the frequency of the desired characteristic is P

if (nargin < 4)
    lower = 0;
end

if (lower)
    p = binocdf(X,N,F);
else
    p = binocdf(X,N,F,'upper') + binopdf(X,N,F);
end

