function p = hygetail(X,M,K,N)
% X = values
% M = the corresponding size of the population
% K = number of items with the desired characteristic in the population
% N = and number of samples drawn 

p = 1 - hygecdf(X,M,K,N) + hygepdf(X,M,K,N);

