function [x,y] = plotCdf(lst)

% This function recieves a list of values.
% Sorts the Data
% fits at normal distribution to the data
% export a x and y coordinate for 100 point along the distribution.

lst_sorted = sort(lst);

fit = fitdist(lst_sorted,'Normal');
mini = min(lst_sorted)-1;
maxi = max(lst_sorted)+1;
x = linspace(mini,maxi,100);

y = cdf(fit,x);