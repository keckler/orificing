function [zones] = logGrouping(Q, Q_sorted, ng, nmethods)

zones = discretize(Q,logspace(log10(min(min(Q))),log10(max(max(Q))),ng));
zones_sorted = sort(zones);
subplot(1,nmethods,2); 
method = 'logarithmically spaced';

plotMethod(ng, Q_sorted, zones_sorted, method)