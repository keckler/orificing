function [zones] = kmedoidsGrouping(Q, Q_sorted, ng, nmethods)

zones = kmedoids(Q,ng);
zones_sorted = sort(zones);
subplot(1,nmethods,4); 
method ='k-medoids';

plotMethod(ng, Q_sorted, zones_sorted, method)