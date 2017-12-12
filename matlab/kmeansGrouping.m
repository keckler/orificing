function [zones] = kmeansGrouping(Q, Q_sorted, ng, nmethods)

zones = kmeans(Q,ng);
zones_sorted = sort(zones);
subplot(1,nmethods,3); 
method = 'k-means';

plotMethod(ng, Q_sorted, zones_sorted, method)