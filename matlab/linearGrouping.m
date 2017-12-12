function [zones] = linearGrouping(Q, Q_sorted, ng, nmethods)

zones = discretize(Q,linspace(min(min(Q)),max(max(Q)),ng));
zones_sorted = sort(zones);
subplot(1,nmethods,1);
method = ('linearly spaced');

plotMethod(ng, Q_sorted, zones_sorted, method)