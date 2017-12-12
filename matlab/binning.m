function zones = binning(Q_ave, Q, nb, groupingMethod)

switch groupingMethod
    case 'lin'
        zones = discretize(Q_ave,linspace(min(min(Q)),max(max(Q)),nb));
    case 'log'
        zones = discretize(Q_ave,logspace(log10(min(min(Q))),log10(max(max(Q))),nb));
    case 'kmeans'
        zones = kmeans(Q_ave,nb);
    case 'kmedoids'
        zones = kmedoids(Q_ave,nb);
end