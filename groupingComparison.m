clear variables;

assemblyPowerThreshold = 0.01E6;
ngroups = [10:40];
nmethods = 4;
powerDetectorFiles = {'~/Downloads/BnB_det0','~/Downloads/BnB_det1','~/Downloads/BnB_det2','~/Downloads/BnB_det3'};

Q = readQ(powerDetectorFiles);
[Q, ~] = formMap(Q, assemblyPowerThreshold);
Q_ave = sum(Q,2)/length(powerDetectorFiles);
Q_ave_sorted = sort(Q_ave);

linearpfs = [];
logpfs = [];
kmeanspfs = [];
kmedpfs = [];

for ng = ngroups
    figure;
    
    linearZones = linearGrouping(Q_ave, Q_ave_sorted, ng, nmethods);
    [linearpf, linearpfgroup] = evaluateGroups(linearZones, Q_ave);
    
    logZones = logGrouping(Q_ave, Q_ave_sorted, ng, nmethods);
    [logpf, logpfgroup] = evaluateGroups(logZones, Q_ave);
    
    kmeansZones = kmeansGrouping(Q_ave, Q_ave_sorted, ng, nmethods);
    [kmeanspf, kmeanspfgroup] = evaluateGroups(kmeansZones, Q_ave);
    
    kmedZones = kmedoidsGrouping(Q_ave, Q_ave_sorted, ng, nmethods);
    [kmedpf, kmedpfgroup] = evaluateGroups(kmedZones, Q_ave);
    
    linearpfs = [linearpfs, linearpf];
    logpfs = [logpfs, logpf];
    kmeanspfs = [kmeanspfs, kmeanspf];
    kmedpfs = [kmedpfs, kmedpf];
    
    drawnow;
    
    disp('----------------------------------------------------')
    disp([ng])
    disp('group peaking factors: linear, log, kmeans, kmedoids')
    disp([linearpf, logpf, kmeanspf, kmedpf])
    disp('peaking factor group: linear, log, kmeans, kmedoids')
    disp([linearpfgroup, logpfgroup, kmeanspfgroup, kmedpfgroup])
    disp('----------------------------------------------------')
end

plotPFs(ngroups, linearpfs, logpfs, kmeanspfs, kmedpfs)