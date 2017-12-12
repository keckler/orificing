function [] = printPFs(ngroups, linearpfs, logpfs, kmeanspfs, kmedspfs)

figure;
semilogy(ngroups, linearpfs, ngroups, logpfs, ngroups, kmeanspfs, ngroups, kmedspfs); 
grid on; xlabel('number of orifice groups'); ylabel('largest maximum/minimum power ratio within an orifice group');
legend('linear spacing', 'log spacing', 'kmeans grouping', 'kmedoids grouping');