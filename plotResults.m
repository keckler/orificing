function [] = plotResults(m, T_out, Q, nbins, j)

figure;
bar(m); hold on; 
plot(T_out); ylim([0,1000]); 
yyaxis right; plot(Q); 
title(sprintf('k = %i',nbins(j))); grid on; 
xlabel('assembly'); yyaxis left; ylabel('T_{out} (C), flow (kg/s)'), yyaxis right; ylabel('power (W)');
drawnow;