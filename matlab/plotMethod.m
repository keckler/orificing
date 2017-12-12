function [] = plotMethod(ng, Q_sorted, zones_sorted, method)

i = 1;
while i < ng+1
    j = 1;
    while j < length(Q_sorted)+1
        if zones_sorted(j) == i
            h(j) = bar(j,Q_sorted(j));
            set(h(j),'FaceColor',[i/ng 1-(i/ng) i/ng], 'EdgeColor',[i/ng 1-(i/ng) i/ng]);
            hold on;
        end
        j = j + 1;
    end
    i = i + 1;
end
grid on; xlabel('assembly'); ylabel('assembly power (W)'); title(method); set(gca, 'XTickMode','Auto'); set(gca,'YScale','log');