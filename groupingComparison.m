figure;

Q_sorted = sort(Q);


zones = discretize(Q,linspace(min(min(Q)),max(max(Q)),20));
zones_sorted = sort(zones);
subplot(1,5,1); 
i = 1;
while i < 20+1
    j = 1;
    while j < length(Q)+1
        if zones_sorted(j) == i
            h(j) = bar(j,Q_sorted(j));
            set(h(j),'FaceColor',[i/20 1-(i/20) i/20], 'EdgeColor',[i/20 1-(i/20) i/20]);
            hold on;
        end
        j = j + 1;
    end
    i = i + 1;
end
grid on; xlabel('assembly'); ylabel('assembly power (W)'); title('linearly spaced'); set(gca, 'XTickMode','Auto'); set(gca,'YScale','log');

zones = discretize(Q,logspace(log10(min(min(Q))),log10(max(max(Q))),20));
zones_sorted = sort(zones);
subplot(1,5,2); 
i = 1;
while i < 20+1
    j = 1;
    while j < length(Q)+1
        if zones_sorted(j) == i
            h(j) = bar(j,Q_sorted(j));
            set(h(j),'FaceColor',[i/20 1-(i/20) i/20], 'EdgeColor',[i/20 1-(i/20) i/20]);
            hold on;
        end
        j = j + 1;
    end
    i = i + 1;
end
grid on; xlabel('assembly'); ylabel('assembly power (W)'); title('logarithmically spaced'); set(gca, 'XTickMode','Auto'); set(gca,'YScale','log');
    
zones = kmeans(Q,20);
zones_sorted = sort(zones);
subplot(1,5,3); 
i = 1;
while i < 20+1
    j = 1;
    while j < length(Q)+1
        if zones_sorted(j) == i
            h(j) = bar(j,Q_sorted(j));
            set(h(j),'FaceColor',[i/20 1-(i/20) i/20], 'EdgeColor',[i/20 1-(i/20) i/20]);
            hold on;
        end
        j = j + 1;
    end
    i = i + 1;
end
grid on; xlabel('assembly'); ylabel('assembly power (W)'); title('k-means clustered'); set(gca, 'XTickMode','Auto'); set(gca,'YScale','log');

zones = kmedoids(Q,20);
zones_sorted = sort(zones);
subplot(1,5,4); 
i = 1;
while i < 20+1
    j = 1;
    while j < length(Q)+1
        if zones_sorted(j) == i
            h(j) = bar(j,Q_sorted(j));
            set(h(j),'FaceColor',[i/20 1-(i/20) i/20], 'EdgeColor',[i/20 1-(i/20) i/20]);
            hold on;
        end
        j = j + 1;
    end
    i = i + 1;
end
grid on; xlabel('assembly'); ylabel('assembly power (W)'); title('k-medoids clustered'); set(gca, 'XTickMode','Auto'); set(gca,'YScale','log');

zones = clusterdata(Q,'maxclust',20);
zones_sorted = sort(zones);
subplot(1,5,5); 
i = 1;
while i < 20+1
    j = 1;
    while j < length(Q)+1
        if zones_sorted(j) == i
            h(j) = bar(j,Q_sorted(j));
            set(h(j),'FaceColor',[i/20 1-(i/20) i/20], 'EdgeColor',[i/20 1-(i/20) i/20]);
            hold on;
        end
        j = j + 1;
    end
    i = i + 1;
end
grid on; xlabel('assembly'); ylabel('assembly power (W)'); title('euclidian distance clustered'); set(gca, 'XTickMode','Auto'); set(gca,'YScale','log');