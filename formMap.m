function [Q, map] = formMap(Q, assemblyPowerThreshold)

%remove zero entries and form map
map = []; %i'th entry corresponds to index of reduced Q vector
Q_new = [];
i = 1;
j = 1;
while i < length(Q(:,1))
    if max(Q(i,:)) > assemblyPowerThreshold
        Q_new(end+1,:) = Q(i,1:end);
        map(i) = j;
        j = j + 1;
    else
        map(i) = NaN;
    end
    i = i + 1;
end
Q = Q_new;