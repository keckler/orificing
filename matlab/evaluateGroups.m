function [pf, pfgroup] = evaluateGroups(zones, Q)

zs = sort(zones);
Qs = sort(Q);

pf = 0.0;
pfgroup = 0;

i = 1;
j = 1;
while i < max(zs)+1
    idx = j;
    while i == zs(j) && j < length(zs)
        idy = j;
        j = j + 1;
    end
    
    pf_tmp = Qs(idy)/Qs(idx);
    if pf_tmp > pf
        pf = pf_tmp;
        pfgroup = i;
    end
    
    i = i + 1;
end