function [adjacentAssemblies] = findAdjacentAssemblies(nass, map, length_Q)

ncols = sqrt(length_Q); %size of the square matrix to be formed by the Q vector before entries are removed
adjacentAssemblies = zeros(nass, 6); %matrix to demarcate which assemblies are adjacent to each other. row corresponds to assembly number, and columns are -- {1:e, 2:ne, 3:nw, 4:w, 5:sw, 6:se}

k = ncols+1; %skip the first row, as it is assumed to be padded
while k < length(map)-ncols %skip the last row, as it is assumed to be padded
    r = floor(k/ncols)+1; %row of current assembly in matrix
    c = mod(k,ncols); %column of current assembly in matrix
    if c == 0 %adjust for special case of last column in each row
        r = r-1; c = ncols;
    end

    [nwr, nwc, ner, nec, wr, wc, er, ec, swr, swc, ser, sec] = getColumnRow(r, c);
    [center, nw, ne, w, e, sw, se] = getAdjacentNumber(r, c, ncols, nwr, ner, wr, er, swr, ser, nwc, nec, wc, ec, swc, sec);
    [center, nw, ne, w, e, sw, se] = convertToReducedMatrix(map, center, nw, ne, w, e, sw, se);

    %constrain
    if isnan(center) %if current assembly is not in reduced matrix, don't do anything
    else
        if isnan(e) %if east assembly is not in reduced matrix, don't do anything
        else
            adjacentAssemblies(center,1) = e;
        end

        if isnan(ne) %if northeast assembly is not in reduced matrix, don't do anything
        else
            adjacentAssemblies(center,2) = ne;
        end
        
        if isnan(nw) %if northwest assembly is not in reduced matrix, don't do anything
        else
            adjacentAssemblies(center,3) = nw;
        end

        if isnan(w) %if west assembly is not in reduced matrix, don't do anything
        else
            adjacentAssemblies(center,4) = w;
        end

        if isnan(sw) %if southwest assembly is not in reduced matrix, don't do anything
        else
            adjacentAssemblies(center,5) = sw;
        end

        if isnan(se) %if southeast assembly is not in reduced matrix, don't do anything
        else
            adjacentAssemblies(center,6) = se;
        end
    end

    k = k + 1;
end