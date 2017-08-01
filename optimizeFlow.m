function [m, cvx_cputime, cvx_optbnd, cvx_optval, cvx_slvitr, cvx_slvtol, cvx_status] = optimizeFlow(alpha, j, T_out_bar, T_out_bar_tol, T_in, dT_max, rho, A_flow, P_w, v_max, f_novendstern, H, dP_max, nbins, zones, map, ncols, adjacent_max_diff)

cvx_precision best
cvx_begin
    variables m(size(alpha)) x(nbins(j)) %inverse flowrate, dummies

    i = 1;
    obj = 0;
    while i < length(alpha(1,:))+1
        obj = obj + var(m(:,i)./alpha(:,i));
        i = i + 1;
    end

    minimize( obj )

    %constraints
    sum(m) >= sum(alpha)./(T_out_bar+T_out_bar_tol-T_in) %weighted outlet temp is less than or equal to T_out_bar plus a tolerance
    sum(m) <= sum(alpha)./(T_out_bar-T_out_bar_tol-T_in) %weighted outlet temp is greater than or equal to T_out_bar minus a tolerance
    m >= 0 %all flows positive
    m./alpha >= 1/dT_max %maximum temperature change over channel
    m/rho/A_flow <= v_max %maximum flow velocity in assembly
    (1.29+2.01+0.363+0.41+0.098+0.79+1.32)*m.^2/2/rho/A_flow^2 + (1.0+1.0+0.022+0.0024+0.000082+0.00035+f_novendstern*H/(4*A_flow/P_w))*m.^2/2/rho/A_flow^2 + rho*9.81*H <= dP_max %maximum pressure drop in assembly, form+friction+gravity, using generic values for form and friction losses from p281 of Waltar

    %constraints so flow in each assembly remains constant over all steps
    i = 2;
    while i < length(alpha(1,:))+1
        m(:,1) == m(:,i)
        i = i + 1;
    end

    %constraints to remain in bins
    k = 1;
    while k < nbins(j)+1
        zp = zones - k;
        l = 1;
        while l < length(zones)+1
            if zp(l) == 0
                m(l,1) == x(k)
            end
            l = l + 1;
        end
        k = k + 1;
    end

    %constraints on adjacent assemblies
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
            if isnan(nw) %if northwest assembly is not in reduced matrix, don't do anything
            else
                m(center,:)./alpha(center,:) >= 1/(1+adjacent_max_diff)*m(nw,:)./alpha(nw,:)
                m(center,:)./alpha(center,:) <= 1/(1-adjacent_max_diff)*m(nw,:)./alpha(nw,:)
            end

            if isnan(ne) %if northeast assembly is not in reduced matrix, don't do anything
            else
                m(center,:)./alpha(center,:) >= 1/(1+adjacent_max_diff)*m(ne,:)./alpha(ne,:)
                m(center,:)./alpha(center,:) <= 1/(1-adjacent_max_diff)*m(ne,:)./alpha(ne,:)
            end

            if isnan(w) %if west assembly is not in reduced matrix, don't do anything
            else
                m(center,:)./alpha(center,:) >= 1/(1+adjacent_max_diff)*m(w,:)./alpha(w,:)
                m(center,:)./alpha(center,:) <= 1/(1-adjacent_max_diff)*m(w,:)./alpha(w,:)
            end

            if isnan(e) %if east assembly is not in reduced matrix, don't do anything
            else
                m(center,:)./alpha(center,:) >= 1/(1+adjacent_max_diff)*m(e,:)./alpha(e,:)
                m(center,:)./alpha(center,:) <= 1/(1-adjacent_max_diff)*m(e,:)./alpha(e,:)
            end

            if isnan(sw) %if southwest assembly is not in reduced matrix, don't do anything
            else
                m(center,:)./alpha(center,:) >= 1/(1+adjacent_max_diff)*m(sw,:)./alpha(sw,:)
                m(center,:)./alpha(center,:) <= 1/(1-adjacent_max_diff)*m(sw,:)./alpha(sw,:)
            end

            if isnan(se) %if southeast assembly is not in reduced matrix, don't do anything
            else
                m(center,:)./alpha(center,:) >= 1/(1+adjacent_max_diff)*m(se,:)./alpha(se,:)
                m(center,:)./alpha(center,:) <= 1/(1-adjacent_max_diff)*m(se,:)./alpha(se,:)
            end
        end

        k = k + 1;
    end

cvx_end