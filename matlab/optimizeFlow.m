function [m, cvx_cputime, cvx_optbnd, cvx_optval, cvx_slvitr, cvx_slvtol, cvx_status] = optimizeFlow(alpha, j, T_out_bar, T_out_bar_tol, T_in, dT_max, rho, A_flow, P_w, v_max, f_novendstern, H, dP_max, nbins, zones, map, ncols, adjacent_max_diff, fm)

cvx_precision best
cvx_begin
    variables m(size(alpha)) x(nbins(j)) %inverse flowrate, dummies

    fprintf(fm, 'var m{i in 1..%i, j in 1..%i};\n', length(alpha(:,1)), length(alpha(1,:)));
    fprintf(fm, 'var x{i in 1..%i} >= 0;\n', nbins(j));
    fprintf(fm, 'var step1ave;\n');
    fprintf(fm, 'var step2ave;\n');
    fprintf(fm, 'var step3ave;\n');
    fprintf(fm, 'var step4ave;\n');
    
    fprintf(fm, 'minimize variance: (sum{i in 1..%i} ( (m[i,1]/alpha[i,1]) - step1ave)^2) + (sum{i in 1..%i} ( (m[i,2]/alpha[i,2]) - step2ave)^2) + (sum{i in 1..%i} ( (m[i,3]/alpha[i,3]) - step3ave)^2) + (sum{i in 1..%i} ( (m[i,4]/alpha[i,4]) - step4ave)^2);\n', length(m(:,1)), length(m(:,1)), length(m(:,1)), length(m(:,1)));
    
    fprintf(fm, 'subject to s1ave: step1ave = (sum{i in 1..%i} m[i,1]/alpha[i,1])/%i;\n', length(m(:,1)), length(m(:,1)));
    fprintf(fm, 'subject to s2ave: step2ave = (sum{i in 1..%i} m[i,2]/alpha[i,2])/%i;\n', length(m(:,2)), length(m(:,2)));
    fprintf(fm, 'subject to s3ave: step3ave = (sum{i in 1..%i} m[i,3]/alpha[i,3])/%i;\n', length(m(:,3)), length(m(:,3)));
    fprintf(fm, 'subject to s4ave: step4ave = (sum{i in 1..%i} m[i,4]/alpha[i,4])/%i;\n', length(m(:,4)), length(m(:,4)));
    
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

    fprintf(fm, 'subject to mixedUpper1: (sum{i in 1..%i} m[i,1]) >= (sum{i in 1..%i} alpha[i,1])/(T_out_bar + T_out_bar_tol - T_in);\n', length(m(:,1)), length(m(:,1)));
    fprintf(fm, 'subject to mixedUpper2: (sum{i in 1..%i} m[i,2]) >= (sum{i in 1..%i} alpha[i,2])/(T_out_bar + T_out_bar_tol - T_in);\n', length(m(:,2)), length(m(:,2)));
    fprintf(fm, 'subject to mixedUpper3: (sum{i in 1..%i} m[i,3]) >= (sum{i in 1..%i} alpha[i,3])/(T_out_bar + T_out_bar_tol - T_in);\n', length(m(:,3)), length(m(:,3)));
    fprintf(fm, 'subject to mixedUpper4: (sum{i in 1..%i} m[i,4]) >= (sum{i in 1..%i} alpha[i,4])/(T_out_bar + T_out_bar_tol - T_in);\n', length(m(:,4)), length(m(:,4)));
    fprintf(fm, 'subject to mixedLower1: (sum{i in 1..%i} m[i,1]) <= (sum{i in 1..%i} alpha[i,1])/(T_out_bar - T_out_bar_tol - T_in);\n', length(m(:,1)), length(m(:,1)));
    fprintf(fm, 'subject to mixedLower2: (sum{i in 1..%i} m[i,2]) <= (sum{i in 1..%i} alpha[i,2])/(T_out_bar - T_out_bar_tol - T_in);\n', length(m(:,2)), length(m(:,2)));
    fprintf(fm, 'subject to mixedLower3: (sum{i in 1..%i} m[i,3]) <= (sum{i in 1..%i} alpha[i,3])/(T_out_bar - T_out_bar_tol - T_in);\n', length(m(:,3)), length(m(:,3)));
    fprintf(fm, 'subject to mixedLower4: (sum{i in 1..%i} m[i,4]) <= (sum{i in 1..%i} alpha[i,4])/(T_out_bar - T_out_bar_tol - T_in);\n', length(m(:,4)), length(m(:,4)));
    fprintf(fm, 'subject to maxincrease1 {i in 1..%i}: m[i,1]/alpha[i,1] >= 1/dT_max;\n', length(m(:,1)));
    fprintf(fm, 'subject to maxincrease2 {i in 1..%i}: m[i,2]/alpha[i,2] >= 1/dT_max;\n', length(m(:,2)));
    fprintf(fm, 'subject to maxincrease3 {i in 1..%i}: m[i,3]/alpha[i,3] >= 1/dT_max;\n', length(m(:,3)));
    fprintf(fm, 'subject to maxincrease4 {i in 1..%i}: m[i,4]/alpha[i,4] >= 1/dT_max;\n', length(m(:,4)));
    fprintf(fm, 'subject to maxflow1 {i in 1..%i}: m[i,1]/(rho*A_flow) <= v_max;\n', length(m(:,1)));
    fprintf(fm, 'subject to maxflow2 {i in 1..%i}: m[i,2]/(rho*A_flow) <= v_max;\n', length(m(:,2)));
    fprintf(fm, 'subject to maxflow3 {i in 1..%i}: m[i,3]/(rho*A_flow) <= v_max;\n', length(m(:,3)));
    fprintf(fm, 'subject to maxflow4 {i in 1..%i}: m[i,4]/(rho*A_flow) <= v_max;\n', length(m(:,4)));
    %fprintf(fm, 'subject to pdrop1 {i in 1..%i}: (1.29+2.01+0.363+0.41+0.098+0.79+1.32)*(m[i,1]^2)/(2*rho*A_flow^2) + (1.0+1.0+0.022+0.0024+0.000082+0.00035+f_novendstern*H/(4*A_flow/P_w))*(m[i,1]^2)/(2*rho*A_flow^2) + rho*9.81*H <= dP_max;\n', length(m(:,1)));
    %fprintf(fm, 'subject to pdrop2 {i in 1..%i}: (1.29+2.01+0.363+0.41+0.098+0.79+1.32)*(m[i,2]^2)/(2*rho*A_flow^2) + (1.0+1.0+0.022+0.0024+0.000082+0.00035+f_novendstern*H/(4*A_flow/P_w))*(m[i,2]^2)/(2*rho*A_flow^2) + rho*9.81*H <= dP_max;\n', length(m(:,2)));
    %fprintf(fm, 'subject to pdrop3 {i in 1..%i}: (1.29+2.01+0.363+0.41+0.098+0.79+1.32)*(m[i,3]^2)/(2*rho*A_flow^2) + (1.0+1.0+0.022+0.0024+0.000082+0.00035+f_novendstern*H/(4*A_flow/P_w))*(m[i,3]^2)/(2*rho*A_flow^2) + rho*9.81*H <= dP_max;\n', length(m(:,3)));
    %fprintf(fm, 'subject to pdrop4 {i in 1..%i}: (1.29+2.01+0.363+0.41+0.098+0.79+1.32)*(m[i,4]^2)/(2*rho*A_flow^2) + (1.0+1.0+0.022+0.0024+0.000082+0.00035+f_novendstern*H/(4*A_flow/P_w))*(m[i,4]^2)/(2*rho*A_flow^2) + rho*9.81*H <= dP_max;\n', length(m(:,4)));
    
    %constraints so flow in each assembly remains constant over all steps
    i = 2;
    while i < length(alpha(1,:))+1
        m(:,1) == m(:,i)
        
        fprintf(fm, 'subject to flowConstantOverSteps%i {i in 1..%i}: m[i,1] = m[i,%i];\n', i, length(m(:,1)), i);
        
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
                
                fprintf(fm, 'subject to chan%iinGroup%i: m[%i,1] = x[%i];\n', l, k, l, k);
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
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: m[%i,i]/alpha[%i,i] >= (1/(1+adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, nw, length(m(1,:)), center, center, nw, nw);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: m[%i,i]/alpha[%i,i] <= (1/(1-adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, nw, length(m(1,:)), center, center, nw, nw);
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, nw, length(m(1,:)), center, center, nw, nw);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, nw, length(m(1,:)), nw, nw, center, center);
                
                fprintf(fm, 'subject to channel%iAdjacent%i_1 {k in 1..%i}: (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) - (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, nw, length(m(1,:)), center, center, nw, nw);
                fprintf(fm, 'subject to channel%iAdjacent%i_2 {k in 1..%i}: -(sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) + (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, nw, length(m(1,:)), center, center, nw, nw);
            end

            if isnan(ne) %if northeast assembly is not in reduced matrix, don't do anything
            else
                m(center,:)./alpha(center,:) >= 1/(1+adjacent_max_diff)*m(ne,:)./alpha(ne,:)
                m(center,:)./alpha(center,:) <= 1/(1-adjacent_max_diff)*m(ne,:)./alpha(ne,:)
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: m[%i,i]/alpha[%i,i] >= (1/(1+adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, ne, length(m(1,:)), center, center, ne, ne);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: m[%i,i]/alpha[%i,i] <= (1/(1-adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, ne, length(m(1,:)), center, center, ne, ne);
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, ne, length(m(1,:)), center, center, ne, ne);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, ne, length(m(1,:)), ne, ne, center, center);
                
                fprintf(fm, 'subject to channel%iAdjacent%i_1 {k in 1..%i}: (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) - (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, ne, length(m(1,:)), center, center, ne, ne);
                fprintf(fm, 'subject to channel%iAdjacent%i_2 {k in 1..%i}: -(sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) + (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, ne, length(m(1,:)), center, center, ne, ne);
            end

            if isnan(w) %if west assembly is not in reduced matrix, don't do anything
            else
                m(center,:)./alpha(center,:) >= 1/(1+adjacent_max_diff)*m(w,:)./alpha(w,:)
                m(center,:)./alpha(center,:) <= 1/(1-adjacent_max_diff)*m(w,:)./alpha(w,:)
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: m[%i,i]/alpha[%i,i] >= (1/(1+adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, w, length(m(1,:)), center, center, w, w);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: m[%i,i]/alpha[%i,i] <= (1/(1-adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, w, length(m(1,:)), center, center, w, w);
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, w, length(m(1,:)), center, center, w, w);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, w, length(m(1,:)), w, w, center, center);
                
                fprintf(fm, 'subject to channel%iAdjacent%i_1 {k in 1..%i}: (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) - (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, w, length(m(1,:)), center, center, w, w);
                fprintf(fm, 'subject to channel%iAdjacent%i_2 {k in 1..%i}: -(sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) + (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, w, length(m(1,:)), center, center, w, w);
            end

            if isnan(e) %if east assembly is not in reduced matrix, don't do anything
            else
                m(center,:)./alpha(center,:) >= 1/(1+adjacent_max_diff)*m(e,:)./alpha(e,:)
                m(center,:)./alpha(center,:) <= 1/(1-adjacent_max_diff)*m(e,:)./alpha(e,:)
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: m[%i,i]/alpha[%i,i] >= (1/(1+adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, e, length(m(1,:)), center, center, e, e);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: m[%i,i]/alpha[%i,i] <= (1/(1-adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, e, length(m(1,:)), center, center, e, e);
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, e, length(m(1,:)), center, center, e, e);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, e, length(m(1,:)), e, e, center, center);
                
                fprintf(fm, 'subject to channel%iAdjacent%i_1 {k in 1..%i}: (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) - (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, e, length(m(1,:)), center, center, e, e);
                fprintf(fm, 'subject to channel%iAdjacent%i_2 {k in 1..%i}: -(sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) + (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, e, length(m(1,:)), center, center, e, e);
            end

            if isnan(sw) %if southwest assembly is not in reduced matrix, don't do anything
            else
                m(center,:)./alpha(center,:) >= 1/(1+adjacent_max_diff)*m(sw,:)./alpha(sw,:)
                m(center,:)./alpha(center,:) <= 1/(1-adjacent_max_diff)*m(sw,:)./alpha(sw,:)
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: m[%i,i]/alpha[%i,i] >= (1/(1+adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, sw, length(m(1,:)), center, center, sw, sw);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: m[%i,i]/alpha[%i,i] <= (1/(1-adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, sw, length(m(1,:)), center, center, sw, sw);
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, sw, length(m(1,:)), center, center, sw, sw);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, sw, length(m(1,:)), sw, sw, center, center);
                
                fprintf(fm, 'subject to channel%iAdjacent%i_1 {k in 1..%i}: (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) - (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, sw, length(m(1,:)), center, center, sw, sw);
                fprintf(fm, 'subject to channel%iAdjacent%i_2 {k in 1..%i}: -(sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) + (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, sw, length(m(1,:)), center, center, sw, sw);
            end

            if isnan(se) %if southeast assembly is not in reduced matrix, don't do anything
            else
                m(center,:)./alpha(center,:) >= 1/(1+adjacent_max_diff)*m(se,:)./alpha(se,:)
                m(center,:)./alpha(center,:) <= 1/(1-adjacent_max_diff)*m(se,:)./alpha(se,:)
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: m[%i,i]/alpha[%i,i] >= (1/(1+adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, se, length(m(1,:)), center, center, se, se);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: m[%i,i]/alpha[%i,i] <= (1/(1-adjacent_max_diff))*(m[%i,i]/alpha[%i,i]);\n', center, se, length(m(1,:)), center, center, se, se);
                
                %fprintf(fm, 'subject to channel%iAdjacent%iGreater {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, se, length(m(1,:)), center, center, se, se);
                %fprintf(fm, 'subject to channel%iAdjacent%iLess {i in 1..%i}: log(alpha[%i,i]) - mu[%i] - log(alpha[%i,i]) + mu[%i] <= log(adjacent_max_diff);\n', center, se, length(m(1,:)), se, se, center, center);
                
                fprintf(fm, 'subject to channel%iAdjacent%i_1 {k in 1..%i}: (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) - (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, se, length(m(1,:)), center, center, se, se);
                fprintf(fm, 'subject to channel%iAdjacent%i_2 {k in 1..%i}: -(sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) + (sum{j in 1..n} delta[%i,j]*Omega[%i,j,k]) <= adjacent_max_diff;\n', center, se, length(m(1,:)), center, center, se, se);
            end
        end

        k = k + 1;
    end

cvx_end