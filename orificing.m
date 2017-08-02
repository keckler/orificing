clear variables;

A_flow = 175.84E-4; %flow area per assembly, m^2
adjacent_max_diff = 0.4; %percentage difference in dT between adjacent assemblies
assemblyPowerThreshold = 0.01E6; %minimum power in an assembly to be considered in the optimization, W
cp = 1272; %average heat capacity of coolant, J/kg/K
dP_max = 1E6; %maximum allowable pressure drop over core, Pa, limit taken from Qvist et al
dT_max = 200; %maximum temp increase allowable, C
f_novendstern = 0.021132; %friction factor in the Novendstern friction loss model, see p282 of Waltar
groupingMethod = 'log'; %choose from 'lin', 'log', 'kmeans', 'kmedoids'
H = 3.18; %height of rods, m
nbins = [28]; %vector of number of orifice zones
P_w = 7.2062; %wetted perimeter per assembly, m
powerDetectorFiles = {'~/Downloads/BnB_det0','~/Downloads/BnB_det1','~/Downloads/BnB_det2','~/Downloads/BnB_det3'}; %paths of Serpent detector files with lattice power detectors
rho = 850; %coolant density at lowest point, kg/m^3
T_in = 355; %coolant inlet temperature, C
T_out_bar = 510; %perfectly mixed coolant outlet plenum temperature, C
T_out_bar_tol = 5; %tolerance on perfectly mixed coolant outlet temperature, C
v_max = 12.0; %maximum coolant velocity allowed in assembly, m/s, limit taken from Qvist et al

Q = readQ(powerDetectorFiles);

ncols = sqrt(length(Q)); %size of the square matrix to be formed by the Q vector before entries are removed

[Q, map] = formMap(Q, assemblyPowerThreshold);

Q_ave = sum(Q,2)/length(powerDetectorFiles); %divide by number of steps to get average assembly power over cycle

alpha = Q/cp; %kg-K/s

%optimize for each number of orifice zones
j = 1;
while j < length(nbins)+1
    nb = nbins(j);
    zones = binning(Q_ave, Q, nb, groupingMethod);
    
    [m, cvx_cputime, cvx_optbnd, cvx_optval, cvx_slvitr, cvx_slvtol, cvx_status] = optimizeFlow(alpha, j, T_out_bar, T_out_bar_tol, T_in, dT_max, rho, A_flow, P_w, v_max, f_novendstern, H, dP_max, nbins, zones, map, ncols, adjacent_max_diff);

    %calculate final parameters
    T_out = Q./m/cp + T_in;
    T_out_mixed = sum(T_out.*m)./sum(m);
    dP = (1.29+2.01+0.363+0.41+0.098+0.79+1.32)*m.^2/2/rho/A_flow^2 + (1.0+1.0+0.022+0.0024+0.000082+0.00035+f_novendstern*H/(4*A_flow/P_w))*m.^2/2/rho/A_flow^2 + rho*9.81*H; %pressure drop without orifices
    v = m/rho/A_flow;
    
    %check if bounds are tight
    if isnan(T_out(1,1)) == 1
            display('problem infeasible')
    else
        tol = 0.001; %tolerance for checking
        i = 1;
        while i < length(powerDetectorFiles) + 1 %cycle through all power profiles
            display('------------------')
            disp([i])
            if sum(T_out(:,i) >= (T_in+dT_max)-tol*(T_in+dT_max))
                display('max temp increase constraint is tight!')
            else
                display('max temp increase constraint is loose')
            end

            if T_out_mixed(i) < T_out_bar(1)-T_out_bar_tol+0.1
                display('mixed outlet temp constraint is lower tight!')
            elseif T_out_mixed > T_out_bar(1)+T_out_bar_tol-0.1
                display('mixed outlet temp constraint is upper tight!')
            else
                display('mixed outlet temp constraint is loose')
            end

            if sum(dP(:,i) >= dP_max-tol*dP_max)
                display('max pressure drop constraint is tight!')
            else
                display('max pressure drop constraint is loose')
            end

            if sum(v(:,i) >= v_max-tol*v_max)
                display('max flow velocity constraint is tight!')
            else
                display('max flow velocity constraint is loose')
            end

            if sum(m(:,i) <= 0.01)
                display('min flowrate constraint is tight!')
            else
                display('min flowrate constraint is loose')
            end
            
            display('------------------')
            i = i + 1;
            
        end

    end
    
    %plot T_out, power, flowrate
    figure;
    bar(m); hold on; 
    plot(T_out); ylim([0,1000]); 
    yyaxis right; plot(Q); 
    title(sprintf('k = %i',nbins(j))); grid on; 
    xlabel('assembly'); yyaxis left; ylabel('T_{out} (C), flow (kg/s)'), yyaxis right; ylabel('power (W)');
    drawnow;
    
    %create hexmap of outlet temps for each power profile
    k = 1;
    while k < length(powerDetectorFiles)+1
        i = 1;
        while i < length(map)+1
            r = floor(i/ncols)+1; %row of current assembly in matrix
            c = mod(i,ncols); %column of current assembly in matrix
            if c == 0 %adjust for special case of last column in each row
                r = r-1; c = ncols;
            end

            if isnan(map(i))
                T_out_hexmap(r,c,k) = 0;
            else
                T_out_hexmap(r,c,k) = T_out(map(i),k);
            end

            i = i + 1;
        end
        k = k + 1;
    end
    
    %create hexmap of flowrates
    i = 1;
    while i < length(map)+1
        r = floor(i/ncols)+1; %row of current assembly in matrix
        c = mod(i,ncols); %column of current assembly in matrix
        if c == 0 %adjust for special case of last column in each row
            r = r-1; c = ncols;
        end
        
        if isnan(map(i))
            m_hexmap(r,c) = 0;
        else
            m_hexmap(r,c) = m(map(i));
        end
        
        i = i + 1;
    end
    
    j = j + 1;
end