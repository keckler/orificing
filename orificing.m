clear all;

A_flow = 175.84E-4; %flow area per assembly, m^2
assemblyPowerThreshold = 0.01E6; %minimum power in an assembly to be considered in the optimization, W
cp = 1272; %average heat capacity of coolant, J/kg/K
powerDetectorFiles = {'~/Downloads/BnB_det0_31','~/Downloads/BnB_det0_32','~/Downloads/BnB_det0_33'}; %paths of Serpent detector files with lattice power detectors
dP_max = 1E6; %maximum allowable pressure drop over core, Pa, limit taken from Qvist et al
dT_max = 300; %maximum temp increase allowable, C
adjacent_max_diff = 0.4; %percentage difference in dT between adjacent assemblies
f_novendstern = 0.021132; %friction factor in the Novendstern friction loss model, see p282 of Waltar
H = 3.18; %height of rods, m
nbins = [20]; %vector of number of orifice zones
P_w = 7.2062; %wetted perimeter per assembly, m
rho = 850; %coolant density at lowest point, kg/m^3
T_in = 355; %coolant inlet temperature, C
T_out_bar = 510; %perfectly mixed coolant outlet plenum temperature, C
T_out_bar_tol = 5; %tolerance on perfectly mixed coolant outlet temperature, C
v_max = 12.0; %maximum coolant velocity allowed in assembly, m/s, limit taken from Qvist et al

Q = [];
for file = powerDetectorFiles
    %read in detector results
    run(file{1})

    %axially integrate assembly powers
    Q_step = DETAssemblyPowerAxial1(:,11)+DETAssemblyPowerAxial2(:,11)+DETAssemblyPowerAxial3(:,11)+DETAssemblyPowerAxial4(:,11)+DETAssemblyPowerAxial5(:,11)+DETAssemblyPowerAxial6(:,11)+DETAssemblyPowerAxial7(:,11)+DETAssemblyPowerAxial8(:,11)+DETAssemblyPowerAxial9(:,11); %fission power in each assembly, W
    Q_gamma_step = DETAssemblyGammaPowerAxial1(:,11)+DETAssemblyGammaPowerAxial2(:,11)+DETAssemblyGammaPowerAxial3(:,11)+DETAssemblyGammaPowerAxial4(:,11)+DETAssemblyGammaPowerAxial5(:,11)+DETAssemblyGammaPowerAxial6(:,11)+DETAssemblyGammaPowerAxial7(:,11)+DETAssemblyGammaPowerAxial8(:,11)+DETAssemblyGammaPowerAxial9(:,11); %gamma power in each assembly, W
    
    %add gamma power into total power
    totalPower = sum(Q_step);
    totalGammaPower = sum(Q_gamma_step);
    fissionFrac = totalPower/(totalPower+totalGammaPower);
    Q_step = fissionFrac*Q_step; %scale total power by the fraction of power from fission, since this is the tally
    Q_step = Q_step + Q_gamma_step; %add in the gamma power
    
    %add assembly powers in current depletion step to running total
    Q(:,end+1) = Q_step;
end

ncols = sqrt(length(Q)); %size of the square matrix to be formed by the Q vector before entries are removed

%remove zero entries and form map
map = []; %i'th entry corresponds to index of reduced Q vector
Q_new = [];
i = 1;
j = 1;
while i < length(Q(:,1))
    if Q(i,1) > assemblyPowerThreshold
        Q_new(end+1,:) = Q(i,1:end);
        map(i) = j;
        j = j + 1;
    else
        map(i) = NaN;
    end
    i = i + 1;
end
Q = Q_new;

Q_ave = sum(Q,2)/length(powerDetectorFiles); %divide by number of steps to get average assembly power over cycle

alpha = Q/cp; %kg-K/s

%optimize for each number of orifice zones
j = 1;
while j < length(nbins)+1
    zones = discretize(Q_ave,logspace(log10(min(min(Q))),log10(max(max(Q))),nbins(j))); %bin it up
    
    cvx_precision best
    cvx_begin
        variables m(size(Q)) x(nbins(j)) %inverse flowrate, dummies
        
        i = 1;
        obj = 0;
        while i < length(Q(1,:))+1
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
        while i < length(Q(1,:))+1
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
            
            %get row and column of each adjacent assembly
            nwr = r-1; nwc = c; %northwest row, column
            ner = r-1; nec = c+1; %northeast row, column
            wr = r; wc = c-1; %west row, column
            er = r; ec = c+1; %east row, column
            swr = r+1; swc = c-1; %southwest row, column
            ser = r+1; sec = c; %southeast row, column
            
            %get assembly number of each adjacent assembly in original matrix
            center = (r-1)*ncols+c;
            nw = (nwr-1)*ncols+nwc;
            ne = (ner-1)*ncols+nec;
            w = (wr-1)*ncols+wc;
            e = (er-1)*ncols+ec;
            sw = (swr-1)*ncols+swc;
            se = (ser-1)*ncols+sec;
            
            %convert assembly number to reduced matrix
            center = map(center);
            nw = map(nw);
            ne = map(ne);
            w = map(w);
            e = map(e);
            sw = map(sw);
            se = map(se);
            
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

    %calculate final parameters
    T_out = Q./m/cp + T_in;
    T_out_mixed = sum(T_out.*m)./sum(m);
    dP = (1.29+2.01+0.363+0.41+0.098+0.79+1.32)*m.^2/2/rho/A_flow^2 + (1.0+1.0+0.022+0.0024+0.000082+0.00035+f_novendstern*H/(4*A_flow/P_w))*m.^2/2/rho/A_flow^2 + rho*9.81*H; %pressure drop without orifices
    v = m/rho/A_flow;
    
    %check if bounds are tight
    tol = 0.001; %tolerance for checking
    
    display('------------------')
    
    if isnan(T_out(1)) == 1
        display('problem infeasible')
    else
        if sum(T_out >= (T_in+dT_max)-tol*(T_in+dT_max))
            display('max temp increase constraint is tight!')
        else
            display('max temp increase constraint is loose')
        end
        
        if T_out_mixed < T_out_bar(1)-T_out_bar_tol+0.1
            display('mixed outlet temp constraint is lower tight!')
        elseif T_out_mixed > T_out_bar(1)+T_out_bar_tol-0.1
            display('mixed outlet temp constraint is upper tight!')
        else
            display('mixed outlet temp constraint is loose')
        end
    
        if sum(dP >= dP_max-tol*dP_max)
            display('max pressure drop constraint is tight!')
        else
            display('max pressure drop constraint is loose')
        end
    
        if sum(v >= v_max-tol*v_max)
            display('max flow velocity constraint is tight!')
        else
            display('max flow velocity constraint is loose')
        end
    
        if sum(m <= 0.01)
            display('min flowrate constraint is tight!')
        else
            display('min flowrate constraint is loose')
        end
    end
    
    display('------------------')
    
    %plot T_out, power, flowrate
    figure;
    bar(m); hold on; 
    plot(T_out); ylim([0,1000]); 
    yyaxis right; plot(Q); 
    legend('flow','T_{out}','power'); title(sprintf('k = %i',nbins(j))); grid on; 
    xlabel('assembly'); yyaxis left; ylabel('T_{out} (C), flow (kg/s)'), yyaxis right; ylabel('power (W)');
    drawnow;
    
    %create hexmap of outlet temps
    i = 1;
    while i < length(map)+1
        r = floor(i/ncols)+1; %row of current assembly in matrix
        c = mod(i,ncols); %column of current assembly in matrix
        if c == 0 %adjust for special case of last column in each row
            r = r-1; c = ncols;
        end
        
        if isnan(map(i))
            T_out_hexmap(r,c) = 0;
        else
            T_out_hexmap(r,c) = T_out(map(i));
        end
        
        i = i + 1;
    end
    
    j = j + 1;
end