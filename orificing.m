clear variables;

A_flow = 123.3184E-4; %flow area per assembly, m^2
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
powerDetectorFiles = {'~/Documents/work/ARC/serpent/BnB/BnB_det0','~/Documents/work/ARC/serpent/BnB/BnB_det1','~/Documents/work/ARC/serpent/BnB/BnB_det2','~/Documents/work/ARC/serpent/BnB/BnB_det3'}; %paths of Serpent detector files with lattice power detectors
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

    [T_out, T_out_mixed, dP, v] = calculateFinalParameters(Q, m, cp, T_in, rho, A_flow, f_novendstern, H, P_w);
    
    checkBoundsTightness(T_out, powerDetectorFiles, T_in, dT_max, T_out_mixed, T_out_bar, T_out_bar_tol, dP, dP_max, v, v_max, m);
    
    plotResults(m, T_out, Q, nbins, j);
    
    T_out_hexmap = generateToutHexmap(powerDetectorFiles, map, ncols, T_out);
    
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