function [T_out, T_out_mixed, dP, v] = calculateFinalParameters(Q, m, cp, T_in, rho, A_flow, f_novendstern, H, P_w)
    
T_out = Q./m/cp + T_in;
T_out_mixed = sum(T_out.*m)./sum(m);
dP = (1.29+2.01+0.363+0.41+0.098+0.79+1.32)*m.^2/2/rho/A_flow^2 + (1.0+1.0+0.022+0.0024+0.000082+0.00035+f_novendstern*H/(4*A_flow/P_w))*m.^2/2/rho/A_flow^2 + rho*9.81*H; %pressure drop without orifices
v = m/rho/A_flow;