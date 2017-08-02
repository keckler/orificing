function [] = checkBoundsTightness(T_out, powerDetectorFiles, T_in, dT_max, T_out_mixed, T_out_bar, T_out_bar_tol, dP, dP_max, v, v_max, m)

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