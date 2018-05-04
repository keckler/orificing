function Q = readQ(powerDetectorFiles)

Q = [];
for file = powerDetectorFiles
    %read in detector results
    run(file{1})

    %axially integrate assembly powers
    Q_step = DETAssemblyPowerAxial1(:,11)+DETAssemblyPowerAxial2(:,11)+DETAssemblyPowerAxial3(:,11)+DETAssemblyPowerAxial4(:,11)+DETAssemblyPowerAxial5(:,11)+DETAssemblyPowerAxial6(:,11)+DETAssemblyPowerAxial7(:,11)+DETAssemblyPowerAxial8(:,11)+DETAssemblyPowerAxial9(:,11)+DETAssemblyPowerAxial10(:,11)+DETAssemblyPowerAxial11(:,11)+DETAssemblyPowerAxial12(:,11)+DETAssemblyPowerAxial13(:,11)+DETAssemblyPowerAxial14(:,11)+DETAssemblyPowerAxial15(:,11)+DETAssemblyPowerAxial16(:,11)+DETAssemblyPowerAxial17(:,11)+DETAssemblyPowerAxial18(:,11)+DETAssemblyPowerAxial19(:,11)+DETAssemblyPowerAxial20(:,11)+DETAssemblyPowerAxial21(:,11)+DETAssemblyPowerAxial22(:,11)+DETAssemblyPowerAxial23(:,11)+DETAssemblyPowerAxial24(:,11); %fission power in each assembly, W
    Q_gamma_step = DETAssemblyGammaPowerAxial1(:,11)+DETAssemblyGammaPowerAxial2(:,11)+DETAssemblyGammaPowerAxial3(:,11)+DETAssemblyGammaPowerAxial4(:,11)+DETAssemblyGammaPowerAxial5(:,11)+DETAssemblyGammaPowerAxial6(:,11)+DETAssemblyGammaPowerAxial7(:,11)+DETAssemblyGammaPowerAxial8(:,11)+DETAssemblyGammaPowerAxial9(:,11)+DETAssemblyGammaPowerAxial10(:,11)+DETAssemblyGammaPowerAxial11(:,11)+DETAssemblyGammaPowerAxial12(:,11)+DETAssemblyGammaPowerAxial13(:,11)+DETAssemblyGammaPowerAxial14(:,11)+DETAssemblyGammaPowerAxial15(:,11)+DETAssemblyGammaPowerAxial16(:,11)+DETAssemblyGammaPowerAxial17(:,11)+DETAssemblyGammaPowerAxial18(:,11)+DETAssemblyGammaPowerAxial19(:,11)+DETAssemblyGammaPowerAxial20(:,11)+DETAssemblyGammaPowerAxial21(:,11)+DETAssemblyGammaPowerAxial22(:,11)+DETAssemblyGammaPowerAxial23(:,11)+DETAssemblyGammaPowerAxial24(:,11); %gamma power in each assembly, W
    
    %add gamma power into total power
    totalPower = sum(Q_step);
    totalGammaPower = sum(Q_gamma_step);
    fissionFrac = totalPower/(totalPower+totalGammaPower);
    Q_step = fissionFrac*Q_step; %scale total power by the fraction of power from fission, since this is the tally
    Q_step = Q_step + Q_gamma_step; %add in the gamma power
    
    %add assembly powers in current depletion step to running total
    Q(:,end+1) = Q_step;
end