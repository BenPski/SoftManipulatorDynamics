function [g,xi,eta,tcaTemps] = initDynamics(n)
    %a better version for initializing the dynamics
    %provides all the possible inits, so for the different dynamics request
    %different different outputs
    
    L = 4e-2; %would like to abstract this away eventually, or get this whole dynamics in its own object
    g = [ones(n,1),zeros(n,3),ones(n,1),zeros(n,3),ones(n,1),zeros(n,2),linspace(0,L,n)'];
    xi = [zeros(n,5),ones(n,1)];
    eta = zeros(n,6);
    tcaTemps = [25;25;25];
    %only works on higher versions of matlab
    %[thermalModel,msh,thermalResults] = generateThermalModel();
end
