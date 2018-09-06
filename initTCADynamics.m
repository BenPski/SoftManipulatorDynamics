function [g,xi,eta,tcaTemps] = initTCADynamics(n)
    [xi0,eta0,g0] = initialDynamics(n);
    xi=transpose(xi0);
    eta=transpose(eta0);
    g=transpose(g0);
    
%     global thermalmodel;
%     global msh;
%     [thermalmodel,msh] = generateThermalModel();
%     global thermalResults;
%     thermalResults = 25;

    %if using neural network for thermal
    
    tcaTemps = [25;25;25];
    
end