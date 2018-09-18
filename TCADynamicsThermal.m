function [g,xi,eta,tcaTemps,thermalResults] = TCADynamicsThermal(q, eta_prev, xi_prev, dt, thermalmodel, msh, thermalResults)
    %the TCA driven manipulator using the actual thermal model
    
    %first advance the thermal model then do the dynamics step
    
    [tcaTemps,centerTemps,thermalResults] = ThermalFromVoltageIterative(thermalmodel,msh,q,dt,thermalResults);
    
    [g,xi,eta] = fastDynamicsStableTCA(tcaTemps,eta_prev,xi_prev,dt);
    
end