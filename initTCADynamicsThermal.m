function [g,xi,eta,tcaTemps,thermalModel,msh,thermalResults] = initTCADynamicsThermal(n)
    %initialize the TCA dynamics and use the original thermal model
    
    [g,xi,eta,tcaTemps] = initTCADynamics(n);
    [thermalModel,msh,thermalResults] = generateThermalModel();
end