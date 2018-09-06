function [g,xi,eta,tcaTemps] = fullTCADynamics(volt,eta,xi,dt,tcaTemps)
    %ties together the thermal and dynamics model for the TCAs
    
    %[tcaTemps,centerTemps] = ThermalFromVoltageIterative(volt,dt);
    %using neural network for the temperature computations
    tcaTemps = interpNet(@(x) thermalNet(x')',reshape(volt,3,1),reshape(tcaTemps,3,1),10,dt);
    
    [g,xi,eta] = fastDynamicsStableTCA(tcaTemps,eta,xi,dt);
end