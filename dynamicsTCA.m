function [g,xi,eta,tcaTemps] = dynamicsTCA(q,eta,xi,dt,tcaTemps)
    tcaTemps = interpNet(@(x) thermalNet(x')',reshape(q,3,1),reshape(tcaTemps,3,1),10,dt);
    [g,xi,eta] = manip_dynamics(1,tcaTemps,eta,xi,dt);
end