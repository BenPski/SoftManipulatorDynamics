function [g,xi,eta,tcaTemps] = dynamicsTCA(q,eta,xi,dt,tcaTemps)
    temps_prev = tcaTemps;
    tcaTemps = interpNet(@(x) thermalNet(x')',reshape(q,3,1),reshape(tcaTemps,3,1),10,dt);
    temps_dot = (tcaTemps-temps_prev)/dt;
    [g,xi,eta] = manip_dynamics(1,tcaTemps,temps_dot,eta,xi,dt);
end