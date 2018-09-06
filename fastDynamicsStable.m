function [g,xi,eta] = fastDynamicsStable(F,eta_prev,xi_prev,dt)
    % have to provide a guess for both the initial xi and most of eta, just
    % pass in eta with a slight modification, the other code should handle
    % it
    options = optimoptions('fsolve','Display','off','UseParallel',false);
%     eta_guess = eta_prev+dt*eta_der;
    eta_guess = eta_prev;
%     eta_guess = (4/3)*eta_prev-(1/3)*eta_prev_prev+(2/3)*dt*eta_der;
    vals = fsolve(@(x) newDynamicsHelper(@(eta) dynamicsStable(F,eta,eta_prev,xi_prev,dt),x),[xi_prev(1,:);eta_guess(2:end,:)],options);
    [y,bcs] = dynamicsStable(F,vals,eta_prev,xi_prev,dt);
%     eta = fsolve(@(x) newDynamicsHelper(@(eta) dynamicsEtaTCA(F,eta,eta_prev,eta_acc,dt),x),eta_prev,options);
%     [y,bcs] = dynamicsEtaTCA(F,eta,eta_prev,eta_acc,dt);
    bcs;
    y;
    g = y(1:12,:)';
    xi = y(13:18,:)';
    eta = y(19:24,:)';
end