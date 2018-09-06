function [xi,eta,g] = initialDynamics(n)
    L = 4e-2;
    g = [ones(1,n);zeros(3,n);ones(1,n);zeros(3,n);ones(1,n);zeros(2,n);linspace(0,L,n)];
    xi = [zeros(5,n);ones(1,n)];
    eta = zeros(6,n);
end