function T = Model_Solve( x , t , D1 , D2 , L1 , L2 )
% Author: Sean McInerney; Last Update: 17/05/2018
% Solves a two layer diffusion problem, given x, t and appropriate
% model parameters. Uses numerical inverse Laplace code, cf.m, on an
% analytical Laplace solution to the problem, as forumlated in
% Laplace_Temperature.m
%
% Solution to Equations (2.11)-(2.18)
%
% FUNCTION INPUTS:
% x - Depth
% t - Time
% D1 - Diffusivity in skin layer
% D2 - Diffusivity in fat layer
% L1 - Depth of layer interface
% L2 - Depth of entire system
% FUNCTION OUTPUT:
% T - Model solution discretised as a matrix [nx by nt]

% Used for inverse Laplace
n = 14;
[z,c] = cf(n);

% Initialise T
N = length(x);
M = length(t);
T = zeros(N,M);

% TBar(S,X) formulates the Laplace solution
TBar = @(S,X) Laplace_Temperature(S,X,D1,D2,L1,L2);
% Loop over space and time
for i = 1:N
    for j = 1:M
        % Check if initial condition needs to be applied
        if t(j) == 0
            T(i,j) = 0;
        else
            % Perform inverse Laplace
            T(i,j) = 0;
            for k = 1:n/2
                I = 2*k-1;
                s = z(I)/t(j);
                T(i,j) = T(i,j) - c(I)*TBar(s,x(i))/t(j);
            end
            T(i,j) = 2*real(T(i,j));
        end
    end
end

end