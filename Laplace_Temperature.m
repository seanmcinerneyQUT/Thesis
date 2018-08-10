function TBar = Laplace_Temperature( s , x , D1 , D2 , L1 , L2 )
% Author: Sean McInerney; Last Update: 17/05/2018
% Given s and x and the appropriate parameters, Laplace_Temperature returns
% the Laplace solution to the following two-layer thermal diffusion problem
%
% x = 0  ___________________________________________ T(0,t) = 1
%         standard diffusion according to D1  
%         T_1(x,0) = 0
%
% x = L1 ___________________________________________ thermal flux conserved
%         standard diffusion according to D2 
%         T_2(x,0) = 0
%
%
% x = L2 ___________________________________________ heat insulated
%
% FUNCTION INPUTS:
% s - Laplace variable
% x - Depth
% D1 - Diffusivity in skin layer
% D2 - Diffusivity in fat layer
% L1 - Depth of layer interface
% L2 - Depth of entire system
% FUNCTION OUTPUT:
% TBar - Laplace solution

% Establish intermediate parameters
xi1 = sqrt(s/D1);
xi2 = sqrt(s/D2);

if x < L1
    % Skin layer solution (Refer to Equation (3.39)_
    TBar = cosh(xi1*(x-L1))/(s*cosh(xi1*L1)) - ...
        (D2*xi2*sinh(xi1*x)*tanh(xi2*(L2-L1))) / ...
        (s*cosh(xi1*L1)^2*(D2*xi2*tanh(xi1*L1)*tanh(xi2*(L2-L1))+D1*xi1));
else
    % Fat layer solution (Refer to Equation (3.40))
   TBar = (D1*xi1*cosh(xi2*(x-L2))*tanh(xi2*(L2-L1))) / ...
       (s*cosh(xi1*L1)*sinh(xi2*(L2-L1))*(D2*xi2*tanh(xi1*L1)*...
       tanh(xi2*(L2-L1))+D1*xi1));
end

end