function TestParameterSpace( p , q )
%% Test Parameter Space
% Author: Sean McInerney; Last Update: 18/06/2018
%
% Used to generate Figures 2.3, 2.4 and 2.6
%
% Is computational costly for high values of Npts
%
% FUNCTION INPUTS:
% p - The first probe location
% q - The second probe location. Input 'none', if no second probe is used

%% Set-up

close all

% Target Parameters
D1_hat = 0.09;
D2_hat = 0.009;

% Spatial Structure
L1 = 1.6;
L2 = 4;

% Set to 'yes' if you would like to save data
saveData = 'yes';

% Define bounds for parameter space
minD1 = D1_hat/20;
maxD1 = D1_hat*5;
minD2 = D2_hat/20;
maxD2 = D2_hat*5;

% Npts determines how fine the mesh is. Mesh is size Npts by Npts
Npts = 21; % Npts = 2001 was selected in the thesis

% Thresholds for indicator function
epsilon = [1.5,1.0,0.5];

% Colour map below. Currently set up for three indicator thresholds
map = [
    0,204/255,102/255; % Green
    1,0,0; % Red
    0,102/255,204/255 % Blue
    ];

%% Background Set-up

if strcmp(q,'none')
    x_data = p;
    nx = 1; % Number of probes
else
    x_data = [p,q];
    nx = 2;
end

% Constructing the mesh in both vector and matrix form
D1_pts = linspace(minD1,maxD1,Npts);
D2_pts = linspace(minD2,maxD2,Npts);
[D1_mesh,D2_mesh] = meshgrid(D1_pts,D2_pts);
D1_test = D1_mesh(:);
D2_test = D2_mesh(:);

%% Generate Data

% Find appropriate point to stop temporal sampling
thresh = 0.99; % Once this temperature is reached, sampling stops
fun = @(t) Model_Solve(L2,t,D1_hat,D2_hat,L1,L2) - thresh;
t_end = fzero(fun,[0,10^10*L2^2/D2_hat]); % Initial guess is scaled accordingly
M = 100; % Number of temporal samples
delta_t = t_end/M;
t_data = delta_t:delta_t:t_end;

T_hat = Model_Solve(x_data,t_data,D1_hat,D2_hat,L1,L2);

%% Evaluate Indicator Function over Parameter Space

% Set up to evaluate indicator function at probes
Nsamples = Npts^2;
Neps = length(epsilon); % Number of different thresholds considered
I_eps = zeros(Nsamples,Neps); % Contains a collection of I for different epsilon

% Plot_heights is used to generate the figure. It keeps track of whether
% the indicator function is satisfied for a given parameter pair. In the
% case that a parameter pairs results in I = 1 for more than one epsilon,
% the strictest such epsilon is used.
Plot_heights = NaN(Nsamples,1);
epsilon = sort(epsilon,'descend');

% Loop over all the different parameter pairs
for ii = 1:Nsamples
    T_test = Model_Solve(x_data,t_data,D1_test(ii),D2_test(ii),L1,L2);
    % Loop over the different epsilon
    for jj = 1:length(epsilon)
        I_eps(ii,jj) = IndicatorFunc(epsilon(jj),T_test,T_hat);
        if I_eps(ii,jj) == 1
            Plot_heights(ii) = epsilon(jj);
        end
    end
end

%% Create a Figure

% Unravel the Plot_heights into a mesh for plotting
Plot_mesh = zeros(Npts,Npts);
for i = 1:Npts
    Plot_mesh(1:Npts,i) = Plot_heights((i-1)*Npts+1:i*Npts);
end

% Plot
figure
hold on
surf(D1_mesh,D2_mesh,Plot_mesh,'EdgeColor','none')

% Plot details
view(2)
colormap(map)
grid off
xlabel('$D_1$ (mm$^2/$s)','Interpreter','Latex')
ylabel('$D_2$ (mm$^2/$s)','Interpreter','Latex')
plot3(D1_hat,D2_hat,epsilon(Neps)+1,'ko','MarkerFaceColor',[0,0,0])
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','Latex')
axis([minD1,maxD1,minD2,maxD2])
box on
set(gca,'Layer','top')

% Save data
if strcmp(saveData,'yes')
    Results = Plot_mesh;
    if nx == 1
        save(['ResultsForSingleProbeAt',num2str(10*p)],'Results')
    else
        save(['ResultsForTwoProbesAt',num2str(10*p),'and',num2str(10*q)],'Results')
    end
end

function I = IndicatorFunc(epsilon,T_test,T_hat)
% Author: Sean McInerney; Last Update: 28/02/2018
% Evaluates the Indicator Function.
%
% Refer to Equations (2.20) and (2.22)
% 
% FUNCTION INPUTS:
% epsilon - Threshold for indicator function
% T_hat - The temperature data, evaluated using target parameters
% T_test - The temperature, evaluated using test parameters
% FUNCTION OUTPUT:
% I - The evaluated indicator function: either 1 or 0

d = sum(abs(T_hat-T_test),2); % Refer to Equaion (2.19)

[nx,~] = size(T_hat); % Number of probes
if nnz(d <= epsilon) == nx % If all discrepancies are below epsilon
    I = 1;
else
    I = 0;
end